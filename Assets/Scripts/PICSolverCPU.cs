using System.Collections;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using UnityEngine;
using UnityEngine.Assertions;

namespace FluidData
{
    [StructLayout(LayoutKind.Sequential, Size = 48)]
    public struct ParticleData
    {
        public BlittableBool active;
        public Vector3 position;

        public Vector3 velocity;
        public float life;

        public Color color;
    }
    public class PICSolverCPU : GPUParticleStructBase<ParticleData>
    {
        protected override void Start()
        {
            this.parameter.numberOfParticles.Value = this.gridSize.x * this.gridSize.y * 4;
            this.ResizeBuffer(this.parameter.numberOfParticles.Value);

            this.InitData();
            this.InitParticle();
        }

        protected override void Update()
        {
            if(Input.GetKeyDown(KeyCode.Space))
            {
                //this.AdvanceSetp(Time.fixedDeltaTime);
                this.AdvanceFrame();
                this.CopyDataToGPU();

                Debug.LogWarning("Done" + Time.frameCount);
            }
        }

        protected float GetCFL()
        {
            var h = cellSpace / this.gridSize.x;
            var norm = FluidHelper.Infnorm(this.velocity);
            float maxv2 = Mathf.Max(h * Gravaty, sqr(norm.x) + sqr(norm.y));
            if (maxv2 < 1e-16) maxv2 = 1e-16f;
            return h / Mathf.Sqrt(maxv2);
        }
        protected void AdvanceFrame()
        {
            var max = 50;
            var count = 0;

            var t = 0f;
            var dt = 0f;
            var finished = false;
            var frametime = Time.fixedDeltaTime;
            while (!finished && count++ < max)
            {
                dt = 2 * GetCFL();
                if (t + dt >= frametime)
                {
                    dt = frametime - t;
                    finished = true;
                }
                else if (t + 1.5f * dt >= frametime)
                    dt = 0.5f * (frametime - t);
                Debug.LogFormat("advancing {0} (to {1} of frame)\n", dt, 100.0 * (t + dt) / frametime);
                this.AdvanceSetp(dt);
                t += dt;
            }
        }

        protected void AdvanceSetp(float delta)
        {
            this.ParticleToGrid();
            this.SaveVelocity();
            this.SolveGravity(delta);
            this.GridToParticle();

            for (int i = 0; i < 5; ++i)
                this.MoveParticle(0.2f * delta);
            return;

            this.BuildSDF();
            this.ExtrapolateVelocityToAir();
            this.AddSolidBoundary();
            this.ProjectPressure();
            this.ExtrapolateVelocityToAir();

            //grid.get_velocity_update();
            //particles.update_from_grid();

            this.GridToParticle();


            for (int i = 0; i < 5; ++i)
                this.MoveParticle(0.2f * delta);
        }

        protected void CopyDataToGPU()
        {
            this.UpdateGPUDataBuffer();
        }

        #region Fuild

        #region Data
        protected const float cellSpace = 1;
        protected const float invSpace = 1.0f / cellSpace;
        protected const float Gravaty = 9.8f;

        protected const int AIR = 0;
        protected const int FLUID = 1;
        protected const int SOLID = 2;

        protected Vector2Int gridSize = new Vector2Int(50, 50);

        protected GridFactory factory = new GridFactory();

        protected FaceCenterdVectorGrid2D velocity;
        protected FaceCenterdVectorGrid2D savedVelocity;
        //protected FaceCenteredGrid2D weightSum;

        protected CellCenteredScalarGrid2D pressure;
        protected CellCenteredScalarGrid2Di marker;
        protected CellCenteredScalarGrid2D phi;

        protected CellCenteredVectorGrid2D poisson;
        protected CellCenteredScalarGrid2D preconditioner;

        protected CellCenteredScalarGrid2D mField;
        protected CellCenteredScalarGrid2D r;
        protected CellCenteredScalarGrid2D z;
        protected CellCenteredScalarGrid2D s;

        #endregion
        protected void InitData()
        {
            Grid2DConfigure config = new Grid2DConfigure();
            config.CellSize = new Vector2(cellSpace, cellSpace);
            config.Resolution = this.gridSize;

            this.velocity = factory.MakeGrid2D(GridFactory.CenterType.FaceCentered, GridFactory.DataType.Vector, config) as FaceCenterdVectorGrid2D;
            this.savedVelocity = factory.MakeGrid2D(GridFactory.CenterType.FaceCentered, GridFactory.DataType.Vector, config) as FaceCenterdVectorGrid2D;
            //this.weightSum = factory.MakeGrid2D(GridFactory.CenterType.FaceCentered, GridFactory.DataType.Vector, config) as FaceCenteredGrid2D;

            this.pressure = factory.MakeGrid2D(GridFactory.CenterType.CellCentered, GridFactory.DataType.Scalar, config) as CellCenteredScalarGrid2D;
            this.marker = factory.MakeGrid2Di(GridFactory.CenterType.CellCentered, GridFactory.DataType.Scalar, config) as CellCenteredScalarGrid2Di;
            this.phi = factory.MakeGrid2D(GridFactory.CenterType.CellCentered, GridFactory.DataType.Scalar, config) as CellCenteredScalarGrid2D;

            this.poisson = factory.MakeGrid2D(GridFactory.CenterType.CellCentered, GridFactory.DataType.Vector, config) as CellCenteredVectorGrid2D;
            this.preconditioner = factory.MakeGrid2D(GridFactory.CenterType.CellCentered, GridFactory.DataType.Scalar, config) as CellCenteredScalarGrid2D;

            this.mField = factory.MakeGrid2D(GridFactory.CenterType.CellCentered, GridFactory.DataType.Scalar, config) as CellCenteredScalarGrid2D;
            this.r = factory.MakeGrid2D(GridFactory.CenterType.CellCentered, GridFactory.DataType.Scalar, config) as CellCenteredScalarGrid2D;
            this.z = factory.MakeGrid2D(GridFactory.CenterType.CellCentered, GridFactory.DataType.Scalar, config) as CellCenteredScalarGrid2D;
            this.s = factory.MakeGrid2D(GridFactory.CenterType.CellCentered, GridFactory.DataType.Scalar, config) as CellCenteredScalarGrid2D;

        }

        protected float FluidPhi(float x, float y, Vector2 domainSize)
        {
            //return a circle with radius r and center(0.5*grid.x,0.5*grid.y)
            //and a bottom water that y = 0.2 * gridSize.y
            var r = 0.3f * domainSize.x;
            var center = new Vector2(0.5f, 0.5f) * domainSize;
            var dx = x - center.x;
            var dy = y - center.y;
            return Mathf.Min(Mathf.Sqrt(dx * dx + dy * dy) - r, y - 0.2f * domainSize.y);
        }

        protected void InitParticle()
        {
            var count = 0;
            for(var i = 0; i < this.gridSize.x; ++i)
            {
                for (var j = 0; j < this.gridSize.y; ++j)
                {
                    for(var ncx = 0; ncx < 2; ++ncx)
                    {
                        for(var ncy = 0; ncy < 2; ++ncy)
                        {
                            var x = (i + (ncx + 0.8f * Random.value)/2 );
                            var y = (j + (ncy + 0.8f * Random.value)/2 );

                            var phi = this.FluidPhi(x, y, this.gridSize);
                            if(phi > 0 || count > this.CPUData.Length-1)
                            {
                                continue;
                            }
                            else
                            {
                                this.CPUData[count].active = true;
                                this.CPUData[count].position.x = x;
                                this.CPUData[count].position.y = y;
                                this.CPUData[count++].color = Color.red;
                            }
                        }
                    }
                }
            }
        }
        protected void ParticleToGrid()
        {
            this.velocity.ResetValueAndWeight();
            this.marker.Reset(AIR);
            foreach (var p in this.CPUData)
            {
                if (p.active == false) continue;

                var pos = p.position;
                var vel = p.velocity;

                this.velocity.AccumulatePoint(pos, vel);             

                Vector2Int posIndex;
                Vector2 posFrac;
                FluidHelper.GetIndexAndFraction(pos, this.marker.Origin, this.marker.CellSize, Vector2Int.zero, this.marker.DataSize, out posIndex, out posFrac);
                this.marker[posIndex.x, posIndex.y] = FLUID;
            }

            this.velocity.NormalizeWeight();
        }

        void SaveVelocity()
        {
            this.velocity.CopyTo(this.savedVelocity);
        }
        protected void Advection()
        {

        }
        protected void BuildSDF()
        {
            this.InitPhi();

            for (var i = 0; i < 2; ++i)
            {
                this.SweepPhi();
            }
        }
        protected void ExtrapolateVelocityToAir()
        {
            for (var i = 0; i < 4; ++i)
            {
                this.SweepVelocity();
            }
        }

        protected void AddSolidBoundary()
        {
            int i, j;
            var markernx = this.gridSize.x;
            var markerny = this.gridSize.y;

            var unx = this.velocity.uDataSize.x;
            var uny = this.velocity.uDataSize.y;
            var vnx = this.velocity.vDataSize.x;
            var vny = this.velocity.vDataSize.y;
            // first mark where solid is
            for (j = 0; j < markerny; ++j)
            {
                this.marker.SetDataToIndex(SOLID, 0, j);
                this.marker.SetDataToIndex(SOLID, markernx - 1, j);
            }
            for (i = 0; i < markernx; ++i)
            {
                this.marker.SetDataToIndex(SOLID, i, 0);
                this.marker.SetDataToIndex(SOLID, i, markerny - 1);
            }
            // now make sure nothing leaves the domain
            for (j = 0; j < uny; ++j)
            {
                this.velocity[MACGrid2DData.DataType.U, 0, j] = 0f;
                this.velocity[MACGrid2DData.DataType.U, 1, j] = 0f;
                this.velocity[MACGrid2DData.DataType.U, unx - 1, j] = 0f;
                this.velocity[MACGrid2DData.DataType.U, unx - 2, j] = 0f;
            }
            for (i = 0; i < vnx; ++i)
            {
                this.velocity[MACGrid2DData.DataType.V, i, 0] = 0f;
                this.velocity[MACGrid2DData.DataType.V, i, 1] = 0f;
                this.velocity[MACGrid2DData.DataType.V, i, vny - 1] = 0f;
                this.velocity[MACGrid2DData.DataType.V, i, vny - 2] = 0f;

            }
        }
        protected void SolveGravity(float timeDelta)
        {
            this.velocity.ForEachvData((ref float value, int i, int j) =>
            {
                value += (-Gravaty * timeDelta);
            });
        }
        protected void SolveViscosity()
        {
            //Skip for now
        }
        protected void ProjectPressure()
        {
            FindDivergence();
            FormPoisson();
            FormPreconditioner();
            SolvePressure(100, 1e-5f);
            AddGradient();
        }
        protected void GridToParticle()
        {
            var sampler = new MacVectorGrid2DSampler();
            for (var i = 0; i < this.CPUData.Length; ++i)
            {
                var pos = this.CPUData[i].position;
                var vel = sampler.Sample(this.velocity, pos);
                this.CPUData[i].velocity = vel;
            }
        }

        protected void MoveParticle(float delta)
        {
            var sampler = new MacVectorGrid2DSampler();
            for (var i = 0; i < this.CPUData.Length; ++i)
            {
                var pos = this.CPUData[i].position;
                //var vel = this.velocity.Sample(pos);
                var vel = this.CPUData[i].velocity;

                var mid = pos + 0.5f * delta * vel;
                mid.x = Mathf.Clamp(mid.x, 0, this.gridSize.x);
                mid.y = Mathf.Clamp(mid.y, 0, this.gridSize.y);

                vel = sampler.Sample(this.velocity, mid);
                pos += delta * vel;
                pos.x = Mathf.Clamp(pos.x, 0, this.gridSize.x);
                pos.y = Mathf.Clamp(pos.y, 0, this.gridSize.y);
                this.CPUData[i].position = pos;
            }
        }

        #region Helper Function
        protected void InitPhi()
        {
            int i, j;
            var size = this.phi.DataSize;
            // start off with indicator inside the fluid and overestimates of distance outside
            float large_distance = size.x * size.y + 2;

            this.phi.ForEachData((ref float value, int[] list) =>
            {
                value = large_distance;
            });

            for (j = 1; j < size.y - 1; ++j)
                for (i = 1; i < size.x - 1; ++i)
                {
                    if (this.marker[i, j] == FLUID)
                    {
                        //negative value means inside the surface
                        this.phi[i, j] = -0.5f;
                    }
                }
        }
        protected void SweepPhi()
        {
            int i, j;
            int nx = this.phi.DataSize.x;
            int ny = this.phi.DataSize.y;

            for (j = 1; j < ny; ++j)
                for (i = 1; i < nx; ++i)
                {
                    var markerVal = Mathf.RoundToInt(this.marker.GetDataFromIndex(i, j));
                    if (markerVal != FLUID)
                    {
                        var phiR = this.phi.GetDataFromIndex(i, j);
                        SolveDistance(this.phi.GetDataFromIndex(i - 1, j), this.phi.GetDataFromIndex(i, j - 1), ref phiR);
                        this.phi.SetDataToIndex(phiR, i, j);
                    }
                }
            for (j = ny - 2; j >= 0; --j)
                for (i = 1; i < nx; ++i)
                {
                    var markerVal = Mathf.RoundToInt(this.marker.GetDataFromIndex(i, j));
                    if (markerVal != FLUID)
                    {
                        var phiR = this.phi.GetDataFromIndex(i, j);
                        SolveDistance(this.phi.GetDataFromIndex(i - 1, j), this.phi.GetDataFromIndex(i, j + 1), ref phiR);
                        this.phi.SetDataToIndex(phiR, i, j);
                    }
                }
            for (j = 1; j < ny; ++j)
                for (i = nx - 2; i >= 0; --i)
                {
                    var markerVal = Mathf.RoundToInt(this.marker.GetDataFromIndex(i, j));
                    if (markerVal != FLUID)
                    {
                        var phiR = this.phi.GetDataFromIndex(i, j);
                        SolveDistance(this.phi.GetDataFromIndex(i + 1, j), this.phi.GetDataFromIndex(i, j - 1), ref phiR);
                        this.phi.SetDataToIndex(phiR, i, j);
                    }
                }
            for (j = ny - 2; j >= 0; --j)
                for (i = nx - 2; i >= 0; --i)
                {
                    var markerVal = Mathf.RoundToInt(this.marker.GetDataFromIndex(i, j));
                    if (markerVal != FLUID)
                    {
                        var phiR = this.phi.GetDataFromIndex(i, j);
                        SolveDistance(this.phi.GetDataFromIndex(i + 1, j), this.phi.GetDataFromIndex(i, j + 1), ref phiR);
                        this.phi.SetDataToIndex(phiR, i, j);
                    }
                }
        }

        protected void SolveDistance(float p, float q, ref float r)
        {
            float d = Mathf.Min(p, q) + 1;
            if (d > Mathf.Max(p, q))
                d = (p + q + Mathf.Sqrt(2 - sqr(p - q))) / 2;
            if (d < r) r = d;
        }

        protected void SweepVelocity()
        {
            int i, j;
            int unx, uny;
            int vnx, vny;

            /*//Debug.LogWarning("Check all size");

            unx = this.velocity.uDataSize.x;
            uny = this.velocity.uDataSize.y;

            vnx = this.velocity.vDataSize.x;
            vny = this.velocity.vDataSize.y;
            // sweep u, only into the air
            SweepU(1, unx - 1, 1, uny - 1);
            SweepU(1, unx - 1, uny - 2, 0);
            SweepU(unx - 2, 0, 1, uny - 1);
            SweepU(unx - 2, 0, uny - 2, 0);
            for (i = 0; i < unx; ++i)
            {
                this.velocity.SetuDataToIndex(this.velocity.GetuDataFromIndex(i, 1), i, 0);
                this.velocity.SetuDataToIndex(this.velocity.GetuDataFromIndex(i, uny - 2), i, uny - 1);
            }
            for (j = 0; j < uny; ++j)
            {
                this.velocity.SetuDataToIndex(this.velocity.GetuDataFromIndex(1, j), 0, j );
                this.velocity.SetuDataToIndex(this.velocity.GetuDataFromIndex(unx - 2, j) ,unx - 1, j);
            }
            // now the same for v
            SweepV(1, vnx - 1, 1, vny - 1);
            SweepV(1, vnx - 1, vny - 2, 0);
            SweepV(vnx - 2, 0, 1, vny - 1);
            SweepV(vnx - 2, 0, vny - 2, 0);
            for (i = 0; i < vnx; ++i)
            {
                this.velocity.SetvDataToIndex(this.velocity.GetvDataFromIndex(i, 1), i, 0);
                this.velocity.SetvDataToIndex(this.velocity.GetvDataFromIndex(i, vny - 2), i, vny - 1);
            }
            for (j = 0; j < vny; ++j)
            {
                this.velocity.SetvDataToIndex(this.velocity.GetvDataFromIndex(1, j), 0, j);
                this.velocity.SetvDataToIndex(this.velocity.GetvDataFromIndex(vnx - 2, j), vnx - 1, j);
            }*/
        }

        protected void SweepU(int i0, int i1, int j0, int j1)
        {
            /*int di = (i0 < i1) ? 1 : -1, dj = (j0 < j1) ? 1 : -1;
            float dp, dq, alpha;
            for (int j = j0; j != j1; j += dj) for (int i = i0; i != i1; i += di)
                    if (Mathf.RoundToInt(this.marker.GetDataFromIndex(i - 1, j)) ==  AIR && Mathf.RoundToInt(this.marker.GetDataFromIndex(i, j)) == AIR)
                    {
                        dp = di * (this.phi.GetDataFromIndex(i, j) - this.phi.GetDataFromIndex(i - 1, j));
                        if (dp < 0) continue; // not useful on this sweep direction
                        dq = 0.5f * (this.phi.GetDataFromIndex(i - 1, j) + this.phi.GetDataFromIndex(i, j) - this.phi.GetDataFromIndex(i - 1, j - dj) - this.phi.GetDataFromIndex(i, j - dj));
                        if (dq < 0) continue; // not useful on this sweep direction
                        if (dp + dq == 0) alpha = 0.5f;
                        else alpha = dp / (dp + dq);

                        var uData = alpha * this.velocity.GetuDataFromIndex(i - di, j) + (1 - alpha) * this.velocity.GetuDataFromIndex(i, j - dj);
                        this.velocity.SetuDataToIndex(uData, i, j);
                    }*/
        }

        protected void SweepV(int i0, int i1, int j0, int j1)
        {
            /*int di = (i0 < i1) ? 1 : -1, dj = (j0 < j1) ? 1 : -1;
            float dp, dq, alpha;
            for (int j = j0; j != j1; j += dj) for (int i = i0; i != i1; i += di)
                    if (Mathf.RoundToInt(this.marker.GetDataFromIndex(i, j - 1)) == AIR && Mathf.RoundToInt(this.marker.GetDataFromIndex(i, j)) == AIR)
                    {
                        dq = dj * (this.phi.GetDataFromIndex(i, j) - this.phi.GetDataFromIndex(i, j - 1));
                        if (dq < 0) continue; // not useful on this sweep direction
                        dp = 0.5f * (this.phi.GetDataFromIndex(i, j - 1) + this.phi.GetDataFromIndex(i, j) - this.phi.GetDataFromIndex(i - di, j - 1) - this.phi.GetDataFromIndex(i - di, j));
                        if (dp < 0) continue; // not useful on this sweep direction
                        if (dp + dq == 0) alpha = 0.5f;
                        else alpha = dp / (dp + dq);
                        var vData = alpha * this.velocity.GetvDataFromIndex(i - di, j) + (1 - alpha) * this.velocity.GetvDataFromIndex(i, j - dj);
                        this.velocity.SetvDataToIndex(vData, i, j);
                    }*/
        }

        void FindDivergence()
        {
            r.Reset(0);
            /*var rnx = this.r.DataSize.x;
            var rny = this.r.DataSize.y;
            for (int j = 0; j < rny; ++j)
                for (int i = 0; i < rnx; ++i)
                {
                    if (Mathf.RoundToInt(this.marker.GetDataFromIndex(i, j)) == FLUID)
                    {
                        var rData = this.velocity.GetuDataFromIndex(i + 1, j) - this.velocity.GetuDataFromIndex(i, j) 
                            + this.velocity.GetvDataFromIndex(i, j + 1) - this.velocity.GetvDataFromIndex(i, j);
                        r.SetDataToIndex(rData, i, j);
                    }
                }*/
        }
        void FormPoisson()
        {
            poisson.Reset(Vector3.zero);
            var poissonnx = this.poisson.DataSize.x;
            var poissonny = this.poisson.DataSize.y;
            for (int j = 1; j < poissonny - 1; ++j) for (int i = 1; i < poissonnx - 1; ++i)
                {
                    if (Mathf.RoundToInt(this.marker.GetDataFromIndex(i, j)) == FLUID)
                    {
                        if (Mathf.RoundToInt(this.marker.GetDataFromIndex(i - 1, j)) != SOLID)
                        {
                            var data = this.poisson.GetDataFromIndex(i, j);
                            data.x += 1;
                            this.poisson.SetDataToIndex(data, i, j);
                        }
                        if (Mathf.RoundToInt(this.marker.GetDataFromIndex(i + 1, j)) != SOLID)
                        {
                            var data = this.poisson.GetDataFromIndex(i, j);
                            data.x += 1;
                            if (Mathf.RoundToInt(this.marker.GetDataFromIndex(i + 1, j)) == FLUID)
                            {
                                data.y = -1;
                            }
                            this.poisson.SetDataToIndex(data, i, j);
                        }
                        if (Mathf.RoundToInt(this.marker.GetDataFromIndex(i, j - 1)) != SOLID)
                        {
                            var data = this.poisson.GetDataFromIndex(i, j);
                            data.x += 1;
                            this.poisson.SetDataToIndex(data, i, j);
                        }
                        if ((this.marker.GetDataFromIndex(i, j + 1)) != SOLID)
                        {
                            var data = this.poisson.GetDataFromIndex(i, j);
                            data.x += 1;
                            if (Mathf.RoundToInt(this.marker.GetDataFromIndex(i, j + 1)) == FLUID)
                            {
                                data.z = -1;
                            }
                            this.poisson.SetDataToIndex(data, i, j);
                        }
                    }
                }
        }
        void ApplyPoisson(CellCenteredScalarGrid2D x, CellCenteredScalarGrid2D y)
        {
            y.Reset(0);
            var poissonnx = this.poisson.DataSize.x;
            var poissonny = this.poisson.DataSize.y;
            for (int j = 1; j < poissonny - 1; ++j)
                for (int i = 1; i < poissonnx - 1; ++i)
                {
                    if (Mathf.RoundToInt(this.marker.GetDataFromIndex(i, j)) == FLUID)
                    {
                        var value = poisson.GetDataFromIndex(i,j).x * x.GetDataFromIndex(i,j)
                            + poisson.GetDataFromIndex(i - 1, j).y * x.GetDataFromIndex(i - 1, j)
                            + poisson.GetDataFromIndex(i, j).y * x.GetDataFromIndex(i + 1, j)
                            + poisson.GetDataFromIndex(i, j - 1).z * x.GetDataFromIndex(i, j - 1)
                            + poisson.GetDataFromIndex(i, j).z * x.GetDataFromIndex(i, j + 1);
                        y.SetDataToIndex(value, i, j);

                    }
                }
        }
        float sqr(float x)
        {
            return x * x;
        }
        void FormPreconditioner()
        {
            const float mic_parameter = 0.99f;
            float d;
            preconditioner.Reset(0);

            var preconditionernx = this.preconditioner.DataSize.x;
            var preconditionerny = this.preconditioner.DataSize.y;

            for (int j = 1; j < preconditionerny - 1; ++j) for (int i = 1; i < preconditionernx - 1; ++i)
                {
                    if (Mathf.RoundToInt(this.marker.GetDataFromIndex(i, j)) == FLUID)
                    {
                        d = poisson.GetDataFromIndex(i, j).x - sqr(poisson.GetDataFromIndex(i - 1, j).y * preconditioner.GetDataFromIndex(i - 1, j))
                                         - sqr(poisson.GetDataFromIndex(i, j - 1).z * preconditioner.GetDataFromIndex(i, j - 1))
                                         - mic_parameter * (poisson.GetDataFromIndex(i - 1, j).y * 
                                                            poisson.GetDataFromIndex(i - 1, j).z * 
                                                            sqr(preconditioner.GetDataFromIndex(i - 1, j))
                                                          + poisson.GetDataFromIndex(i, j - 1).z * 
                                                            poisson.GetDataFromIndex(i, j - 1).y * 
                                                            sqr(preconditioner.GetDataFromIndex(i, j - 1)));
                        preconditioner.SetDataToIndex(1 / Mathf.Sqrt(d + 1e-6f), i, j);
                    }
                }
        }
        void apply_preconditioner(CellCenteredScalarGrid2D x, CellCenteredScalarGrid2D y, CellCenteredScalarGrid2D m)
        {
            //Maybe Incomplete Cholesky??
            int i, j;
            float d;
            m.Reset(0);

            var xnx = x.DataSize.x;
            var xny = x.DataSize.y;
            // solve L*m=x
            for (j = 1; j < xny - 1; ++j)
                for (i = 1; i < xnx - 1; ++i)
                    if (Mathf.RoundToInt(marker.GetDataFromIndex(i, j)) == FLUID)
                    {
                        d = x.GetDataFromIndex(i, j) - poisson.GetDataFromIndex(i - 1, j).y * preconditioner.GetDataFromIndex(i - 1, j) * m.GetDataFromIndex(i - 1, j)
                                 - poisson.GetDataFromIndex(i, j - 1).z * preconditioner.GetDataFromIndex(i, j - 1) * m.GetDataFromIndex(i, j - 1);
                        var value = preconditioner.GetDataFromIndex(i, j) * d;
                        m.SetDataToIndex(value, i, j); 
                    }

            // solve L'*y=m
            y.Reset(0);
            for (j = xny - 2; j > 0; --j)
                for (i = xnx - 2; i > 0; --i)
                    if (Mathf.RoundToInt(marker.GetDataFromIndex(i, j)) == FLUID)
                    {
                        d = m.GetDataFromIndex(i, j) - poisson.GetDataFromIndex(i, j).y * preconditioner.GetDataFromIndex(i, j) * y.GetDataFromIndex(i + 1, j)
                                 - poisson.GetDataFromIndex(i, j).z * preconditioner.GetDataFromIndex(i, j) * y.GetDataFromIndex(i, j + 1);
                        var value = preconditioner.GetDataFromIndex(i, j) * d;
                        y.SetDataToIndex(value, i, j);
                    }
        }

        void SolvePressure(int maxits, float tolerance)
        {
            /*int its;
            var informR = FluidHelper.Infnorm(r);
            float tol = tolerance * informR;
            pressure.Reset(0);
            if (informR == 0)
                return;
            apply_preconditioner(r, z, mField);
            z.CopyTo(s);

            float rho = FluidHelper.Dot(z,r);
            if (rho == 0)
                return;
            for (its = 0; its < maxits; ++its)
            {
                ApplyPoisson(s, z);
                float alpha = rho / FluidHelper.Dot(s,z);
                FluidHelper.Increment(pressure, s, alpha);
                FluidHelper.Increment(r, z, -alpha);
                if (FluidHelper.Infnorm(r) <= tol)
                {
                    Debug.LogFormat("pressure converged to {0} in {1} iterations\n", FluidHelper.Infnorm(r), its);
                    return;
                }
                apply_preconditioner(r, z, mField);
                float rhonew = FluidHelper.Dot(z, r);
                float beta = rhonew / rho;
                FluidHelper.SacleAndIncrement(s, z, beta);
                rho = rhonew;
            }
            Debug.LogFormat("Didn't converge in pressure solve (its={0}, tol={1}, |r|={2})\n", its, tol, FluidHelper.Infnorm(r));*/
        }
        void AddGradient()
        {
            int i, j;

            //Assert.IsTrue(false);
            var unx = this.velocity.uDataSize.x;
            var uny = this.velocity.uDataSize.y;


            var vnx = this.velocity.vDataSize.x;
            var vny = this.velocity.vDataSize.y;

            /*for (j = 1; j < uny - 1; ++j)
                for (i = 2; i < unx - 2; ++i)
                {
                    if (!(Mathf.RoundToInt(marker.GetDataFromIndex(i - 1, j)) != FLUID && Mathf.RoundToInt(marker.GetDataFromIndex(i, j)) != FLUID))
                    { // if at least one is FLUID, neither is SOLID
                        var value = this.velocity.GetuDataFromIndex(i, j);
                        value += pressure.GetDataFromIndex(i, j) - pressure.GetDataFromIndex(i - 1, j);
                        this.velocity.SetuDataToIndex(value, i, j); 
                    }
                }
            for (j = 2; j < vny - 2; ++j)
                for (i = 1; i < vnx - 1; ++i)
                {
                    if (!(Mathf.RoundToInt(marker.GetDataFromIndex(i, j - 1)) != FLUID && Mathf.RoundToInt(marker.GetDataFromIndex(i, j)) != FLUID))
                    { // if at least one is FLUID, neither is SOLID
                        var value = this.velocity.GetvDataFromIndex(i, j);
                        value += pressure.GetDataFromIndex(i, j) - pressure.GetDataFromIndex(i, j - 1);
                        this.velocity.SetvDataToIndex(value, i, j);
                    }
                }*/
        }

        #endregion
        #endregion

        #region Debug
        protected void OnDrawGizmos()
        {
            this.DrawMarker();
            this.DrawVelocity();
        }

        protected void DrawMarker()
        {
            var old = Gizmos.color;
            this.marker?.ForEachData((ref int value, int[] list) =>
            {
                var dataOrg = this.marker.DataOrigin;
                var c = Color.white;
                if (value == FLUID) c = Color.blue;
                else if (value == SOLID) c = Color.red;
                c.a = 0.5f;
                Gizmos.color = c;
                Gizmos.DrawCube(new Vector3(list[0], list[1], 0) + new Vector3(dataOrg.x, dataOrg.y, 0), new Vector3(cellSpace,cellSpace,0.1f) * 0.9f);
            });
            Gizmos.color = old;
        }

        protected void DrawVelocity()
        {
            if (this.velocity == null) return;

            var norm = FluidHelper.Infnorm(this.velocity);
            var scale = Vector2.one / norm * 0.4f * new Vector2(cellSpace, cellSpace);

            var old = Gizmos.color;
            this.velocity?.ForEachuData((ref float value, int i, int j) =>
            {
                var dataOrg = this.velocity.uDataOrigin;
                var from = new Vector3(i, j, 0) + new Vector3(dataOrg.x, dataOrg.y, 0);
                var to = from + new Vector3(value, 0, 0) * scale.x ;

                Gizmos.DrawLine(from, to);
            });

            this.velocity?.ForEachvData((ref float value, int i, int j) =>
            {
                var dataOrg = this.velocity.vDataOrigin;
                var from = new Vector3(i, j, 0) + new Vector3(dataOrg.x, dataOrg.y, 0);
                var to = from + new Vector3(0, value, 0) * scale.y;

                Gizmos.DrawLine(from, to);
            });
            Gizmos.color = old;
        }
        #endregion


    }
}
