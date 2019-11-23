using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace UnityFluid
{

    public class PICMain : GPUParticleStructBase<ParticleData>
    {

        [SerializeField] ComputeShader solver;
        [SerializeField] ComputeShader dataUpdate;
        public enum SimType
        {
            PIC,
            APIC,
        }
        public enum LinearSolver
        {
            Jacobi,
            PCG,
        }
        public enum Method
        {
            CPU,
            GPU
        }

        public SimType type = SimType.PIC;
        public LinearSolver linearSolver = LinearSolver.PCG;
        public Method method = Method.CPU;

        PICGrid grid;
        PICParticles particles;

        Vector2Int gridSize = new Vector2Int(50, 50);
        protected override void Start()
        {
            this.grid = new PICGrid(-9.8f, this.gridSize.x, this.gridSize.y, 1);
            this.particles = new PICParticles(grid, this.type);

            this.AddParticle();

            this.parameter.numberOfParticles.Value = this.gridSize.x * this.gridSize.y * 4;
            this.ResizeBuffer(this.parameter.numberOfParticles.Value);

            this.grid.solver = this.solver;
            this.grid.dataUpdate = this.dataUpdate;
        }

        protected override void Update()
        {
            if (Input.GetKey(KeyCode.Space))
            {
                this.AdvanceOneFrame(this.grid, this.particles, 1 / 30f);
                this.CopyDataToGPU();
            }
            if(Input.GetKeyDown(KeyCode.I))
            {
                this.AddDrop();
            }
        }

        void AdvanceOneFrame(PICGrid grid, PICParticles particles, float frametime)
        {
            float t = 0;
            float dt;
            bool finished = false;
            while (!finished)
            {
                dt = 2 * grid.GetCFL();
                if (t + dt >= frametime)
                {
                    dt = frametime - t;
                    finished = true;
                }
                else if (t + 1.5 * dt >= frametime)
                    dt = 0.5f * (frametime - t);
                Debug.LogFormat("advancing {0} (to {1}% of frame)\n", dt, 100.0 * (t + dt) / frametime);
                AdvanceOneStep(grid, particles, dt);
                t += dt;
            }
        }

        void AdvanceOneStep(PICGrid grid, PICParticles particles, float dt)
        {
            particles.TransferToGrid();
            //grid.save_velocities();
            grid.SolveGravity(dt);
            grid.BuildSDF();
            grid.ExternVelocity();
            grid.ApplySolidSDF();
            grid.SolveIncompressible(this.method == Method.CPU, this.linearSolver == LinearSolver.Jacobi);
            grid.ExternVelocity();
            //grid.get_velocity_update();
            particles.GridToParticle();

            for (int i = 0; i < 5; ++i)
                particles.MoveParticle(0.2f * dt);
        }
        protected float FluidPhi(float x, float y)
        {
            //return a circle with radius r and center(0.5*grid.x,0.5*grid.y)
            //and a bottom water that y = 0.2 * gridSize.y
            var r = 0.3f;
            var center = new Vector2(0.5f, 0.7f);
            var dx = x - center.x;
            var dy = y - center.y;
            return Mathf.Min(Mathf.Sqrt(dx * dx + dy * dy) - r, y - 0.2f);
        }
        protected float FluidDrop(float x, float y)
        {
            //return a circle with radius r and center(0.5*grid.x,0.5*grid.y)
            //and a bottom water that y = 0.2 * gridSize.y
            var r = 0.1f;
            var center = new Vector2(0.7f, 0.7f);
            var dx = x - center.x;
            var dy = y - center.y;
            return Mathf.Sqrt(dx * dx + dy * dy) - r;
        }
        protected void AddDrop(int na = 2, int nb = 2)
        {
            int i, j, a, b;
            float x, y, phi;
            float vx, vy;
            vx = 0;
            vy = 0;

            for (i = 1; i < grid.marker.nx - 1; ++i)
            {
                for (j = 1; j < grid.marker.ny - 1; ++j)
                {
                    for (a = 0; a < na; ++a)
                    {
                        for (b = 0; b < nb; ++b)
                        {
                            //x and y is random pos in a every cell
                            //they are distributed on 4 sub-cell of main cell
                            x = (i + (a + 0.1f + 0.8f * Random.value) / na) * grid.h;  //=>grid.h will map grid coord to domain coord
                            y = (j + (b + 0.1f + 0.8f * Random.value) / nb) * grid.h;
                            phi = FluidDrop(x, y);
                            if (phi > -0.25 * grid.h / na)
                                continue;
                            /*else if (phi > -1.5 * grid.h / na)
                            {
                                //this will make a smooth edge of phi broader
                                //should not effect to velocity field
                                FluidPhi(grid, x, y, phi, -0.75 * grid.h / na);
                                phi = FluidPhi(grid, x, y);
                                project(grid, x, y, phi, -0.75 * grid.h / na);
                                phi = fluidphi(grid, x, y);
                            }*/
                            if (particles.x.Count < this.parameter.numberOfParticles.Value - 1)
                            {
                                particles.AddParticle(new Vector2(x, y), new Vector2(vx, vy));
                            }
                        }
                    }
                }
            }
        }
        protected void AddParticle(int na = 2, int nb = 2)
        {
            int i, j, a, b;
            float x, y, phi;
            float vx, vy;
            vx = 0;
            vy = 0;

            for (i = 1; i < grid.marker.nx - 1; ++i)
            {
                for (j = 1; j < grid.marker.ny - 1; ++j)
                {
                    for (a = 0; a < na; ++a)
                    {
                        for (b = 0; b < nb; ++b)
                        {
                            //x and y is random pos in a every cell
                            //they are distributed on 4 sub-cell of main cell
                            x = (i + (a + 0.1f + 0.8f * Random.value) / na) * grid.h;  //=>grid.h will map grid coord to domain coord
                            y = (j + (b + 0.1f + 0.8f * Random.value) / nb) * grid.h;
                            phi = FluidPhi(x, y);
                            if (phi > -0.25 * grid.h / na)
                                continue;
                            /*else if (phi > -1.5 * grid.h / na)
                            {
                                //this will make a smooth edge of phi broader
                                //should not effect to velocity field
                                FluidPhi(grid, x, y, phi, -0.75 * grid.h / na);
                                phi = FluidPhi(grid, x, y);
                                project(grid, x, y, phi, -0.75 * grid.h / na);
                                phi = fluidphi(grid, x, y);
                            }*/
                            particles.AddParticle(new Vector2(x, y), new Vector2(vx, vy));
                        }
                    }
                }
            }
        }
        protected void CopyDataToGPU()
        {
            for(var i = 0; i < this.CPUData.Length; ++i)
            {
                if (i < this.particles.x.Count)
                {
                    var p = this.particles.x[i] * this.gridSize;
                    this.CPUData[i].position = this.transform.position + new Vector3(p.x, p.y, 0);
                    this.CPUData[i].color = this.type == SimType.PIC ? Color.green : Color.red;
                }
                else
                {
                    this.CPUData[i].color = Color.clear;
                }
            }
            this.UpdateGPUDataBuffer();
        }

        protected void OnDrawGizmos()
        {
            if (this.grid == null) return;
            var gridDimension = new Vector2Int(this.grid.pressure.nx, this.grid.pressure.ny);
            if (this.grid.b != null) Gizmos.DrawGUITexture(new Rect(0, 0, gridDimension.x, gridDimension.y), this.grid.b);
            if (this.grid.x != null) Gizmos.DrawGUITexture(new Rect(0, gridDimension.y, gridDimension.x, gridDimension.y), this.grid.x);
            if (this.grid.markerTex != null) Gizmos.DrawGUITexture(new Rect(0, gridDimension.y * 2, gridDimension.x, gridDimension.y), this.grid.markerTex);
            //if (this.difference != null) Gizmos.DrawGUITexture(new Rect(0, gridDimension.y * 2, gridDimension.x, gridDimension.y), this.difference);
        }
    }
}