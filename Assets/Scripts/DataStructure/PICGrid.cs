using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace UnityFluid
{
    public class PICGrid
    {
        public const int AIR = 0;
        public const int FLUID = 1;
        public const int SOLID = 2;

        public float gravity;
        public float cellx, celly;
        public float h, invH;

        public Array2Df u, v;
        public Array2Df du, dv;
        public Array2Di marker;
        public Array2Df phi;

        public Array2Df pressure;
        public Array2Dx3f poisson;

        public Array2Df preconditioner;
        public Array2Df m;
        public Array2Df r, z, s;
        float sqr(float x)
        {
            return x * x;
        }

        public void bary_x(float x, ref int i, ref float fx)
        {
            //i is the cell index
            //fx is the fractional part of sx
            float sx = x * invH;
            i = (int)sx;
            fx = sx - Mathf.Floor(sx);
        }

        public void bary_x_centre(float x, ref int i, ref float fx)
        {
            float sx = x * invH - 0.5f;
            i = (int)sx;
            if (i < 0) { i = 0; fx = 0.0f; }
            else if (i > pressure.nx - 2) { i = pressure.nx - 2; fx = 1.0f; }
            else { fx = sx - Mathf.Floor(sx); }
        }

        public void bary_y(float y, ref int j, ref float fy)
        {
            float sy = y * invH;
            j = (int)sy;
            fy = sy - Mathf.Floor(sy);
        }

        public void bary_y_centre(float y, ref int j, ref float fy)
        {
            float sy = y * invH - 0.5f;//this overh will scale pos in [0,1]to grid size [0,50]
            j = (int)sy;
            if (j < 0) { j = 0; fy = 0.0f; }
            else if (j > pressure.ny - 2) { j = pressure.ny - 2; fy = 1.0f; }
            else { fy = sy - Mathf.Floor(sy); }
        }

        public void bilerp_uv(float px, float py, ref float pu, ref float pv)
        {
            int i = 0, j = 0;
            float fx = 0, fy = 0;
            bary_x(px, ref i, ref fx);
            bary_y_centre(py, ref j, ref fy);
            pu = u.BiLerp(i, j, fx, fy);
            bary_x_centre(px, ref i, ref fx);
            bary_y(py, ref j, ref fy);
            pv = v.BiLerp(i, j, fx, fy);
        }

        public PICGrid(float gravity, int nx, int ny, float cellsize)
        {
            this.gravity = gravity;

            cellx = cellsize;
            celly = ny * cellx / nx;
            h = cellx / nx;
            invH = 1 / h;

            u = new Array2Df(nx + 1, ny);
            v = new Array2Df(nx, ny + 1);
            pressure = new Array2Df(nx, ny);
            marker = new Array2Di(nx, ny);
            phi = new Array2Df(nx, ny);
            du = new Array2Df(nx, ny);
            dv = new Array2Df(nx, ny);
            poisson = new Array2Dx3f(nx, ny);
            preconditioner = new Array2Df(nx, ny);
            m = new Array2Df(nx, ny);
            r = new Array2Df(nx, ny);
            z = new Array2Df(nx, ny);
            s = new Array2Df(nx, ny);
        }

        public float GetCFL()
        {
            float maxv2 = Mathf.Max(h * gravity, sqr(u.InfNorm()) + sqr(v.InfNorm()));
            if (maxv2 < 1e-16) maxv2 = 1e-16f;
            return h / Mathf.Sqrt(maxv2);
        }

        public void SolveGravity(float timeDelta)
        {
            for (int i = 0; i < v.size; ++i)
            {
                v.data[i] -= (-gravity * timeDelta);
            }
        }

        public void BuildSDF()
        {
            this.InitPhi();

            for (var i = 0; i < 2; ++i)
            {
                this.SweepPhi();
            }
        }

        protected void InitPhi()
        {
            int i, j;
            // start off with indicator inside the fluid and overestimates of distance outside
            float large_distance = phi.nx + phi.ny + 2;
            for (i = 0; i < phi.size; ++i)
                phi.data[i] = large_distance;
            for (j = 1; j < phi.ny - 1; ++j) for (i = 1; i < phi.nx - 1; ++i)
                {
                    if (marker[i, j] == FLUID)
                    {
                        //negative value means inside the surface
                        phi[i, j] = -0.5f;
                    }
                }
        }
        public void SolveDistance(float p, float q, ref float r)
        {
            float d = Mathf.Min(p, q) + 1;
            if (d > Mathf.Max(p, q))
                d = (p + q + Mathf.Sqrt(2 - sqr(p - q))) / 2;
            if (d < r) r = d;
        }

        protected void SweepPhi()
        {
            // fast sweeping outside the fluid in all four sweep directions
            //Fluid Simulation for Computer Graphics P51 for more info
            int i, j;
            for (j = 1; j < phi.ny; ++j) for (i = 1; i < phi.nx; ++i)
                    if (marker[i, j] != FLUID)
                    {
                        var phiValue = phi[i, j];
                        SolveDistance(phi[i - 1, j], phi[i, j - 1], ref phiValue);
                        phi[i, j] = phiValue;
                    }
            for (j = phi.ny - 2; j >= 0; --j) for (i = 1; i < phi.nx; ++i)
                    if (marker[i, j] != FLUID)
                    {
                        var phiValue = phi[i, j];
                        SolveDistance(phi[i - 1, j], phi[i, j + 1], ref phiValue);
                        phi[i, j] = phiValue;
                    }
            for (j = 1; j < phi.ny; ++j) for (i = phi.nx - 2; i >= 0; --i)
                    if (marker[i, j] != FLUID)
                    {
                        var phiValue = phi[i, j];
                        SolveDistance(phi[i + 1, j], phi[i, j - 1], ref phiValue);
                        phi[i, j] = phiValue;
                    }
            for (j = phi.ny - 2; j >= 0; --j) for (i = phi.nx - 2; i >= 0; --i)
                    if (marker[i, j] != FLUID)
                    {
                        var phiValue = phi[i, j];
                        SolveDistance(phi[i + 1, j], phi[i, j + 1], ref phiValue);
                        phi[i, j] = phiValue;
                    }
        }

        protected void SweepU(int i0, int i1, int j0, int j1)
        {
            int di = (i0 < i1) ? 1 : -1, dj = (j0 < j1) ? 1 : -1;
            float dp, dq, alpha;
            for (int j = j0; j != j1; j += dj) for (int i = i0; i != i1; i += di)
                    if (marker[i - 1, j] == AIR && marker[i, j] == AIR)
                    {
                        dp = di * (phi[i, j] - phi[i - 1, j]);
                        if (dp < 0) continue; // not useful on this sweep direction
                        dq = 0.5f * (phi[i - 1, j] + phi[i, j] - phi[i - 1, j - dj] - phi[i, j - dj]);
                        if (dq < 0) continue; // not useful on this sweep direction
                        if (dp + dq == 0) alpha = 0.5f;
                        else alpha = dp / (dp + dq);
                        u[i, j] = alpha * u[i - di, j] + (1 - alpha) * u[i, j - dj];
                    }
        }

        protected void SweepV(int i0, int i1, int j0, int j1)
        {
            int di = (i0 < i1) ? 1 : -1, dj = (j0 < j1) ? 1 : -1;
            float dp, dq, alpha;
            for (int j = j0; j != j1; j += dj) for (int i = i0; i != i1; i += di)
                    if (marker[i, j - 1] == AIR && marker[i, j] == AIR)
                    {
                        dq = dj * (phi[i, j] - phi[i, j - 1]);
                        if (dq < 0) continue; // not useful on this sweep direction
                        dp = 0.5f * (phi[i, j - 1] + phi[i, j] - phi[i - di, j - 1] - phi[i - di, j]);
                        if (dp < 0) continue; // not useful on this sweep direction
                        if (dp + dq == 0) alpha = 0.5f;
                        else alpha = dp / (dp + dq);
                        v[i, j] = alpha * v[i - di, j] + (1 - alpha) * v[i, j - dj];
                    }
        }

        public void ExternVelocity()
        {
            for (int i = 0; i < 4; ++i)
                SweepVelocity();
        }

        public void SweepVelocity()
        {
            int i, j;
            // sweep u, only into the air
            SweepU(1, u.nx - 1, 1, u.ny - 1);
            SweepU(1, u.nx - 1, u.ny - 2, 0);
            SweepU(u.nx - 2, 0, 1, u.ny - 1);
            SweepU(u.nx - 2, 0, u.ny - 2, 0);
            for (i = 0; i < u.nx; ++i)
            {
                u[i, 0] = u[i, 1]; u[i, u.ny - 1] = u[i, u.ny - 2];
            }
            for (j = 0; j < u.ny; ++j)
            {
                u[0, j] = u[1, j]; u[u.nx - 1, j] = u[u.nx - 2, j];
            }
            // now the same for v
            SweepV(1, v.nx - 1, 1, v.ny - 1);
            SweepV(1, v.nx - 1, v.ny - 2, 0);
            SweepV(v.nx - 2, 0, 1, v.ny - 1);
            SweepV(v.nx - 2, 0, v.ny - 2, 0);
            for (i = 0; i < v.nx; ++i)
            {
                v[i, 0] = v[i, 1]; v[i, v.ny - 1] = v[i, v.ny - 2];
            }
            for (j = 0; j < v.ny; ++j)
            {
                v[0, j] = v[1, j]; v[v.nx - 1, j] = v[v.nx - 2, j];
            }
        }

        public void FindDivergence()
        {
            r.Reset();
            for (int j = 0; j < r.ny; ++j) for (int i = 0; i < r.nx; ++i)
                {
                    if (marker[i, j] == FLUID)
                        r[i, j] = u[i + 1, j] - u[i, j] + v[i, j + 1] - v[i, j];
                }
        }
        public void FormPoisson()
        {
            poisson.Reset();
            for (int j = 1; j < poisson.ny - 1; ++j) for (int i = 1; i < poisson.nx - 1; ++i)
                {
                    if (marker[i, j] == FLUID)
                    {
                        if (marker[i - 1, j] != SOLID)
                            poisson[i, j, 0] += 1;
                        if (marker[i + 1, j] != SOLID)
                        {
                            poisson[i, j, 0] += 1;
                            if (marker[i + 1, j] == FLUID)
                                poisson[i, j, 1] = -1;
                        }
                        if (marker[i, j - 1] != SOLID)
                            poisson[i, j, 0] += 1;
                        if (marker[i, j + 1] != SOLID)
                        {
                            poisson[i, j, 0] += 1;
                            if (marker[i, j + 1] == FLUID)
                                poisson[i, j, 2] = -1;
                        }
                    }
                }
        }
        public void ApplyPoisson(Array2Df x, Array2Df y)
        {
            y.Reset();
            for (int j = 1; j < poisson.ny - 1; ++j) for (int i = 1; i < poisson.nx - 1; ++i)
                {
                    if (marker[i, j] == FLUID)
                    {
                        y[i, j] = poisson[i, j, 0] * x[i, j] + poisson[i - 1, j, 1] * x[i - 1, j]
                                                     + poisson[i, j, 1] * x[i + 1, j]
                                                     + poisson[i, j - 1, 2] * x[i, j - 1]
                                                     + poisson[i, j, 2] * x[i, j + 1];
                    }
                }
        }
        void FormPreconditioner()
        {
            const float mic_parameter = 0.99f;
            float d;
            preconditioner.Reset();
            for (int j = 1; j < preconditioner.ny - 1; ++j) for (int i = 1; i < preconditioner.nx - 1; ++i)
                {
                    if (marker[i, j] == FLUID)
                    {
                        d = poisson[i, j, 0] - sqr(poisson[i - 1, j, 1] * preconditioner[i - 1, j])
                                         - sqr(poisson[i, j - 1, 2] * preconditioner[i, j - 1])
                                         - mic_parameter * (poisson[i - 1, j, 1] * poisson[i - 1, j, 2] * sqr(preconditioner[i - 1, j])
                                                          + poisson[i, j - 1, 2] * poisson[i, j - 1, 1] * sqr(preconditioner[i, j - 1]));
                        preconditioner[i, j] = 1 / Mathf.Sqrt(d + 1e-6f);
                    }
                }
        }
        void ApplyPreconditioner(Array2Df x, Array2Df y, Array2Df m)
        {
            //Maybe Incomplete Cholesky??
            int i, j;
            float d;
            m.Reset();
            // solve L*m=x
            for (j = 1; j < x.ny - 1; ++j) for (i = 1; i < x.nx - 1; ++i)
                    if (marker[i, j] == FLUID)
                    {
                        d = x[i, j] - poisson[i - 1, j, 1] * preconditioner[i - 1, j] * m[i - 1, j]
                                 - poisson[i, j - 1, 2] * preconditioner[i, j - 1] * m[i, j - 1];
                        m[i, j] = preconditioner[i, j] * d;
                    }
            // solve L'*y=m
            y.Reset();
            for (j = x.ny - 2; j > 0; --j) for (i = x.nx - 2; i > 0; --i)
                    if (marker[i, j] == FLUID)
                    {
                        d = m[i, j] - poisson[i, j, 1] * preconditioner[i, j] * y[i + 1, j]
                                 - poisson[i, j, 2] * preconditioner[i, j] * y[i, j + 1];
                        y[i, j] = preconditioner[i, j] * d;
                    }
        }
        public void ApplySolidSDF()
        {
            int i, j;
            // first mark where solid is
            for (j = 0; j < marker.ny; ++j)
                marker[0, j] = marker[marker.nx - 1, j] = SOLID;
            for (i = 0; i < marker.nx; ++i)
            {
                marker[i, 0] = marker[i, marker.ny - 1] = SOLID;
               /* if (i > marker.nx / 2)
                {
                    //for (var k = i; k < marker.nx; ++k)
                    {
                        marker[i, marker.ny / 2] = SOLID;
                    }
                }*/
            }
            // now make sure nothing leaves the domain
            for (j = 0; j < u.ny; ++j)
            {
                u[0, j] = u[1, j] = u[u.nx - 1, j] = u[u.nx - 2, j] = 0;

                /*if (j == u.ny / 2 || j == u.ny / 2 -1)
                {
                    for (var k = j; k < u.nx; ++k)
                    {
                        u[k, u.ny / 2] = 0;
                    }
                }*/
            }
            for (i = 0; i < v.nx; ++i)
            {
                v[i, 0] = v[i, 1] = v[i, v.ny - 1] = v[i, v.ny - 2] = 0;

                /*if (i >  v.nx / 2 - 1)
                {
                    //for (var k = i; k < v.ny; ++k)
                    {
                        //v[i, v.ny / 2 - 1] = 0;
                        v[i, v.ny / 2] = 0;
                    }
                }*/
            }
        }
        public void SolveIncompressible()
        {
            FindDivergence();
            FormPoisson();
            FormPreconditioner();
            SolvePressure(100, 1e-5f);
            AddGradient();
        }
        public void SolvePressure(int maxits, float tolerance)
        {
            int its;
            float tol = tolerance * r.InfNorm();
            pressure.Reset();
            if (r.InfNorm() == 0)
                return;
            ApplyPreconditioner(r, z, m);
            z.CopyTo(s);
            float rho = z.Dot(r);
            if (rho == 0)
                return;
            for (its = 0; its < maxits; ++its)
            {
                ApplyPoisson(s, z);
                float alpha = rho / s.Dot(z);
                pressure.Increment(alpha, s);
                r.Increment(-alpha, z);
                if (r.InfNorm() <= tol)
                {
                    Debug.LogFormat("pressure converged to {0} in {1} iterations\n", r.InfNorm(), its);
                    return;
                }
                ApplyPreconditioner(r, z, m);
                float rhonew = z.Dot(r);
                float beta = rhonew / rho;
                s.ScaleAndIncrement(beta, z);
                rho = rhonew;
            }
            Debug.LogFormat("Didn't converge in pressure solve (its={0}, tol={1}, |r|={2})\n", its, tol, r.InfNorm());
        }
        public void AddGradient()
        {
            int i, j;
            for (j = 1; j < u.ny - 1; ++j) for (i = 2; i < u.nx - 2; ++i)
                {
                    if (!(marker[i - 1, j] != FLUID && marker[i, j] != FLUID))
                    { // if at least one is FLUID, neither is SOLID
                        u[i, j] += pressure[i, j] - pressure[i - 1, j];
                    }
                }
            for (j = 2; j < v.ny - 2; ++j) for (i = 1; i < v.nx - 1; ++i)
                {
                    if (!(marker[i, j - 1] != FLUID && marker[i, j] != FLUID))
                    { // if at least one is FLUID, neither is SOLID
                        v[i, j] += pressure[i, j] - pressure[i, j - 1];
                    }
                }
        }
    }

    public class Array2D<T>
    {
        public int nx, ny;
        public int size;

        public T[] data;

        public Array2D(int nx, int ny)
        {
            this.nx = nx;
            this.ny = ny;

            this.size = this.nx * this.ny;

            this.data = new T[this.size];
        }

        public void Reset(T value = default)
        {
            for (var i = 0; i < this.size; ++i)
            {
                this.data[i] = value;
            }
        }
        public T this[int i, int j]
        {
            get { return this.data[i + nx * j]; }
            set { this.data[i + nx * j] = value; }
        }

        public void CopyTo(Array2D<T> target)
        {
            target.data = (T[])this.data.Clone();
        }
        public virtual T BiLerp(int i, int j, float fx, float fy)
        {
            return default;
        }
        public virtual T InfNorm()
        {
            return default;
        }

        public virtual T Dot(Array2D<T> other)
        {
            return default;
        }

        public virtual void Increment(float scale, Array2D<T> other)
        {
        }
        public virtual void ScaleAndIncrement(float scale, Array2D<T> other)
        {
        }
    }
    public class Array2Di : Array2D<int>
    {
        public Array2Di(int nx, int ny) : base(nx, ny)
        {
        }
    }

    public class Array2Df : Array2D<float>
    {
        public Array2Df(int nx, int ny) : base(nx, ny)
        {
        }

        public override float BiLerp(int i, int j, float fx, float fy)
        {
            var f00 = this[i, j];
            var f10 = this[i + 1, j];
            var f01 = this[i, j + 1];
            var f11 = this[i + 1, j + 1];

            return FluidHelper.BiLerp(f00, f10, f01, f11, fx, fy);
        }

        public override float InfNorm()
        {
            float r = 0;
            for (int i = 0; i < size; ++i)
                if (!(Mathf.Abs(data[i]) <= r)) r = Mathf.Abs(data[i]);
            return r;
        }

        public override void Increment(float scale, Array2D<float> other)
        {
            for (int i = 0; i < size; ++i) data[i] += scale * other.data[i];
        }

        public override void ScaleAndIncrement(float scale, Array2D<float> other)
        {
            for (int i = 0; i < size; ++i) data[i] = scale * data[i] + other.data[i];
        }

        public override float Dot(Array2D<float> other)
        {
            float r = 0;
            for (int i = 0; i < size; ++i)
                r += data[i] * other.data[i];
            return r;
        }
    }

    public class Array2Dx3f
    {
        public int nx, ny, nz;
        public int size;

        public float[] data;

        public Array2Dx3f(int nx, int ny)
        {
            this.nx = nx;
            this.ny = ny;

            this.size = this.nx * this.ny * 3;

            this.data = new float[this.size];
        }

        public void Reset(float value = default)
        {
            for (var i = 0; i < this.size; ++i)
            {
                this.data[i] = value;
            }
        }
        public float this[int i, int j, int k]
        {
            get { return this.data[(i + nx * j) * 3 + k]; }
            set { this.data[(i + nx * j) * 3 + k] = value; }
        }

        public float InfNorm()
        {
            float r = 0;
            for (int i = 0; i < size; ++i)
                if (!(Mathf.Abs(data[i]) <= r)) r = Mathf.Abs(data[i]);
            return r;
        }
    }
}