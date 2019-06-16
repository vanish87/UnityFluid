using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace UnityFluid
{
    public class PICParticles
    {
        PICGrid grid;
        int np; // number of particles
        public List<Vector2> x = new List<Vector2>();
        public List<Vector2> u = new List<Vector2>(); // positions and velocities
                                 /* TODO: add helper variables */
                                 //std::vector<Vec2f> cx, cy; // c vectors stored, times h

        // transfer stuff
        Array2Df sum;

        public PICParticles(PICGrid grid)
        {
            this.grid = grid;
            np = 0;
            sum = new Array2Df(grid.pressure.nx + 1, grid.pressure.ny + 1);
        }

        public void AddParticle(Vector2 px, Vector2 pu)
        {
            x.Add(px);
            u.Add(pu);

            ++np;
        }
        void accumulate(Array2Df accum, float q, int i, int j, float fx, float fy)
        {
            float weight;

            weight = (1 - fx) * (1 - fy);
            accum[i, j] += weight * q;
            sum[i, j] += weight;

            weight = fx * (1 - fy);
            accum[i + 1, j] += weight * q;
            sum[i + 1, j] += weight;

            weight = (1 - fx) * fy;
            accum[i, j + 1] += weight * q;
            sum[i, j + 1] += weight;

            weight = fx * fy;
            accum[i + 1, j + 1] += weight * q;
            sum[i + 1, j + 1] += weight;
        }

        public void TransferToGrid()
        {
            int p = 0, i = 0, ui = 0, j = 0, vj = 0;
            float fx = 0, ufx = 0, fy = 0, vfy = 0;

            grid.u.Reset();
            sum.Reset();
            for (p = 0; p < np; ++p)
            {
                grid.bary_x(x[p][0], ref ui, ref ufx);
                grid.bary_y_centre(x[p][1], ref j, ref fy);
                accumulate(grid.u, u[p][0], ui, j, ufx, fy);
                /* TODO: call affineFix to incorporate c_px^n into the grid.u update */
                // affineFix([FILL THIS IN]);
            }
            for (j = 0; j < grid.u.ny; ++j) for (i = 0; i < grid.u.nx; ++i)
                {
                    if (sum[i, j] != 0) grid.u[i, j] /= sum[i, j];
                }

            grid.v.Reset();
            sum.Reset();
            for (p = 0; p < np; ++p)
            {
                grid.bary_x_centre(x[p][0], ref i, ref fx);
                grid.bary_y(x[p][1], ref vj, ref vfy);
                accumulate(grid.v, u[p][1], i, vj, fx, vfy);
                /* TODO: call affineFix to incorporate c_py^n into the grid.v update */
                // affineFix([FILL THIS IN]);
            }
            for (j = 0; j < grid.v.ny; ++j) for (i = 0; i < grid.v.nx; ++i)
                {
                    if (sum[i, j] != 0) grid.v[i, j] /= sum[i, j];
                }

            // identify where particles are in grid
            grid.marker.Reset();
            for (p = 0; p < np; ++p)
            {
                grid.bary_x(x[p][0], ref i, ref fx);
                grid.bary_y(x[p][1], ref j, ref fy);
                grid.marker[i, j] = PICGrid.FLUID;
            }
        }
        public void update_from_grid()
        {
            int p;
            int i = 0, ui = 0, j = 0, vj = 0;
            float fx = 0, ufx = 0, fy = 0, vfy = 0;
            for (p = 0; p < np; ++p)
            {
                grid.bary_x(x[p][0], ref ui, ref ufx);
                grid.bary_x_centre(x[p][0], ref i, ref fx);
                grid.bary_y(x[p][1], ref vj, ref vfy);
                grid.bary_y_centre(x[p][1], ref j, ref fy);
                u[p] = new Vector2(grid.u.BiLerp(ui, j, ufx, fy), grid.v.BiLerp(i, vj, fx, vfy)); // PIC and APIC

            }
        }
        public void move_particles_in_grid(float dt)
        {
            Vector2 midx, gu = default;
            float xmin = 1.001f * grid.h, xmax = grid.cellx - 1.001f * grid.h;
            float ymin = 1.001f * grid.h, ymax = grid.celly - 1.001f * grid.h;
            for (int p = 0; p < np; ++p)
            {
                // first stage of Runge-Kutta 2 (do a half Euler step)
                var gu0 = gu.x;
                var gu1 = gu.y;
                grid.bilerp_uv(x[p][0], x[p][1], ref gu0, ref gu1);
                gu = new Vector2(gu0, gu1);
                midx = x[p] + 0.5f * dt * gu;
                midx[0] = Mathf.Clamp(midx[0], xmin, xmax);
                midx[1] = Mathf.Clamp(midx[1], ymin, ymax);
                // second stage of Runge-Kutta 2
                gu0 = gu.x;
                gu1 = gu.y;
                grid.bilerp_uv(midx[0], midx[1], ref gu0, ref gu1);
                gu = new Vector2(gu0, gu1);
                x[p] += dt * gu;

                var xp0 = Mathf.Clamp(x[p][0], xmin, xmax);
                var xp1 = Mathf.Clamp(x[p][1], ymin, ymax);

                x[p] = new Vector2(xp0, xp1);
            }
        }
    }
}