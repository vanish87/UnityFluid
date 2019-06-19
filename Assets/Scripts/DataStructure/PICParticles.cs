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
        public List<Vector2> cx = new List<Vector2>();
        public List<Vector2> cy = new List<Vector2>();

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

            /* TODO: initialize the variables you created in particles.h */
            cx.Add(Vector2.zero);
            cy.Add(Vector2.zero);

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
                affineFix(grid.u, cx[p], ui, j, ufx, fy);
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
                affineFix(grid.v, cy[p], i, vj, fx, vfy);
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
        public void GridToParticle()
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

                /* TODO: call computeC with the right indices to compute c_px^n and c_py^n */
                cx[p] = computeC(grid.u, ui, j, ufx, fy); // APIC
                cy[p] = computeC(grid.v, i, vj, fx, vfy); // APIC
            }
        }

        /* this function computes c from the gradient of w and the velocity field from the grid. */
        Vector2 computeC(Array2Df ufield, int i, int j, float fx, float fy)
        {
            var f00 = ufield[i, j];
            var f10 = ufield[i + 1, j];
            var f01 = ufield[i, j + 1];
            var f11 = ufield[i + 1, j + 1];

            var w00 = new Vector2((1 - fx), (1 - fy));
            var w10 = new Vector2(fx, (1 - fy));
            var w01 = new Vector2((1 - fx), fy);
            var w11 = new Vector2(fx, fy);

            var gw00 = new Vector2(w10.x - w00.x, w01.y - w00.y) * grid.h;
            var gw10 = new Vector2(w10.x - w00.x, w11.y - w10.y) * grid.h;
            var gw01 = new Vector2(w11.x - w01.x, w01.y - w00.y) * grid.h;
            var gw11 = new Vector2(w11.x - w01.x, w11.y - w10.y) * grid.h;
            /* TODO: fill this in */
            return gw00 * f00 * w00.x * w00.y + gw10 * f10 * w10.x * w10.y + 
                   gw01 * f01 * w01.x * w01.y + gw11 * f11 * w11.x * w11.y;
        }

        /* call this function to incorporate c[] when transfering particles to grid */
        /* This function should take the c_pa^n values from c, and update them, with proper weighting, */
        /*  into the correct grid velocity values in accum */
        void affineFix(Array2Df accum, Vector2 c, int i, int j, float fx, float fy)
        {

            accum[i, j] += c.x * (1 - fx) + c.y * (1 - fy);

            accum[i + 1, j] += c.x * fx + c.y * (1 - fy);

            accum[i, j + 1] += c.x * (1 - fx) + c.y * fy;

            accum[i + 1, j + 1] += c.x * fx + c.y * fy;
        }

        public void MoveParticle(float dt)
        {
            Vector2 midx, gu = default;
            float xmin = 1.001f * grid.h, xmax = grid.domainx - 1.001f * grid.h;
            float ymin = 1.001f * grid.h, ymax = grid.domainy - 1.001f * grid.h;
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