using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace UnityFluid
{
    public class GridSolver : Animation
    {
        GridSystemData data = new GridSystemData();
        protected override void onAdvanceTimeStep(float delta)
        {
            //begin

            //force
            //time integration
            //resolve collision

            ComputeExternalForces(delta);
            ComputeViscosity(delta);
            ComputePressure(delta);
            ComputeAdvection(delta);

            //end
        }

        protected virtual void ComputeExternalForces(float delta)
        {
            //
            //velocity(i, j, k) += timeIntervalInSeconds * _gravity

        }
        protected virtual void ComputeViscosity(float delta)
        {

        }
        protected virtual void ComputePressure(float delta)
        {

        }
        protected virtual void ComputeAdvection(float delta)
        {
            var solver = new SemiLagrangian();
            var input = this.data.GetData(GridSystemData.DataType.Pressure) as VectorField2D;

            solver.Advect(input, input, delta, out input);
        }
        protected void Start()
        {
            //Init data
        }

        protected void InitData()
        {

        }
    }

    public class PICSolver : GridSolver
    {
        protected override void onAdvanceTimeStep(float delta)
        {
            TransferParticleToGrid();
            //buildMarkers();
            //extrapolateVelocityToAir();
            //applyBoundaryCondition();


            ComputeExternalForces(delta);
            ComputeViscosity(delta);
            ComputePressure(delta);

            //advection step
            //extrapolateVelocityToAir();
            //applyBoundaryCondition();
            //transferFromGridsToParticles();
            MoveParticle(delta);
        }
        protected void TransferParticleToGrid()
        {

        }

        protected void TransferGridToParticle()
        {

        }

        protected void MoveParticle(float delta)
        {

        }
    }
}