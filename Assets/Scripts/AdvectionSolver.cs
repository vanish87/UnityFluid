using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace UnityFluid
{
    public class AdvectionSolver
    {
        /*public virtual void Advect(VectorField2D input, VectorField2D flow, float delta, out VectorField2D output)
        {
            output = default;
        }*/
    }

    public class SemiLagrangian : AdvectionSolver
    {


        /*protected Vector2 BackTrace(Vector2 pos, float deltaTime, VectorField2D flow)
        {
            var samplePos0 = flow.Sample(pos);
            var mid = pos - 0.5f * (new Vector2(samplePos0.x, samplePos0.y) * deltaTime);

            var samplePos1 = flow.Sample(mid);
            var pos1 = pos - (new Vector2(samplePos1.x, samplePos1.y) * deltaTime);

            //TODO clamp boundary here

            return pos1;
        }

        public override void Advect(VectorField2D input, VectorField2D flow, float delta, out VectorField2D output)
        {
            //TODO
            output = default;
            //foreach(var d in input)
            {
                Vector2 d = Vector2.zero;
                var backTracedPos = this.BackTrace(d, delta, flow);
                var newVal = input.Sample(backTracedPos);
                //output.set(i, j, newVal)
            }
        }*/
    }
}