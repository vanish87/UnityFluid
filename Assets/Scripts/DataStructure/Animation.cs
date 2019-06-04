using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace UnityFluid
{
    public class Animation : MonoBehaviour
    {
        //SimulationData data;

        protected virtual void onAdvanceTimeStep(float delta)
        {
            //begin

            //force
            //time integration
            //resolve collision

            //end
        }

        protected void FixedUpdate()
        {
            var delta = Time.fixedDeltaTime;
            this.onAdvanceTimeStep(delta);
        }
    }
}