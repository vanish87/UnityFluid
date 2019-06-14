using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace UnityFluid
{
    public class PICGrid
    {
        public Vector2Int Size { get; }
        protected float[] uVelocity;
        protected float[] vVelocity;

        public PICGrid(Vector2Int size)
        {
            this.Size = size;

            var dataSize = this.Size.x * this.Size.y;

            this.uVelocity = new float[dataSize + 1];
            this.vVelocity = new float[dataSize + 1];
        }


    }
}