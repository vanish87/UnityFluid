using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace UnityFluid
{
    public class Matrix2D<T>
    {
        protected T[,] data;
        protected Texture2D gpuData;

        protected Vector2Int dimension;
        public Matrix2D(Vector2Int dim)
        {
            this.dimension = dim;
            this.data = new T[dimension.x, dimension.y];
            this.gpuData = new Texture2D(dimension.x, dimension.y, TextureFormat.RFloat, false);
        }

        void test()
        {

        }
    }
}
