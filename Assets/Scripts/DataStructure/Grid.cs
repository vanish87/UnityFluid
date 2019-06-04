using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Assertions;

namespace UnityFluid
{
    public interface Grid2D
    {
        Vector2 Origin { get; set; }
        Vector2 CellSize { get; set; }//spacing size of one cell
    }
    public struct Grid2DConfig : Grid2D
    {
        public Vector2 Resolution { get; set; }
        public Vector2 Origin { get; set; }
        public Vector2 CellSize { get; set; }
    }
    public abstract class ScalarGrid2D : ScaleField2D, Grid2D
    {
        public Vector2 Origin { get; set; }
        public Vector2 CellSize { get; set; }
        public ScalarGrid2D(Grid2DConfig config): base(config.Resolution)
        {
            this.Origin = config.Origin;
            this.CellSize = config.CellSize;
        }

        protected abstract Vector2 dataSize();
        protected abstract Vector2 dataOrigin();

        public override Vector2 SampleGrandient(Vector2 position)
        {
            Vector2 frac = Vector2.zero;
            Vector2Int index = Vector2Int.zero;

            Helper.GetIndexAndFraction(position, this.Origin, this.CellSize, Vector2Int.zero, new Vector2Int(data.Length - 1, data.Length - 1), out index, out frac);

            var f00 = this.GetGrandientFromIndex(index.x, index.y);
            var f10 = this.GetGrandientFromIndex(index.x + 1, index.y);
            var f01 = this.GetGrandientFromIndex(index.x, index.y + 1);
            var f11 = this.GetGrandientFromIndex(index.x + 1, index.y + 1);
            return Helper.BiLerp(f00, f10, f01, f11, frac.x, frac.y);
        }

        public override Vector2 GetGrandientFromIndex(params int[] list)
        {
            Assert.IsTrue(list.Length == 2);
            var index = new Vector2Int(list[0], list[1]);
            var left    = this.GetDataFromIndex(index.x - 1 , index.y);
            var right   = this.GetDataFromIndex(index.x + 1 , index.y);
            var up      = this.GetDataFromIndex(index.x     , index.y - 1);
            var down    = this.GetDataFromIndex(index.x + 1 , index.y + 1);

            return 0.5f * new Vector2(right - left, up - down) / this.CellSize;
        }
    }

    public abstract class VectorGrid2D : VectorField2D, Grid2D
    {
        public Vector2 Origin { get; set; }
        public Vector2 CellSize { get; set; }
        public VectorGrid2D(Grid2DConfig config) : base(config.Resolution)
        {
            this.Origin = config.Origin;
            this.CellSize = config.CellSize;
        }

    }

    public abstract class CollocatedVectorGrid2D : VectorGrid2D
    {
        public CollocatedVectorGrid2D(Grid2DConfig config) : base(config)
        {
        }

        protected abstract Vector2 dataSize();
        protected abstract Vector2 dataOrigin();
    }

    public class CellCenteredScalarGrid2D : ScalarGrid2D
    {
        public CellCenteredScalarGrid2D(Grid2DConfig config) : base(config)
        {
        }

        protected override Vector2 dataOrigin()
        {
            return this.Origin + (0.5f * this.CellSize);
        }

        protected override Vector2 dataSize()
        {
            return this.Resolution;
        }
    }

    public class VertexCenteredScalarGrid2D : ScalarGrid2D
    {
        public VertexCenteredScalarGrid2D(Grid2DConfig config) : base(config)
        {
        }

        protected override Vector2 dataOrigin()
        {
            return this.Origin;
        }

        protected override Vector2 dataSize()
        {
            return this.Resolution + new Vector2(1, 1);
        }
    }

    public class CellCenteredCollocatedVectorGrid2D : CollocatedVectorGrid2D
    {
        public CellCenteredCollocatedVectorGrid2D(Grid2DConfig config) : base(config)
        {
        }

        protected override Vector2 dataOrigin()
        {
            return this.Origin + (0.5f * this.CellSize);
        }

        protected override Vector2 dataSize()
        {
            return this.Resolution;
        }
    }

    public class VertexCenteredCollocatedVectorGrid2D : CollocatedVectorGrid2D
    {
        public VertexCenteredCollocatedVectorGrid2D(Grid2DConfig config) : base(config)
        {
        }

        protected override Vector2 dataOrigin()
        {
            return this.Origin + (0.5f * this.CellSize);
        }

        protected override Vector2 dataSize()
        {
            return this.Resolution;
        }
    }


    /// <summary>
    /// MAC Grid
    /// Note that this.data variable is not used in this kind of grid
    /// all data is stored at uData, vData, wData
    /// </summary>
    public class FaceCenteredGrid2D : VectorGrid2D
    {
        protected Vector2[] uData, vData;
        protected Vector2 uDataOrigin, vDataOrigin;

        public FaceCenteredGrid2D(Grid2DConfig config): base(config)
        {
            this.data = null;

            //TODO not correct
            var totalCount = Mathf.CeilToInt(this.Resolution.x * this.Resolution.y);
            this.uData = new Vector2[totalCount];
            this.vData = new Vector2[totalCount];

            this.uDataOrigin = this.Origin;
            this.vDataOrigin = this.Origin;
        }
    }

    public class NullGrid : Grid2D
    {
        static NullGrid grid = new NullGrid();
        public static NullGrid NullInstance() { return grid; }
        public Vector2 Origin { get => throw new System.NotImplementedException(); set => throw new System.NotImplementedException(); }
        public Vector2 CellSize { get => throw new System.NotImplementedException(); set => throw new System.NotImplementedException(); }
    }
}