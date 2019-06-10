using System;
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

        protected abstract Vector2 dataOrigin();
        

        public override Vector2 GetGrandientFromIndex(params int[] list)
        {
            Assert.IsTrue(list.Length == 2);
            var index = new Vector2Int(list[0], list[1]);
            var left    = this.GetDataFromIndex(index.x - 1 , index.y);
            var right   = this.GetDataFromIndex(index.x + 1 , index.y);
            var up      = this.GetDataFromIndex(index.x     , index.y - 1);
            var down    = this.GetDataFromIndex(index.x     , index.y + 1);

            return 0.5f * new Vector2(right - left, up - down) / this.CellSize;
        }

        public override Vector2 SampleGrandient(Vector2 position)
        {
            Vector2 frac = Vector2.zero;
            Vector2Int index = Vector2Int.zero;

            Helper.GetIndexAndFraction(position, this.Origin, this.CellSize, Vector2Int.zero, new Vector2Int(data.Length - 1, data.Length - 1), out index, out frac);

            var f00 = this.GetGrandientFromIndex(index.x    , index.y);
            var f10 = this.GetGrandientFromIndex(index.x + 1, index.y);
            var f01 = this.GetGrandientFromIndex(index.x    , index.y + 1);
            var f11 = this.GetGrandientFromIndex(index.x + 1, index.y + 1);

            Debug.LogWarning("Verify this");

            return Helper.BiLerp(f00, f10, f01, f11, frac.x, frac.y);
        }

         public override float GetLaplacianFromIndex(params int[] list)
        {
            Assert.IsTrue(list.Length == 2);
            var index = new Vector2Int(list[0], list[1]);
            var center  = this.GetDataFromIndex(index.x     , index.y);
            var left    = this.GetDataFromIndex(index.x - 1 , index.y);
            var right   = this.GetDataFromIndex(index.x + 1 , index.y);
            var up      = this.GetDataFromIndex(index.x     , index.y - 1);
            var down    = this.GetDataFromIndex(index.x     , index.y + 1);

            var spaceSqure = this.CellSize * this.CellSize;

            return (right - 2 * center + left) / spaceSqure.x
                +  (up    - 2 * center + down) / spaceSqure.y;
        }

        public override float SampleLaplacian(Vector2 position)
        {
            Vector2 frac = Vector2.zero;
            Vector2Int index = Vector2Int.zero;

            Helper.GetIndexAndFraction(position, this.Origin, this.CellSize, Vector2Int.zero, new Vector2Int(data.Length - 1, data.Length - 1), out index, out frac);

            var f00 = this.GetLaplacianFromIndex(index.x    , index.y);
            var f10 = this.GetLaplacianFromIndex(index.x + 1, index.y);
            var f01 = this.GetLaplacianFromIndex(index.x    , index.y + 1);
            var f11 = this.GetLaplacianFromIndex(index.x + 1, index.y + 1);

            Debug.LogWarning("Verify this");

            return Helper.BiLerp(f00, f10, f01, f11, frac.x, frac.y);
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
        public override float GetDivergenceFromIndex(params int[] list)
        {
            Assert.IsTrue(list.Length == 2);
            var index = new Vector2Int(list[0], list[1]);
            var left    = this.GetDataFromIndex(index.x - 1 , index.y);
            var right   = this.GetDataFromIndex(index.x + 1 , index.y);
            var up      = this.GetDataFromIndex(index.x     , index.y - 1);
            var down    = this.GetDataFromIndex(index.x     , index.y + 1);

            return 0.5f * (right.x - left.x) / this.CellSize.x
                 + 0.5f * (up.y - down.y)    / this.CellSize.y;
        }
        public override float SampleDivergence(Vector2 position)
        {
            Vector2 frac = Vector2.zero;
            Vector2Int index = Vector2Int.zero;

            Helper.GetIndexAndFraction(position, this.Origin, this.CellSize, Vector2Int.zero, new Vector2Int(data.Length - 1, data.Length - 1), out index, out frac);

            var f00 = this.GetDivergenceFromIndex(index.x    , index.y);
            var f10 = this.GetDivergenceFromIndex(index.x + 1, index.y);
            var f01 = this.GetDivergenceFromIndex(index.x    , index.y + 1);
            var f11 = this.GetDivergenceFromIndex(index.x + 1, index.y + 1);

            Debug.LogWarning("Verify this");

            return Helper.BiLerp(f00, f10, f01, f11, frac.x, frac.y);
        }
        public override Vector3 SampleCurl(Vector2 position)
        {
            throw new NotImplementedException();
        }
    }

    public abstract class CollocatedVectorGrid2D : VectorGrid2D
    {
        public CollocatedVectorGrid2D(Grid2DConfig config) : base(config)
        {
        }
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

        protected override Vector2 DataSize 
        {
            get { return this.Resolution; }
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

        protected override Vector2 DataSize
        {
            get{ return this.Resolution + new Vector2(1, 1); }
        }

        public override void InitData()
        {
            var totalCount = this.DataSize.x * this.DataSize.y;
            this.data = new float[(int)totalCount];
        }
    }

    public class CellCenteredVectorGrid2D : CollocatedVectorGrid2D
    {
        public CellCenteredVectorGrid2D(Grid2DConfig config) : base(config)
        {
        }

        protected override Vector2 dataOrigin()
        {
            return this.Origin + (0.5f * this.CellSize);
        }

        protected override Vector2 DataSize
        {
            get
            {
                return this.Resolution;
            }
        }
    }

    public class VertexCenteredVectorGrid2D : CollocatedVectorGrid2D
    {
        public VertexCenteredVectorGrid2D(Grid2DConfig config) : base(config)
        {
        }

        protected override Vector2 dataOrigin()
        {
            return this.Origin;
        }

        protected override Vector2 DataSize
        {
            get
            {
                return this.Resolution + new Vector2(1, 1);
            }
        }
        public override void InitData()
        {
            var totalCount = this.DataSize.x * this.DataSize.y;
            this.data = new Vector3[(int)totalCount];
        }
    }


    /// <summary>
    /// MAC Grid
    /// Note that this.data variable is not used in this kind of grid
    /// all data is stored at uData, vData, wData
    /// </summary>
    public class FaceCenteredGrid2D : VectorGrid2D
    {
        protected float[] uData, vData;
        protected Vector2 uDataOrigin, vDataOrigin;

        public FaceCenteredGrid2D(Grid2DConfig config): base(config)
        {
        }

        public override void InitData()
        {
            //this.data is not used
            this.data = null;

            var totalU = Mathf.CeilToInt((this.Resolution.x + 1) * this.Resolution.y);
            var totalV = Mathf.CeilToInt(this.Resolution.x * (this.Resolution.y + 1));
            this.uData = new float[totalU];
            this.vData = new float[totalV];

            this.uDataOrigin = this.Origin + new Vector2(this.CellSize.x * 0.5f, 0);
            this.vDataOrigin = this.Origin + new Vector2(0, this.CellSize.y * 0.5f);
        }

        public override Vector3 GetDataFromIndex(params int[] list)
        {
            Assert.IsTrue(list.Length == 2);

            //do not use this to access data
            Assert.IsTrue(false);
            return default;
        }

        public float GetuDataFromIndex(Vector2Int uIndex)
        {
            Vector2Int index = new Vector2Int(Mathf.Clamp(uIndex.x, 0, (int)this.DataSize.x - 1), Mathf.Clamp(uIndex.y, 0, (int)this.DataSize.y - 1));

            var dataIndex = index.x * (int)this.DataSize.y + index.y;
            return this.uData[dataIndex];
        }

        public float GetvDataFromIndex(Vector2Int vIndex)
        {
            Vector2Int index = new Vector2Int(Mathf.Clamp(vIndex.x, 0, (int)this.DataSize.x - 1), Mathf.Clamp(vIndex.y, 0, (int)this.DataSize.y - 1));

            var dataIndex = index.x * (int)this.DataSize.y + index.y;
            return this.vData[dataIndex];
        }

        public override float GetDivergenceFromIndex(params int[] list)
        {
            Assert.IsTrue(list.Length == 2);
            var index = new Vector2Int(list[0], list[1]);
            var left    = this.GetuDataFromIndex(new Vector2Int(index.x     , index.y));
            var right   = this.GetuDataFromIndex(new Vector2Int(index.x + 1 , index.y));
            var up      = this.GetvDataFromIndex(new Vector2Int(index.x     , index.y));
            var down    = this.GetvDataFromIndex(new Vector2Int(index.x     , index.y + 1));

            return 0.5f * (right - left) / this.CellSize.x
                 + 0.5f * (up - down)    / this.CellSize.y;
        }
        public override float SampleDivergence(Vector2 position)
        {
            Vector2 frac = Vector2.zero;
            Vector2Int index = Vector2Int.zero;

            Helper.GetIndexAndFraction(position, this.Origin, this.CellSize, Vector2Int.zero, new Vector2Int(data.Length - 1, data.Length - 1), out index, out frac);

            var f00 = this.GetDivergenceFromIndex(index.x    , index.y);
            var f10 = this.GetDivergenceFromIndex(index.x + 1, index.y);
            var f01 = this.GetDivergenceFromIndex(index.x    , index.y + 1);
            var f11 = this.GetDivergenceFromIndex(index.x + 1, index.y + 1);

            Debug.LogWarning("Verify this");

            return Helper.BiLerp(f00, f10, f01, f11, frac.x, frac.y);
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