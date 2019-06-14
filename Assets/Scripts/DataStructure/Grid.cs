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
        Vector2 DataOrigin { get;}
        Vector2 CellSize { get; set; }//spacing size of one cell
    }
    public struct Grid2DConfig : Grid2D
    {
        public Vector2Int Resolution { get; set; }
        public Vector2 Origin { get; set; }
        public Vector2 CellSize { get; set; }
        public Vector2 DataOrigin { get;}
    }
    public abstract class ScalarGrid2D : ScaleField2D, Grid2D
    {
        public Vector2 Origin { get; set; }
        public Vector2 CellSize { get; set; }
        public abstract Vector2 DataOrigin { get;}
        public ScalarGrid2D(Grid2DConfig config): base(config.Resolution)
        {
            this.Origin = config.Origin;
            this.CellSize = config.CellSize;
        }        

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

        public Vector2 SampleGrandient(Vector2 position)
        {
            Vector2 frac = Vector2.zero;
            Vector2Int index = Vector2Int.zero;

            FluidHelper.GetIndexAndFraction(position, this.Origin, this.CellSize, Vector2Int.zero, new Vector2Int(data.Length - 1, data.Length - 1), out index, out frac);

            var f00 = this.GetGrandientFromIndex(index.x    , index.y);
            var f10 = this.GetGrandientFromIndex(index.x + 1, index.y);
            var f01 = this.GetGrandientFromIndex(index.x    , index.y + 1);
            var f11 = this.GetGrandientFromIndex(index.x + 1, index.y + 1);

            Debug.LogWarning("Verify this");

            return FluidHelper.BiLerp(f00, f10, f01, f11, frac.x, frac.y);
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

        public float SampleLaplacian(Vector2 position)
        {
            Vector2 frac = Vector2.zero;
            Vector2Int index = Vector2Int.zero;

            FluidHelper.GetIndexAndFraction(position, this.Origin, this.CellSize, Vector2Int.zero, new Vector2Int(data.Length - 1, data.Length - 1), out index, out frac);

            var f00 = this.GetLaplacianFromIndex(index.x    , index.y);
            var f10 = this.GetLaplacianFromIndex(index.x + 1, index.y);
            var f01 = this.GetLaplacianFromIndex(index.x    , index.y + 1);
            var f11 = this.GetLaplacianFromIndex(index.x + 1, index.y + 1);

            Debug.LogWarning("Verify this");

            return FluidHelper.BiLerp(f00, f10, f01, f11, frac.x, frac.y);
        }
    }

    public abstract class VectorGrid2D : VectorField2D, Grid2D
    {
        public Vector2 Origin { get; set; }
        public Vector2 CellSize { get; set; }
        public abstract Vector2 DataOrigin { get;}
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
        public float SampleDivergence(Vector2 position)
        {
            Vector2 frac = Vector2.zero;
            Vector2Int index = Vector2Int.zero;

            FluidHelper.GetIndexAndFraction(position, this.Origin, this.CellSize, Vector2Int.zero, new Vector2Int(data.Length - 1, data.Length - 1), out index, out frac);

            var f00 = this.GetDivergenceFromIndex(index.x    , index.y);
            var f10 = this.GetDivergenceFromIndex(index.x + 1, index.y);
            var f01 = this.GetDivergenceFromIndex(index.x    , index.y + 1);
            var f11 = this.GetDivergenceFromIndex(index.x + 1, index.y + 1);

            Debug.LogWarning("Verify this");

            return FluidHelper.BiLerp(f00, f10, f01, f11, frac.x, frac.y);
        }
        public Vector3 SampleCurl(Vector2 position)
        {
            throw new NotImplementedException();
        }
    }

    public abstract class CollocatedVectorGrid2D : VectorGrid2D
    {
        public CollocatedVectorGrid2D(Grid2DConfig config) : base(config)
        {
            this.InitData();
        }
    }

    public class CellCenteredScalarGrid2D : ScalarGrid2D
    {
        public CellCenteredScalarGrid2D(Grid2DConfig config) : base(config)
        {
            this.InitData();
        }

        public override Vector2 DataOrigin
        {
            get { return this.Origin + (0.5f * this.CellSize); }
        }
    }

    public class VertexCenteredScalarGrid2D : ScalarGrid2D
    {
        public VertexCenteredScalarGrid2D(Grid2DConfig config) : base(config)
        {
            this.InitData();
        }

        public override Vector2 DataOrigin
        {
            get { return this.Origin; }
        }

        public override Vector2Int DataSize
        {
            get{ return this.Resolution + new Vector2Int(1, 1); }
        }
    }

    public class CellCenteredVectorGrid2D : CollocatedVectorGrid2D
    {
        public CellCenteredVectorGrid2D(Grid2DConfig config) : base(config)
        {
            this.InitData();
        }

        public override Vector2 DataOrigin
        {
            get { return this.Origin + (0.5f * this.CellSize); }
        }
    }

    public class VertexCenteredVectorGrid2D : CollocatedVectorGrid2D
    {
        public VertexCenteredVectorGrid2D(Grid2DConfig config) : base(config)
        {
            this.InitData();
        }

        public override Vector2 DataOrigin
        {
            get { return this.Origin; }
        }

        public override Vector2Int DataSize
        {
            get { return this.Resolution + new Vector2Int(1, 1); }
        }
    }


    /// <summary>
    /// MAC Grid
    /// Note that this.data variable is not used in this kind of grid
    /// all data is stored at uData, vData, wData
    /// </summary>
    public class FaceCenteredGrid2D: Grid2D
    {
        protected VertexCenteredScalarGrid2D uData;
        protected VertexCenteredScalarGrid2D vData;

        public Vector2 uDataOrigin { get; set; }
        public Vector2 vDataOrigin { get; set; }
        public Vector2Int uDataSize { get; set; }
        public Vector2Int vDataSize { get; set; }

        public FaceCenteredGrid2D(Grid2DConfig config)
        {
            this.Origin = config.Origin;
            this.CellSize = config.CellSize;

            this.uDataSize = config.Resolution + new Vector2Int(1, 0);
            this.uDataOrigin = config.Origin + new Vector2(0, config.CellSize.y * 0.5f);

            this.vDataSize = config.Resolution + new Vector2Int(0, 1);
            this.vDataOrigin = config.Origin + new Vector2(this.CellSize.x * 0.5f, 0);

            var factory = new GridFactory();

            var uConfig = new Grid2DConfig()
            {
                Resolution = uDataSize,
                Origin = uDataOrigin,
                CellSize = config.CellSize,
            };
            this.uData = factory.MakeGrid2D(GridFactory.CenterType.VertexCentered, GridFactory.DataType.Scalar, uConfig) as VertexCenteredScalarGrid2D;

            var vConfig = new Grid2DConfig()
            {
                Resolution = vDataSize,
                Origin = vDataOrigin,
                CellSize = config.CellSize,
            };
            this.vData = factory.MakeGrid2D(GridFactory.CenterType.VertexCentered, GridFactory.DataType.Scalar, vConfig) as VertexCenteredScalarGrid2D;

        }

        public Vector2 Origin { get; set; }

        public Vector2 DataOrigin { get=>throw new NotImplementedException(); }

        public Vector2 CellSize { get; set; }

/*
        public void InitData()
        {
            //this.data is not used
            this.data = null;

            var totalU = Mathf.CeilToInt((this.Resolution.x + 1) * this.Resolution.y);
            var totalV = Mathf.CeilToInt(this.Resolution.x * (this.Resolution.y + 1));
            this.uData = new float[totalU];
            this.vData = new float[totalV];

            //note origin is inversed, u is in cell y, v is in cell x
            this.uDataOrigin = this.Origin + new Vector2(0, this.CellSize.y * 0.5f);
            this.vDataOrigin = this.Origin + new Vector2(this.CellSize.x * 0.5f, 0);
        }*/
        public void Reset(Vector3 value = default)
        {
            this.uData.Reset(value.x);
            this.vData.Reset(value.y);
        }        
                
        public float GetuDataFromIndex(int x, int y)
        {
            return GetuDataFromIndex(new Vector2Int(x, y));
        }
        public void SetuDataToIndex(float value, int x, int y)
        {
            SetuDataToIndex(value, new Vector2Int(x, y));
        }

        public float GetvDataFromIndex(int x, int y)
        {
            return GetvDataFromIndex(new Vector2Int(x, y));
        }
        public void SetvDataToIndex(float value, int x, int y)
        {
            SetvDataToIndex(value, new Vector2Int(x, y));
        }

        public float GetuDataFromIndex(Vector2Int uIndex)
        {
            return this.uData.GetDataFromIndex(uIndex.x, uIndex.y);
        }        

        public float GetvDataFromIndex(Vector2Int vIndex)
        {
            return this.vData.GetDataFromIndex(vIndex.x, vIndex.y);
        }

        public void SetuDataToIndex(float value, Vector2Int uIndex)
        {
            this.uData.SetDataToIndex(value, uIndex.x, uIndex.y);
        }

        public void SetvDataToIndex(float value, Vector2Int vIndex)
        {
            this.vData.SetDataToIndex(value, vIndex.x, vIndex.y);

        }

        public void ForEachuData(Func<Vector2Int, float, float> func)
        {
            for (var i = 0; i < this.uDataSize.x; ++i)
                for (var j = 0; j < this.uDataSize.y; ++j)
                {
                    var index = new Vector2Int(i, j);
                    this.uData.SetDataToIndex(func(index, this.GetuDataFromIndex(index)), i, j);
                }
        }

        public void ForEachvData(Func<Vector2Int, float, float> func)
        {
            for (var i = 0; i < this.vDataSize.x; ++i)
                for (var j = 0; j < this.vDataSize.y; ++j)
                {
                    var index = new Vector2Int(i, j);
                    this.uData.SetDataToIndex(func(index, this.GetvDataFromIndex(index)), i, j);
                }
        }

        public void AccuDataToIndexWithWeight(float value, Vector2Int[] index, Vector2[] weights)
        {
            for (var i = 0; i < index.Length; ++i)
            {
                var accValue = this.GetuDataFromIndex(index[i]);
                accValue += value * weights[i].x * weights[i].y;
                this.SetuDataToIndex(accValue, index[i]);
            }
        }
        public void AccvDataToIndexWithWeight(float value, Vector2Int[] index, Vector2[] weights)
        {
            for (var i = 0; i < index.Length; ++i)
            {
                var accValue = this.GetvDataFromIndex(index[i]);
                accValue += value * weights[i].x * weights[i].y;
                this.SetvDataToIndex(accValue, index[i]);
            }
        }

        public void CopyTo(FieldData<Vector3, Vector2> target)
        {/*
            var t = target as FaceCenteredGrid2D;
            Assert.IsNotNull(t);

            t.uData = (float[])this.uData.Clone();
            t.vData = (float[])this.vData.Clone();*/
        }

        public float GetDivergenceFromIndex(params int[] list)
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
        public float SampleDivergence(Vector2 position)
        {
            return 0;
            /*Vector2 frac = Vector2.zero;
            Vector2Int index = Vector2Int.zero;

            FluidHelper.GetIndexAndFraction(position, this.Origin, this.CellSize, Vector2Int.zero, new Vector2Int(data.Length - 1, data.Length - 1), out index, out frac);

            var f00 = this.GetDivergenceFromIndex(index.x    , index.y);
            var f10 = this.GetDivergenceFromIndex(index.x + 1, index.y);
            var f01 = this.GetDivergenceFromIndex(index.x    , index.y + 1);
            var f11 = this.GetDivergenceFromIndex(index.x + 1, index.y + 1);

            Debug.LogWarning("Verify this");

            return FluidHelper.BiLerp(f00, f10, f01, f11, frac.x, frac.y);*/
        }

        public Vector3 Sample(Vector2 position)
        {
            Vector2 frac = Vector2.zero;
            Vector2Int index = Vector2Int.zero;

            var low = Vector2Int.zero;
            var high = this.uData.DataSize - new Vector2Int(0, 1);

            FluidHelper.GetIndexAndFraction(position, this.uDataOrigin, this.CellSize, low, high, out index, out frac);
            FluidHelper.ClampIndexAndWeight(low, high, ref index, ref frac);
            var f00 = this.GetuDataFromIndex(index.x, index.y);
            var f10 = this.GetuDataFromIndex(index.x + 1, index.y);
            var f01 = this.GetuDataFromIndex(index.x, index.y + 1);
            var f11 = this.GetuDataFromIndex(index.x + 1, index.y + 1);
            var uValue = FluidHelper.BiLerp(f00, f10, f01, f11, frac.x, frac.y);

            frac = Vector2.zero;
            index = Vector2Int.zero;
            high = this.vData.DataSize - new Vector2Int(1, 0);

            FluidHelper.GetIndexAndFraction(position, this.vDataOrigin, this.CellSize, low, high, out index, out frac);
            FluidHelper.ClampIndexAndWeight(low, high, ref index, ref frac);
            f00 = this.GetvDataFromIndex(index.x, index.y);
            f10 = this.GetvDataFromIndex(index.x + 1, index.y);
            f01 = this.GetvDataFromIndex(index.x, index.y + 1);
            f11 = this.GetvDataFromIndex(index.x + 1, index.y + 1);
            var vValue = FluidHelper.BiLerp(f00, f10, f01, f11, frac.x, frac.y);

            //Debug.LogWarning("Verify this");

            return new Vector3(uValue, vValue, 0);
        }

        public void OnDebugDraw()
        {
            return;
            var size = this.uDataSize;
            var u = new float[size.x];
            var v = new float[size.y];
            this.ForEachuData((index, value) =>
            {
                u[index.x] = value;
                return value;
            });

            this.ForEachvData((index, value) =>
            {
                v[index.y] = value;
                return value;
            });

            for (var i = 0; i < size.x; ++i)
            {
                for(var j = 0; j < size.y; ++j)
                {
                    var org = new Vector3(i, j, 0);
                    var to = org + new Vector3(u[i], v[j], 0);
                    Gizmos.DrawLine(org, to);
                }
            }
        }

        
    }
    public class ParticleFaceCenteredGrid2D
    {
        public FaceCenteredGrid2D Velocity { get => this.uvVelocity; }

        protected FaceCenteredGrid2D uvVelocity;
        protected FaceCenteredGrid2D weights;
        public ParticleFaceCenteredGrid2D(Grid2DConfig config)
        {
            var factory = new GridFactory();
            this.uvVelocity = factory.MakeGrid2D(GridFactory.CenterType.FaceCentered, GridFactory.DataType.Vector, config) as FaceCenteredGrid2D;
            this.weights = factory.MakeGrid2D(GridFactory.CenterType.FaceCentered, GridFactory.DataType.Vector, config) as FaceCenteredGrid2D;
        }

        public void Reset()
        {
            this.uvVelocity.Reset();
            this.weights.Reset();
        }

        public void TransferValueToGrid(Vector2 pos, Vector2 value)
        {
            Vector2Int[] index;
            Vector2[] weights;

            var cellSpace = this.uvVelocity.CellSize;

            var org = this.uvVelocity.uDataOrigin;
            var dataSize = this.uvVelocity.uDataSize - new Vector2Int(0, 1);//take minus 1 to make sure value weight is 1 for upper position

            //note index max should be size-1, clamp did this check
            FluidHelper.GetIndexAndWeight(pos, org, cellSpace,
                                                Vector2Int.zero, dataSize,
                                                out index, out weights);

            this.uvVelocity.AccuDataToIndexWithWeight(value.x, index, weights);
            this.weights.AccuDataToIndexWithWeight(1, index, weights);


            org = this.uvVelocity.vDataOrigin;
            dataSize = this.uvVelocity.vDataSize - new Vector2Int(1, 0);//take minus 1 to make sure value weight is 1 for right position

            FluidHelper.GetIndexAndWeight(pos, org, cellSpace,
                                                Vector2Int.zero, dataSize,
                                                out index, out weights);

            this.uvVelocity.AccvDataToIndexWithWeight(value.x, index, weights);
            this.weights.AccvDataToIndexWithWeight(1, index, weights);
        }

        public void NormalizeWeight()
        {
            this.uvVelocity.ForEachuData((index, value) =>
            {
                var w = this.weights.GetuDataFromIndex(index.x, index.y);
                return value / (w > 0 ? w : 1);
            });
            this.uvVelocity.ForEachvData((index, value) =>
            {
                var w = this.weights.GetvDataFromIndex(index.x, index.y);
                return value / (w > 0 ? w : 1);
            });
        }

        public Vector2 GetInfNorm()
        {
            return FluidHelper.Infnorm(this.uvVelocity);
        }

    }
    public class NullGrid : Grid2D
    {
        static NullGrid grid = new NullGrid();
        public static NullGrid NullInstance() { return grid; }
        public Vector2 Origin { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }
        public Vector2 CellSize { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }
        public Vector2 DataOrigin => throw new NotImplementedException();
    }
}