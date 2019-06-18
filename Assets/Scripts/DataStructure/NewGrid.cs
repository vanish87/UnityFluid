using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Assertions;

namespace FluidData
{
    public interface ScalarGridOperation<DataType, GrandientDataType>
    {
        GrandientDataType GetGrandientFromIndex(params int[] list);
        DataType GetLaplacianFromIndex(params int[] list);
    }
    public interface VectorGridOperation<DataType, DivergenceType>
    {
        DivergenceType GetDivergenceFromIndex(params int[] list);
    }
    public interface GridDataOperation<DataType>
    {
        // a = a + value;
        DataType Add(DataType value, params int[] list);
        DataType GetDataFromIndex(params int[] list);
        void SetDataToIndex(DataType value, params int[] list);
    }
    public abstract class Grid<Dimension, SizeType>
    {
        protected Dimension Resolution { get; }
        protected SizeType Origin { get; set; }
        protected SizeType CellSize { get; set; }//spacing size of one cell

        public Grid(Dimension resolution, SizeType cellSize, SizeType origin)
        {
            this.Resolution = resolution;
            this.Origin = origin;
            this.CellSize = cellSize;
        }
        ~Grid()
        {

        }

        public class NullGrid : Grid<int, int>
        {
            static protected NullGrid instance = new NullGrid(default);

            public NullGrid(int resolution) : base(resolution, default, default)
            {
            }

            static public NullGrid Instance { get { return instance; } }
        }
    }

    public abstract class MACGrid2DData : Grid<Vector2Int, Vector2>, GridDataOperation<Vector2>
    {
        protected float[] uData, vData;

        public Vector2 uDataOrigin { get; set; }
        public Vector2 vDataOrigin { get; set; }
        public Vector2Int uDataSize { get; set; }
        public Vector2Int vDataSize { get; set; }

        public MACGrid2DData(Vector2Int resolution, Vector2 cellSize, Vector2 origin) 
            : base(resolution, cellSize, origin)
        {
            this.uDataSize = resolution + new Vector2Int(1, 0);
            this.uDataOrigin = origin + new Vector2(0, cellSize.y * 0.5f);

            this.vDataSize = resolution + new Vector2Int(0, 1);
            this.vDataOrigin = origin + new Vector2(cellSize.x * 0.5f, 0);


            this.uData = new float[this.uDataSize.x * this.uDataSize.y];
            this.vData = new float[this.vDataSize.x * this.vDataSize.y];
        }
        protected int AbsoluteUIndex(int i, int j)
        {
            return this.AbsoluteUIndex(new Vector2Int(i, j));
        }
        protected int AbsoluteVIndex(int i, int j)
        {
            return this.AbsoluteVIndex(new Vector2Int(i, j));
        }
        protected int AbsoluteUIndex(Vector2Int index)
        {
            return index.x + (this.uDataSize.x * index.y);
        }
        protected int AbsoluteVIndex(Vector2Int index)
        {
            return index.x + (this.vDataSize.x * index.y);
        }
        public Vector2 Add(Vector2 value, params int[] list)
        {
            throw new System.NotImplementedException();
        }

        public Vector2 GetDataFromIndex(params int[] list)
        {
            Assert.IsTrue(list.Length == 4);

            var ret = new Vector2();

            Assert.IsTrue(0 <= list[0] && list[0] < this.uDataSize.x);
            Assert.IsTrue(0 <= list[1] && list[1] < this.uDataSize.y);
            Vector2Int index = FluidHelper.ClampIndex(new Vector2Int(list[0], list[1]), Vector2Int.zero, this.uDataSize);
            ret.x = this.uData[index.x + (this.uDataSize.x * index.y)];

            Assert.IsTrue(0 <= list[2] && list[2] < this.vDataSize.x);
            Assert.IsTrue(0 <= list[3] && list[3] < this.vDataSize.y);
            index = FluidHelper.ClampIndex(new Vector2Int(list[2], list[3]), Vector2Int.zero, this.vDataSize);
            ret.y = this.vData[index.x + (this.vDataSize.x * index.y)];

            return ret;
        }

        public void SetDataToIndex(Vector2 value, params int[] list)
        {
            throw new System.NotImplementedException();
        }

        protected void Reset(float value = default)
        {
            for (var i = 0; i < this.uData.Length; ++i) this.uData[i] = value;
            for (var i = 0; i < this.vData.Length; ++i) this.vData[i] = value;
        }


        public delegate void DataFunction(ref float value, int i, int j);
        public void ForEachuData(DataFunction func)
        {
            for (var i = 0; i < this.uDataSize.x; ++i)
                for (var j = 0; j < this.uDataSize.y; ++j)
                {
                    var dataIndex = this.AbsoluteUIndex(i, j);
                    func(ref this.uData[dataIndex], i, j);
                }
        }
        public void ForEachvData(DataFunction func)
        {
            for (var i = 0; i < this.vDataSize.x; ++i)
                for (var j = 0; j < this.vDataSize.y; ++j)
                {
                    var dataIndex = this.AbsoluteUIndex(i, j);
                    func(ref this.vData[dataIndex], i, j);
                }
        }
    }

    public abstract class GridData<DataType, Dimension, SizeType> : Grid<Dimension, SizeType>, GridDataOperation<DataType>
    {
        //all data is one dimensional data
        //so DataSize returns dimension data size
        //but AbsoluteDataIndex will convert index to 1d data index
        protected DataType[] data;
        protected Dimension DataSize { get; set; }
        protected SizeType DataOrigin { get; set; }

        protected abstract int AbsoluteDataIndex(params int[] list);
        protected abstract int AbsoluteDataLength();

        public virtual void InitData()
        {
            this.data = new DataType[this.AbsoluteDataLength()];
        }
        public virtual void Reset(DataType value = default)
        {
            for (var i = 0; i < this.data.Length; ++i) this.data[i] = value;
        }

        public delegate void DataFunction(ref DataType value, params int[] list);
        public abstract void ForEachData(DataFunction func);

        public GridData(Dimension resolution, SizeType cellSize, SizeType origin, SizeType dataOrigin, Dimension dataSize) 
            : base(resolution, cellSize, origin)
        {
            this.DataOrigin = dataOrigin;
            this.DataSize = dataSize;

            this.InitData();
        }

        public virtual void CopyTo(GridData<DataType, Dimension, SizeType> target)
        {
            target.data = (DataType[])this.data.Clone();
        }
        public DataType this[params int[] list]
        {
            get { return this.GetDataFromIndex(list); }
            set { this.SetDataToIndex(value, list); }
        }

        public abstract DataType Add(DataType value, params int[] list);
        public DataType GetDataFromIndex(params int[] list)
        {
            var dataIndex = this.AbsoluteDataIndex(list);
            return this.data[dataIndex];
        }
        public void SetDataToIndex(DataType value, params int[] list)
        {
            var dataIndex = this.AbsoluteDataIndex(list);
            this.data[dataIndex] = value;
        }
    }
    public abstract class ScalarGrid2D<DataType> : GridData<DataType, Vector2Int, Vector2>, ScalarGridOperation<DataType, Vector2>
    {
        public ScalarGrid2D(Vector2Int resolution, Vector2 cellSize, Vector2 origin, Vector2 dataOrigin, Vector2Int dataSize) 
            : base(resolution, cellSize, origin, dataOrigin, dataSize)
        {
        }
        protected override int AbsoluteDataLength()
        {
            return this.DataSize.x * this.DataSize.y;
        }
        protected override int AbsoluteDataIndex(params int[] list)
        {
            Assert.IsTrue(list.Length == 2);
            Assert.IsTrue(0 <= list[0] && list[0] < this.DataSize.x);
            Assert.IsTrue(0 <= list[1] && list[1] < this.DataSize.y);
            Vector2Int index = FluidHelper.ClampIndex(new Vector2Int(list[0], list[1]), Vector2Int.zero, this.DataSize);

            return index.x + (this.DataSize.x * index.y);
        }

        /*public float Sample(Vector2 position)
        {
            Vector2 frac = Vector2.zero;
            Vector2Int index = Vector2Int.zero;

            //get index from position
            FluidHelper.GetIndexAndFraction(position, Vector2.zero, Vector2.one, Vector2Int.zero, new Vector2Int(data.Length - 1, data.Length - 1), out index, out frac);

            //then get 4 corner point in this cell
            var f00 = this.GetDataFromIndex(index.x, index.y);
            var f10 = this.GetDataFromIndex(index.x + 1, index.y);
            var f01 = this.GetDataFromIndex(index.x, index.y + 1);
            var f11 = this.GetDataFromIndex(index.x + 1, index.y + 1);

            //do lerp for position
            return FluidHelper.BiLerp(f00, f10, f01, f11, frac.x, frac.y);
        }*/

        //DOTO implement Grid operations
        public virtual Vector2 GetGrandientFromIndex(params int[] list)
        {
            throw new System.NotImplementedException();
        }
        public virtual DataType GetLaplacianFromIndex(params int[] list)
        {
            throw new System.NotImplementedException();
        }

        public override void ForEachData(DataFunction func)
        {
            for (var i = 0; i < this.DataSize.x; ++i)
                for (var j = 0; j < this.DataSize.y; ++j)
                {
                    var dataIndex = this.AbsoluteDataIndex(i, j);
                    func(ref this.data[dataIndex], i, j);
                }
        }
    }
    public class ScalarGrid2Df : ScalarGrid2D<float>
    {
        public ScalarGrid2Df(Vector2Int resolution, Vector2 cellSize, Vector2 origin = default, Vector2 dataOrigin = default, Vector2Int dataSize = default)
            : base(resolution, cellSize, origin, dataOrigin, dataSize)
        {
        }

        public override float Add(float value, params int[] list)
        {
            return this[list] + value;
        }
        public override Vector2 GetGrandientFromIndex(params int[] list)
        {
            Assert.IsTrue(list.Length == 2);
            var index = new Vector2Int(list[0], list[1]);
            var left = this.GetDataFromIndex(index.x - 1, index.y);
            var right = this.GetDataFromIndex(index.x + 1, index.y);
            var up = this.GetDataFromIndex(index.x, index.y - 1);
            var down = this.GetDataFromIndex(index.x, index.y + 1);

            return 0.5f * new Vector2(right - left, up - down) / this.CellSize;
        }

    }
    public class ScalarGrid2Di : ScalarGrid2D<int>
    {
        public ScalarGrid2Di(Vector2Int resolution, Vector2 cellSize, Vector2 origin = default, Vector2 dataOrigin = default, Vector2Int dataSize = default)
            : base(resolution, cellSize, origin, dataOrigin, dataSize)
        {
        }

        public override int Add(int value, params int[] list)
        {
            return this[list] + value;
        }
        public override Vector2 GetGrandientFromIndex(params int[] list)
        {
            Assert.IsTrue(list.Length == 2);
            var index = new Vector2Int(list[0], list[1]);
            var left = this.GetDataFromIndex(index.x - 1, index.y);
            var right = this.GetDataFromIndex(index.x + 1, index.y);
            var up = this.GetDataFromIndex(index.x, index.y - 1);
            var down = this.GetDataFromIndex(index.x, index.y + 1);

            return 0.5f * new Vector2(right - left, up - down) / this.CellSize;
        }
    }

    public class ScalarField2Df : ScalarGrid2Df
    {
        public ScalarField2Df(Vector2Int resolution) : base(resolution, Vector2.one)
        {
            //Field is a special type of grid
            //which has DataOrigin 0 and Grid size 1
            this.DataOrigin = Vector2.zero;
            this.CellSize = Vector2.one;
        }
    }
    public class VectorGrid2D : GridData<Vector3, Vector2Int, Vector2>, VectorGridOperation<Vector3, float>
    {
        public VectorGrid2D(Vector2Int resolution, Vector2 cellSize, Vector2 origin, Vector2 dataOrigin, Vector2Int dataSize)
            : base(resolution, cellSize, origin, dataOrigin, dataSize)
        {
        }
        protected override int AbsoluteDataLength()
        {
            return this.DataSize.x * this.DataSize.y;
        }
        protected override int AbsoluteDataIndex(params int[] list)
        {
            Assert.IsTrue(list.Length == 2);
            Assert.IsTrue(0 <= list[0] && list[0] < this.DataSize.x);
            Assert.IsTrue(0 <= list[1] && list[1] < this.DataSize.y);
            Vector2Int index = FluidHelper.ClampIndex(new Vector2Int(list[0], list[1]), Vector2Int.zero, this.DataSize);

            return index.x + (this.DataSize.x * index.y);
        }

        /*public Vector3 Sample(Vector2 position)
        {
            Vector2 frac = Vector2.zero;
            Vector2Int index = Vector2Int.zero;

            FluidHelper.GetIndexAndFraction(position, Vector2.zero, Vector2.one, Vector2Int.zero, new Vector2Int(data.Length - 1, data.Length - 1), out index, out frac);

            var f00 = this.GetDataFromIndex(index.x, index.y);
            var f10 = this.GetDataFromIndex(index.x + 1, index.y);
            var f01 = this.GetDataFromIndex(index.x, index.y + 1);
            var f11 = this.GetDataFromIndex(index.x + 1, index.y + 1);

            return FluidHelper.BiLerp(f00, f10, f01, f11, frac.x, frac.y);
        }*/


        //DOTO implement Grid operations
        //public virtual Vector3 SampleCurl(Vector2 position) { return default; }
        //public virtual float SampleDivergence(Vector2 position) { return default; }
        public virtual float GetDivergenceFromIndex(params int[] list)        
        {
            Assert.IsTrue(list.Length == 2);
            var index = new Vector2Int(list[0], list[1]);
            var left = this.GetDataFromIndex(index.x - 1, index.y);
            var right = this.GetDataFromIndex(index.x + 1, index.y);
            var up = this.GetDataFromIndex(index.x, index.y - 1);
            var down = this.GetDataFromIndex(index.x, index.y + 1);

            return 0.5f * (right.x - left.x) / this.CellSize.x
                 + 0.5f * (up.y - down.y) / this.CellSize.y;
        }

        public override void ForEachData(DataFunction func)
        {
            for (var i = 0; i < this.DataSize.x; ++i)
                for (var j = 0; j < this.DataSize.y; ++j)
                {
                    var dataIndex = this.AbsoluteDataIndex(i, j);
                    func(ref this.data[dataIndex], i, j);
                }
        }

        public override Vector3 Add(Vector3 value, params int[] list)
        {
            return this[list] + value;
        }
    }


    public class CellCenteredVectorGrid2D : VectorGrid2D
    {
        public CellCenteredVectorGrid2D(Vector2Int resolution, Vector2 cellSize, Vector2 origin = default) 
            : base(resolution, cellSize, origin, origin + (0.5f * cellSize), resolution)
        {
        }
    }

    public class VertexCenteredVectorGrid2D : VectorGrid2D
    {
        public VertexCenteredVectorGrid2D(Vector2Int resolution, Vector2 cellSize, Vector2 origin = default) 
            : base(resolution, cellSize, origin, origin, resolution + new Vector2Int(1,1))
        {
        }
    }

    public class FaceCenterdVectorGrid2D : MACGrid2DData
    {
        protected ScalarField2Df uWeightSum;
        protected ScalarField2Df vWeightSum;

        public FaceCenterdVectorGrid2D(Vector2Int resolution, Vector2 cellSize, Vector2 origin) 
            : base(resolution, cellSize, origin)
        {
            this.uWeightSum = new ScalarField2Df(this.uDataSize);
            this.vWeightSum = new ScalarField2Df(this.vDataSize);
        }
        public void ResetValueAndWeight()
        {
            this.Reset();
            this.uWeightSum.Reset();
            this.vWeightSum.Reset();
        }
        public void AccumulatePoint(Vector2 pos, Vector2 value)
        {
            Vector2Int[] index;
            Vector2[] weights;

            var cellSpace = this.CellSize;

            var org = this.uDataOrigin;
            var dataSize = this.uDataSize - new Vector2Int(0, 1);//take minus 1 to make sure value weight is 1 for upper position

            //note index max should be size-1, clamp did this check
            FluidHelper.GetIndexAndWeight(pos, org, cellSpace,
                                                Vector2Int.zero, dataSize,
                                                out index, out weights);

            for (var i = 0; i < index.Length; ++i)
            {
                var id = this.AbsoluteUIndex(index[i]);
                var weight = weights[i].x * weights[i].y;
                this.uData[id] += value.x * weight;
                this.uWeightSum[index[i].x, index[i].y] += weight;
            }

            org = this.vDataOrigin;
            dataSize = this.vDataSize - new Vector2Int(1, 0);//take minus 1 to make sure value weight is 1 for upper position

            //note index max should be size-1, clamp did this check
            FluidHelper.GetIndexAndWeight(pos, org, cellSpace,
                                                Vector2Int.zero, dataSize,
                                                out index, out weights);

            for (var i = 0; i < index.Length; ++i)
            {
                var id = this.AbsoluteVIndex(index[i]);
                var weight = weights[i].x * weights[i].y;
                this.vData[id] += value.y * weight;
                this.vWeightSum[index[i].x, index[i].y] += weight;
            }
        }

        public void NormalizeWeight()
        {
            this.ForEachuData((ref float value, int i , int j) =>
            {
                var w = this.uWeightSum[i, j];
                value /= w != 0 ? w : 1;
            });

            this.ForEachvData((ref float value, int i, int j) =>
            {
                var w = this.vWeightSum[i, j];
                value /= w != 0 ? w : 1;
            });
        }
    }
}
