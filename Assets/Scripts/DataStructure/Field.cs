using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Assertions;

namespace UnityFluid
{
    public interface ScalarFieldOperation<DataType, GrandientDataType>
    {
        //Field<DataType> Grandient();
        //Field<DataType> Laplacian();

        //DataType SampleGrandient(Dimension position);
        GrandientDataType GetGrandientFromIndex(params int[] list);
        //DataType SampleLaplacian(SamplerType position);
        DataType GetLaplacianFromIndex(params int[] list);
    }
    public interface VectorFieldOperation<DataType, DivergenceType>
    {
        //Field<DataType, SamplerType> Divergence();
        //Field<DataType, SamplerType> Curl();

        //DivergenceType SampleDivergence(SamplerType position);
        DivergenceType GetDivergenceFromIndex(params int[] list);
        //DataType SampleCurl(SamplerType position);
    }
    public interface FieldDataOperation<DataType>
    {
        // a = a + value;
        DataType Add(DataType value, params int[] list);
        DataType GetDataFromIndex(params int[] list);
        void SetDataToIndex(DataType value, params int[] list);
    }
    public abstract class Field<Dimension, SizeType>
    {
        protected Dimension Resolution { get; }
        protected SizeType Origin { get; set; }
        protected SizeType CellSize { get; set; }//spacing size of one cell

        public Field(Dimension resolution)
        {
            this.Resolution = resolution;
        }
        ~Field()
        {

        }

        public class NullField : Field<int, int>
        {
            static protected NullField instance = new NullField(default);

            public NullField(int resolution) : base(resolution)
            {
            }

            static public NullField Instance { get { return instance; } }
        }
    }
    public abstract class FieldData<DataType, Dimension, SizeType> : Field<Dimension, SizeType>, FieldDataOperation<DataType>
    {
        //all data is one dimensional data
        //so DataSize returns dimension data size
        //but AbsoluteDataIndex will convert index to 1d data index
        protected DataType[] data;
        public virtual Dimension DataSize { get { return this.Resolution; } }
        public SizeType DataOrigin { get; }

        protected abstract int AbsoluteDataIndex(params int[] list);
        protected abstract int AbsoluteDataLength();

        public virtual void InitData()
        {
            this.data = new DataType[this.AbsoluteDataLength()];
        }
        public virtual void Reset(DataType value = default)
        {
            for (var i = 0; i < this.data.Length; ++i)this.data[i] = value;
        }

        public delegate void DataFunction(ref DataType value, params int[] list);
        public abstract void ForEachData(DataFunction func);

        public FieldData(Dimension resolution) : base(resolution)
        {
        }

        public virtual void CopyTo(FieldData<DataType, Dimension, SizeType> target)
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

    /*public class ScaleField1D : FieldData<float, float>, FieldOperation<float, float>
    {
        public Field<float, float> Grandient()
        {
            throw new System.NotImplementedException();
        }

        public Field<float, float> Laplacian()
        {
            throw new System.NotImplementedException();
        }

        public override float Sample(float position)
        {
            throw new System.NotImplementedException();
        }

        public float SampleGrandient(float position)
        {
            throw new System.NotImplementedException();
        }
    }
    public class VectorField1D : FieldData<Vector3, float>, FieldOperation<Vector3, float>, VectorFieldOperation<Vector3, float>
    {
        public Field<Vector3, float> Curl()
        {
            throw new System.NotImplementedException();
        }

        public Field<Vector3, float> Divergence()
        {
            throw new System.NotImplementedException();
        }

        public Field<Vector3, float> Grandient()
        {
            throw new System.NotImplementedException();
        }

        public Field<Vector3, float> Laplacian()
        {
            throw new System.NotImplementedException();
        }

        public override Vector3 Sample(float position)
        {
            throw new System.NotImplementedException();
        }

        public Vector3 SampleGrandient(float position)
        {
            throw new System.NotImplementedException();
        }
    }*/

    public abstract class ScalarField2D<DataType> : FieldData<DataType, Vector2Int, Vector2>, ScalarFieldOperation<DataType, Vector2>
    {
        public ScalarField2D(Vector2Int resolution) : base(resolution)
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

        //DOTO implement field operations
        public virtual Vector2 GetGrandientFromIndex(params int[] list) { return default; }
        public virtual DataType GetLaplacianFromIndex(params int[] list) { return default; }
        
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

    public class ScalarField2Df : ScalarField2D<float>
    {
        public ScalarField2Df(Vector2Int resolution) : base(resolution)
        {
        }

        public override float Add(float value, params int[] list)
        {
            return this[list] + value;
        }

    }
    public class ScalarField2Di : ScalarField2D<int>
    {
        public ScalarField2Di(Vector2Int resolution) : base(resolution)
        {
        }

        public override int Add(int value, params int[] list)
        {
            return this[list] + value;
        }
    }

    public class VectorField2D : FieldData<Vector3, Vector2Int, Vector2>, VectorFieldOperation<Vector3, float>
    {
        public VectorField2D(Vector2Int resolution) : base(resolution)
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

        public Vector3 Sample(Vector2 position)
        {
            Vector2 frac = Vector2.zero;
            Vector2Int index = Vector2Int.zero;

            FluidHelper.GetIndexAndFraction(position, Vector2.zero, Vector2.one, Vector2Int.zero, new Vector2Int(data.Length - 1, data.Length - 1), out index, out frac);

            var f00 = this.GetDataFromIndex(index.x     , index.y);
            var f10 = this.GetDataFromIndex(index.x + 1 , index.y);
            var f01 = this.GetDataFromIndex(index.x     , index.y + 1);
            var f11 = this.GetDataFromIndex(index.x + 1 , index.y + 1);

            return FluidHelper.BiLerp(f00, f10, f01, f11, frac.x, frac.y);
        }


        //DOTO implement field operations
        //public virtual Vector3 SampleCurl(Vector2 position) { return default; }
        //public virtual float SampleDivergence(Vector2 position) { return default; }
        public virtual float GetDivergenceFromIndex(params int[] list) { return default; }
        
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
}