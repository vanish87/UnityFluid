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
        void Acc(DataType value, params int[] list);
        DataType GetDataFromIndex(params int[] list);
        void SetDataToIndex(DataType value, params int[] list);
    }
    public abstract class Field<Dimension>
    {
        protected Dimension Resolution { get;}

        public Field(Dimension resolution)
        {
            this.Resolution = resolution;
        }
        ~Field()
        {

        }

        public class NullField : Field<int>
        {
            static protected NullField instance = new NullField(default);

            public NullField(int resolution) : base(resolution)
            {
            }

            static public NullField Instance { get { return instance; } }
        }
    }
    public abstract class FieldData<DataType, Dimension> : Field<Dimension>, FieldDataOperation<DataType>
    {
        protected DataType[] data;
        public virtual Dimension DataSize { get { return this.Resolution; } }
        protected abstract int AbsoluteIndex(params int[] list);

        public abstract void InitData();
        public virtual void Reset(DataType value = default)
        {
            for (var i = 0; i < this.data.Length; ++i)
            {
                this.data[i] = value;
            }
        }

        public delegate DataType DataFunction(DataType value, params int[] list);
        public abstract void ForEachData(DataFunction func);

        public FieldData(Dimension resolution) : base(resolution)
        {
        }

        public virtual void CopyTo(FieldData<DataType, Dimension> target)
        {
            target.data = (DataType[])this.data.Clone();
        }

        public abstract void Acc(DataType value, params int[] list);
        public abstract DataType GetDataFromIndex(params int[] list);
        public abstract void SetDataToIndex(DataType value, params int[] list);
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

    public class ScaleField2D : FieldData<float, Vector2Int>, ScalarFieldOperation<float, Vector2>
    {
        public ScaleField2D(Vector2Int resolution) : base(resolution)
        {
        }
        protected override int AbsoluteIndex(params int[] list)
        {
            Assert.IsTrue(list.Length == 2);
            Assert.IsTrue(0 <= list[0] && list[0] < this.DataSize.x);
            Assert.IsTrue(0 <= list[1] && list[1] < this.DataSize.y);
            Vector2Int index = FluidHelper.ClampIndex(new Vector2Int(list[0], list[1]), Vector2Int.zero, this.DataSize);

            return index.x + (this.DataSize.x * index.y);
        }

        public float Sample(Vector2 position)
        {
            Vector2 frac = Vector2.zero;
            Vector2Int index = Vector2Int.zero;

            //get index from position
            FluidHelper.GetIndexAndFraction(position, Vector2.zero, Vector2.one, Vector2Int.zero, new Vector2Int(data.Length-1, data.Length - 1), out index, out frac);

            //then get 4 corner point in this cell
            var f00 = this.GetDataFromIndex(index.x    , index.y    );
            var f10 = this.GetDataFromIndex(index.x + 1, index.y    );
            var f01 = this.GetDataFromIndex(index.x    , index.y + 1);
            var f11 = this.GetDataFromIndex(index.x + 1, index.y + 1);

            //do lerp for position
            return FluidHelper.BiLerp(f00, f10, f01, f11, frac.x, frac.y);
        }

        //DOTO implement field operations
        public virtual Vector2 GetGrandientFromIndex(params int[] list) { return default; }
        public virtual float GetLaplacianFromIndex(params int[] list) { return default; }

        public override float GetDataFromIndex(params int[] list)
        {
            var dataIndex = this.AbsoluteIndex(list);
            return this.data[dataIndex];
        }

        public override void SetDataToIndex(float value, params int[] list)
        {
            var dataIndex = this.AbsoluteIndex(list);
            this.data[dataIndex] = value;
        }

        public override void InitData()
        {
            var totalCount = this.DataSize.x * this.DataSize.y;
            this.data = new float[totalCount];
        }

        public override void ForEachData(DataFunction func)
        {
            for (var i = 0; i < this.DataSize.x; ++i)
                for (var j = 0; j < this.DataSize.y; ++j)
                {
                    var dataIndex = this.AbsoluteIndex(i, j);
                    this.data[dataIndex] = func(this.data[dataIndex], i, j);
                }
        }

        public override void Acc(float value, params int[] list)
        {
            this.ForEachData((data, index) => data + value);
        }
    }
    public class VectorField2D : FieldData<Vector3, Vector2Int>, VectorFieldOperation<Vector3, float>
    {
        public VectorField2D(Vector2Int resolution) : base(resolution)
        {
        }
        protected override int AbsoluteIndex(params int[] list)
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

        public override Vector3 GetDataFromIndex(params int[] list)
        {
            var dataIndex = this.AbsoluteIndex(list);
            return this.data[dataIndex];
        }
        public override void SetDataToIndex(Vector3 value, params int[] list)
        {
            var dataIndex = this.AbsoluteIndex(list);
            this.data[dataIndex] = value;
        }

        public override void InitData()
        {
            var totalCount = this.DataSize.x * this.DataSize.y;
            this.data = new Vector3[totalCount];
        }
        public override void ForEachData(DataFunction func)
        {
            for (var i = 0; i < this.DataSize.x; ++i)
                for (var j = 0; j < this.DataSize.y; ++j)
                {
                    var dataIndex = this.AbsoluteIndex(i, j);
                    this.data[dataIndex] = func(this.data[dataIndex], i, j);
                }
        }

        public override void Acc(Vector3 value, params int[] list)
        {
            this.ForEachData((data, index) => data + value);
        }
    }
}