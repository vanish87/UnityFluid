using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Assertions;

namespace UnityFluid
{
    public interface ScalarFieldOperation<DataType, SamplerType>
    {
        Field<DataType, SamplerType> Grandient();
        Field<DataType, SamplerType> Laplacian();

        SamplerType SampleGrandient(SamplerType position);
        SamplerType GetGrandientFromIndex(params int[] list);
        DataType SampleLaplacian(SamplerType position);
        DataType GetLaplacianFromIndex(params int[] list);
    }
    public interface VectorFieldOperation<DataType, SamplerType, DivergenceType>
    {
        Field<DataType, SamplerType> Divergence();
        Field<DataType, SamplerType> Curl();

        DivergenceType SampleDivergence(SamplerType position);
        DivergenceType GetDivergenceFromIndex(params int[] list);
        DataType SampleCurl(SamplerType position);
    }
    public abstract class Field<DataType, SamplerType>
    {
        protected SamplerType Resolution { get;}

        public Field(SamplerType resolution)
        {
            this.Resolution = resolution;
        }
        ~Field()
        {

        }

        public abstract DataType Sample(SamplerType position);

        public class NullField : Field<DataType, SamplerType>
        {
            static protected NullField instance = new NullField(default);

            public NullField(SamplerType resolution) : base(resolution)
            {
            }

            static public NullField Instance { get { return instance; } }

            public override DataType Sample(SamplerType position)
            {
                throw new System.NotImplementedException();
            }
        }
    }

    public abstract class FieldData<T, S> : Field<T, S>
    {
        protected T[] data;
        protected virtual S DataSize { get { return this.Resolution; } }
        public abstract T GetDataFromIndex(params int[] list);
        public abstract void SetDataToIndex(T value, params int[] list);
        public abstract void InitData();
        public virtual void Reset(T value = default)
        {
            for (var i = 0; i < this.data.Length; ++i)
            {
                this.data[i] = value;
            }
        }

        public delegate T DataFunction(T value, params int[] list);
        public abstract void ForEachData(DataFunction func);

        public FieldData(S resolution) : base(resolution)
        {
            this.InitData();
        }

        public virtual void CopyTo(FieldData<T, S> target)
        {
            target.data = (T[])this.data.Clone();
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

    public class ScaleField2D : FieldData<float, Vector2>, ScalarFieldOperation<float, Vector2>
    {
        public ScaleField2D(Vector2 resolution) : base(resolution)
        {
        }

        public override float Sample(Vector2 position)
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
        public virtual Field<float, Vector2> Grandient() { return NullField.Instance; }
        public virtual Field<float, Vector2> Laplacian() { return NullField.Instance; }
        public virtual Vector2 SampleGrandient(Vector2 position) { return default; }
        public virtual Vector2 GetGrandientFromIndex(params int[] list) { return default; }
        public virtual float SampleLaplacian(Vector2 position) { return default; }
        public virtual float GetLaplacianFromIndex(params int[] list) { return default; }

        public override float GetDataFromIndex(params int[] list)
        {
            Assert.IsTrue(list.Length == 2);

            Vector2Int index = new Vector2Int(Mathf.Clamp(list[0], 0, (int)this.Resolution.x - 1), Mathf.Clamp(list[1], 0, (int)this.Resolution.y - 1));

            var dataIndex = index.x * (int)this.Resolution.y + index.y;
            return this.data[dataIndex];
        }

        public override void SetDataToIndex(float value, params int[] list)
        {
            Assert.IsTrue(list.Length == 2);

            Vector2Int index = new Vector2Int(Mathf.Clamp(list[0], 0, (int)this.Resolution.x - 1), Mathf.Clamp(list[1], 0, (int)this.Resolution.y- 1));

            var dataIndex = index.x * (int)this.Resolution.y + index.y;
            this.data[dataIndex] = value;
        }

        public override void InitData()
        {
            var totalCount = this.Resolution.x * this.Resolution.y;
            this.data = new float[(int)totalCount];
        }

        public override void ForEachData(DataFunction func)
        {
            for (var i = 0; i < this.DataSize.x; ++i)
                for (var j = 0; j < this.DataSize.y; ++j)
                {
                    Vector2Int index = new Vector2Int(i,j);
                    var dataIndex = index.x * (int)this.DataSize.y + index.y;

                    this.data[dataIndex] = func(this.data[dataIndex], i, j);
                }
        }
    }
    public class VectorField2D : FieldData<Vector3, Vector2>, VectorFieldOperation<Vector3, Vector2, float>
    {
        public VectorField2D(Vector2 resolution) : base(resolution)
        {
        }

        public override Vector3 Sample(Vector2 position)
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
        public virtual Field<Vector3, Vector2> Curl() { return NullField.Instance; }
        public virtual Field<Vector3, Vector2> Divergence() { return NullField.Instance; }
        public virtual Vector3 SampleCurl(Vector2 position) { return default; }
        public virtual float SampleDivergence(Vector2 position) { return default; }
        public virtual float GetDivergenceFromIndex(params int[] list) { return default; }

        public override Vector3 GetDataFromIndex(params int[] list)
        {
            Assert.IsTrue(list.Length == 2);
            Vector2Int index = new Vector2Int(Mathf.Clamp(list[0], 0, (int)this.Resolution.x - 1), Mathf.Clamp(list[1], 0, (int)this.Resolution.y - 1));

            var dataIndex = index.x * (int)this.Resolution.y + index.y;
            return this.data[dataIndex];
        }
        public override void SetDataToIndex(Vector3 value, params int[] list)
        {
            Assert.IsTrue(list.Length == 2);

            Vector2Int index = new Vector2Int(Mathf.Clamp(list[0], 0, (int)this.Resolution.x - 1), Mathf.Clamp(list[1], 0, (int)this.Resolution.y - 1));

            var dataIndex = index.x * (int)this.Resolution.y + index.y;
            this.data[dataIndex] = value;
        }

        public override void InitData()
        {
            var totalCount = this.Resolution.x * this.Resolution.y;
            this.data = new Vector3[(int)totalCount];
        }
        public override void ForEachData(DataFunction func)
        {
            for (var i = 0; i < this.DataSize.x; ++i)
                for (var j = 0; j < this.DataSize.y; ++j)
                {
                    Vector2Int index = new Vector2Int(i, j);
                    var dataIndex = index.x * (int)this.DataSize.y + index.y;

                    this.data[dataIndex] = func(this.data[dataIndex], i, j);
                }
        }
    }
}