using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Assertions;

namespace UnityFluid
{
    public interface FieldOperation<DataType, SamplerType>
    {
        Field<DataType, SamplerType> Grandient();
        Field<DataType, SamplerType> Laplacian();

        DataType Sample(SamplerType position);
        SamplerType SampleGrandient(SamplerType position);
        SamplerType GetGrandientFromIndex(params int[] list);
        DataType SampleLaplacian(SamplerType position);
    }
    public interface VectorFieldOperation<DataType, SamplerType>
    {
        Field<DataType, SamplerType> Divergence();
        Field<DataType, SamplerType> Curl();

        DataType SampleDivergence(SamplerType position);
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

        public class NullField : Field<DataType, SamplerType>
        {
            static protected NullField instance = new NullField(default);

            public NullField(SamplerType resolution) : base(resolution)
            {
            }

            static public NullField Instance { get { return instance; } }
        }
    }

    public abstract class FieldData<T, S> : Field<T, S>
    {
        protected T[] data;
        public abstract T GetDataFromIndex(params int[] list);
        public abstract void InitData();

        public FieldData(S resolution) : base(resolution)
        {
            this.InitData();
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

    public class ScaleField2D : FieldData<float, Vector2>, FieldOperation<float, Vector2>
    {
        public ScaleField2D(Vector2 resolution) : base(resolution)
        {
        }

        public virtual float Sample(Vector2 position)
        {
            Vector2 frac = Vector2.zero;
            Vector2Int index = Vector2Int.zero;

            Helper.GetIndexAndFraction(position, Vector2.zero, Vector2.one, Vector2Int.zero, new Vector2Int(data.Length-1, data.Length - 1), out index, out frac);

            var f00 = this.GetDataFromIndex(index.x    , index.y    );
            var f10 = this.GetDataFromIndex(index.x + 1, index.y    );
            var f01 = this.GetDataFromIndex(index.x    , index.y + 1);
            var f11 = this.GetDataFromIndex(index.x + 1, index.y + 1);

            return Helper.BiLerp(f00,f10, f01, f11, frac.x, frac.y);
        }
        
        Field<float, Vector2> FieldOperation<float, Vector2>.Grandient() { return NullField.Instance; }

        Field<float, Vector2> FieldOperation<float, Vector2>.Laplacian() { return NullField.Instance; }

        public virtual Vector2 SampleGrandient(Vector2 position) { return default; }
        public virtual Vector2 GetGrandientFromIndex(params int[] list) { return default; }
        public virtual float SampleLaplacian(Vector2 position) { return default; }

        public override float GetDataFromIndex(params int[] list)
        {
            Assert.IsTrue(list.Length == 2);

            Vector2Int index = new Vector2Int(Mathf.Clamp(list[0], 0, (int)this.Resolution.x - 1), Mathf.Clamp(list[1], 0, (int)this.Resolution.y- 1));

            var dataIndex = index.x * (int)this.Resolution.y + index.y;
            return this.data[dataIndex];
        }
        public override void InitData()
        {
            var totalCount = this.Resolution.x * this.Resolution.y;
            this.data = new float[(int)totalCount];
        }

    }
    public class VectorField2D : FieldData<Vector3, Vector2>, FieldOperation<Vector3, Vector2>, VectorFieldOperation<Vector3, Vector2>
    {
        public VectorField2D(Vector2 resolution) : base(resolution)
        {
        }

        public Field<Vector3, Vector2> Curl() { return NullField.Instance; }
        public Field<Vector3, Vector2> Divergence() { return NullField.Instance; }

        public Vector3 SampleCurl(Vector2 position) { return default; }

        public Vector3 SampleDivergence(Vector2 position) { return default; }

        Field<Vector3, Vector2> FieldOperation<Vector3, Vector2>.Grandient() { return NullField.Instance; }
        Field<Vector3, Vector2> FieldOperation<Vector3, Vector2>.Laplacian() { return NullField.Instance; }

        Vector3 FieldOperation<Vector3, Vector2>.Sample(Vector2 position)
        {
            Vector2 frac = Vector2.zero;
            Vector2Int index = Vector2Int.zero;

            Helper.GetIndexAndFraction(position, Vector2.zero, Vector2.one, Vector2Int.zero, new Vector2Int(data.Length - 1, data.Length - 1), out index, out frac);

            var f00 = this.GetDataFromIndex(index.x     , index.y);
            var f10 = this.GetDataFromIndex(index.x + 1 , index.y);
            var f01 = this.GetDataFromIndex(index.x     , index.y + 1);
            var f11 = this.GetDataFromIndex(index.x + 1 , index.y + 1);

            return Helper.BiLerp(f00, f10, f01, f11, frac.x, frac.y);
        }

        public virtual Vector2 GetGrandientFromIndex(params int[] list) { return default; }
        Vector2 FieldOperation<Vector3, Vector2>.SampleGrandient(Vector2 position) { return default; }
        Vector3 FieldOperation<Vector3, Vector2>.SampleLaplacian(Vector2 position) { return default; }

        public override Vector3 GetDataFromIndex(params int[] list)
        {
            Assert.IsTrue(list.Length == 2);
            Vector2Int index = new Vector2Int(Mathf.Clamp(list[0], 0, (int)this.Resolution.x - 1), Mathf.Clamp(list[1], 0, (int)this.Resolution.y - 1));

            var dataIndex = index.x * (int)this.Resolution.y + index.y;
            return this.data[dataIndex];
        }
        public override void InitData()
        {
            var totalCount = this.Resolution.x * this.Resolution.y;
            this.data = new Vector3[(int)totalCount];
        }
    }
}