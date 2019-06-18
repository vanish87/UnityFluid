using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace UnityFluid
{
    public abstract class FieldSampler<DataType, Dimension, SamplerType>
    {
        public abstract DataType Sample(FieldData<DataType, Dimension, SamplerType> field, SamplerType input);
    }

    public class ScalarField2DSampler: FieldSampler<float, Vector2Int, Vector2>
    {
        public override float Sample(FieldData<float, Vector2Int, Vector2> field, Vector2 input)
        {
            Vector2 frac = Vector2.zero;
            Vector2Int index = Vector2Int.zero;

            var size = field.DataSize;

            //get index from position
            FluidHelper.GetIndexAndFraction(input, Vector2.zero, Vector2.one, Vector2Int.zero, size - Vector2Int.one, out index, out frac);

            //then get 4 corner point in this cell
            var f00 = field.GetDataFromIndex(index.x, index.y);
            var f10 = field.GetDataFromIndex(index.x + 1, index.y);
            var f01 = field.GetDataFromIndex(index.x, index.y + 1);
            var f11 = field.GetDataFromIndex(index.x + 1, index.y + 1);

            //do lerp for position
            return FluidHelper.BiLerp(f00, f10, f01, f11, frac.x, frac.y);
        }
    }

    public class VectorField2DSampler : FieldSampler<Vector3, Vector2Int, Vector2>
    {
        public override Vector3 Sample(FieldData<Vector3, Vector2Int, Vector2> field, Vector2 input)
        {
            Vector2 frac = Vector2.zero;
            Vector2Int index = Vector2Int.zero;

            var size = field.DataSize;

            //get index from position
            FluidHelper.GetIndexAndFraction(input, Vector2.zero, Vector2.one, Vector2Int.zero, size - Vector2Int.one, out index, out frac);

            //then get 4 corner point in this cell
            var f00 = field.GetDataFromIndex(index.x, index.y);
            var f10 = field.GetDataFromIndex(index.x + 1, index.y);
            var f01 = field.GetDataFromIndex(index.x, index.y + 1);
            var f11 = field.GetDataFromIndex(index.x + 1, index.y + 1);

            //do lerp for position
            return FluidHelper.BiLerp(f00, f10, f01, f11, frac.x, frac.y);
        }
    }

    public class ScalarGrid2DSampler : FieldSampler<float, Vector2Int, Vector2>
    {
        public override float Sample(FieldData<float, Vector2Int, Vector2> field, Vector2 input)
        {
            Vector2 frac = Vector2.zero;
            Vector2Int index = Vector2Int.zero;

            var dataOrg = (field as Grid2D).DataOrigin;
            var dataSize = field.DataSize;

            var cellSize = (field as Grid2D).CellSize;

            //get index from position
            FluidHelper.GetIndexAndFraction(input, dataOrg, cellSize, Vector2Int.zero, dataSize - Vector2Int.one, out index, out frac);

            //then get 4 corner point in this cell
            var f00 = field.GetDataFromIndex(index.x, index.y);
            var f10 = field.GetDataFromIndex(index.x + 1, index.y);
            var f01 = field.GetDataFromIndex(index.x, index.y + 1);
            var f11 = field.GetDataFromIndex(index.x + 1, index.y + 1);

            //do lerp for position
            return FluidHelper.BiLerp(f00, f10, f01, f11, frac.x, frac.y);
        }
    }

    public class VectorGrid2DSampler : FieldSampler<Vector3, Vector2Int, Vector2>
    {
        public override Vector3 Sample(FieldData<Vector3, Vector2Int, Vector2> field, Vector2 input)
        {
            Vector2 frac = Vector2.zero;
            Vector2Int index = Vector2Int.zero;

            var dataOrg = (field as Grid2D).DataOrigin;
            var dataSize = field.DataSize;

            var cellSize = (field as Grid2D).CellSize;

            //get index from position
            FluidHelper.GetIndexAndFraction(input, dataOrg, cellSize, Vector2Int.zero, dataSize - Vector2Int.one, out index, out frac);

            //then get 4 corner point in this cell
            var f00 = field.GetDataFromIndex(index.x, index.y);
            var f10 = field.GetDataFromIndex(index.x + 1, index.y);
            var f01 = field.GetDataFromIndex(index.x, index.y + 1);
            var f11 = field.GetDataFromIndex(index.x + 1, index.y + 1);

            //do lerp for position
            return FluidHelper.BiLerp(f00, f10, f01, f11, frac.x, frac.y);
        }
    }
}