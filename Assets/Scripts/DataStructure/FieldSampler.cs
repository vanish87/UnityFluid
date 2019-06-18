using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace FluidData
{
    public abstract class GridSampler<DataType, Dimension, SizeType, SamplerType>
    {
        public abstract DataType Sample(GridInterface<DataType, Dimension, SizeType> grid, SamplerType input);
    }

    public class ScalarGrid2DSampler : GridSampler<float, Vector2Int, Vector2, Vector2>
    {
        public override float Sample(GridInterface<float, Vector2Int, Vector2> grid, Vector2 input)
        {
            Vector2 frac = Vector2.zero;
            Vector2Int index = Vector2Int.zero;

            var g = grid as GridData<float, Vector2Int, Vector2>;
            var size = g.DataSize;

            //get index from position
            FluidHelper.GetIndexAndFraction(input, g.DataOrigin, g.CellSize, Vector2Int.zero, size - Vector2Int.one, out index, out frac);

            //then get 4 corner point in this cell
            var f00 = grid.GetDataFromIndex(index.x, index.y);
            var f10 = grid.GetDataFromIndex(index.x + 1, index.y);
            var f01 = grid.GetDataFromIndex(index.x, index.y + 1);
            var f11 = grid.GetDataFromIndex(index.x + 1, index.y + 1);

            //do lerp for position
            return FluidHelper.BiLerp(f00, f10, f01, f11, frac.x, frac.y);
        }
    }

    public class VectorGrid2DSampler : GridSampler<Vector3, Vector2Int, Vector2, Vector2>
    {

        public override Vector3 Sample(GridInterface<Vector3, Vector2Int, Vector2> grid, Vector2 input)
        {
            Vector2 frac = Vector2.zero;
            Vector2Int index = Vector2Int.zero;

            var g = grid as GridData<Vector3, Vector2Int, Vector2>;
            var size = g.DataSize;

            //get index from position
            FluidHelper.GetIndexAndFraction(input, g.DataOrigin, g.CellSize, Vector2Int.zero, size - Vector2Int.one, out index, out frac);

            //then get 4 corner point in this cell
            var f00 = grid.GetDataFromIndex(index.x, index.y);
            var f10 = grid.GetDataFromIndex(index.x + 1, index.y);
            var f01 = grid.GetDataFromIndex(index.x, index.y + 1);
            var f11 = grid.GetDataFromIndex(index.x + 1, index.y + 1);

            //do lerp for position
            return FluidHelper.BiLerp(f00, f10, f01, f11, frac.x, frac.y);
        }
    }

    public class MacVectorGrid2DSampler : GridSampler<Vector2, Vector2Int, Vector2, Vector2>
    {
        public override Vector2 Sample(GridInterface<Vector2, Vector2Int, Vector2> grid, Vector2 input)
        {
            Vector2 frac = Vector2.zero;
            Vector2Int index = Vector2Int.zero;

            var g = grid as MACGrid2DData;

            var low = Vector2Int.zero;
            var high = g.uDataSize - new Vector2Int(0, 1);

            FluidHelper.GetIndexAndFraction(input, g.uDataOrigin, g.CellSize, low, high, out index, out frac);
            //FluidHelper.ClampIndexAndWeight(low, high, ref index, ref frac);
            var f00 = g[MACGrid2DData.DataType.U, index.x, index.y];
            var f10 = g[MACGrid2DData.DataType.U, index.x + 1, index.y];
            var f01 = g[MACGrid2DData.DataType.U, index.x, index.y + 1];
            var f11 = g[MACGrid2DData.DataType.U, index.x + 1, index.y + 1];
            var uValue = FluidHelper.BiLerp(f00, f10, f01, f11, frac.x, frac.y);

            frac = Vector2.zero;
            index = Vector2Int.zero;
            high = g.vDataSize - new Vector2Int(1, 0);

            FluidHelper.GetIndexAndFraction(input, g.vDataOrigin, g.CellSize, low, high, out index, out frac);
            //FluidHelper.ClampIndexAndWeight(low, high, ref index, ref frac);
            f00 = g[MACGrid2DData.DataType.V, index.x, index.y];
            f10 = g[MACGrid2DData.DataType.V, index.x + 1, index.y];
            f01 = g[MACGrid2DData.DataType.V, index.x, index.y + 1];
            f11 = g[MACGrid2DData.DataType.V, index.x + 1, index.y + 1];
            var vValue = FluidHelper.BiLerp(f00, f10, f01, f11, frac.x, frac.y);

            //Debug.LogWarning("Verify this");

            return new Vector3(uValue, vValue, 0);
        }
    }


}