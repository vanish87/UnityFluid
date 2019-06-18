using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Assertions;

namespace FluidData
{
    public class GridFactory
    {
        public enum CenterType
        {
            CellCentered,
            VertexCentered,
            FaceCentered,
        }

        public enum DataType
        {
            Scalar,
            Vector
        }
        public Grid<Vector2Int, Vector2> MakeGrid2D(CenterType type, DataType dataType, Grid2DConfigure config)
        {
            Grid<Vector2Int, Vector2> ret = Grid<Vector2Int, Vector2>.NullGrid.Instance;
            switch (dataType)
            {
                case DataType.Scalar: ret = this.MakeScalar(type, config);break;
                case DataType.Vector: ret = this.MakeVector(type, config);break;
                default: Assert.IsFalse(true); break;
            }

            return ret;
        }
        public Grid<Vector2Int, Vector2> MakeGrid2Di(CenterType type, DataType dataType, Grid2DConfigure config)
        {
            Grid<Vector2Int, Vector2> ret = Grid<Vector2Int, Vector2>.NullGrid.Instance;
            switch (dataType)
            {
                case DataType.Scalar:
                    {
                        if (type == CenterType.CellCentered)
                            ret = new CellCenteredScalarGrid2Di(config.Resolution, config.CellSize, config.Origin);
                        else
                            Assert.IsFalse(true);
                    }
                    break;
                case DataType.Vector: //no vector int grid for now 
                default: Assert.IsFalse(true); break;
            }

            return ret;
        }

        protected Grid<Vector2Int, Vector2> MakeScalar(CenterType center, Grid2DConfigure config)
        {
            switch (center)
            {
                case CenterType.CellCentered:
                    return new CellCenteredScalarGrid2D(config.Resolution, config.CellSize, config.Origin);
                case CenterType.VertexCentered:
                    return new VertexCenteredScalarGrid2D(config.Resolution, config.CellSize, config.Origin);
                case CenterType.FaceCentered:
                default:
                    Assert.IsFalse(true);
                    return Grid<Vector2Int, Vector2>.NullGrid.Instance;
            }
        }
        protected Grid<Vector2Int, Vector2> MakeVector(CenterType center, Grid2DConfigure config)
        {
            switch (center)
            {
                case CenterType.CellCentered:
                    return new CellCenteredVectorGrid2D(config.Resolution, config.CellSize, config.Origin);
                case CenterType.VertexCentered:
                    return new VertexCenteredVectorGrid2D(config.Resolution, config.CellSize, config.Origin);
                case CenterType.FaceCentered:
                    return new FaceCenterdVectorGrid2D(config.Resolution, config.CellSize, config.Origin);
                default:
                    Assert.IsFalse(true);
                    return Grid<Vector2Int, Vector2>.NullGrid.Instance;
            }
        }
    }



}