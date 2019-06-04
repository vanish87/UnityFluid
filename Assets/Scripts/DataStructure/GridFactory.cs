using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Assertions;

namespace UnityFluid
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
        public Grid2D MakeGrid2D(CenterType type, DataType dataType, Grid2DConfig config)
        {
            Grid2D ret = NullGrid.NullInstance();
            switch (dataType)
            {
                case DataType.Scalar: ret = this.MakeScalar(type, config);break;
                case DataType.Vector: ret = this.MakeVector(type, config);break;
                default: Assert.IsFalse(true); break;
            }
            return ret;
        }

        protected Grid2D MakeScalar(CenterType center, Grid2DConfig config)
        {
            switch (center)
            {
                case CenterType.CellCentered:
                    return new CellCenteredScalarGrid2D(config);
                case CenterType.VertexCentered:
                    return new VertexCenteredScalarGrid2D(config);
                case CenterType.FaceCentered:
                default:
                    Assert.IsFalse(true);
                    return NullGrid.NullInstance();
            }
        }
        protected Grid2D MakeVector(CenterType center, Grid2DConfig config)
        {
            switch (center)
            {
                case CenterType.CellCentered:
                    return new CellCenteredCollocatedVectorGrid2D(config);
                case CenterType.VertexCentered:
                    return new VertexCenteredCollocatedVectorGrid2D(config);
                case CenterType.FaceCentered:
                    return new FaceCenteredGrid2D(config);
                default:
                    Assert.IsFalse(true);
                    return NullGrid.NullInstance();
            }
        }
    }



}