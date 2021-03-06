﻿using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Assertions;

namespace UnityFluid
{
    public class DataContainer<DataType>
    {
        protected List<DataType> dataList = new List<DataType>();

        public DataType GetData(int index)
        {
            Assert.IsTrue(0 <= index && index < this.dataList.Count);
            return this.dataList[index];
        }
        public void AddData(DataType newData)
        {
            if (this.dataList.Contains(newData)) return;

            this.dataList.Add(newData);
        }
        public void Remove(int index)
        {
            Assert.IsTrue(0 <= index && index < this.dataList.Count);
            this.dataList.RemoveAt(index);
        }
    }
    public class GridSystemData
    {
        public enum DataType
        {
            Pressure = 0,
        }
        protected FaceCenterdVectorGrid2D velocity;
        protected DataContainer<ScalarGrid2D<float>> scalarData = new DataContainer<ScalarGrid2D<float>>();
        protected DataContainer<VectorGrid2D> vectorData = new DataContainer<VectorGrid2D>();

        public VectorGrid2D GetData(DataType type)
        {
            return default;
        }

        void Test()
        {
            var d = this.vectorData.GetData((int)DataType.Pressure);

            var config = new Grid2DConfigure();
            config.Resolution = new Vector2Int(10, 10);
            config.CellSize = new Vector2(1, 1);
            config.Origin = Vector2.zero;

            /*var fac = new GridFactory();
            var grid = fac.MakeGrid2D(GridFactory.CenterType.CellCentered, GridFactory.DataType.Scalar, config);

            velocity = fac.MakeGrid2D(GridFactory.CenterType.FaceCentered, GridFactory.DataType.Vector, config) as FaceCenteredGrid2D;

            this.scalarData.AddData(grid);*/
        }
    }
}