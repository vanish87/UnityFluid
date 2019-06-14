using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Assertions;

public class FluidHelper
{
    static public Vector2Int Clamp(Vector2Int value, Vector2Int min, Vector2Int max)
    {
        return new Vector2Int(Mathf.Clamp(value.x, min.x, max.x), Mathf.Clamp(value.y, min.y, max.y));
    }
    static public Vector2Int ClampIndex(Vector2Int value, Vector2Int min, Vector2Int max)
    {
        return new Vector2Int(Mathf.Clamp(value.x, min.x, max.x-1), Mathf.Clamp(value.y, min.y, max.y-1));
    }
    static public float Lerp(float f0, float f1, float t)
    {
        return ((1 - t) * f0) + (t * f1);
    }

    static public Vector3 Lerp(Vector3 f0, Vector3 f1, float t)
    {
        return ((1 - t) * f0) + (t * f1);
    }

    static public float BiLerp(float f00, float f10, float f01, float f11, float tx, float ty)
    {
        return Lerp(Lerp(f00, f10, tx), Lerp(f01, f11, tx), ty);
    }

    static public Vector3 BiLerp(Vector3 f00, Vector3 f10, Vector3 f01, Vector3 f11, float tx, float ty)
    {
        return Lerp(Lerp(f00, f10, tx), Lerp(f01, f11, tx), ty);
    }

    static public void GetIndexAndFraction(float pos, out int index, out float frac)
    {
        //Note we use another function to Clamp index in range of low and high
        index = Mathf.FloorToInt(pos);
        frac = pos - index;
    }

    static public void GetIndexAndFraction(
        Vector2 position, Vector2 origin, Vector2 spacing, 
        Vector2Int low  , Vector2Int high, 
        out Vector2Int index, out Vector2 frac)
    {
        var x = 0;
        var y = 0;
        var fx = 0f;
        var fy = 0f;

        Assert.IsTrue(spacing.x > 0 && spacing.y > 0);
        var normalizedPos = (position - origin) / spacing;

        GetIndexAndFraction(normalizedPos.x, out x, out fx);
        GetIndexAndFraction(normalizedPos.y, out y, out fy);

        index = new Vector2Int(x, y);
        frac =  new Vector2(fx, fy);
        
        ClampIndexAndWeight(low, high, ref index, ref frac);

        //Assert.IsTrue(low.x <= index.x && index.x <= high.x);
        //Assert.IsTrue(low.y <= index.y && index.y <= high.y);
    }

    static public void ClampIndexAndFrac(int low, int high, ref int index, ref float frac)
    {
        if (index < low)
        {
            index = low;
            frac = 0;
        }
        else
        if (index > high - 1)
        {
            index = high - 1;
            frac = 1;
        }
    }
    static public void ClampIndexAndWeight(Vector2Int low, Vector2Int high, ref Vector2Int index, ref Vector2 frac)
    {
        var x = index.x;
        var y = index.y;
        ClampIndexAndFrac(low.x, high.x, ref x, ref frac.x);
        ClampIndexAndFrac(low.y, high.y, ref y, ref frac.y);

        index.x = x;
        index.y = y;
    }

    static public void GetIndexAndWeight(Vector2 position, Vector2 origin, Vector2 spacing,
        Vector2Int low, Vector2Int high, out Vector2Int[] index, out Vector2[] weights)
    {
        var x = 0;
        var y = 0;
        var fx = 0f;
        var fy = 0f;

        Assert.IsTrue(spacing.x > 0 && spacing.y > 0);
        var normalizedPos = (position - origin) / spacing;

        GetIndexAndFraction(normalizedPos.x, out x, out fx);
        GetIndexAndFraction(normalizedPos.y, out y, out fy);

        ClampIndexAndFrac(low.x, high.x, ref x, ref fx);
        ClampIndexAndFrac(low.y, high.y, ref y, ref fy);

        index = new Vector2Int[]
            {
                new Vector2Int(x    ,y),
                new Vector2Int(x +1 ,y),
                new Vector2Int(x    ,y+1),
                new Vector2Int(x +1 ,y+1),
            };

        weights = new Vector2[]
            {
                new Vector2(1-fx,1-fy),
                new Vector2(fx  ,1-fy),
                new Vector2(1-fx,fy),
                new Vector2(fx  ,fy),
            };
    }

   static public float Infnorm(UnityFluid.CellCenteredScalarGrid2D field) 
   { 
      float r = 0;
        field.ForEachData((value, index) =>
        {
            if (!(Mathf.Abs(value) <= r))
                r = Mathf.Abs(value);
            return value;
        });
      return r;
   }
    static public Vector2 Infnorm(UnityFluid.FaceCenteredGrid2D field)
    {
        float ru = 0, rv = 0;
        field.ForEachuData((index, value) =>
        {
            if (!(Mathf.Abs(value) <= ru))
                ru = Mathf.Abs(value);
            return value;
        });

        field.ForEachvData((index, value) =>
        {
            if (!(Mathf.Abs(value) <= rv))
                rv = Mathf.Abs(value);
            return value;
        });
        Debug.Log(ru + " " + rv);
        return new Vector2(ru, rv);
    }

    static public float Dot(UnityFluid.CellCenteredScalarGrid2D lhs, UnityFluid.CellCenteredScalarGrid2D rhs)
    {
        float r = 0;
        lhs.ForEachData((value, index) => {r += value * rhs.GetDataFromIndex(index[0], index[1]); return value; });
        return r;
    }

    static public void Increment(UnityFluid.CellCenteredScalarGrid2D lhs, UnityFluid.CellCenteredScalarGrid2D rhs, float scale)
    {
        lhs.ForEachData((value, index) => { value += scale * rhs.GetDataFromIndex(index[0], index[1]); return value; });
    }
    static public void SacleAndIncrement(UnityFluid.CellCenteredScalarGrid2D lhs, UnityFluid.CellCenteredScalarGrid2D rhs, float scale)
    {
        lhs.ForEachData((value, index) => { value = value * scale + rhs.GetDataFromIndex(index[0], index[1]); return value; });
    }
}
