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

    static public void GetIndexAndFraction(float pos, int min, int max, out int index, out float frac)
    {
        if (pos > max)
        {
            //here will return start index for bilerp
            index = max - 1;
            frac = 1;
        }
        else
        if(pos < min)
        {
            index = min;
            frac = 0;
        }
        else
        {
            index = Mathf.FloorToInt(pos);
            frac = pos - index;
        }
    }

    static public void GetIndexAndFraction(
        Vector2 position, Vector2 origin, Vector2 spacing,  //coordinate in space
        Vector2Int low  , Vector2Int high,                  //coordinate in DATA SPACE, so make sure high is n-1 as input
        out Vector2Int index, out Vector2 frac)
    {
        var x = 0;
        var y = 0;
        var fx = 0f;
        var fy = 0f;

        Assert.IsTrue(spacing.x > 0 && spacing.y > 0);
        var normalizedPos = (position - origin) / spacing;

        GetIndexAndFraction(normalizedPos.x, low.x, high.x, out x, out fx);
        GetIndexAndFraction(normalizedPos.y, low.y, high.y, out y, out fy);

        index = new Vector2Int(x, y);
        frac =  new Vector2(fx, fy);
        
        Assert.IsTrue(low.x <= index.x && index.x < high.x);
        Assert.IsTrue(low.y <= index.y && index.y < high.y);
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

        //get and clamp pos in space with [min, max]=[low, high]
        //index low/high is also space min/max
        //here use low/high as space min/max
        GetIndexAndFraction(normalizedPos.x, low.x, high.x, out x, out fx);
        GetIndexAndFraction(normalizedPos.y, low.y, high.y, out y, out fy);

        //clamp index for  [0, high-1]
        //here use high as index high
        //var x1 = Mathf.Min(high.x - 1, x + 1);
        //var y1 = Mathf.Min(high.y - 1, y + 1);

        var x1 = x + 1;
        var y1 = y + 1;

        Assert.IsTrue(x1 <= high.x);
        Assert.IsTrue(y1 <= high.y);

        index = new Vector2Int[]
            {
                new Vector2Int(x  ,y),
                new Vector2Int(x1 ,y),
                new Vector2Int(x  ,y1),
                new Vector2Int(x1 ,y1),
            };

        weights = new Vector2[]
            {
                new Vector2(1-fx,1-fy),
                new Vector2(fx  ,1-fy),
                new Vector2(1-fx,fy),
                new Vector2(fx  ,fy),
            };
    }

   static public float Infnorm(FluidData.CellCenteredScalarGrid2D field)
    {
        float r = 0;
        field.ForEachData((ref float value, int[] list) =>
        {
            if (!(Mathf.Abs(value) <= r))
                r = Mathf.Abs(value);
        });
        return r;
    }
    static public Vector2 Infnorm(FluidData.FaceCenterdVectorGrid2D grid)
    {
        float ru = 0, rv = 0;
        grid.ForEachuData((ref float value, int i, int j) =>
        {
            if (!(Mathf.Abs(value) <= ru))
                ru = Mathf.Abs(value);
        });

        grid.ForEachvData((ref float value, int i, int j) =>
        {
            if (!(Mathf.Abs(value) <= rv))
                rv = Mathf.Abs(value);
        });
        //Debug.Log(ru + " " + rv);
        return new Vector2(ru, rv);
    }

    static public float Dot(FluidData.CellCenteredScalarGrid2D lhs, FluidData.CellCenteredScalarGrid2D rhs)
    {
        float r = 0;
        lhs.ForEachData((ref float value, int[] index) => {r += value * rhs[index[0], index[1]];});
        return r;
    }

    static public void Increment(FluidData.CellCenteredScalarGrid2D lhs, FluidData.CellCenteredScalarGrid2D rhs, float scale)
    {
        lhs.ForEachData((ref float value, int[] index) => { value += scale * rhs[index[0], index[1]];});
    }
    static public void ScaleAndIncrement(FluidData.CellCenteredScalarGrid2D lhs, FluidData.CellCenteredScalarGrid2D rhs, float scale)
    {
        lhs.ForEachData((ref float value, int[] index) => { value = value * scale + rhs[index[0], index[1]];});
    }
}
