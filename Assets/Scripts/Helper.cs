using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Assertions;

public class Helper
{
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

    static public void GetIndexAndFraction(float pos, int low, int high, out int index, out float frac)
    {
        //TODO: Clamp index in range of low and high
        index = Mathf.FloorToInt(pos);
        frac = pos - index;

        Assert.IsTrue(low <= index && index <= high);
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

        GetIndexAndFraction(normalizedPos.x, low.x, high.x, out x, out fx);
        GetIndexAndFraction(normalizedPos.y, low.y, high.y, out y, out fy);

        index = new Vector2Int(x, y);
        frac =  new Vector2(fx, fy);
    }
}
