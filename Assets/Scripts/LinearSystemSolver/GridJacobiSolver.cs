using UnityEngine;
using UnityEngine.Assertions;

public class GridJacobiSolver : MonoBehaviour
{
    //here  x b are 2D so that we can access them by grid index
    //we will do coefficient handling in compute shader 
    //because this kind of system has at most 4 coefficients in one cell
    //it is easy to combine them
    //A is not required because we build them in compute shader

    [SerializeField] protected RenderTexture x;
    [SerializeField] protected RenderTexture xn;
    [SerializeField] protected RenderTexture b;
    [SerializeField] protected RenderTexture difference;

    [SerializeField] protected Texture2D marker;

    [SerializeField] protected Vector2Int gridDimension;

    [SerializeField] protected ComputeShader solver;

    protected Color[] xCPU;
    protected Color[] bCPU;
    protected Color[] markerCPU;

    protected int kernal = -1;

    protected void Start()
    {
        this.x = new RenderTexture(gridDimension.x, gridDimension.y, 24, RenderTextureFormat.ARGBFloat);
        this.xn = new RenderTexture(gridDimension.x, gridDimension.y, 24, RenderTextureFormat.ARGBFloat);
        this.b = new RenderTexture(gridDimension.x, gridDimension.y, 24, RenderTextureFormat.ARGBFloat);
        this.difference = new RenderTexture(gridDimension.x, gridDimension.y, 24, RenderTextureFormat.ARGBFloat);

        this.x.enableRandomWrite = true;
        this.xn.enableRandomWrite = true;
        this.b.enableRandomWrite = true;
        this.difference.enableRandomWrite = true;

        this.x.Create();
        this.xn.Create();
        this.b.Create();
        this.difference.Create();

        this.kernal = solver.FindKernel("Solver");

        this.InitB();
    }

    protected void Solve()
    {
        var delta = Time.fixedDeltaTime;
        var density = 1;
        var cellSpace = 1;

        var scale = delta / (density * cellSpace * cellSpace);

        for (var i = 0; i < 50; ++i)
        {
            this.solver.SetTexture(this.kernal, "x", this.x);
            this.solver.SetTexture(this.kernal, "xn", this.xn);
            this.solver.SetTexture(this.kernal, "b", this.b);

            this.solver.SetTexture(this.kernal, "marker", this.marker);
            this.solver.SetVector("gridDimension", new Vector4(this.gridDimension.x, this.gridDimension.y));
            this.solver.SetVector("scale", new Vector4(scale, scale, scale));

            this.solver.Dispatch(this.kernal, this.gridDimension.x, this.gridDimension.y, 1);

            var t = this.x;
            this.x = this.xn;
            this.xn = t;
        }
    }

    protected void Verify()
    {
        var delta = Time.fixedDeltaTime;
        var density = 1;
        var cellSpace = 1;

        var scale = delta / (density * cellSpace * cellSpace);
        var k = solver.FindKernel("Verify");

        this.solver.SetTexture(k, "x", this.x);
        this.solver.SetTexture(k, "b", this.b);
        this.solver.SetTexture(k, "difference", this.difference);

        this.solver.SetTexture(k, "marker", this.marker);
        this.solver.SetVector("gridDimension", new Vector4(this.gridDimension.x, this.gridDimension.y));
        this.solver.SetVector("scale", new Vector4(scale, scale, scale));

        this.solver.Dispatch(k, this.gridDimension.x, this.gridDimension.y, 1);
    }

    protected void InitB()
    {
        var k = solver.FindKernel("Init");
        this.solver.SetTexture(k, "b", this.b);
        this.solver.Dispatch(k, this.gridDimension.x, this.gridDimension.y, 1);

    }

    protected void Update()
    {
        if(Input.anyKeyDown)
        {
            this.Solve();
            this.Verify();
        }
    }

    protected void OnDrawGizmos()
    {
        if(this.b != null) Gizmos.DrawGUITexture(new Rect(0, 0, gridDimension.x, gridDimension.y), this.b);
        if(this.x != null) Gizmos.DrawGUITexture(new Rect(0, gridDimension.y, gridDimension.x, gridDimension.y), this.x);
        if (this.difference != null) Gizmos.DrawGUITexture(new Rect(0, gridDimension.y*2, gridDimension.x, gridDimension.y), this.difference);
    }
}
