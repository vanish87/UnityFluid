using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Assertions;

namespace UnityFluid
{
    public class PCGSolver : MonoBehaviour
    {
        [SerializeField] protected RenderTexture A;
        [SerializeField] protected RenderTexture x;
        [SerializeField] protected RenderTexture b;

        [SerializeField] protected RenderTexture z;
        [SerializeField] protected RenderTexture s;


        [SerializeField] protected Vector2Int matrixDimension;

        [SerializeField] protected ComputeShader solver;

        protected int kernal = -1;

        protected void Start()
        {
            Assert.IsTrue(matrixDimension.x == matrixDimension.y);
            this.A = new RenderTexture(matrixDimension.x, matrixDimension.y, 24, RenderTextureFormat.ARGBFloat);
            this.x = new RenderTexture(matrixDimension.x, 1 , 24, RenderTextureFormat.ARGBFloat);
            this.b = new RenderTexture(matrixDimension.x, 1, 24, RenderTextureFormat.ARGBFloat);


            this.A.enableRandomWrite = true;
            this.x.enableRandomWrite = true;
            this.b.enableRandomWrite = true;

            this.z = new RenderTexture(matrixDimension.x, 1, 24, RenderTextureFormat.ARGBFloat);
            this.s = new RenderTexture(matrixDimension.x, 1, 24, RenderTextureFormat.ARGBFloat);

            this.z.enableRandomWrite = true;
            this.s.enableRandomWrite = true;


            this.kernal = solver.FindKernel("PCGSolver");
        }

        protected void Update()
        {
            this.solver.SetTexture(this.kernal, "A", this.A);
            this.solver.SetTexture(this.kernal, "x", this.x);
            this.solver.SetTexture(this.kernal, "b", this.b);

            this.solver.Dispatch(this.kernal, this.matrixDimension.x, this.matrixDimension.y, 1);
        }
    }
}