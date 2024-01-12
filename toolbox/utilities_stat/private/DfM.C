/* DfM.C Deviations from the Mean.
 * Syntax: DfM(X,Mu,XO,n,d); */

#include "mex.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{  
    int i,j;
    
    double *X = mxGetPr(prhs[0]);           /* Get the pointer to the data of X */
    double *Mu = mxGetPr(prhs[1]);          /* Get the pointer to the data of Mu */
    double *XO = mxGetPr(prhs[2]);          /* Get ... */
    int n =((int)(mxGetScalar(prhs[3])));   /* Get ... */
    int d =((int)(mxGetScalar(prhs[4])));   /* Get ... */
    

    for (i=0;i<d;i++)
        for (j=0;j<n;j++)
            *XO++ =-Mu[i] + *X++;
    return;
}