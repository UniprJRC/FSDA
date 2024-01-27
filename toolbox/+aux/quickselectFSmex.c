#include "mex.h"
// #include "matrix.h"

// The gateway
void
        mexFunction (int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    
    // Input variables.
    double* A;
    int n,k;
    
    // Output variable.
    int kE=0;
    
    
    // Local variable declaration and initialization.
    int left=0;
    int right;
    
    int position=-1;
    
    double pivot;
    double buffer;
    
    int i;
    
    
    // Assign input variables to local variable.
    A=mxGetPr(prhs[0]);
    n=(int)mxGetScalar(prhs[1]);
    k=(int)mxGetScalar(prhs[2]);
    
    // Initialize the rigth sntinel.
    right=n-1;
    
    while ( position != k ) //((left < right) && ( position != k ))
    {
        pivot=A[k];
        A[k]=A[right];
        A[right]=pivot;
        
        position=left;
        
        for (i = left; i <= right; i++)
        {         
            if(A[i]<pivot)
            {      
                buffer=A[i];
                A[i]=A[position];
                A[position]=buffer;
                
                position=position+1;
            }
        }
        
        A[right]=A[position];
        A[position]=pivot;
        
        if  (position < k)
        {
            left  = position + 1;
        }
        
        else
        {
            right = position - 1;
        }
    }
    
    plhs[0]=mxCreateDoubleScalar(A[k]);
    return;
}
