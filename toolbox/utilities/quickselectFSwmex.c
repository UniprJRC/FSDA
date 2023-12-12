#include "mex.h"
#include <math.h>

// The gateway
void
        mexFunction (int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    
    // Input variables.
    double* D;
    double* W;
    double p;
    int n;
    
    // Assign input variables to local variable.
    D = mxGetPr(prhs[0]);
    W = mxGetPr(prhs[1]);
    p = (double)mxGetScalar(prhs[2]);
    n = (int)mxGetScalar(prhs[3]);
    
    // Output variable.
    double kD;
    double kW;
    int kstar;
    
    // Local variable declaration and initialization.
    int left            = 0;
    int right           = n-1;
    int position        = -1;
    int k               = ceil(n*p)-1;
    int BleichOverton   = 1;
    
    double pivotD;
    double pivotW;
    double buffer;
    double Le;
    
    while ( BleichOverton )
    {
        position=-1;
        
        //The internal loop is like in quickselectFS %%  
        while ( position != k )
        {
            
            pivotD   = D[k];
            D[k]     = D[right];
            D[right] = pivotD;
            
            pivotW   = W[k];
            W[k]     = W[right];
            W[right] = pivotW;
            
            position = left;
            
            for (int i = left; i <= right; i++)
            {
                if(D[i]<pivotD)
                {
                    buffer = D[i];
                    D[i]   = D[position];
                    D[position] = buffer;
                    
                    buffer = W[i];
                    W[i]   = W[position];
                    W[position] = buffer;
                    
                    position=position+1;
                }
            }
            
            D[right]    = D[position];
            D[position] = pivotD;
            
            W[right]    = W[position];
            W[position] = pivotW;
            
            if  (position < k)
            {
                left  = position + 1;
            }
            else
            {
                right = position - 1;
            }
            
        }
        
        // Checks on weights extends Bleich-Overton %%
        double Le = 0;          // Le = sum(W(1:k-1));
        for(int j = 0; j < k; j++)
        {
            Le += W[j];
        }
        
        if ((Le-p<=0) && (p-Le-W[k]<=0))
        {
            // The condition is met: stop computation.
            kD    = D[k];
            kW    = W[k];
            kstar = k+1;
            BleichOverton = 0;
        }
        else
        //    The conditions not met: go back to quickselectFS with new
        //    sentinels - (k,n) or (1,k) - and new order statistics - 
        //    k+1 or k-1 (see point 8 and 9.).
        {
            if  (W[k] < 2*(p-Le))
            {
                // Need to add weight to reach the condition in point 5. 
                // Add an element (weight and data) to the left part.
                k     = k+1;
                left  = k;
                right = n-1;
            }
            else
            {
                // Here, we have that D(k,2)> 2*(p-Le). Need to remove an 
                // element (weight and data) from the right part.
                k     = k-1;
                left  = 0;
                right = k;
            }
        }
        
    }
    
    // Return the values in output
    plhs[0]=mxCreateDoubleScalar(kD);
    plhs[1]=mxCreateDoubleScalar(kW);
    plhs[2]=mxCreateDoubleScalar(kstar);
    return;
}

