#include "mex.h"
#include <math.h>
#include <stdbool.h>

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
    
    // Pre-processing to apply Arce (1998)
    
    // Arce(1998) step 1a: take the sign of weights
    double* sW = (double*)mxMalloc(n * sizeof(double));
    for(int i = 0; i < n; i++)
    {
        if (W[i] > 0.0) sW[i] = 1.0;
        else if (W[i] < 0.0) sW[i] = -1.0;
        else sW[i] = 0.0;
    }
    
    // Boolean array representing the position of the positive weights
    bool* sWb = (bool*)mxMalloc(n * sizeof(bool));
    for(int i = 0; i < n; i++)
    {
        sWb[i] = (sW[i] == 1.0);
    }
    
    // Arce(1998) step 1b: take the absolute value of the weights
    // Arce(1998) step 1c: take the "W-signed" observations
    for(int i = 0; i < n; i++)
    {
        D[i] = sW[i] * D[i];  // W-signed observations
        W[i] = fabs(W[i]);     // Absolute weights
    }
    
    // Local variable declaration and initialization.
    int left            = 0;
    int right           = n-1;
    int position        = -1;
    int k               = ceil(n*p)-1;
    int BleichOverton   = 1;
    
    double pivotD;
    double pivotW;
    bool pivotSWb;
    double buffer;
    bool bufferBool;
    double Le;
    
    // Compute sum of weights for tolerance calculation
    double sumW = 0.0;
    for(int i = 0; i < n; i++)
    {
        sumW += W[i];
    }
    
    // Scale p by sum of weights
    p = p * sumW;
    
    // Numerical tolerance for Bleich-Overton condition
    double tol = fabs(sumW) * 2.2204460492503131e-16 * 100.0;  // 100 * eps(sumW)
    
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
            
            // Swap sWb if signs are different (XOR logic)
            if (sWb[k] != sWb[right])
            {
                sWb[k] = !sWb[k];
                sWb[right] = !sWb[right];
            }
            
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
                    
                    // Swap sWb if signs are different
                    if (sWb[i] != sWb[position])
                    {
                        sWb[i] = !sWb[i];
                        sWb[position] = !sWb[position];
                    }
                    
                    position=position+1;
                }
            }
            
            D[right]    = D[position];
            D[position] = pivotD;
            
            W[right]    = W[position];
            W[position] = pivotW;
            
            // Swap sWb if signs are different
            if (sWb[right] != sWb[position])
            {
                sWb[right] = !sWb[right];
                sWb[position] = !sWb[position];
            }
            
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
        Le = 0.0;          // Le = sum(W(1:k-1));
        for(int j = 0; j < k; j++)
        {
            Le += W[j];
        }
        
        // Apply tolerance to Bleich-Overton condition
        if ((Le-p<=tol) && (p-Le-W[k]<=tol))
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
    
    // Post-processing, following Arce (1998)
    // Reconstruct original sign at position kstar
    // Maps: true->1, false->-1
    double original_sign = sWb[kstar-1] ? 1.0 : -1.0;
    
    kD = kD * original_sign;
    kW = kW * original_sign;
    
    // Free allocated memory
    mxFree(sW);
    mxFree(sWb);
    
    // Return the values in output
    plhs[0]=mxCreateDoubleScalar(kD);
    plhs[1]=mxCreateDoubleScalar(kW);
    plhs[2]=mxCreateDoubleScalar(kstar);
    return;
}
