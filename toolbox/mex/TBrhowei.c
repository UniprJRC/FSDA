/*
 * TBrhowei.c — Fused Tukey biweight rho + scale update + weights
 *
 * [scale, weights] = TBrhowei(residuals, scale, c, kc)
 *
 * Computes in a single pass:
 *   1. rho values (Tukey bisquare)
 *   2. Updated scale = scale * sqrt(mean(rho) / kc)
 *   3. Weights = (1 - (u/c)^2)^2 for |u| <= c, 0 otherwise
 *      where u = residuals / updated_scale
 *
 * This eliminates two separate passes over the data and all temporary
 * array allocations that the M-code version requires.
 */

#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    double *res, scale, c, kc;
    double *weights, *pScale;
    mwSize n, i;
    double c2, c2_6, u, u2, uc2, rhosum, meanrho, t;

    /* Input validation */
    if (nrhs != 4)
        mexErrMsgIdAndTxt("FSDA:TBrhowei:nrhs", "Four inputs required.");
    if (nlhs != 2)
        mexErrMsgIdAndTxt("FSDA:TBrhowei:nlhs", "Two outputs required.");

    res   = mxGetPr(prhs[0]);
    scale = mxGetScalar(prhs[1]);
    c     = mxGetScalar(prhs[2]);
    kc    = mxGetScalar(prhs[3]);

    n = mxGetNumberOfElements(prhs[0]);

    c2   = c * c;
    c2_6 = c2 / 6.0;

    /* Pass 1: compute rho values and sum them */
    rhosum = 0.0;
    for (i = 0; i < n; i++) {
        u = res[i] / scale;
        if (fabs(u) <= c) {
            u2  = u * u;
            uc2 = u2 / c2;
            rhosum += (u2 / 2.0) * (1.0 - uc2 + uc2 * uc2 / 3.0);
        } else {
            rhosum += c2_6;
        }
    }

    /* Update scale */
    meanrho = rhosum / (double)n;
    scale = scale * sqrt(meanrho / kc);

    /* Pass 2: compute weights with updated scale */
    plhs[0] = mxCreateDoubleScalar(scale);
    plhs[1] = mxCreateDoubleMatrix(n, 1, mxREAL);
    weights  = mxGetPr(plhs[1]);

    for (i = 0; i < n; i++) {
        u = res[i] / scale;
        if (fabs(u) <= c) {
            t = 1.0 - (u / c) * (u / c);
            weights[i] = t * t;
        } else {
            weights[i] = 0.0;
        }
    }
}
