/* all ideas and most code from Marsaglia & Tsang 2002 */
#include <math.h>
#include <mex.h>
#include <stdint.h>

#define SHR3 (jz=jsr, jsr^=(jsr<<13), jsr^=(jsr>>17), jsr^=(jsr<<5),jz+jsr)
#define UNI (.5 + (signed) SHR3*.2328306e-9)

void mexFunction(
        int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    
    uint32_t jz, jsr, n, j;
    float *w;    
    
    if (nrhs != 2 || nlhs != 1) {
        mexErrMsgTxt("must have two inputs and one output");
    }    
    if (mxGetClassID(prhs[0]) != mxSINGLE_CLASS) {
        mexErrMsgTxt("first input must be single");
    }
    if ((mxGetNumberOfElements(prhs[1]) != 1) | (mxGetClassID(prhs[1]) != mxUINT32_CLASS)) {
        mexErrMsgTxt("second input must be a uint32 scalar");
    }
    plhs[0] = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
    n = (uint32_t) mxGetNumberOfElements(prhs[0]);
    jsr = *((uint32_t *) mxGetData(prhs[1]));    
    
    w = (float *) mxGetData(prhs[0]);
    for (j = 0; j < n; j++) {
        w[j] = UNI;
    }
    *((uint32_t *) mxGetData(plhs[0])) = jsr; /* set the second input as the new seed */
    
}