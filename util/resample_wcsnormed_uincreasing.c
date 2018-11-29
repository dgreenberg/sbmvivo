#include <mex.h>
#include <stdint.h>
//Matlab Syntax:
//i = resample_wcsnormed_uincreasing(wcs, u)
//u should be between 0 and 1 and sorted
//wcs(end) should be 1

void resamplef_d(double *wcs, double *u, uint32_t *i, uint32_t n);
void resamplef_s(float *wcs, float *u, uint32_t *i, uint32_t n);

void mexFunction(
        int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    uint32_t n;
    int outdims[2];    
    uint32_t *i;
    if (nrhs != 2) {
        mexErrMsgTxt("must have two inputs");
    }
    if (mxGetClassID(prhs[0]) != mxGetClassID(prhs[1])) {
        mexErrMsgTxt("wcs and u must be the same numeric type");
    }
    if ((mxDOUBLE_CLASS != mxGetClassID(prhs[0])) & (mxSINGLE_CLASS != mxGetClassID(prhs[0]))) {
        mexErrMsgTxt("inputs of type of single or double only");
    }
    n = (uint32_t) mxGetNumberOfElements(prhs[0]);
    if ((double) n > pow(2.0, 32.0) - 2) {
        mexErrMsgTxt("max input size is 2^32 - 2 elements");
    }
    if ((n != mxGetNumberOfElements(prhs[1])) | (n == 0)) {
        mexErrMsgTxt("w and u must be the same size and non-empty");
    }
    if (nlhs == 0) {
        return;
    }
    outdims[0] = n;
    outdims[1] = 1;
    plhs[0] = mxCreateNumericArray(2, (const int *) outdims,
            mxUINT32_CLASS, mxREAL);
    i = (uint32_t *) mxGetData(plhs[0]);
    if (mxGetClassID(prhs[0]) == mxDOUBLE_CLASS) {
        resamplef_d((double *)mxGetData(prhs[0]), (double *)mxGetData(prhs[1]), i, n);
    } else { //single precision        
        resamplef_s((float *)mxGetData(prhs[0]), (float *)mxGetData(prhs[1]), i, n);
    }
}

void resamplef_d(double *wcs, double *u, uint32_t *i, uint32_t n) {
    uint32_t j = 0, k = 0; //j is an index into the output categories, k into the input weights    
    double nextu;
    for (j = 0; j < n; j++) {
        nextu = u[j];
        while ((wcs[k] < nextu) && (k < n - 1)) { //if the weights sum to 1 the second condition shouldn't be necessary but is kept just in case to prevent an infinite loop for bad input            
            k++;
        }
        i[j] = k + 1; //note that if the last weight added is e.g. the 3rd weight w[2], then i[j] will end up being 3. thus i is a matlab index, not a c index.
    }
}

void resamplef_s(float *wcs, float *u, uint32_t *i, uint32_t n) {
    uint32_t j = 0, k = 0; //j is an index into the output categories, k into the input weights    
    float nextu;
    for (j = 0; j < n; j++) {
        nextu = u[j];
        while ((wcs[k] < nextu) && (k < n - 1)) { //if the weights sum to 1 the second condition shouldn't be necessary but is kept just in case to prevent an infinite loop for bad input            
            k++;
        }
        i[j] = k + 1; //note that if the last weight added is e.g. the 3rd weight w[2], then i[j] will end up being 3. thus i is a matlab index, not a c index.
    }
}