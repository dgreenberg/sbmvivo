#include "mex.h"
#include "matrix.h"

int RetrieveDoubleArrayCols(const mxArray *arrayvar, double **doublevars, int ncols) {
    
    int s, nrows, ncols_present;
    double *d;
    
    if (mxGetClassID(arrayvar) != mxDOUBLE_CLASS) {
        
        mexPrintf("not a double array");
        return -1;
        
    }
    
    nrows = (int) mxGetM(arrayvar);
    ncols_present = (int) mxGetN(arrayvar);
    if (ncols_present < ncols) {
        
        mexPrintf("need at least %d columns, found %d\n", ncols, ncols_present);
        return -2;
        
    }
    
    
    d = mxGetData(arrayvar);
    
    for (s = 0; s < ncols; s++) {
        
        if (nrows == 0) {
            
            doublevars[s] = NULL;
            
        } else {
            
            doublevars[s] = d + s * nrows;
            
        }
        
    }
    
    return 1;
    
}


double *RetrieveDoubleScalars(const mxArray *vectorvar, int n) {
    
    if (mxGetClassID(vectorvar) != mxDOUBLE_CLASS) {
        
        mexPrintf("not a double array\n");
        return NULL;
        
    }
    if (!((mxGetM(vectorvar) == 1 || mxGetN(vectorvar) == 1))) {
        
        mexPrintf("not a vector\n");
        return NULL;
        
    }
    if (mxGetNumberOfElements(vectorvar) < n) {
        
        mexPrintf("array is too small\n");
        return NULL;
        
    }
    
    return mxGetData(vectorvar);
    
}


int RetrieveDoubleArrays(const mxArray *cellvar, double **doublevars, int n) {
    
    int s;
    mxArray *a;
    
    if (mxGetClassID(cellvar) != mxCELL_CLASS) {
        
        mexPrintf("not a cell array");
        return -1;
        
    }
    
    if (mxGetNumberOfElements(cellvar) < n) {
        
        mexPrintf("array is too small");
        return -2;
        
    }
    
    for (s = 0; s < n; s++) {
        
        a = mxGetCell(cellvar, s);
        if (a == NULL) {
            
            mexPrintf("invalid data");
            return -3;
            
        }
        if (mxGetClassID(a) != mxDOUBLE_CLASS) {
            
            mexPrintf("cell array must contain double arrays");
            return -4;
            
        }
        
        doublevars[s] = mxGetData(a);
        
    }
    
    return 1;
    
}


int round2int(double x) {
    
    int y;
    y = (int) x;
    if ((double) y + 0.5 < x) {
        y++;
    }
    return y;
    
}