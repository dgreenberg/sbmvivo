#include <mex.h>
#include <matrix.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#define max(A,B) ( (A) < (B) ? (B):(A))
#define min(A,B) ( (A) > (B) ? (B):(A))
//Matlab Syntax:
//baseline = kinetic2dff_windowed_mainloop(r,p,w,init_perm,init_sort);

double calc_mean(double *d,int n);

int roundup(double x);

void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
//baseline = kinetic2dff_windowed_mainloop(r,p,w,init_perm,init_sort);
    // initialize
    double *data, *init_perm, *init_sorted_vals, *baseline;
    double p;
    int m,n,w,ws,i,j,c,prev,next;
    int lowest_val_ind;
    int oldest_val;
    int nvals_available;
    int nvals_to_use;
    double new_val,next_bl_val;    
    int *next_val_ind;
    int *prev_val_ind;
    double *win_vals;
    data = mxGetPr(prhs[0]);
    init_perm = mxGetPr(prhs[3]);
    init_sorted_vals = mxGetPr(prhs[4]);
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    //mexPrintf("%d\n",m);
    //mexPrintf("%d\n",n); 
    p = mxGetScalar(prhs[1]) / 100;
    w = (int) mxGetScalar(prhs[2]);
    ws = (w-1)/2; //w is odd
   
    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
    
    baseline = mxGetPr(plhs[0]);
    
    win_vals = (double *)mxMalloc(sizeof(double) * w);
    prev_val_ind = (int *)mxMalloc(sizeof(int) * w);
    next_val_ind = (int *)mxMalloc(sizeof(int) * w);
    for(c = 0; c < n ; c++) { //c is current column
        //initialize "linked list" via arrays
        lowest_val_ind = (int)init_perm[0 + c*ws] - 1;
        prev_val_ind[lowest_val_ind] = - 1; //min value from initial sort
        next = (int)init_perm[0 + c*ws] - 1; //matlab index -> C index
        for(i = 0; i < ws - 1; i++) {
            //mexPrintf("%d\n",i); 
            prev = next;
            next = (int)init_perm[i + 1 + c*ws] - 1; //matlab index -> C index
            next_val_ind[prev] = next;
            prev_val_ind[next] = prev;
            win_vals[prev] = init_sorted_vals[i + c*ws];
        }
        //max value from initial sort
        next_val_ind[next] = -1;
        win_vals[next] = init_sorted_vals[ws - 1 + c*ws];
        
        //main loop
        oldest_val = ws; //would be ws + 1 in matlab array indexing style
        for(i = 0; i < m; i++) {
            //mexPrintf("%d\n", i);
            nvals_available = 1 + min(i+ws,m - 1) - max(i - ws,0);
            nvals_to_use = roundup(p * (double) nvals_available);
            //mexPrintf("%d\n", nvals_to_use);
            if (i + ws < m) {
                new_val = data[i + ws + c*m];
                win_vals[oldest_val] = new_val;
            }
            if (i < ws + 1) {                   //add a value to the sorted list
                if (new_val <= win_vals[lowest_val_ind]) {
                    prev_val_ind[oldest_val] = -1;
                    next_val_ind[oldest_val] = lowest_val_ind;
                    prev_val_ind[lowest_val_ind] = oldest_val;
                    lowest_val_ind = oldest_val;
                }
                else {
                    next = lowest_val_ind;
                    while((next != -1) && (win_vals[next] < new_val)) {
                        prev = next;
                        next = next_val_ind[next];
                    }
                    if (next != -1)
                        prev_val_ind[next] = oldest_val;
                    next_val_ind[prev] = oldest_val;
                    prev_val_ind[oldest_val] = prev;
                    next_val_ind[oldest_val] = next;
                }                
            }            
            else if (i + ws < m) {          //replace a value in the sorted list
                next = next_val_ind[oldest_val];
                prev = prev_val_ind[oldest_val];
                if ((next != -1) && (win_vals[next] < new_val)) {
                    prev_val_ind[next] = prev; //possibly -1
                    if (prev != -1)
                        next_val_ind[prev] = next;
                    else
                       lowest_val_ind = next;
                    while((next != -1) && (win_vals[next] < new_val)) {
                        prev = next;
                        next = next_val_ind[next];
                    }
                    if (next != -1)
                        prev_val_ind[next] = oldest_val;
                    next_val_ind[prev] = oldest_val;
                    prev_val_ind[oldest_val] = prev;
                    next_val_ind[oldest_val] = next;
                }
                else if ((prev != -1) && (win_vals[prev] > new_val)) {
                    next_val_ind[prev] = next; //possibly -1
                    if (next != -1)
                        prev_val_ind[next] = prev;
                    while((prev != -1) && (win_vals[prev] > new_val)) {
                        next = prev;
                        prev = prev_val_ind[prev];
                    }
                    if (prev != -1)
                        next_val_ind[prev] = oldest_val;
                    else
                        lowest_val_ind = oldest_val;
                    prev_val_ind[oldest_val] = prev;
                    prev_val_ind[next] = oldest_val;
                    next_val_ind[oldest_val] = next;
                }
            }
            else  {                        //remove a value from the sorted list
                next = next_val_ind[oldest_val];
                prev = prev_val_ind[oldest_val];
                if (next != -1)
                    prev_val_ind[next] = prev;
                if (prev != -1)
                    next_val_ind[prev] = next;
                else
                    lowest_val_ind = next;
            }
            next_bl_val = 0;
            next = lowest_val_ind;
            for (j=0; j<nvals_to_use; j++ ) {
                next_bl_val += win_vals[next];
                next = next_val_ind[next];
            }
            next_bl_val /= nvals_to_use;
            baseline[i + c*m] = next_bl_val;
            oldest_val++;
            oldest_val = oldest_val % w; //loop around
        }
    }
    mxFree((void *) win_vals);
    mxFree((void *) next_val_ind);
    mxFree((void *) prev_val_ind);
}

int roundup(double x) {
    int ii = (int)x;
    if ((double)ii < x)
        ii++;
    return ii;
}
