/* calling syntax:
 * fspikecountsfit_5state_odesolver_euler_mainloop_5_2( ...
 *       spikecounts, A, states_eq, kon', koff', kon_B', koff_B', Btot', g, c0, niter, states, subseglist, stepsize);
 * all array parameters inputs (kon etc.) need to have each column corresponding to a data segment and each row corresponding to a parameter dimension
 *
 * compilation with gcc on linux: 
 * mex fspikecountsfit_5state_odesolver_backwardeuler_mainloop_5_2.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
 *
 * compilation with MSVC on windows:
 * mex fspikecountsfit_5state_odesolver_backwardeuler_mainloop_5_2.c COMPFLAGS="/openmp $COMPFLAGS"
 */
#include <math.h>
#include <string.h> /* for memcpy */
#include <omp.h>
#include "sbm_fit_mextools.h"  /* includes mex.h, matrix.h */

#define MAX(A,B) ( (A) < (B) ? (B):(A))

#define VAR_spikecounts prhs[0]
#define VAR_A prhs[1]
#define VAR_initial_state prhs[2]
#define VAR_kon prhs[3]
#define VAR_koff prhs[4]
#define VAR_kon_B prhs[5]
#define VAR_koff_B prhs[6]
#define VAR_Btot prhs[7]
#define VAR_g prhs[8]
#define VAR_c0 prhs[9]
#define VAR_niter prhs[10]
#define VAR_states prhs[11]
#define VAR_subseglist prhs[12]
#define VAR_stepsize prhs[13]

#define NINPUTS_MIN 14
#define NINPUTS_MAX 14
#define NOUTPUTS 0

/* this is a specialized ode solver for 4 binding steps (5 states) and 2 buffers */
void odesolve_5sb_backward_5_2(int T, double stepsize,
        double *states, double *spikecounts,
        double *initial_state, int niter,
        double g, double c0,
        double kon0, double kon1, double kon2, double kon3,
        double koff0, double koff1, double koff2, double koff3,
        double kon_B0, double kon_B1,
        double koff_B0, double koff_B1,
        double Btot0, double Btot1,
        double A
        )
{
    int i, k, offset;
    double c, cprev, M_c_c, z_c, a, denom;
    
    double rate0, rate1, rate2, rate3;
    double brate0, brate1;
    
    double s0, s1, s2, s3, s4;
    double sprev0, sprev1, sprev2, sprev3, sprev4;    
    double M_c_s0, M_c_s1, M_c_s2, M_c_s3, M_c_s4;
    double M_s_c0, M_s_c1, M_s_c2, M_s_c3, M_s_c4;
    double M_s_s_diag0, M_s_s_diag1, M_s_s_diag2, M_s_s_diag3, M_s_s_diag4;    
    double M_s_s_diagp10, M_s_s_diagp11, M_s_s_diagp12, M_s_s_diagp13;
    double z_s0, z_s1, z_s2, z_s3, z_s4;
    
    double b0, b1;
    double bprev0, bprev1;    
    double M_c_b0, M_c_b1;
    double M_b_c0, M_b_c1;
    double M_b_b0, M_b_b1;
    double z_b0, z_b1;
    
    /* initialize states */
    cprev = initial_state[0];
    sprev0 = initial_state[1];
    sprev1 = initial_state[2];
    sprev2 = initial_state[3];
    sprev3 = initial_state[4];
    sprev4 = initial_state[5];
    bprev0 = initial_state[6];
    bprev1 = initial_state[7];
    
    for (k = 0; k < T; k++) {
        
        cprev += spikecounts[k] * A; /* spiking */
        
        /* initialize based on previous states (could we do better some other way)? */
        c = cprev;
        s0 = sprev0; s1 = sprev1; s2 = sprev2; s3 = sprev3; s4 = sprev4;
        b0 = bprev0; b1 = bprev1;
        
        for (i = 0; i < niter; i++) { /* Newton iterations */
            
            M_c_c = 1.0 + stepsize * (g + kon0 * s0 + kon1 * s1 + kon2 * s2 + kon3 * s3 + kon_B0 * (Btot0 - b0) + kon_B1 * (Btot1 - b1));
            
            rate0 = kon0 * c * s0 - koff0 * s1;
            rate1 = kon1 * c * s1 - koff1 * s2;
            rate2 = kon2 * c * s2 - koff2 * s3;
            rate3 = kon3 * c * s3 - koff3 * s4;
            
            z_s0 = sprev0 - s0 + stepsize * (-rate0);
            z_s1 = sprev1 - s1 + stepsize * (rate0 - rate1);
            z_s2 = sprev2 - s2 + stepsize * (rate1 - rate2);
            z_s3 = sprev3 - s3 + stepsize * (rate2 - rate3);
            z_s4 = sprev4 - s4 + stepsize * (rate3);
            
            brate0 = kon_B0 * c * (Btot0 - b0) - koff_B0 * b0;
            brate1 = kon_B1 * c * (Btot1 - b1) - koff_B1 * b1;
            
            z_b0 = bprev0 - b0 + stepsize * (brate0);
            z_b1 = bprev1 - b1 + stepsize * (brate1);
            
            z_c = cprev - c - stepsize * ((c - c0) * g + rate0 + rate1 + rate2 + rate3 + brate0 + brate1);
            
            M_s_c0 = stepsize * kon0 * s0;
            M_s_c1 = stepsize * (kon1 * s1 - kon0 * s0);
            M_s_c2 = stepsize * (kon2 * s2 - kon1 * s1);
            M_s_c3 = stepsize * (kon3 * s3 - kon2 * s2);
            M_s_c4 = stepsize * (-kon3 * s3);
            
            M_c_s0 = stepsize * kon0 * c;
            M_c_s1 = stepsize * (kon1 * c - koff0);
            M_c_s2 = stepsize * (kon2 * c - koff1);
            M_c_s3 = stepsize * (kon3 * c - koff2);
            M_c_s4 = stepsize * (-koff3);
            
            M_s_s_diag0 = 1.0 + stepsize * kon0 * c;
            M_s_s_diag1 = 1.0 + stepsize * (kon1 * c + koff0);
            M_s_s_diag2 = 1.0 + stepsize * (kon2 * c + koff1);
            M_s_s_diag3 = 1.0 + stepsize * (kon3 * c + koff2);
            M_s_s_diag4 = 1.0 + stepsize * (koff3);
            
            M_s_s_diagp10 = -stepsize * koff0;
            M_s_s_diagp11 = -stepsize * koff1;
            M_s_s_diagp12 = -stepsize * koff2;
            M_s_s_diagp13 = -stepsize * koff3;
            
            M_c_b0 = -stepsize * (kon_B0 * c + koff_B0);
            M_c_b1 = -stepsize * (kon_B1 * c + koff_B1);
            
            M_b_c0 = -stepsize * kon_B0 * (Btot0 - b0);
            M_b_c1 = -stepsize * kon_B1 * (Btot1 - b1);
            
            M_b_b0 = 1.0 - M_c_b0;
            M_b_b1 = 1.0 - M_c_b1;
            
            /* solve linear system by Gaussian elimination */
            
            /* forward pass sets below-diagonal to zero and diagonal to one, except first row and column which change but remain dense */
            M_s_s_diagp10 /= M_s_s_diag0;
            z_s0          /= M_s_s_diag0;
            M_s_c0        /= M_s_s_diag0;
            
            a = -stepsize * kon0 * c; /* first below diagonal element of M */
            denom = M_s_s_diag1 - a * M_s_s_diagp10;
            z_s1   = (z_s1    - a *   z_s0) / denom;
            M_s_c1 = (M_s_c1  - a * M_s_c0) / denom;
            M_s_s_diagp11 /= denom;
            
            a = -stepsize * kon1 * c; /* second below diagonal element of M */
            denom = M_s_s_diag2 - a * M_s_s_diagp11;
            z_s2   = (z_s2    - a *   z_s1) / denom;
            M_s_c2 = (M_s_c2  - a * M_s_c1) / denom;
            M_s_s_diagp12 /= denom;
            
            a = -stepsize * kon2 * c; /* third below diagonal element of M */
            denom = M_s_s_diag3 - a * M_s_s_diagp12;
            z_s3   = (z_s3    - a *   z_s2) / denom;
            M_s_c3 = (M_s_c3  - a * M_s_c2) / denom;
            M_s_s_diagp13 /= denom;
                        
            a = -stepsize * kon3 * c; /* fourth below diagonal element of M */
            denom = M_s_s_diag4 - a * M_s_s_diagp13;
            z_s4   = (z_s4    - a *   z_s3) / denom;
            M_s_c4 = (M_s_c4  - a * M_s_c3) / denom;
            
            /* backward pass sets above-diagonal and M_c_s to zero */
            z_s3   -= M_s_s_diagp13 * z_s4;
            M_s_c3 -= M_s_s_diagp13 * M_s_c4;
            z_c    -= M_c_s4 * z_s4;
            M_c_c  -= M_c_s4 * M_s_c4;
            
            z_s2   -= M_s_s_diagp12 * z_s3;
            M_s_c2 -= M_s_s_diagp12 * M_s_c3;
            z_c    -= M_c_s3 * z_s3;
            M_c_c  -= M_c_s3 * M_s_c3;
            
            z_s1   -= M_s_s_diagp11 * z_s2;
            M_s_c1 -= M_s_s_diagp11 * M_s_c2;
            z_c    -= M_c_s2 * z_s2;
            M_c_c  -= M_c_s2 * M_s_c2;
            
            z_s0   -= M_s_s_diagp10 * z_s1;
            M_s_c0 -= M_s_s_diagp10 * M_s_c1;
            z_c    -= M_c_s1 * z_s1;
            M_c_c  -= M_c_s1 * M_s_c1;
            
            z_c    -= M_c_s0 * z_s0;
            M_c_c  -= M_c_s0 * M_s_c0;
            
            /* forward pass only for buffers, sets M_b_b to one M_c_b to zero */
            z_b0   /= M_b_b0;
            M_b_c0 /= M_b_b0;
            z_c   -= M_c_b0 *   z_b0;
            M_c_c -= M_c_b0 * M_b_c0;
            
            z_b1   /= M_b_b1;
            M_b_c1 /= M_b_b1;
            z_c   -= M_c_b1 *   z_b1;
            M_c_c -= M_c_b1 * M_b_c1;
            
            z_c /= M_c_c; /* set M_c_c to one */
            
            /* forward pass eliminates M_s_c */
            z_s0 -= M_s_c0 * z_c;
            z_s1 -= M_s_c1 * z_c;
            z_s2 -= M_s_c2 * z_c;
            z_s3 -= M_s_c3 * z_c;
            z_s4 -= M_s_c4 * z_c;
            
            /* forward pass eliminates M_b_c */
            z_b0 -= M_b_c0 * z_c;
            z_b1 -= M_b_c1 * z_c;
            
            /* we've now reduced M to the identity matrix, so \Delta[c; s; b] = [z_c; z_s; z_b] */
            c  += z_c;
            s0 += z_s0;
            s1 += z_s1;
            s2 += z_s2;
            s3 += z_s3;
            s4 += z_s4;
            b0 += z_b0;
            b1 += z_b1;
            
        }
        
        cprev = c;
        sprev0 = s0; sprev1 = s1; sprev2 = s2; sprev3 = s3; sprev4 = s4;
        bprev0 = b0; bprev1 = b1;
        
        /* store results */
        offset = k * 8;
        states[offset] = cprev;
        states[offset + 1] = sprev0;
        states[offset + 2] = sprev1;
        states[offset + 3] = sprev2;
        states[offset + 4] = sprev3;
        states[offset + 5] = sprev4;
        states[offset + 6] = bprev0;
        states[offset + 7] = bprev1;
        
    }
    
}


int calc_temparray_size(ns, nbuffers, nseg) {
    int nbytes_needed;
    nbytes_needed =            
            nseg * sizeof(int) + /* T */
            8 * sizeof(double *) * nseg; /* states, spikecounts, initial_state, kon, koff, Btot, kon_B, koff_B */
    return nbytes_needed;
}


void mexFunction(
        int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    int nbindingsteps, nbuffers, nsegments_used, nsegments_available, r, ii, nbytes_needed, ns, offset;
    double *subsegs;
    char *tmpdata;
    double **states, **spikecounts; /* outputs passed as inputs from matlab */
    double **initial_state, **kon, **koff, **Btot, **kon_B, **koff_B; /* vector parameters */
    double *g, *c0, *A; /* scalar parameters */    
    int *T; /* temporary arrays that we will allocate and deallocate */
    int niter;  /* number of newton iterations */
    double *stepsize;
    
    if (nrhs < NINPUTS_MIN || nrhs > NINPUTS_MAX || nlhs > NOUTPUTS) {
        
        mexErrMsgTxt("incorrect number of inputs or outputs");
        
    }
    
    niter = round2int(mxGetScalar(VAR_niter));
    nsegments_available = (int) mxGetNumberOfElements(VAR_states);
    nsegments_used = (int) mxGetNumberOfElements(VAR_subseglist);
    subsegs = RetrieveDoubleScalars(VAR_subseglist, nsegments_used);
    nbuffers  = (int) mxGetM(VAR_Btot);
    nbindingsteps = (int) mxGetM(VAR_kon);
    ns = nbindingsteps + 1;  /* number of binding states */
    
    stepsize = RetrieveDoubleScalars(VAR_stepsize, nsegments_available);
    
    if (subsegs == NULL || stepsize == NULL) {
        
        mexErrMsgTxt("Invalid input");
        
    }
    
    /* initialize temporary array */
    nbytes_needed = calc_temparray_size(ns, nbuffers, nsegments_available);
    tmpdata = (char *) mxMalloc(nbytes_needed);
    
    offset = 0;
    
    T = (int * ) (tmpdata + offset);
    offset += nsegments_available * sizeof(int);
    
    states = (double **) (tmpdata + offset);
    offset += nsegments_available * sizeof(double *);
    
    spikecounts = (double **) (tmpdata + offset);
    offset += nsegments_available * sizeof(double *);
    
    initial_state = (double **) (tmpdata + offset);
    offset += nsegments_available * sizeof(double *);
    
    kon = (double **) (tmpdata + offset);
    offset += nsegments_available * sizeof(double *);
    
    koff = (double **) (tmpdata + offset);
    offset += nsegments_available * sizeof(double *);
    
    Btot = (double **) (tmpdata + offset);
    offset += nsegments_available * sizeof(double *);
    
    kon_B = (double **) (tmpdata + offset);
    offset += nsegments_available * sizeof(double *);
    
    koff_B = (double **) (tmpdata + offset);
    offset += nsegments_available * sizeof(double *);
    
    if (
            RetrieveDoubleArrays(VAR_states, states, nsegments_available) < 0 ||
            RetrieveDoubleArrays(VAR_spikecounts, spikecounts, nsegments_available) < 0
            ) {
        
        mxFree((void *) tmpdata);
        mexErrMsgTxt("Failed to retrieve double arrays");
        
    }
    
    for (ii = 0; ii < nsegments_available; ii++) {
        
        /* no need for input checking, as that's already been done above */
        T[ii] = mxGetNumberOfElements(mxGetCell(VAR_spikecounts, ii));
        
    }
    
    g = RetrieveDoubleScalars(VAR_g, nsegments_available);
    c0 = RetrieveDoubleScalars(VAR_c0, nsegments_available);
    A = RetrieveDoubleScalars(VAR_A, nsegments_available);
    
    if (g == NULL || c0 == NULL) {
        
        mxFree((void *) tmpdata);
        mexErrMsgTxt("Failed to retrieve double scalar parameters");
        
    }
    
    if (
            RetrieveDoubleArrayCols(VAR_kon, kon, nsegments_available) < 0 ||
            RetrieveDoubleArrayCols(VAR_koff, koff, nsegments_available) < 0 ||
            RetrieveDoubleArrayCols(VAR_initial_state, initial_state, nsegments_available) < 0 ||
            RetrieveDoubleArrayCols(VAR_Btot, Btot, nsegments_available) < 0 ||
            RetrieveDoubleArrayCols(VAR_kon_B, kon_B, nsegments_available) < 0 ||
            RetrieveDoubleArrayCols(VAR_koff_B, koff_B, nsegments_available) < 0
            ) {
        
        mxFree((void *) tmpdata);
        mexErrMsgTxt("Failed to retrieve array columns");
        
    }
    
    omp_set_num_threads(32);
    
#pragma omp parallel for
    for (r = 0; r < nsegments_used; r++) { /* solve the ODEs */
        
        ii = round2int(subsegs[r]) - 1; /* convert from matlab indexing to C indexing */
        
        odesolve_5sb_backward_5_2(T[ii], stepsize[ii],
                states[ii], spikecounts[ii],
                initial_state[ii], niter,                
                g[ii], c0[ii],
                kon[ii][0], kon[ii][1], kon[ii][2], kon[ii][3],
                koff[ii][0], koff[ii][1], koff[ii][2], koff[ii][3],
                kon_B[ii][0], kon_B[ii][1],
                koff_B[ii][0], koff_B[ii][1],
                Btot[ii][0], Btot[ii][1],
                A[ii]);
        
    }
    
    /* free temporary array */
    mxFree((void *) tmpdata);
    
}