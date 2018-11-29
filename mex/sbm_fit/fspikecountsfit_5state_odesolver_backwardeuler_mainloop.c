/* calling syntax:
 * fspikecountsfit_5state_odesolver_euler_mainloop( ...
 *       spikecounts, A, states_eq, kon', koff', kon_B', koff_B', Btot', g, c0, niter, states, subseglist, stepsize);
 * all array parameters inputs (kon etc.) need to have each column corresponding to a data segment and each row corresponding to a parameter dimension
 *
 * compilation with gcc on linux: 
 * mex fspikecountsfit_5state_odesolver_backwardeuler_mainloop.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
 *
 * compilation with MSVC on windows:
 * mex fspikecountsfit_5state_odesolver_backwardeuler_mainloop.c COMPFLAGS="/openmp $COMPFLAGS"
  */
#include <math.h>
#include <string.h> /* for memcpy */
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


/*this is the general ode solver. note that it does not do any memory (de)allocation,
 * all arrays must be pre-allocated*/
void odesolve_5sb_backward(int T, double stepsize,
        double *states, double *spikecounts,
        double *initial_state, int niter,
        int nbindingsteps, double *s, double *sprev, double *ds_dt,
        int nbuffers, double *b, double *bprev, double *db_dt,
        double *M_c_s, double *M_s_c,
        double *M_s_s_diag, double *M_s_s_diagm1, double *M_s_s_diagp1,
        double *M_c_b, double *M_b_c, double *M_b_b,
        double *z_s, double *z_b,
        double g, double c0,
        double *kon, double *koff,
        double *kon_B, double *koff_B, double *Btot,
        double A
        )
{
    
    int nstates, i, j, k, offset;
    double c, cprev, rate, dc_dt, M_c_c, z_c, a, denom;
    
    nstates = 1 + (nbindingsteps + 1) + nbuffers; /*c, s, b*/
    
    /* initialize states */
    cprev = initial_state[0];
    offset = 1;
    for (j = 0; j < nbindingsteps + 1; j++) {
        
        sprev[j] = initial_state[offset + j];
        
    }
    offset = 1 + (nbindingsteps + 1);
    for (j = 0; j < nbuffers; j++) {
        
        bprev[j] = initial_state[offset + j];
        
    }
    
    for (k = 0; k < T; k++) {
        
        cprev += spikecounts[k] * A; /* spiking */
        
        c = cprev;
        memcpy(s, sprev, (nbindingsteps + 1) * sizeof(double));
        memcpy(b, bprev, nbuffers * sizeof(double));
        
        for (i = 0; i < niter; i++) { /* Newton iterations */
            
            /* calculate derivatives and Jacobian */
            
            /* ligand */
            dc_dt = -(c - c0) * g;
            M_c_c = stepsize * g + 1.0;
            
            /* macromolecule */
            ds_dt[0] = 0.0;
            M_c_s[0] = 0.0;
            M_s_c[0] = 0.0;
            M_s_s_diag[0] = 1.0;
            for (j = 0; j < nbindingsteps; j++) {
                
                rate          = kon[j] * c * s[j] - koff[j] * s[j + 1];
                ds_dt[j]     -= rate;
                ds_dt[j + 1]  = rate;
                dc_dt        -= rate;
                
                M_c_c             +=  stepsize * kon[j] * s[j];
                M_c_s[j]          +=  stepsize * kon[j] * c;
                M_c_s[j + 1]       = -stepsize * koff[j];
                M_s_c[j]          +=  stepsize * kon[j] * s[j];
                M_s_c[j + 1]       = -stepsize * kon[j] * s[j];
                M_s_s_diag[j]     +=  stepsize * kon[j] * c;
                M_s_s_diag[j + 1]  =  stepsize * koff[j] + 1.0;
                M_s_s_diagm1[j]    = -stepsize * kon[j] * c;  /* 1 below diagonal */
                M_s_s_diagp1[j]    = -stepsize * koff[j];     /* 1 above diagonal */
                
            }
            
            /* buffers */
            for (j = 0; j < nbuffers; j++) {
                
                rate      = kon_B[j] * c * (Btot[j] - b[j]) - koff_B[j] * b[j];
                db_dt[j]  = rate;
                dc_dt    -= rate;
                
                M_c_c    +=  stepsize * kon_B[j] * (Btot[j] - b[j]);
                M_c_b[j]  = -stepsize * (kon_B[j] * c + koff_B[j]);
                M_b_c[j]  = -stepsize * kon_B[j] * (Btot[j] - b[j]);
                M_b_b[j]  =  1.0 - M_c_b[j];
                
            }
            
            /* set up RHS of linear system */
            z_c = cprev - c + stepsize * dc_dt;
            for (j = 0; j < nbindingsteps + 1; j++) {
                
                z_s[j] = sprev[j] - s[j] + stepsize * ds_dt[j];
                
            }
            for (j = 0; j < nbuffers; j++){
                
                z_b[j] = bprev[j] - b[j] + stepsize * db_dt[j];
                
            }
            
            /* for debugging purposes:
             *
             * if (k == 4 && i == 0) {
             * mexPrintf("\nM_c_c\t%f\tz_c\t%f\n",M_c_c,z_c);
             * for (j = 0; j < nbindingsteps + 1; j++) {
             * mexPrintf("\n%d\tM_s_s_diag\t%f\tM_s_c\t%f\tM_c_s\t%f", j, M_s_s_diag[j], M_s_c[j], M_c_s[j]);
             * if (j < nbindingsteps) {
             * mexPrintf("\tM_s_s_diagm1\t%f\tM_s_s_diagp1\t%f",M_s_s_diagm1[j],M_s_s_diagp1[j]);
             * }
             * }
             * mexPrintf("\n");
             * for (j = 0; j < nbuffers; j++) {
             * mexPrintf("\n%d\tM_b_b\t%f\tM_b_c\t%f\tM_c_b\t%f\tz_b\t%f",j, M_b_b[j], M_b_c[j], M_c_b[j], z_b[j]);
             * }
             * mexPrintf("\n");
             *
             * }
             */
            
            /* solve linear system by Gaussian elimination */
            
            /* forward pass sets below-diagonal to zero and diagonal to one, except first row and column which change but remain dense */
            M_s_s_diagp1[0] /= M_s_s_diag[0];
            z_s[0]          /= M_s_s_diag[0];
            M_s_c[0]        /= M_s_s_diag[0];
            for (j = 1; j < nbindingsteps; j++) {
                
                a = M_s_s_diagm1[j - 1];
                denom = M_s_s_diag[j] - a * M_s_s_diagp1[j - 1];
                z_s[j]   = (z_s[j]    - a *   z_s[j - 1]) / denom;
                M_s_c[j] = (M_s_c[j]  - a * M_s_c[j - 1]) / denom;
                M_s_s_diagp1[j] /= denom;
                
            }
            a = M_s_s_diagm1[nbindingsteps - 1];
            denom = M_s_s_diag[nbindingsteps] - a * M_s_s_diagp1[nbindingsteps - 1];
            z_s[nbindingsteps]   = (z_s[nbindingsteps]   - a *   z_s[nbindingsteps - 1]) / denom;
            M_s_c[nbindingsteps] = (M_s_c[nbindingsteps] - a * M_s_c[nbindingsteps - 1]) / denom;
            
            /* backward pass sets above-diagonal and M_c_s to zero */
            for (j = nbindingsteps; j != 0; j--) {
                
                /* subtract M_s_s_diagp1[j - 1] * (row j of M_*_s) from (row j -1 of M_*_s) */
                z_s[j - 1]   -= M_s_s_diagp1[j - 1] * z_s[j];
                M_s_c[j - 1] -= M_s_s_diagp1[j - 1] * M_s_c[j];
                
                /* subtract M_c_s[j] * (row j of M_*_s) from (M_c_s) */
                z_c   -= M_c_s[j] * z_s[j];
                M_c_c -= M_c_s[j] * M_s_c[j];
                
            }
            z_c   -= M_c_s[0] * z_s[0];
            M_c_c -= M_c_s[0] * M_s_c[0];
            
            /* forward pass only for buffers, sets M_b_b to one M_c_b to zero */
            for (j = 0; j < nbuffers; j++) {
                
                z_b[j]   /= M_b_b[j];
                M_b_c[j] /= M_b_b[j];
                
                z_c   -= M_c_b[j] *   z_b[j];
                M_c_c -= M_c_b[j] * M_b_c[j];
                
            }
            
            z_c /= M_c_c; /* set M_c_c to one */
            
            /* forward pass eliminates M_s_c */
            for (j = 0; j < nbindingsteps + 1; j++) {
                
                z_s[j] -= M_s_c[j] * z_c;
                
            }
            
            /* forward pass eliminates M_b_c */
            for (j = 0; j < nbuffers; j++) {
                
                z_b[j] -= M_b_c[j] * z_c;
                
            }
            
            /* we've now reduced M to the identity matrix, so \Delta[c; s; b] = [z_c; z_s; z_b] */
            c += z_c;
            for (j = 0; j < nbindingsteps + 1; j++) {
                
                s[j] += z_s[j];
                
            }
            for (j = 0; j < nbuffers; j++) {
                
                b[j] += z_b[j];
                
            }
            
        }
        
        cprev = c;
        memcpy(sprev, s, (nbindingsteps + 1) * sizeof(double));
        memcpy(bprev, b, nbuffers * sizeof(double));
        
        /* store results */
        offset = k * nstates;
        states[offset] = cprev;
        offset++;
        for (j = 0; j < nbindingsteps + 1; j++) {
            
            states[offset + j] = s[j];
            
        }
        offset += (nbindingsteps + 1);
        for (j = 0; j < nbuffers; j++) {
            
            states[offset + j] = b[j];
            
        }
        
    }
    
}


int calc_temparray_size(ns, nbuffers, nseg) {
    int nbytes_needed;
    nbytes_needed =
            7 * sizeof(double) * nseg * ns + /* s, sprev, ds_dt, M_c_s, M_s_c, M_s_s_diag, z_s */
            2 * sizeof(double) * nseg * (ns - 1) + /* M_s_s_diagm1, M_s_s_diagp1 */
            7 * sizeof(double) * nseg * nbuffers + /* b, bprev, db_dt, M_c_b, M_b_c, M_b_b, z_b */
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
    double *s, *sprev, *ds_dt, *b, *bprev, *db_dt, *M_c_s, *M_s_c, *M_c_b, *M_b_c, *M_b_b, *M_s_s_diag, *M_s_s_diagm1, *M_s_s_diagp1, *z_s, *z_b; /* temporary arrays that we will allocate and deallocate */
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
    
    s = (double *) (tmpdata + offset);
    offset += nsegments_available * ns * sizeof(double);
    
    sprev = (double *) (tmpdata + offset);
    offset += nsegments_available * ns * sizeof(double);
    
    ds_dt = (double *) (tmpdata + offset);
    offset += nsegments_available * ns * sizeof(double);
    
    M_c_s = (double *) (tmpdata + offset);
    offset += nsegments_available * ns * sizeof(double);
    
    M_s_c = (double *) (tmpdata + offset);
    offset += nsegments_available * ns * sizeof(double);
    
    M_s_s_diag = (double *) (tmpdata + offset);
    offset += nsegments_available * ns * sizeof(double);
    
    z_s = (double *) (tmpdata + offset);
    offset += nsegments_available * ns * sizeof(double);
    
    M_s_s_diagm1 = (double *) (tmpdata + offset);
    offset += nsegments_available * (ns - 1) * sizeof(double);
    
    M_s_s_diagp1 = (double *) (tmpdata + offset);
    offset += nsegments_available * (ns - 1) * sizeof(double);
    
    b = (double *) (tmpdata + offset);
    offset += nsegments_available * nbuffers * sizeof(double);
    
    bprev = (double *) (tmpdata + offset);
    offset += nsegments_available * nbuffers * sizeof(double);
    
    db_dt = (double *) (tmpdata + offset);
    offset += nsegments_available * nbuffers * sizeof(double);
    
    M_c_b = (double *) (tmpdata + offset);
    offset += nsegments_available * nbuffers * sizeof(double);
    
    M_b_c = (double *) (tmpdata + offset);
    offset += nsegments_available * nbuffers * sizeof(double);
    
    M_b_b = (double *) (tmpdata + offset);
    offset += nsegments_available * nbuffers * sizeof(double);
    
    z_b = (double *) (tmpdata + offset);
    offset += nsegments_available * nbuffers * sizeof(double);
    
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
    
#pragma omp parallel for
    for (r = 0; r < nsegments_used; r++) { /* solve the ODEs */
        
        ii = round2int(subsegs[r]) - 1; /* convert from matlab indexing to C indexing */
        
        odesolve_5sb_backward(T[ii], stepsize[ii],
                states[ii], spikecounts[ii],
                initial_state[ii], niter,
                nbindingsteps, s + ii * ns, sprev + ii * ns, ds_dt + ii * ns,
                nbuffers, b + ii * nbuffers, bprev + ii * nbuffers, db_dt + ii * nbuffers,
                M_c_s + ii * ns, M_s_c + ii * ns,
                M_s_s_diag + ii * ns, M_s_s_diagm1 + ii * (ns - 1), M_s_s_diagp1 + ii * (ns - 1),
                M_c_b + ii * nbuffers, M_b_c + ii * nbuffers, M_b_b + ii * nbuffers,
                z_s + ii * ns, z_b + ii * nbuffers,
                g[ii], c0[ii],
                kon[ii], koff[ii],
                kon_B[ii], koff_B[ii], Btot[ii],
                A[ii]);
        
    }
    
    /* free temporary array */
    mxFree((void *) tmpdata);
    
}