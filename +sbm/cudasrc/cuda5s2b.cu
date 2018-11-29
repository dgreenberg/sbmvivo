#include <iostream>
#include <curand_kernel.h> //need to include this before cuda5s.h so we don't redefine device prop. struct
#include "cuda5s2b.h"
#include <sys/time.h>
#include <cfloat>
#include <cub/cub.cuh>

#define RTYPE curandStateXORWOW
#define SHAREDMEMRESAMPLE true
#define GUARDFAC 1.0204
#define R_INIT_SD 4.0
#define GUARDMODE true
#define TIMECALC false
//#define logfunc __logf
#define logfunc logf
//#define expfunc __expf
#define expfunc expf

#ifndef INT_MIN
#define INT_MIN -32767
#endif

void cufreeifnotnull(void *&);
timespec tdiff(timespec start, timespec end);
void process_weights(gpu5s_problem *g, int t, int offset);
__global__ void resample_shared(const int * __restrict__ O, int np, int prevoffset, int offset,
		float4 * states_all, float4 * cbr_all, float *log_w, float *log_w_corrected, int *ancestor);
__global__ void advancepf(float fcorr, int nsteps, float stepsize, int n_newton_iterations,
		float4 * states_all, float4 * cbr_all, float * __restrict__ gp,
		bool * __restrict__ ns,
		float * __restrict__ log_wall, float * __restrict__ log_w_correctedall,
		const float * __restrict__ q_spike, RTYPE *rngstates, int np, bool resamplenow, float equal_log_w);

bool resamplenow;
float h_neff;

__device__ float max_log_w, sum_w, max_log_w_corrected, sum_w_corrected, sum_wsq, log_total_sum_weights, log_total_sum_weights_c, d_neff;
__device__ int mlw_total, mlw_c_total;
__device__ unsigned int retirementCount; //for use with threadfence reduction
__constant__ float FBGp1, basevar, gain, p_spike, c0, maxex, kd_ex, sigr_sqrtdt, S;
__constant__ float kon0, kon1, kon2, kon3, koff0, koff1, koff2, koff3;
__constant__ float kon_B0, koff_B0, kon_B1, koff_B1, Btot0, Btot1;
__constant__ float db1, db2, db3, db4;

/* #define CUCALL(val)    cucheck( (val), #val, __FILE__, __LINE__ )
void cucheck(cudaError_t  e, char const *const func, const char *const file, int const line) {
	if (e != cudaSuccess) {
		std::cerr << "CUDA error at " << file << ":" << line << " code=" << static_cast<unsigned int>(e) << " \"" << func << "\" \n";
		exit(EXIT_FAILURE);
	}
} */

//FIXME decouple reset / list functions from cuda5s code, compile seperately
int cudaDeviceList_wrapper(cudaDeviceProp *pDeviceList, int MaxDevices) {
	//pDeviceList is a pointer into a pre-initialized array of cudaDeviceProp structs of size MaxDevices
	//returns number of available devices, or a negative error code
	//for errors during calls to cudaGetDeviceProperties, the negative error code is decremented by 1000
	int count = 0;
	cudaError_t DeviceCountResult = cudaGetDeviceCount ( &count );

	if (DeviceCountResult != cudaSuccess) {
		return -((int) DeviceCountResult);
	}

	int n = MIN(count, MaxDevices); //number of devices to query
	for (int j = 0; j < n; j++) {

		cudaError_t DevicePropertiesResult = cudaGetDeviceProperties (pDeviceList + j, j);
		if (DevicePropertiesResult != cudaSuccess) {
			return -1000 - ((int) DevicePropertiesResult);
		}
	}
	return count;
}

void cudaDeviceReset_wrapper() {
	cudaDeviceReset();
}

void cudaMemcpy_d2h_wrapper(void *h, void *d, int nbytes) {
	cudaMemcpy(h, d, nbytes, cudaMemcpyDeviceToHost); //fixme: check for success here
}

void cudaMemcpy_h2d_wrapper(void *d, void *h, int nbytes) {
	cudaMemcpy(d, h, nbytes, cudaMemcpyHostToDevice); //fixme: check for success here
}

__global__ void fillarray_kernel(float *x, float v, int np) {
	int ii = threadIdx.x + blockIdx.x * BLOCKSIZE;
	while (ii < np) {
		x[ii] = v;
		ii += BLOCKSIZE * gridDim.x; //grid strides
	}
}

void fillarray(float *x, float v, int n) {
	int nblocks = (n % NV) ? n / NV + 1: n / NV;
	fillarray_kernel<<<nblocks, NT>>>(x, v, n);
}

__global__ void initancestors_noresample(int *ancestor, int np) {
	int ii = threadIdx.x + blockIdx.x * BLOCKSIZE;
	while (ii < np) {
		ancestor[ii] = ii; //note that the next time step is the same as K time steps back. it's ok to overwrite this since we've already copied out the relevant values as a_gs
		ii += BLOCKSIZE * gridDim.x;
	}
}

__global__ void W2O(float u, int np, const float * __restrict__ W, int * __restrict__ O) {
	//fixme should combine this kernel with a max-int-scan
	int ii = threadIdx.x + blockIdx.x * BLOCKSIZE;
	float Wfac = ((float) np) / W[np - 1];
	while (ii < np) {
		O[ii] = MIN(np, (int) (W[ii] * Wfac + u));
		ii += BLOCKSIZE * gridDim.x;
	}
}

__global__ void subtract_and_exponentiate_correctedonly(const float* __restrict__ log_w_corrected, float *w_corrected, int np) {
	//fixme should combine this kernel with an add-float-scan
	int ii = threadIdx.x + blockIdx.x * BLOCKSIZE;
	while (ii < np) {
		w_corrected[ii] = expfunc(log_w_corrected[ii] - max_log_w_corrected); //at present, w_corrected is used here only to store linear domain values for a float-add-reduction
		ii += BLOCKSIZE * gridDim.x;
	}
}

float gpu5s_marglik(gpu5s_problem *g) {
	float *ll = (float *) malloc(sizeof(float) * g->T);
	cudaMemcpy(ll, g->d.log_sum_raw_w, g->T * sizeof(float), cudaMemcpyDeviceToHost); //fixme check for success
	double s = 0.0;
	for (int ii = 0; ii < g->T; ii++) //fixme do this on the gpu
		s += (double) ll[ii];
	free(ll);
	return (float) s;
}

void FreeCubTemporaryArrays(gpu5s_problem *g) {
	cufreeifnotnull(g->d.d_temp_storage_float_max_reduction);
	cufreeifnotnull(g->d.d_temp_storage_float_add_reduction);
	cufreeifnotnull(g->d.d_temp_storage_float_add_scan);
	cufreeifnotnull(g->d.d_temp_storage_uint_max_scan);
}

int InitializeCubTemporaryArrays(gpu5s_problem *g) { //initialize temporary arrays for cub scans and reductions

	//dummy variables for the calls to determine temporary storage requirements
	float *pd_mlw_corrected = NULL;
	float *pd_sw_corrected = NULL;

	//make calls to reduction and scan functions with temp storage pointers set to NULL, to determine storage requirements.
	//no actual scan/reduction work is done by these calls
	//FIXME these temporary arrays should be allocated along with other device memory, not here in the main GPU-PF function
	if (	cub::DeviceReduce::Max(NULL, 		 g->d.temp_storage_bytes_float_max_reduction, g->d.log_w_corrected, pd_mlw_corrected, g->options.nparticles) != cudaSuccess ||
			cub::DeviceReduce::Sum(NULL, 		 g->d.temp_storage_bytes_float_add_reduction, g->d.w_corrected,     pd_sw_corrected,  g->options.nparticles) != cudaSuccess ||
			cub::DeviceScan::InclusiveSum(NULL,  g->d.temp_storage_bytes_float_add_scan,      g->d.w, 				g->d.W, 	      g->options.nparticles) != cudaSuccess ||
			cub::DeviceScan::InclusiveScan(NULL, g->d.temp_storage_bytes_uint_max_scan,       g->d.O, 				g->d.Oi, cub::Max(), g->options.nparticles)	!= cudaSuccess) {
		std::cerr << "failed to determine temporary array requirements for cub\n";
		return -1;
	}

	//allocate temporary arrays of the needed sizes
	if (	cudaMalloc((void **) &g->d.d_temp_storage_float_max_reduction, g->d.temp_storage_bytes_float_max_reduction) != cudaSuccess ||
			cudaMalloc((void **) &g->d.d_temp_storage_float_add_reduction, g->d.temp_storage_bytes_float_add_reduction) != cudaSuccess ||
			cudaMalloc((void **) &g->d.d_temp_storage_float_add_scan, 	   g->d.temp_storage_bytes_float_add_scan)  	!= cudaSuccess ||
			cudaMalloc((void **) &g->d.d_temp_storage_uint_max_scan, 	   g->d.temp_storage_bytes_uint_max_scan) 	    != cudaSuccess) {
		std::cerr << "failed to initialize one or more arrays requirements for cub\n";
		FreeCubTemporaryArrays(g);
		return -2;
	}
	return 1;
}

int runpfmainloop(gpu5s_problem *g){
	timespec ts, ts2, td, tsum, t0, t1, trand;
	if (TIMECALC)
		clock_gettime(CLOCK_REALTIME, &trand);
	int np = g->options.nparticles;
	int nresamples = 0;
	//get pointers to reduction outputs, which are global device memory variables
	float *pd_mlw_corrected, *pd_sw_corrected;
	cudaGetSymbolAddress((void **) &pd_mlw_corrected, max_log_w_corrected);
	cudaGetSymbolAddress((void **) &pd_sw_corrected, sum_w_corrected);

	//initialize reduction values
	float z = 0.0; float o = 1.0;
	cudaMemcpyToSymbol(max_log_w_corrected, &z,  sizeof(float));
	cudaMemcpyToSymbol(sum_w_corrected,     &o,  sizeof(float));
	unsigned int uz = 0;
	cudaMemcpyToSymbol(retirementCount,     &uz, sizeof(unsigned int));
	float foffset = g->params.fdc;

	int offset_Kbuffers; //offset into array that stores nparticles values buffered K times
	int offset_ancestors; //offset into array that stores ancestors, maxK + 1 columns (a power of 2)
	int offset_ns;  //offset into spiking, K * nsteps columns

	int prevoffset, offset = 0; //offset for the previous and current time steps into buffers with 2 * nparticles capacity
	float stepsize    = g->dt / ((float) g->options.nsteps);
	float equal_log_w = -logf((float) np);
	int n_newton_iterations = g->options.n_newton_iterations;

	cudaMemcpy(&h_neff, g->d.neff, sizeof(float), cudaMemcpyDeviceToHost); //neff for t = 0 should already be set
	resamplenow = h_neff < g->options.resamplethreshold; //do we need resampling now?
	if (!resamplenow && g->options.computenmean)
		initancestors_noresample<<<g->options.nblocks, NT>>>(g->d.ancestor + np, np);
	if (TIMECALC) {
		clock_gettime(CLOCK_REALTIME, &t0);
		tsum.tv_sec = 0; tsum.tv_nsec = 0;
	}
	for (int t = 1; t < g->T; t++) {

		offset_ancestors = (t & g->d.maxK) * np; //ancestors are always circularly buffered with K+1 columns so we can do fast modular arithmetic in calc_moments
		offset_Kbuffers = (t % (g->options.K + 1)) * np; //gp doesn't need to have a power of 2 number of columns since we don't need modular arithmetic on the GPU for it
		offset_ns = offset_Kbuffers * g->options.nsteps; //ns is circularly buffered with a number of columns that may not be a power of 2 to save memory since we don't need fast modular arithmetic to index its columns inside a GPU kernel.

		if (resamplenow) {
			nresamples++;
			prevoffset = offset;
			offset = np - offset;
			cub::DeviceScan::InclusiveSum(g->d.d_temp_storage_float_add_scan, g->d.temp_storage_bytes_float_add_scan, g->d.w, g->d.W, np);
			W2O<<<g->options.nblocks, BLOCKSIZE>>>(g->h.u[t], np, g->d.W, g->d.O);
			cub::DeviceScan::InclusiveScan(g->d.d_temp_storage_uint_max_scan, g->d.temp_storage_bytes_uint_max_scan,  g->d.O, g->d.Oi, cub::Max(), np);

			resample_shared<<<g->options.nblocks_rs, BLOCKSIZE_RS>>>(g->d.Oi, np, prevoffset, offset, (float4 *) g->d.states, (float4 *) g->d.cbr,
					g->d.log_w, g->d.log_w_corrected, g->d.ancestor + offset_ancestors);

			//need to normalize corrected weights to get back to a discrete probability distribution:
			cub::DeviceReduce::Max(g->d.d_temp_storage_float_max_reduction, g->d.temp_storage_bytes_float_max_reduction, g->d.log_w_corrected + offset, pd_mlw_corrected, np); //calculate max log weight
			subtract_and_exponentiate_correctedonly<<<g->options.nblocks, BLOCKSIZE>>>(g->d.log_w_corrected + offset, g->d.w_corrected, np); //subtract log of max weight and exponentiate
			cub::DeviceReduce::Sum(g->d.d_temp_storage_float_add_reduction, g->d.temp_storage_bytes_float_add_reduction, g->d.w_corrected, pd_sw_corrected,  np); //starting from raw corrected weights after resampling, we've now calculated their sum divided by their maximum
			//at this point, we haven't actually normalized the weights but we have the normalization constant which we'll use in later kernel calls
		}

		if (TIMECALC)
			clock_gettime(CLOCK_REALTIME, &ts);

		//advance the particle filter by one time step (i.e. nsteps substeps)
		advancepf<<<g->options.nblocks_pf, BLOCKSIZE_PF>>>(g->h.fobs[t] - foffset, g->options.nsteps, stepsize, n_newton_iterations,
				((float4 *) g->d.states) + offset, ((float4 *) g->d.cbr) + offset, g->d.gp + offset_Kbuffers,
				g->d.ns + offset_ns, //offset into spiking matrix that will be used to compute moments
				g->d.log_w, g->d.log_w_corrected + offset,
				g->d.q_spike + (t + g->options.ntimepoints_pre) * g->options.nsteps, (RTYPE *) g->d.rngstates, np, resamplenow, equal_log_w);

		//normalize weights, calculate their sums, etc.
		process_weights(g, t, offset);
		cudaMemcpyFromSymbol(&h_neff, d_neff, sizeof(float));

		if (TIMECALC) {
			clock_gettime(CLOCK_REALTIME, &ts2);
			td = tdiff(ts, ts2);
			tsum.tv_sec += td.tv_sec; tsum.tv_nsec += td.tv_nsec;
		}

	}

	if (TIMECALC) {
		clock_gettime(CLOCK_REALTIME, &t1);
		td = tdiff(t0, t1);
		std::cerr << "resampled on " << nresamples << " / " << g->T << " time steps\n";
		std::cerr << "total time: " << ((float) td.tv_sec + ((float) td.tv_nsec ) / 1000000000.0) << " s, ";
		std::cerr << "excluding resampling: " << ((float) tsum.tv_sec + ((float) tsum.tv_nsec ) / 1000000000.0) << " s\n";
	}

	return 1;
}

//this function numerically approximates the rate equation while randomly generating APs according to the particle filter's proposal distribution
//a backward (implicit) Euler method is used
__device__ float SimulateKinetics(float4 * __restrict__ pstates, float4 * __restrict__ pcbr,
		bool * __restrict__ ns,
		const float * __restrict__ q_spike,
		int n_newton_iterations,
		float stepsize,
		int nsteps, int np, int ii, //need np for the stride when assigning ns
		RTYPE *prngstate) {

	// load data:
	float s1prev = pstates->w;
	float s2prev = pstates->x;
	float s3prev = pstates->y;
	float s4prev = pstates->z;
	float s0prev = S - (s1prev + s2prev + s3prev + s4prev);  // calcium-free state is not explicitly stored
	float cprev = pcbr->w;
	float b0prev = pcbr->x;
	float b1prev = pcbr->y;

	float pq_spiking = 1.0;

	//local variables:
	int iter, jj;
	float c, a, denom, M_c_c, z_c;
	float rate0, rate1, rate2, rate3;
	float brate0, brate1;

	float s0, s1, s2, s3, s4;
	float M_c_s0, M_c_s1, M_c_s2, M_c_s3, M_c_s4;
    float M_s_c0, M_s_c1, M_s_c2, M_s_c3, M_s_c4;
    float M_s_s_diag0, M_s_s_diag1, M_s_s_diag2, M_s_s_diag3, M_s_s_diag4;
    float M_s_s_diagp10, M_s_s_diagp11, M_s_s_diagp12, M_s_s_diagp13;
    float z_s0, z_s1, z_s2, z_s3, z_s4;

    float b0, b1;
    float M_c_b0, M_c_b1;
    float M_b_c0, M_b_c1;
    float M_b_b0, M_b_b1;
    float z_b0, z_b1;
    float c_leak;

    if (USEKDEX) {
    	c_leak = maxex * c0 / (c0 + kd_ex); /* inward leak current */
    }

	for (jj = 0; jj < nsteps; jj++) {

		//sample spikes
		float next_q = q_spike[jj];
		bool spikenow = curand_uniform(prngstate) < next_q;
		if (spikenow) {

			pq_spiking *= (p_spike / next_q);
			cprev += 1.0;

		} else {

			pq_spiking *= ((1.0 - p_spike) / (1.0 - next_q));

		}
		ns[jj * np + ii] = spikenow;

		// initialize based on previous states (could we do better some other way)?
		c = cprev;
		s0 = s0prev; s1 = s1prev; s2 = s2prev; s3 = s3prev; s4 = s4prev;
		b0 = b0prev; b1 = b1prev;

		for (iter = 0; iter < n_newton_iterations; iter++) {

			if (USEKDEX) {
				M_c_c = 1.0 + stepsize * (maxex * kd_ex / ((c + kd_ex) * (c + kd_ex)) + kon0 * s0 + kon1 * s1 + kon2 * s2 + kon3 * s3 + kon_B0 * (Btot0 - b0) + kon_B1 * (Btot1 - b1));
			} else {
				M_c_c = 1.0 + stepsize * (maxex + kon0 * s0 + kon1 * s1 + kon2 * s2 + kon3 * s3 + kon_B0 * (Btot0 - b0) + kon_B1 * (Btot1 - b1));
			}

			rate0 = kon0 * c * s0 - koff0 * s1;
			rate1 = kon1 * c * s1 - koff1 * s2;
			rate2 = kon2 * c * s2 - koff2 * s3;
			rate3 = kon3 * c * s3 - koff3 * s4;

			z_s0 = s0prev - s0 + stepsize * (-rate0);
			z_s1 = s1prev - s1 + stepsize * (rate0 - rate1);
			z_s2 = s2prev - s2 + stepsize * (rate1 - rate2);
			z_s3 = s3prev - s3 + stepsize * (rate2 - rate3);
			z_s4 = s4prev - s4 + stepsize * (rate3);

			brate0 = kon_B0 * c * (Btot0 - b0) - koff_B0 * b0;
			brate1 = kon_B1 * c * (Btot1 - b1) - koff_B1 * b1;

			z_b0 = b0prev - b0 + stepsize * (brate0);
			z_b1 = b1prev - b1 + stepsize * (brate1);

			if (USEKDEX) {
				z_c = cprev - c - stepsize * (maxex * c / (c + kd_ex) - c_leak + rate0 + rate1 + rate2 + rate3 + brate0 + brate1);
			} else {
				z_c = cprev - c - stepsize * ((c - c0) * maxex + rate0 + rate1 + rate2 + rate3 + brate0 + brate1);
			}

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

			a = -stepsize * kon0 * c; // first below diagonal element of M
			denom = M_s_s_diag1 - a * M_s_s_diagp10;
			z_s1   = (z_s1    - a *   z_s0) / denom;
			M_s_c1 = (M_s_c1  - a * M_s_c0) / denom;
			M_s_s_diagp11 /= denom;

			a = -stepsize * kon1 * c; // second below diagonal element of M
			denom = M_s_s_diag2 - a * M_s_s_diagp11;
			z_s2   = (z_s2    - a *   z_s1) / denom;
			M_s_c2 = (M_s_c2  - a * M_s_c1) / denom;
			M_s_s_diagp12 /= denom;

			a = -stepsize * kon2 * c; // third below diagonal element of M
			denom = M_s_s_diag3 - a * M_s_s_diagp12;
			z_s3   = (z_s3    - a *   z_s2) / denom;
			M_s_c3 = (M_s_c3  - a * M_s_c2) / denom;
			M_s_s_diagp13 /= denom;

			a = -stepsize * kon3 * c; // fourth below diagonal element of M
			denom = M_s_s_diag4 - a * M_s_s_diagp13;
			z_s4   = (z_s4    - a *   z_s3) / denom;
			M_s_c4 = (M_s_c4  - a * M_s_c3) / denom;

			/* backward pass sets above-diagonal and M_c_s to zero */
			z_s3   -= M_s_s_diagp13 * z_s4;
			M_s_c3 -= M_s_s_diagp13 * M_s_c4;
			z_c   -= M_c_s4 * z_s4;
			M_c_c -= M_c_s4 * M_s_c4;

			z_s2   -= M_s_s_diagp12 * z_s3;
			M_s_c2 -= M_s_s_diagp12 * M_s_c3;
			z_c   -= M_c_s3 * z_s3;
			M_c_c -= M_c_s3 * M_s_c3;

			z_s1   -= M_s_s_diagp11 * z_s2;
			M_s_c1 -= M_s_s_diagp11 * M_s_c2;
			z_c   -= M_c_s2 * z_s2;
			M_c_c -= M_c_s2 * M_s_c2;

			z_s0   -= M_s_s_diagp10 * z_s1;
			M_s_c0 -= M_s_s_diagp10 * M_s_c1;
			z_c   -= M_c_s1 * z_s1;
			M_c_c -= M_c_s1 * M_s_c1;

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

		// update states
		cprev = c;
		s0prev = s0; s1prev = s1; s2prev = s2; s3prev = s3; s4prev = s4;
		b0prev = b0; b1prev = b1;

	}

	pstates->w = s1prev;
	pstates->x = s2prev;
	pstates->y = s3prev;
	pstates->z = s4prev;

	pcbr->w = cprev;
	pcbr->x = b0prev;
	pcbr->y = b1prev;

	return pq_spiking;

}


//this kernel advances states and weights one time step in the particle filter
__global__ void advancepf(float fcorr, int nsteps, float stepsize, int n_newton_iterations,
		float4 *states_all, float4 *cbr_all, float * __restrict__ gp,
		bool * __restrict__ ns,
		float * __restrict__ log_wall, float * __restrict__ log_w_correctedall,
		const float * __restrict__ q_spike, RTYPE *rngstates, int np, bool resamplenow, float equal_log_w) {
	int tid = threadIdx.x;
	int ii_base = tid + blockIdx.x * BLOCKSIZE_PF;
	float4 states, cbr;
	float pq_spiking, vFtotal, log_pobs, nextgpval;
	RTYPE rngstate = rngstates[ii_base];
	int ii = ii_base;
	__shared__ float log_w_corrected_offset;

	if (!threadIdx.x) {
		if (resamplenow) //to get normalized weights, we divide by max(w) * sum(w / max(w))
			log_w_corrected_offset = max_log_w_corrected + logfunc(sum_w_corrected);
		else
			log_w_corrected_offset = 0.0;
	}

	__syncthreads();

	while (ii < np) {

		states = states_all[ii];
		cbr = cbr_all[ii];

		pq_spiking = SimulateKinetics(&states, &cbr, ns, q_spike, n_newton_iterations, stepsize, nsteps, np, ii, &rngstate);
		cbr.z += curand_normal(&rngstate) * sigr_sqrtdt; // baseline drift

		states_all[ii] = states;
		cbr_all[ii] = cbr;

		nextgpval = expfunc(cbr.z) * (FBGp1 + (states.w * db1 + states.x * db2 + states.y * db3 + states.z * db4) / S);  // gain * photon flux

		gp[ii] = nextgpval;

		vFtotal = basevar + gain * nextgpval;
		log_pobs = -0.5 * ((nextgpval - fcorr) * (nextgpval - fcorr) / vFtotal + logfunc(((float) (2.0 * M_PI)) * vFtotal));

		if (resamplenow)
			log_wall[ii] = equal_log_w + log_pobs;
		else
			log_wall[ii] += log_pobs;
		log_w_correctedall[ii] += (log_pobs + logfunc(pq_spiking) - log_w_corrected_offset);
		ii += BLOCKSIZE_PF * gridDim.x;

	}
	rngstates[ii_base] = rngstate;
	if (!ii_base) { //reset weight maxima to save two cudamemcpytosymbol calls
		mlw_total = INT_MIN;
		mlw_c_total = INT_MIN;
	}
}


__global__ void InitializePF(
		float fcorr, float fcorr_approx, int nsteps, float stepsize, int n_newton_itrations, //fcorr is observed fluorescence - fdc
		float4 *states_all, float4 *cbr_all, float * gp,
		bool *ns_pre, bool *ns,
		float *log_wall, float *log_w_correctedall,
		const float * __restrict__ q_spike, RTYPE *rngstates,
		int np, int ntimepoints_pre,
		float4 states0, float4 cbr0,
		float mu_r_init, float var_r_init, float vFtotal_init) {

	int ii_base = threadIdx.x + blockIdx.x * BLOCKSIZE_PF;
	RTYPE rngstate = rngstates[ii_base];
	int ii = ii_base;
	float pq_spiking, b, rmaxexp, rmax, dfac, likvarr, likmeanr, v_r, m_r, r, log_priorp_r, log_samplingp_r, expr, vFtotal, log_pobs, log_w_total;

	float4 states, cbr;

	while (ii < np) {

		states = states0;
		cbr = cbr0;

		pq_spiking  = SimulateKinetics(&states, &cbr, ns_pre, q_spike, 							  n_newton_itrations, stepsize, nsteps * ntimepoints_pre, np, ii, &rngstate); //simulate up to one full dt before first measurement
		pq_spiking *= SimulateKinetics(&states, &cbr, ns,     q_spike + nsteps * ntimepoints_pre, n_newton_itrations, stepsize, nsteps,                   np, ii, &rngstate); //simulate up to time of first measurement
		states_all[ii] = states;

		//we now do conditional sampling of r given s and f
		b = FBGp1 + (states.w * db1 + states.x * db2 + states.y * db3 + states.z * db4) / S;

		//do a linear approximation of the exponential function at rmax:
		rmaxexp = fcorr_approx / b; //always positive
		rmax = logfunc(rmaxexp); //r with maximum likelihood given s, F - fdc =  fcorr_approx. this is where we'll do the linearization of exp(r)
		dfac = rmaxexp * b; //denominator factor for calculating a Gaussian approximation to P[F | r]. since rmax has been chosen to match F(1), this is usually F(1) - P.fdc
		likvarr = vFtotal_init / (dfac * dfac); //variance of a Gaussian approximation of likelihood P[F | r]
		likmeanr = fcorr / dfac + rmax - 1.0; //mean of a Gaussian approximation of likelihood P[F | r]. note that unlike fcorr_approx, fcorr can be negative

		v_r = 1.0 / (1.0 / likvarr + 1.0 / var_r_init);
		m_r = v_r * (likmeanr / likvarr + mu_r_init / var_r_init);

		r = m_r + sqrt(v_r) * curand_normal(&rngstate);
		log_priorp_r    = -0.5 * ((r - mu_r_init) * (r - mu_r_init) / var_r_init + logfunc(((float) (2.0 * M_PI)) * var_r_init));
		log_samplingp_r = -0.5 * ((r - m_r)       * (r - m_r)       / v_r        + logfunc(((float) (2.0 * M_PI)) * v_r));
		expr = expfunc(r);
		gp[ii] = expr * b; //gain * photon flux
		vFtotal = basevar + gain * gp[ii];

		log_pobs = -0.5 * ((gp[ii] - fcorr) * (gp[ii] - fcorr) / vFtotal + logfunc(((float) (2.0 * M_PI)) * vFtotal));
		log_w_total = log_pobs + log_priorp_r - log_samplingp_r;
		cbr.z = r;
		cbr_all[ii] = cbr;

		log_wall[ii] = log_w_total;
		log_w_correctedall[ii] = log_w_total + logfunc(pq_spiking);
		ii += BLOCKSIZE_PF * gridDim.x;

	}
	rngstates[ii_base] = rngstate;
	if (!ii_base) { //reset weight maxima to save two cudamemcpytosymbol calls
		mlw_total = INT_MIN;
		mlw_c_total = INT_MIN;
	}

}


void gpu5s_initialstates(gpu5s_problem *g) {
	//calculate the equilibrium solution to the rate equations at resting calcium:
	float Ka;
	float4 states0, cbr0;
	float c0 = g->params.c0;
	float S = g->params.S;
	float Btot0 = g->params.Btot0;
	float Btot1 = g->params.Btot1;

	cbr0.w = c0;
	if (c0 == 0.0) {

		states0.w = 0.0; states0.x = 0.0; states0.y = 0.0; states0.z = 0.0;
		cbr0.x = 0.0; cbr0.y = 0.0;

	} else {

		float v1 = 	    c0 * (g->params.kon0 / g->params.koff0);
		float v2 = v1 * c0 * (g->params.kon1 / g->params.koff1);
		float v3 = v2 * c0 * (g->params.kon2 / g->params.koff2);
		float v4 = v3 * c0 * (g->params.kon3 / g->params.koff3);
		float P = 1.0 + v1 + v2 + v3 + v4;
		states0.w = S * v1 / P;
		states0.x = S * v2 / P;
		states0.y = S * v3 / P;
		states0.z = S * v4 / P;

		Ka = g->params.kon_B0 / g->params.koff_B0;
		cbr0.x = Btot0 * c0 * Ka / (1.0 + c0 * Ka);
		Ka = g->params.kon_B1 / g->params.koff_B1;
		cbr0.y = Btot1 * c0 * Ka / (1.0 + c0 * Ka);

	}
	//std::cerr << "x0 " << x0 << " s0 " << s0 << " c0 " <<  g->params.c0 << "\n";

	float fmin = fabs(g->h.fobs[0] * 0.05); //minimum value of f from which we would initialize r etc., regardless of the value of fdc. FIXME kind of a hack
	float fcorr = g->h.fobs[0] - g->params.fdc; //can be negative due to noise, wrong parameters, etc.
	float fcorr_approx; //value of corrected (by fdc) fluorescence we'll use to determine the prior on r, linearize the exponential function, and calculate the laplace approximation of the observation function
	if (fcorr < fmin) {
		std::cerr << "Warning: initial corrected fluorescence is too low for standard initialization technique, using 5% of first uncorrected fluorescence value\n";
		fcorr_approx = fmin;
	} else {
		fcorr_approx = 0.95 * fcorr;
	}

	float beq = g->params.FBGp1 + (states0.w * g->params.db1 + states0.x * g->params.db2 + states0.y * g->params.db3 + states0.z * g->params.db4) / g->params.S;
	float mu_r_init = logf(fcorr_approx / beq);
	float var_r_init = log(R_INIT_SD) * log(R_INIT_SD); //prior variance on r, FIXME defs and units. also for 3s? why is this passed to the kernel? should be pure preprocessor defs.

	float vFtotal_init = g->params.vF + g->params.gain * fcorr_approx;
	float stepsize = g->dt / ((float) g->options.nsteps);
	unsigned int uz = 0;
	cudaMemcpyToSymbol(retirementCount,     &uz, sizeof(unsigned int));

	InitializePF<<<g->options.nblocks_pf, BLOCKSIZE_PF>>>(
			fcorr, fcorr_approx, g->options.nsteps, stepsize, g->options.n_newton_iterations,
			(float4 *) g->d.states, (float4 *) g->d.cbr, g->d.gp,
			g->d.ns_pre, g->d.ns,
			g->d.log_w, g->d.log_w_corrected,
			g->d.q_spike, (RTYPE *) g->d.rngstates,
			g->options.nparticles, g->options.ntimepoints_pre,
			states0, cbr0,
			mu_r_init, var_r_init, vFtotal_init);

	process_weights(g, 0, 0);
}


__global__ void resample_shared(const int * __restrict__ O, int np, int prevoffset, int offset,
		float4 *states_all, float4 *cbr_all, float *log_w, float *log_w_corrected, int *ancestor) {
	//this version of the resampling kernel uses shared memory to spread work out equally among threads, so we don't have to wait a long time for a single thread to do a lot of work.

	//note that we don't reorder the elements ns or or gp, but instead keep track of ancestors for them so we can calculate moments later
	int tid = threadIdx.x;
	int ii = threadIdx.x + blockIdx.x * BLOCKSIZE_RS;
	int prevO, jj, threadO, parent;
	float4 states_new, cbr_new;
	float lwcnew;
	__shared__ int sO[BLOCKSIZE_RS];
	__shared__ int blockwise_ancestor_index[BLOCKSIZE_RS];
	__shared__ int blockprevO, blockO;
	while (ii < np) {
		if (!ii) {
			prevO = 0;
		} else {
			prevO = O[ii - 1];
		}
		threadO = O[ii];
		sO[tid] = threadO; //copy particle counts to shared memory
		if (!tid) {
			blockprevO = prevO;
		} else if (tid == BLOCKSIZE_RS - 1) {
			blockO = threadO;
		}
		__syncthreads(); //broadcast updates to shared values to other threads
		int nvals_per_thread = (blockO - blockprevO + BLOCKSIZE_RS - 1) / BLOCKSIZE_RS; //divide blockwide number of values by number of threads in the block, and round up
		//we are now going to build a list of initial ancestors indices for each thread in the block.
		//each thread will write up to nvals_per_thread ancestors total.
		for (jj = (prevO - blockprevO + nvals_per_thread - 1) / nvals_per_thread; jj * nvals_per_thread < threadO - blockprevO; jj++)
			blockwise_ancestor_index[jj] = tid; //identify all threads which should start with a given blockwise ancestor index. some ancestors with low weights may not have any threads to give them child particles
		__syncthreads(); //broadcast updates to shared values to other threads
		int thread_offset = blockprevO + tid * nvals_per_thread; //initial offset into OUTPUT particles for this thread
		int thread_offset_max = MIN(blockprevO + (tid + 1) * nvals_per_thread, blockO); //final offset into OUTPUT particles for this thread, plus one
		int thread_blockwise_ancestor = blockwise_ancestor_index[tid]; //initial offset into this block's INPUT particles for this thread
		while (thread_offset < thread_offset_max) {
			parent = blockIdx.x * BLOCKSIZE_RS + thread_blockwise_ancestor;
			states_new = states_all[parent + prevoffset];
			cbr_new = cbr_all[parent + prevoffset];

			lwcnew = log_w_corrected[parent + prevoffset] - log_w[parent]; //log_w is not double buffered
			while (thread_offset < thread_offset_max && thread_offset < sO[thread_blockwise_ancestor]) {

				states_all[thread_offset + offset] = states_new;
				cbr_all[thread_offset + offset] = cbr_new;

				log_w_corrected[thread_offset + offset] = lwcnew; //log_w is not double buffered
				ancestor[thread_offset] = parent; //FIXME don't do this if we're not calculating any moments
				thread_offset++; //move on to the next value

			}
			while (thread_offset < thread_offset_max && thread_offset >= sO[thread_blockwise_ancestor]) //find the next non empty bin if we've already exhausted the present bin and haven't yet assigned all this thread's values.
				thread_blockwise_ancestor++;
		}
		ii += BLOCKSIZE_RS * gridDim.x; //grid stride loops
	}
	//FIXME now we need to do a reduction on the corrected weights, currently this is in another kernel
}


__device__ float floataddblockreduce_fromscalar(float x, float *sdata, int tid) {
	const int SecSize = NT / WARPSIZE;
	int lane = (SecSize - 1) & tid;    //lane within a section, not within a warp
	int sec = tid / SecSize;
	#pragma unroll
	for (int offset = 1; offset < SecSize; offset *= 2)  //values are further reduced within a section
		x += __shfl_down(x, offset, SecSize);
	if (!lane)
		sdata[sec] = x; //write each section's reduction into shared memory
	__syncthreads();
	if(tid < WARPSIZE) { //we now reduce the remaining values into a single block reduction, using a single warp
		x = sdata[tid];
	#pragma unroll
		for(int offset = 1; offset < WARPSIZE; offset *= 2) {
			x += __shfl_down(x, offset);
		}
		sdata[tid] = x;
	}
	__syncthreads();
	float blockreduction = sdata[0];
	//__syncthreads(); this line is present in MGPU, not really clear on why
	return blockreduction;
}

__device__ float floataddblockreduce(float *rdata, float *sdata, int tid) {
	float x = rdata[0];
	#pragma unroll
	for (int i = 1; i < VT; i++)
		x += rdata[i]; //each thread reduces VT values sequentially
	return floataddblockreduce_fromscalar(x, sdata, tid);
}

__device__ float floatmaxblockreduce_fromscalar(float x, float *sdata, int tid) {
	float xn;
	const int SecSize = NT / WARPSIZE;
	int lane = (SecSize - 1) & tid;    //lane within a section, not within a warp
	int sec = tid / SecSize;
	#pragma unroll
	for (int offset = 1; offset < SecSize; offset *= 2) { //values are further reduced within a section
		xn = __shfl_down(x, offset, SecSize);
		x = MAX(x, xn);
	}
	if (!lane)
		sdata[sec] = x; //write each section's reduction into shared memory
	__syncthreads();
	if(tid < WARPSIZE) { //we now reduce the remaining values into a single block reduction, using a single warp
		x = sdata[tid];
	#pragma unroll
		for(int offset = 1; offset < WARPSIZE; offset *= 2) {
			xn = __shfl_down(x, offset);
			x = MAX(x, xn);
		}
		sdata[tid] = x;
	}
	__syncthreads();
	float blockreduction = sdata[0];
	//__syncthreads(); this line is present in MGPU, not really clear on why
	return blockreduction;
}

__device__ float floatmaxblockreduce(float *rdata, float *sdata, int tid) {
	float x = rdata[0];
	#pragma unroll
	for (int i = 1; i < VT; i++)
		x = MAX(x, rdata[i]); //each thread reduces VT values sequentially
	x = floatmaxblockreduce_fromscalar(x, sdata, tid);
	return x;
}

__device__ void global2reg(float *g, float *r, int tid, int count, float initval) { //copies data from device global memory to per-thread register array
	if (count >= NT * VT) { //FIXME >= NT * VT - 1 ????
	#pragma unroll
		for (int i = 0; i < VT; i++)
			r[i] = g[NT * i + tid];
	} else {
	#pragma unroll
		for (int i = 0; i < VT; i++) {
			int index = NT * i + tid;
			r[i] = initval;
			if (index < count) r[i] = g[index];
		}
	}
}

//FIXME should template bool args
__device__ float log2linblocksum(int tid, float * sdata, float x0, float *r, float *g_x, float &sumsq, int count, bool writelin, bool reducesquare) {
	float sum = 0.0;
	for (int i = 0; i < VT; i++) {
		if (NT * i + tid < count) {
			float v = expfunc(r[i] - x0);
			sum += v;
			if (writelin)
				g_x[NT * i + tid] = v;
			if (reducesquare)
				sumsq += v * v;
		}
	}
	__syncthreads();
	sum = floataddblockreduce_fromscalar(sum, sdata, tid);
	if (reducesquare) {
		__syncthreads();
		sumsq = floataddblockreduce_fromscalar(sumsq, sdata, tid);
	}
	return sum;
}


__inline__ __device__ void findlastblock(bool &amLast) { //amLast should be a shared bool variable
	__threadfence(); //make sure all previous atomic operations (e.g. from reductions) are flushed before we take a ticket
	if (!threadIdx.x)
	{
		unsigned int ticket = atomicInc(&retirementCount, gridDim.x);
		// If the ticket ID is equal to the number of blocks, we are the last block!
		amLast = (ticket == gridDim.x - 1);
		if (amLast)
			retirementCount = 0;
	}
	__syncthreads(); //so that all threads get the updated value of the share bool amLast
}


__global__ void reduce_weights(float *lw, float *lw_c, float *w, float *w_c, int np,
		float *block_lw_max, float *block_sw_over_max, float *block_swsq_over_maxsq,
		float *block_lw_c_max, float *block_sw_c_over_max,
		int t, float *log_sum_raw_w, float *neff, bool calcmoments) {
	//this function and its subfunctions accomplish several tasks:
	//1) calculates logs of sums of weights (corrected and uncorrected)
	//2) calculates neff = 1 / sum(w * w), where w has been normalized (uncorrected weights only)
	//3) populates array of linear uncorrected weights divided by max within their block (log weights are NOT adjusted by this function from their raw values!)
	//4) stores max of log corrected/uncorrected weights for each block (sums per block of weights divide by max within the block are also stored, but not later used)
	//based in part on MGPU's reduction algorithm (https://github.com/NVlabs/moderngpu) and the cuda sample "threadfencereduction"
	int tid = threadIdx.x;
	__shared__ float sdata[WARPSIZE];
	int block_rw_offset = blockIdx.x * NT * VT; //offset for this block when reading from or writing to an array of size np in global device memory
	float swsq = 0.0, swsq_c; //hopefully swsq_c is optimized out as nothing is done to it
	float rdata[VT];

	global2reg(lw + block_rw_offset, rdata, tid, np - block_rw_offset, -FLT_MAX);
	__syncthreads();
	float mlw = floatmaxblockreduce(rdata, sdata, tid); //max of log weights for this block
	float sw   = log2linblocksum(tid, sdata, mlw,   rdata, w   + block_rw_offset, swsq,   np - block_rw_offset, true,        true); //we need sum of uncorrected weights and their sum of squares for resampling

	global2reg(lw_c + block_rw_offset, rdata, tid, np - block_rw_offset, -FLT_MAX);
	__syncthreads();
	float mlw_c = floatmaxblockreduce(rdata, sdata, tid); //max of log weights for this block
	float sw_c = log2linblocksum(tid, sdata, mlw_c, rdata, w_c + block_rw_offset, swsq_c, np - block_rw_offset, calcmoments, false); //we need corrected weights for moment calculation. we pass swsq_c but nothing should be done to it

	if (!tid) { //first thread of each block stores the results and calls atomic max functions
		block_lw_max[blockIdx.x] = mlw;
		block_sw_over_max[blockIdx.x] = sw;
		block_swsq_over_maxsq[blockIdx.x] = swsq;
		block_lw_c_max[blockIdx.x] = mlw_c;
		block_sw_c_over_max[blockIdx.x] = sw_c;
		//do an atomic max to store the truncated max logarithms of corrected/uncorrected weights weight
		atomicMax(&mlw_total, (int) mlw);
		atomicMax(&mlw_c_total, (int) mlw_c);
	}
	__shared__ bool amLast;
	findlastblock(amLast);
	if (amLast) { //the last block completes the reductions
		__shared__ float log_totalsw;
		float mlw_r = (float) mlw_total;
		float mlw_c_r = (float) mlw_c_total;
		float threadsw = 0.0;
		float threadsw_c = 0.0;
		for (int i = tid; i < gridDim.x; i += NT) { //whereas normally each thread reads VT values, here each thread reads nblocks / NT values
			threadsw   += expfunc(block_lw_max[i]   - mlw_r  ) * block_sw_over_max[i]; //each term is corrected so that it contributed the block's sum of weight divided by TOTAL max weight
			threadsw_c += expfunc(block_lw_c_max[i] - mlw_c_r) * block_sw_c_over_max[i];
		}
		__syncthreads();
		threadsw = floataddblockreduce_fromscalar(threadsw, sdata, tid); //ok if some threads haven't done anything in the above for loop since threadsw was zeroed.
		threadsw_c = floataddblockreduce_fromscalar(threadsw_c, sdata, tid); //ok if some threads haven't done anything in the above for loop since threadsw was zeroed.
		if (!tid)
			log_totalsw = logfunc(threadsw) + mlw_r; //log of sum of ALL weights
		__syncthreads();
		float log_totalsw_r = log_totalsw;
		float threadswsq = 0.0;
		for (int i = tid; i < gridDim.x; i += NT) {
			float qfac = expfunc(block_lw_max[i] - log_totalsw_r);
			threadswsq += qfac * qfac * block_swsq_over_maxsq[i];
		}
		__syncthreads();
		threadswsq = floataddblockreduce_fromscalar(threadswsq, sdata, tid); //ok if some threads haven't done anything in the above for loop since threadsw was zeroed.
		if(!tid) {
			log_total_sum_weights = log_totalsw_r;
			log_total_sum_weights_c = logfunc(threadsw_c) + mlw_c_r;
			log_sum_raw_w[t] = log_total_sum_weights_c;
			neff[t] = 1.0 / threadswsq;
			d_neff  = 1.0 / threadswsq;
		}
	}
}

//fixme template bool params
__global__ void normalize_weights(float *w, float *lw, float *w_c, float *lw_c, float *block_lw_max, float *block_lw_c_max, int np, bool resamplenow, bool calcmoments) {
	//we need to adjust weights, log weights and log corrected weights
	int tid = threadIdx.x;
	int block_rw_offset = blockIdx.x * NT * VT; //offset for this block when reading from or writing to an array of size np in global device memory
	float dlog = log_total_sum_weights; //calculate the necessary offset for log weights in this block
	float dlog_c = log_total_sum_weights_c; //calculate the necessary offset for log corrected weights in this block
	float dlin, dlin_c;
	if (calcmoments)
		dlin_c = expfunc(log_total_sum_weights_c - block_lw_c_max[blockIdx.x]); //we've already divided each weight by the block maximum
	if (resamplenow)
		dlin   = expfunc(log_total_sum_weights   - block_lw_max[blockIdx.x]); //we've already divided each weight by the block maximum
	if (np >= NT * VT + block_rw_offset) {
		for (int i = 0; i < VT; i++) {
			int index = block_rw_offset + NT * i + tid;
			lw[index] -= dlog;
			lw_c[index] -= dlog_c;
			if (calcmoments)
				w_c[index] /= dlin_c;
			if (resamplenow)
				w[index] /= dlin;
		}
	} else {
		for (int i = 0; i < VT; i++) {
			int index = block_rw_offset + NT * i + tid;
			if (index < np) {
				lw[index] -= dlog;
				lw_c[index] -= dlog_c;
				if (calcmoments)
					w_c[index] /= dlin_c;
				if (resamplenow)
					w[index] /= dlin;
			}
		}
	}
	if (threadIdx.x + blockIdx.x * BLOCKSIZE == 0) {
		//set max to 0 and sum to 1, to indicate that if we don't resample then no more normalization needs to be performed
		max_log_w_corrected = 0.0;
		sum_w_corrected = 1.0;
	}
}

//Fixme __restrict__ might help here?
__global__ void calc_moments(float *w_c, bool *ns, float *gp, int *ancestor,
		float *nmean, float *gpmean, float *gpsqmean,
		float *block_sw_spike, float *block_wgp, float * block_wgpsq,
		int t, int np, int maxK, int K, int nsteps, bool initnextancestors) {
	int ii_base = threadIdx.x + blockIdx.x * BLOCKSIZE;
	int ii = ii_base;
	int a;
	int a_gs[VT];
	int ancestor_offset = 0; //this lets the below code work for K = 0
	int tid = threadIdx.x;
	float w_gs[VT];
	float spikew;
	__shared__ float sdata[WARPSIZE];
	//store the ancestor K steps back for (up to) VT particles, going by grid strides
	for (int gs = 0; gs < VT; gs++) { //grid strides
		if (ii < np) {
			a = ii;
			for (int s = 0; s < K; s++) {
				ancestor_offset = ((t - s) & maxK) * np; //x & maxK == x % (maxK + 1) since maxK + 1 is a power of 2
				a = ancestor[ancestor_offset + a];
			}
			a_gs[gs] = a; //note this is a register array
			w_gs[gs] = w_c[ii]; //note this is a register array. we could in theory do the memory read only if spikes occured, not sure if that would speed things up.
		}
		ii += BLOCKSIZE * gridDim.x; //grid strides
	}

	//if the next time step doesn't require resampling we set ancestors in advance now and skip the resampling kernel entirely
	if (initnextancestors) {
		ancestor_offset = ((t + 1) & maxK) * np; //x & maxK == x % (maxK + 1) since maxK + 1 is a power of 2
		ii = ii_base;
		while (ii < np) {
			ancestor[ancestor_offset + ii] = ii; //note that when K==maxK, the next time step is the same as K time steps back. but that's ok since this function only reads spikes K time steps back, not ancestors.
			ii += BLOCKSIZE * gridDim.x;
		}
	}

	//now that we've determined ancestors K time steps back, we calculate the moments
	//first, calculate reductions for moments with each thread block

	//spiking moments:
	for (int jj = 0; jj < nsteps; jj++) { //loop over substeps within the time point K full steps back
		ii = ii_base;
		spikew = 0.0;
		for (int gs = 0; gs < VT; gs++) { //grid strides
			if (ii < np &&
					ns[jj * np + a_gs[gs]]) { //this second condition indicates a spike occurred
				spikew += w_gs[gs]; //add up weights of spiking particles over grid strides
			}
			ii += BLOCKSIZE * gridDim.x; //grid strides
		}
		__syncthreads(); //necessary?
		spikew = floataddblockreduce_fromscalar(spikew, sdata, tid); //add up the weights of all particles in this block that spiked for this substep
		if (!tid) //even if there we no spikes for this block, we still need to perform the global memory write to get the right answer from the reduction
			block_sw_spike[blockIdx.x + jj * gridDim.x] = spikew; //write the sum of weights with spikes for this block to global memory
	}

	if (gp != NULL) {
		//moments of gain * photon flux
		float wgp = 0.0, wgpsq = 0.0;
		ii = ii_base;
		for (int gs = 0; gs < VT; gs++) { //grid strides
			if (ii < np) {
				wgp   += w_gs[gs] * gp[a_gs[gs]];
				wgpsq += w_gs[gs] * gp[a_gs[gs]] * gp[a_gs[gs]];
			}
			ii += BLOCKSIZE * gridDim.x; //grid strides
		}
		__syncthreads(); //necessary?
		wgp   = floataddblockreduce_fromscalar(wgp, sdata, tid);
		//no __syncthreads() here, as there's one near the end of floataddblockreduce_fromscalar
		wgpsq = floataddblockreduce_fromscalar(wgpsq, sdata, tid);
		if (!tid) {
			block_wgp[blockIdx.x]   = wgp; //write the sum of weights with spikes for this block to global memory
			block_wgpsq[blockIdx.x] = wgpsq; //write the sum of weights with spikes for this block to global memory
		}
	}

	//if we're the last block, reduce over blocks
	__shared__ bool amLast;
	findlastblock(amLast);
	if (amLast) { //last block completes reduction over weights to calculate moments
		//reduction over blocks for spiking moments:
		for (int jj = 0; jj < nsteps; jj++) {
			float thread_sw_spike = 0.0;
			for (int kk = tid; kk < gridDim.x; kk += NT) //whereas normally each thread reads up to VT values, here each thread reads up to nblocks / NT values
				thread_sw_spike += block_sw_spike[kk + jj * gridDim.x];
			__syncthreads(); //necessary?
			thread_sw_spike = floataddblockreduce_fromscalar(thread_sw_spike, sdata, tid); //ok if some threads haven't done anything in the above for loop since thread_sw_spike was zeroed
			if (!tid)  //first thread of last block performs global device memory write
				nmean[jj] = thread_sw_spike; //sum of all weights for particles with spikes on this substep
		}

		if (gp != NULL) {
			//reduction over blocks for gain * photon flux
			float thread_wgp = 0.0, thread_wgpsq = 0.0;
			for (int kk = tid; kk < gridDim.x; kk += NT) { //whereas normally each thread reads up to VT values, here each thread reads up to nblocks / NT values
				thread_wgp   += block_wgp[kk];
				thread_wgpsq += block_wgpsq[kk];
			}
			__syncthreads(); //necessary?
			thread_wgp = floataddblockreduce_fromscalar(thread_wgp, sdata, tid);
			//no __syncthreads() here, as there's one near the end of floataddblockreduce_fromscalar
			thread_wgpsq = floataddblockreduce_fromscalar(thread_wgpsq, sdata, tid);
			if (!tid) { //first thread of last block performs global device memory write
				gpmean[0]   = thread_wgp; //sum of all weights for particles with spikes on this substep
				gpsqmean[0] = thread_wgpsq; //sum of all weights for particles with spikes on this substep
			}
		}
	}
}

void calc_moments_lastK(gpu5s_problem * g) {
	int t = g->T - 1;
	bool initnextancestors = false;
	for (int k = MIN(g->options.K - 1, g->T - 1); k >= 0; k--) {
		int np = g->options.nparticles;

		//get offsets into the arrays we're going to pass to calc_moments
		int gp_offset     = ((t - k) % (g->options.K + 1)) * np;
		int ns_offset     = ((t - k) % (g->options.K + 1)) * np  * g->options.nsteps;
		int gpmean_offset = t - k;
		int nmean_offset  = (t + g->options.ntimepoints_pre - k) * g->options.nsteps;

		calc_moments<<<g->options.nblocks, NT>>>(g->d.w_corrected, g->d.ns + ns_offset, g->d.gp + gp_offset, g->d.ancestor,
				g->d.nmean + nmean_offset, g->d.gpmean + gpmean_offset, g->d.gpsqmean + gpmean_offset,
				g->d.block_sw_spike, g->d.block_wgp, g->d.block_wgpsq,
				t, np, g->d.maxK, k, g->options.nsteps, initnextancestors);

	}
}


void process_weights(gpu5s_problem *g, int t, int offset) {
	int np = g->options.nparticles;
	reduce_weights<<<g->options.nblocks, NT>>>(g->d.log_w, g->d.log_w_corrected + offset, g->d.w, g->d.w_corrected, np,
			g->d.block_lw_max, g->d.block_sw_over_max, g->d.block_swsq_over_maxsq, g->d.block_lw_c_max, g->d.block_sw_c_over_max,
			t, g->d.log_sum_raw_w, g->d.neff, g->options.computenmean);
	cudaMemcpyFromSymbol(&h_neff, d_neff, sizeof(float));

	bool laststep = t == g->T - 1;
	resamplenow = (h_neff < g->options.resamplethreshold) && !laststep; //do we need resampling now?
	bool initnextancestors = !resamplenow && !laststep;
	normalize_weights<<<g->options.nblocks, NT>>>(g->d.w, g->d.log_w, g->d.w_corrected, g->d.log_w_corrected + offset, g->d.block_lw_max, g->d.block_lw_c_max, np, resamplenow, g->options.computenmean);
	if (g->options.computenmean) {
		if (t >= g->options.K) { //standard filter smoother K time points back

			//get offsets into the arrays we're going to pass to calc_moments
			int gp_offset     = ((t - g->options.K) % (g->options.K + 1)) * np;
			int gpmean_offset = t - g->options.K;
			int nmean_offset  = (t + g->options.ntimepoints_pre - g->options.K) * g->options.nsteps;
			int ns_offset     = ((t - g->options.K) % (g->options.K + 1)) * np * g->options.nsteps;

			calc_moments<<<g->options.nblocks, NT>>>(g->d.w_corrected, g->d.ns + ns_offset, g->d.gp + gp_offset, g->d.ancestor,
					g->d.nmean + nmean_offset, g->d.gpmean + gpmean_offset, g->d.gpsqmean + gpmean_offset,
					g->d.block_sw_spike, g->d.block_wgp, g->d.block_wgpsq,
					t, np, g->d.maxK, g->options.K, g->options.nsteps, initnextancestors);

			if (t == g->options.K) //filter smoother for states before the first observation's (full) time step
				calc_moments<<<g->options.nblocks, NT>>>(g->d.w_corrected, g->d.ns_pre,     NULL, g->d.ancestor,
						g->d.nmean, NULL, NULL,
						g->d.block_sw_spike_pre, NULL, NULL,
						t, np, g->d.maxK, g->options.K, g->options.nsteps * g->options.ntimepoints_pre, false);

			if (laststep)
				calc_moments_lastK(g); //calculate moments for the final K time steps, for which the filter smoother window will be less than K full time steps since we've run out of data

		} else if (initnextancestors)
			initancestors_noresample<<<g->options.nblocks, NT>>>(g->d.ancestor + ((t + 1) % (g->options.K + 1)) * np, np); //for t >= K this is done inside of calc_moments instead
	}
}


void pushparameterstodevice(gpu5s_problem *gpr) {
	//put parameter variables on device

	float prob_spike = gpr->params.lambda / ((float) gpr->options.nsteps);
	float sigr_sqrtdt_h = gpr->params.sigma_r * sqrtf(gpr->dt);
	cudaMemcpyToSymbol(sigr_sqrtdt, &sigr_sqrtdt_h, sizeof(float));
	cudaMemcpyToSymbol(p_spike, &prob_spike, sizeof(float));

	cudaMemcpyToSymbol(db1, &(gpr->params.db1), sizeof(float));
	cudaMemcpyToSymbol(db2, &(gpr->params.db2), sizeof(float));
	cudaMemcpyToSymbol(db3, &(gpr->params.db3), sizeof(float));
	cudaMemcpyToSymbol(db4, &(gpr->params.db4), sizeof(float));

	cudaMemcpyToSymbol(kon0, &(gpr->params.kon0), sizeof(float));
	cudaMemcpyToSymbol(kon1, &(gpr->params.kon1), sizeof(float));
	cudaMemcpyToSymbol(kon2, &(gpr->params.kon2), sizeof(float));
	cudaMemcpyToSymbol(kon3, &(gpr->params.kon3), sizeof(float));

	cudaMemcpyToSymbol(koff0, &(gpr->params.koff0), sizeof(float));
	cudaMemcpyToSymbol(koff1, &(gpr->params.koff1), sizeof(float));
	cudaMemcpyToSymbol(koff2, &(gpr->params.koff2), sizeof(float));
	cudaMemcpyToSymbol(koff3, &(gpr->params.koff3), sizeof(float));

	cudaMemcpyToSymbol(kon_B0,  &(gpr->params.kon_B0), sizeof(float));
	cudaMemcpyToSymbol(koff_B0, &(gpr->params.koff_B0), sizeof(float));
	cudaMemcpyToSymbol(kon_B1,  &(gpr->params.kon_B1), sizeof(float));
	cudaMemcpyToSymbol(koff_B1, &(gpr->params.koff_B1), sizeof(float));

	cudaMemcpyToSymbol(FBGp1,   &(gpr->params.FBGp1), sizeof(float));
	cudaMemcpyToSymbol(basevar, &(gpr->params.vF), sizeof(float));
	cudaMemcpyToSymbol(gain,    &(gpr->params.gain), sizeof(float));

	cudaMemcpyToSymbol(c0,      &(gpr->params.c0), sizeof(float));
	cudaMemcpyToSymbol(maxex,   &(gpr->params.maxex), sizeof(float));
	cudaMemcpyToSymbol(kd_ex,   &(gpr->params.kd_ex), sizeof(float));

	cudaMemcpyToSymbol(S,     &(gpr->params.S),     sizeof(float));
	cudaMemcpyToSymbol(Btot0, &(gpr->params.Btot0), sizeof(float));
	cudaMemcpyToSymbol(Btot1, &(gpr->params.Btot1), sizeof(float));

}

__global__ void setup_kernel(RTYPE *rngstates, int nrng, unsigned long long seedval)
{
	int id = threadIdx.x + blockIdx.x * BLOCKSIZE_PF;
	/* Each thread gets same seed, a different sequence number, no offset */
	if (id < nrng)
		curand_init(seedval, id, 0, &rngstates[id]);
}

void cudaseedrng(int nblocks_pf, void *rngstates, int nrng, unsigned long long seedval){
	setup_kernel<<<nblocks_pf, BLOCKSIZE_PF>>>((RTYPE *) rngstates, nrng, seedval);
}

int cudainitrng(int nblocks_pf, void *&rngstates, unsigned long long seedval) {
	int nrng = nblocks_pf * BLOCKSIZE_PF; //one rng per thread, regardless of how much data each thread processes on various kernels
	if (cudaMalloc(&rngstates, nrng * sizeof(RTYPE)) != cudaSuccess)
		return -100;
	cudaseedrng(nblocks_pf, rngstates, nrng, seedval);
	return 1;
}

timespec tdiff(timespec start, timespec end)
{
	timespec temp;
	if ((end.tv_nsec-start.tv_nsec)<0) {
		temp.tv_sec = end.tv_sec-start.tv_sec-1;
		temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
	} else {
		temp.tv_sec = end.tv_sec-start.tv_sec;
		temp.tv_nsec = end.tv_nsec-start.tv_nsec;
	}
	return temp;
}

int setgridsize(gpu5s_problem *g, int np) {
	if (np > g->d.maxnparticles)
		return -1;
	int nblocks = (np % NV) ? np / NV + 1 : np / NV;
	if (nblocks > g->d.maxnblocks)
		return -2;
	g->options.nparticles = np;
	g->options.nblocks = nblocks;
	g->options.nblocks_pf = (np % NV_PF) ? np / NV_PF + 1 : np / NV_PF;
	g->options.nblocks_rs = (np % NV_RS) ? np / NV_RS + 1 : np / NV_RS;
	return 1;
}

int allocate_gpupfdata(gpu5s_problem *g, int maxT, int np, int maxnsteps, int maxK, int maxKp1substeps,
		int maxtotalsubsteps, int maxpresubsteps, unsigned long long seedval) {
	if (g == NULL)
		return -1;
	if (maxK <= 0 || ((maxK + 1) & maxK ))
		return -2; //maxK + 1 must be positive and a power of 2 so we can do fast modular arithmetic using bitwise AND
	g->d.maxT = maxT;
	g->d.maxK = maxK;
	g->d.maxnsteps = maxnsteps;
	g->d.maxKp1substeps = maxKp1substeps;
	g->d.maxtotalsubsteps = maxtotalsubsteps;
	g->d.maxpresubsteps = maxpresubsteps;
	g->d.maxnparticles = np;
	g->d.maxnblocks = (np % NV) ? np / NV + 1 : np / NV;

	g->options.vt = VT; g->options.vt_pf = 1; g->options.vt_rs = VT_RS;
	g->options.nt = BLOCKSIZE; g->options.nt_pf = BLOCKSIZE_PF; g->options.nt_rs = BLOCKSIZE_RS;
	//initialize options to max array sizes
	if (setgridsize(g, np) < 0) {
		free_gpupfdata(g);
		return -3;
	}
	g->options.K = maxK;
	//set all CPU/GPU pointers to NULL
	g->h.fobs = NULL; g->h.u = NULL;
	g->d.rngstates = NULL;
	g->d.d_temp_storage_float_max_reduction = NULL; g->d.d_temp_storage_float_add_reduction = NULL;
	g->d.d_temp_storage_float_add_scan = NULL; g->d.d_temp_storage_uint_max_scan = NULL;
	g->d.states = NULL; g->d.cbr = NULL;
	g->d.ns = NULL; g->d.ns_pre = NULL;
	g->d.log_w = NULL; g->d.w = NULL; g->d.log_w_corrected = NULL; g->d.w_corrected = NULL; g->d.W = NULL;
	g->d.block_lw_max = NULL; g->d.block_sw_over_max = NULL; g->d.block_swsq_over_maxsq = NULL; g->d.block_lw_c_max = NULL; g->d.block_sw_c_over_max = NULL;
	g->d.block_sw_spike = NULL; g->d.block_sw_spike_pre = NULL; g->d.block_wgp = NULL; g->d.block_wgpsq = NULL;
	g->d.log_sum_raw_w = NULL; g->d.neff = NULL;
	g->d.q_spike = NULL; g->d.log_pq_spike = NULL; g->d.log_pq_nospike = NULL;
	g->d.nmean = NULL; g->d.gpmean = NULL; g->d.gpsqmean = NULL;
	g->d.O = NULL; g->d.Oi = NULL; g->d.ancestor = NULL;

	int npf  = np  			    * sizeof(float);
	int npi  = np   		    * sizeof(int);
	int Tf   = maxT 			* sizeof(float);
	int nsf  = maxtotalsubsteps * sizeof(float);
	int nbf  = g->d.maxnblocks  * sizeof(float);

	if (
			cudaMalloc((void **) &(g->d.states), 2 * np * sizeof(float4))				!= cudaSuccess ||
			cudaMalloc((void **) &(g->d.cbr),    2 * np * sizeof(float4))				!= cudaSuccess ||
			cudaMalloc((void **) &(g->d.gp), 	npf * (maxK + 1)) 						!= cudaSuccess ||
			cudaMalloc((void **) &(g->d.ns),     maxKp1substeps * np * sizeof(bool)) 	!= cudaSuccess ||
			cudaMalloc((void **) &(g->d.ns_pre), maxpresubsteps * np * sizeof(bool))	!= cudaSuccess ||
			cudaMalloc((void **) &(g->d.log_w), npf) 									!= cudaSuccess || //weights (only raw corrected weights are double buffered)
			cudaMalloc((void **) &(g->d.w),     npf) 									!= cudaSuccess ||
			cudaMalloc((void **) &(g->d.log_w_corrected), 2 * npf) 						!= cudaSuccess ||
			cudaMalloc((void **) &(g->d.w_corrected),     npf) 							!= cudaSuccess ||
			cudaMalloc((void **) &(g->d.W), npf) 					 					!= cudaSuccess ||
			cudaMalloc((void **) &(g->d.block_lw_max), 			nbf) 					!= cudaSuccess ||
			cudaMalloc((void **) &(g->d.block_sw_over_max),		nbf) 					!= cudaSuccess ||
			cudaMalloc((void **) &(g->d.block_swsq_over_maxsq), nbf) 					!= cudaSuccess ||
			cudaMalloc((void **) &(g->d.block_lw_c_max), 		nbf) 					!= cudaSuccess || //block reduction arrays
			cudaMalloc((void **) &(g->d.block_sw_c_over_max), 	nbf) 					!= cudaSuccess ||
			cudaMalloc((void **) &(g->d.block_sw_spike),     	nbf * maxnsteps) 		!= cudaSuccess ||
			cudaMalloc((void **) &(g->d.block_sw_spike_pre), 	nbf * maxpresubsteps) 	!= cudaSuccess ||
			cudaMalloc((void **) &(g->d.block_wgp), 		 	nbf)			 		!= cudaSuccess ||
			cudaMalloc((void **) &(g->d.block_wgpsq), 		 	nbf)			 		!= cudaSuccess ||
			cudaMalloc((void **) &(g->d.log_sum_raw_w), Tf)  							!= cudaSuccess || //T length arrays
			cudaMalloc((void **) &(g->d.neff), Tf)           							!= cudaSuccess ||
			cudaMalloc((void **) &(g->d.q_spike), nsf) 									!= cudaSuccess || //spiking proposal
			cudaMalloc((void **) &(g->d.log_pq_spike), nsf) 							!= cudaSuccess ||
			cudaMalloc((void **) &(g->d.log_pq_nospike), nsf) 							!= cudaSuccess ||
			cudaMalloc((void **) &(g->d.nmean),   nsf) 									!= cudaSuccess || //moments
			cudaMalloc((void **) &(g->d.gpmean),   Tf) 									!= cudaSuccess || //moments
			cudaMalloc((void **) &(g->d.gpsqmean), Tf) 									!= cudaSuccess || //moments
			cudaMalloc((void **) &(g->d.O), npi)					   	  				!= cudaSuccess ||
			cudaMalloc((void **) &(g->d.Oi), npi)					     				!= cudaSuccess ||
			cudaMalloc((void **) &(g->d.ancestor), npi * (maxK + 1)) 					!= cudaSuccess
	) {
		free_gpupfdata(g);
		return -4;
	}

	//FIXME have an option for a default seed where the rngstates are precomputed and compiled into the program
	if (cudainitrng(g->options.nblocks_pf, g->d.rngstates, seedval) < 0) { //always initialize with the same seed
		free_gpupfdata(g);
		return -5;
	}

	if (InitializeCubTemporaryArrays(g) < 0) {
		free_gpupfdata(g);
		return -6;
	}

	return 1;
}

void destroy_gpu5s_problem(gpu5s_problem *g) {
	free_gpupfdata(g);

	if (g->h.fobs != NULL)
		free(g->h.fobs);
	if (g->h.u != NULL)
		free(g->h.u);

	g->h.fobs = NULL; g->h.u = NULL;
	free(g);
}

//FIXME template these freeifnotnull functions
void cufreeifnotnull(float *&p) {
	if (p != NULL) {
		cudaFree(p);
		p = NULL;
	}
}

void cufreeifnotnull(void *&p) {
	if (p != NULL) {
		cudaFree(p);
		p = NULL;
	}
}

void cufreeifnotnull(bool *&p) {
	if (p != NULL) {
		cudaFree(p);
		p = NULL;
	}
}

void cufreeifnotnull(int *&p) {
	if (p != NULL) {
		cudaFree(p);
		p = NULL;
	}
}

void free_gpupfdata(gpu5s_problem *g) {
	if (g == NULL)
		return;

	cufreeifnotnull(g->d.rngstates);
	cufreeifnotnull(g->d.cbr);
	cufreeifnotnull(g->d.states);
	cufreeifnotnull(g->d.gp);
	cufreeifnotnull(g->d.ns);
	cufreeifnotnull(g->d.ns_pre);
	cufreeifnotnull(g->d.log_w);
	cufreeifnotnull(g->d.w);
	cufreeifnotnull(g->d.log_w_corrected);
	cufreeifnotnull(g->d.w_corrected);
	cufreeifnotnull(g->d.W);
	cufreeifnotnull(g->d.block_lw_max);
	cufreeifnotnull(g->d.block_sw_over_max);
	cufreeifnotnull(g->d.block_swsq_over_maxsq);
	cufreeifnotnull(g->d.block_lw_c_max);
	cufreeifnotnull(g->d.block_sw_c_over_max);
	cufreeifnotnull(g->d.block_sw_spike);
	cufreeifnotnull(g->d.block_sw_spike_pre);
	cufreeifnotnull(g->d.block_wgp);
	cufreeifnotnull(g->d.block_wgpsq);
	cufreeifnotnull(g->d.log_sum_raw_w);
	cufreeifnotnull(g->d.neff);
	cufreeifnotnull(g->d.q_spike);
	cufreeifnotnull(g->d.log_pq_spike);
	cufreeifnotnull(g->d.log_pq_nospike);
	cufreeifnotnull(g->d.nmean);
	cufreeifnotnull(g->d.gpmean);
	cufreeifnotnull(g->d.gpsqmean);
	cufreeifnotnull(g->d.O);
	cufreeifnotnull(g->d.Oi);
	cufreeifnotnull(g->d.ancestor);

	FreeCubTemporaryArrays(g); //calls cufreeifnotnull
}

int file2gpu(void * darray, FILE *pFile, int nbytes) {
	void * harray = malloc(nbytes);
	if (harray == NULL) //failed to allocate
			return 0;
	int nread = fread(harray, 1, nbytes, pFile);
	if (nread != nbytes) {
		free(harray);
		return -nread;
	}
	cudaMemcpy(darray, harray, nbytes, cudaMemcpyHostToDevice); //fixme: check for success here
	free(harray);
	return nread;
}

int save_gpu5sresults_toraw(FILE *pFile, gpu5s_problem * g) {
	//fixme save moments too
	float *ll = (float *) malloc(sizeof(float) * g->T);
	if (ll == NULL) {
		printf("Failed allocate memory!\n");
		return -2;
	}
	float *neff = (float *) malloc(sizeof(float) * g->T);
	if (neff == NULL) {
		printf("Failed allocate memory!\n");
		free(ll);
		return -3;
	}
	cudaMemcpy(ll,   g->d.log_sum_raw_w, g->T * sizeof(float), cudaMemcpyDeviceToHost); //fixme check for success
	cudaMemcpy(neff, g->d.neff,          g->T * sizeof(float), cudaMemcpyDeviceToHost); //fixme check for success
	if (fwrite(ll, sizeof(float),   g->T, pFile) != g->T) {
		printf("Failed to write sum of log weights!\n");
	}
	if (fwrite(neff, sizeof(float), g->T, pFile) != g->T) {
		printf("Failed to write neff!\n");
	}
	free(ll);
	free(neff);
	return 0;
}

gpu5s_problem * init_gpu5sproblem_fromraw(FILE *pFile) {
	int T, np, nsteps, K, ntimepoints_pre;
	if (fread((void *) &T, sizeof(int), 1, pFile) != 1) {
		return NULL;
	}
	if (fread((void *) &np, sizeof(int), 1, pFile) != 1) {
		return NULL;
	}
	if (fread((void *) &nsteps, sizeof(int), 1, pFile) != 1) {
		return NULL;
	}
	if (fread((void *) &K, sizeof(int), 1, pFile) != 1) {
		return NULL;
	}
	if (fread((void *) &ntimepoints_pre, sizeof(int), 1, pFile) != 1) {
		return NULL;
	}
	gpu5s_problem *g = (gpu5s_problem *) malloc(sizeof(gpu5s_problem));
	g->h.fobs = NULL;
	g->h.u = NULL;
	int maxKp1 = 1;
	while(maxKp1 < K + 1)
		maxKp1 *= 2; //set maxK so that maxK+1 is a power of 2 and maxK is sufficiently large
	int maxK = maxKp1 - 1;
	if (allocate_gpupfdata(g, T, np, nsteps, maxK, (K + 1) * nsteps, (T + ntimepoints_pre) * nsteps, ntimepoints_pre * nsteps) < 0) {
		std::cerr << "Allocation failed!\n";
		return NULL;
	}
	g->options.K = K;
	g->options.nsteps = nsteps;
	g->options.ntimepoints_pre = ntimepoints_pre;
	g->T = T;

	g->options.ntimepoints_pre = ntimepoints_pre;
	g->options.computenmean = true;
	if (fread((void *) &(g->dt), sizeof(float), 1, pFile) != 1) {
		destroy_gpu5s_problem(g);
		return NULL;
	}
	if (fread((void *) &(g->options.resamplethreshold), sizeof(float), 1, pFile) != 1) {
		destroy_gpu5s_problem(g);
		return NULL;
	}
	if (file2gpu(g->d.neff, pFile, sizeof(float)) != sizeof(float)) { //neff0
		destroy_gpu5s_problem(g);
		return NULL;
	}
	if (file2gpu(g->d.log_sum_raw_w, pFile, sizeof(float)) != sizeof(float)) { //log_sum_raw_w0
		destroy_gpu5s_problem(g);
		return NULL;
	}
	g->h.fobs   = (float *) malloc(sizeof(float) * g->T);
	g->h.u      = (float *) malloc(sizeof(float) * g->T);
	if (g->h.fobs == NULL || g->h.u == NULL) {
		destroy_gpu5s_problem(g);
		return NULL;
	}
	if (fread((void *) (g->h.fobs), sizeof(float), g->T, pFile) != g->T) {
		destroy_gpu5s_problem(g);
		return NULL;
	}
	if (fread((void *) (g->h.u), sizeof(float), g->T, pFile) != g->T) {
		destroy_gpu5s_problem(g);
		return NULL;
	}
	int nbytes_q = sizeof(float) * g->options.nsteps * (g->T + g->options.ntimepoints_pre);
	if (file2gpu(g->d.q_spike, pFile, nbytes_q) != nbytes_q) {
		destroy_gpu5s_problem(g);
		return NULL;
	}
	if (file2gpu(g->d.log_pq_spike, pFile, nbytes_q) != nbytes_q) {
		destroy_gpu5s_problem(g);
		return NULL;
	}
	if (file2gpu(g->d.log_pq_nospike, pFile, nbytes_q) != nbytes_q) {
		destroy_gpu5s_problem(g);
		return NULL;
	}

	if (fread((void *) &(g->params.db1), sizeof(float), 1, pFile) != 1) { destroy_gpu5s_problem(g); return NULL; }
	if (fread((void *) &(g->params.db2), sizeof(float), 1, pFile) != 1) { destroy_gpu5s_problem(g); return NULL; }
	if (fread((void *) &(g->params.db3), sizeof(float), 1, pFile) != 1) { destroy_gpu5s_problem(g); return NULL; }
	if (fread((void *) &(g->params.db4), sizeof(float), 1, pFile) != 1) { destroy_gpu5s_problem(g); return NULL; }
	if (fread((void *) &(g->params.vF), sizeof(float), 1, pFile) != 1) { destroy_gpu5s_problem(g); return NULL; }
	if (fread((void *) &(g->params.sigma_r), sizeof(float), 1, pFile) != 1) { destroy_gpu5s_problem(g); return NULL; }
	if (fread((void *) &(g->params.lambda), sizeof(float), 1, pFile) != 1) { destroy_gpu5s_problem(g); return NULL; }
	if (fread((void *) &(g->params.S), sizeof(float), 1, pFile) != 1) { destroy_gpu5s_problem(g); return NULL; }
	if (fread((void *) &(g->params.FBGp1), sizeof(float), 1, pFile) != 1) { destroy_gpu5s_problem(g); return NULL; }
	if (fread((void *) &(g->params.gain), sizeof(float), 1, pFile) != 1) { destroy_gpu5s_problem(g); return NULL; }
	if (fread((void *) &(g->params.kd_ex), sizeof(float), 1, pFile) != 1) { destroy_gpu5s_problem(g); return NULL; }
	if (fread((void *) &(g->params.maxex), sizeof(float), 1, pFile) != 1) { destroy_gpu5s_problem(g); return NULL; }
	if (fread((void *) &(g->params.c0), sizeof(float), 1, pFile) != 1) { destroy_gpu5s_problem(g); return NULL; }
	if (fread((void *) &(g->params.fdc), sizeof(float), 1, pFile) != 1) { destroy_gpu5s_problem(g); return NULL; }
	if (fread((void *) &(g->params.Btot0), sizeof(float), 1, pFile) != 1) { destroy_gpu5s_problem(g); return NULL; }
	if (fread((void *) &(g->params.Btot1), sizeof(float), 1, pFile) != 1) { destroy_gpu5s_problem(g); return NULL; }
	if (fread((void *) &(g->params.kon0), sizeof(float), 1, pFile) != 1) { destroy_gpu5s_problem(g); return NULL; }
	if (fread((void *) &(g->params.kon1), sizeof(float), 1, pFile) != 1) { destroy_gpu5s_problem(g); return NULL; }
	if (fread((void *) &(g->params.kon2), sizeof(float), 1, pFile) != 1) { destroy_gpu5s_problem(g); return NULL; }
	if (fread((void *) &(g->params.kon3), sizeof(float), 1, pFile) != 1) { destroy_gpu5s_problem(g); return NULL; }
	if (fread((void *) &(g->params.koff0), sizeof(float), 1, pFile) != 1) { destroy_gpu5s_problem(g); return NULL; }
	if (fread((void *) &(g->params.koff1), sizeof(float), 1, pFile) != 1) { destroy_gpu5s_problem(g); return NULL; }
	if (fread((void *) &(g->params.koff2), sizeof(float), 1, pFile) != 1) { destroy_gpu5s_problem(g); return NULL; }
	if (fread((void *) &(g->params.koff3), sizeof(float), 1, pFile) != 1) { destroy_gpu5s_problem(g); return NULL; }
	if (fread((void *) &(g->params.kon_B0), sizeof(float), 1, pFile) != 1) { destroy_gpu5s_problem(g); return NULL; }
	if (fread((void *) &(g->params.koff_B0), sizeof(float), 1, pFile) != 1) { destroy_gpu5s_problem(g); return NULL; }
	if (fread((void *) &(g->params.kon_B1), sizeof(float), 1, pFile) != 1) { destroy_gpu5s_problem(g); return NULL; }
	if (fread((void *) &(g->params.koff_B1), sizeof(float), 1, pFile) != 1) { destroy_gpu5s_problem(g); return NULL; }

	int state_nbytes = sizeof(float4) * g->options.nparticles;
	if (file2gpu(g->d.states, pFile, state_nbytes) != state_nbytes) {
		destroy_gpu5s_problem(g);
		return NULL;
	}
	if (file2gpu(g->d.cbr, pFile, state_nbytes) != state_nbytes) {
		destroy_gpu5s_problem(g);
		return NULL;
	}

	int npf_nbytes = sizeof(float) * g->options.nparticles;
	if (file2gpu(g->d.w, pFile, npf_nbytes) != npf_nbytes) {
		destroy_gpu5s_problem(g);
		return NULL;
	}
	if (file2gpu(g->d.log_w, pFile, npf_nbytes) != npf_nbytes) {
		destroy_gpu5s_problem(g);
		return NULL;
	}
	if (file2gpu(g->d.log_w_corrected, pFile, npf_nbytes) != npf_nbytes) {
		destroy_gpu5s_problem(g);
		return NULL;
	}
	int ns_pre_nbytes = g->options.nparticles * g->options.nsteps * g->options.ntimepoints_pre;
	if (file2gpu(g->d.ns_pre, pFile, ns_pre_nbytes) != ns_pre_nbytes) {
		destroy_gpu5s_problem(g);
		return NULL;
	}
	if (file2gpu(g->d.gp, pFile, npf_nbytes) != npf_nbytes) {
		destroy_gpu5s_problem(g);
		return NULL;
	}
	return g;
}
