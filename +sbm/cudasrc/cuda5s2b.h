/*
 * cuda5s.h
 *
 *  Created on: Jul 24, 2014
 *      Author: greenberg
 */
#include <stdio.h>

# ifndef CUDA5S_H_
# define CUDA5S_H_


# ifndef M_PI
# define M_PI (3.14159265358979323846)
# endif
# ifndef MAX
# define MAX(a,b) (((a)>(b))?(a):(b))
# endif
# ifndef MIN
# define MIN(a,b) (((a)<(b))?(a):(b))
# endif

//block size, values processed per thread, etc.
# ifndef CUDA5S_PARALLELOPTS_H_
# include "cuda5s_parallelopts.h"
# endif

#define USEKDEX false

//define the struct here for files including cuda5s.h that will call cudaDeviceList_wrapper
/**
 * CUDA device properties
 */
#if !defined(__DRIVER_TYPES_H__)
struct cudaDeviceProp
{
    char   name[256];                  /**< ASCII string identifying device */
    size_t totalGlobalMem;             /**< Global memory available on device in bytes */
    size_t sharedMemPerBlock;          /**< Shared memory available per block in bytes */
    int    regsPerBlock;               /**< 32-bit registers available per block */
    int    warpSize;                   /**< Warp size in threads */
    size_t memPitch;                   /**< Maximum pitch in bytes allowed by memory copies */
    int    maxThreadsPerBlock;         /**< Maximum number of threads per block */
    int    maxThreadsDim[3];           /**< Maximum size of each dimension of a block */
    int    maxGridSize[3];             /**< Maximum size of each dimension of a grid */
    int    clockRate;                  /**< Clock frequency in kilohertz */
    size_t totalConstMem;              /**< Constant memory available on device in bytes */
    int    major;                      /**< Major compute capability */
    int    minor;                      /**< Minor compute capability */
    size_t textureAlignment;           /**< Alignment requirement for textures */
    size_t texturePitchAlignment;      /**< Pitch alignment requirement for texture references bound to pitched memory */
    int    deviceOverlap;              /**< Device can concurrently copy memory and execute a kernel. Deprecated. Use instead asyncEngineCount. */
    int    multiProcessorCount;        /**< Number of multiprocessors on device */
    int    kernelExecTimeoutEnabled;   /**< Specified whether there is a run time limit on kernels */
    int    integrated;                 /**< Device is integrated as opposed to discrete */
    int    canMapHostMemory;           /**< Device can map host memory with cudaHostAlloc/cudaHostGetDevicePointer */
    int    computeMode;                /**< Compute mode (See ::cudaComputeMode) */
    int    maxTexture1D;               /**< Maximum 1D texture size */
    int    maxTexture1DMipmap;         /**< Maximum 1D mipmapped texture size */
    int    maxTexture1DLinear;         /**< Maximum size for 1D textures bound to linear memory */
    int    maxTexture2D[2];            /**< Maximum 2D texture dimensions */
    int    maxTexture2DMipmap[2];      /**< Maximum 2D mipmapped texture dimensions */
    int    maxTexture2DLinear[3];      /**< Maximum dimensions (width, height, pitch) for 2D textures bound to pitched memory */
    int    maxTexture2DGather[2];      /**< Maximum 2D texture dimensions if texture gather operations have to be performed */
    int    maxTexture3D[3];            /**< Maximum 3D texture dimensions */
    int    maxTexture3DAlt[3];         /**< Maximum alternate 3D texture dimensions */
    int    maxTextureCubemap;          /**< Maximum Cubemap texture dimensions */
    int    maxTexture1DLayered[2];     /**< Maximum 1D layered texture dimensions */
    int    maxTexture2DLayered[3];     /**< Maximum 2D layered texture dimensions */
    int    maxTextureCubemapLayered[2];/**< Maximum Cubemap layered texture dimensions */
    int    maxSurface1D;               /**< Maximum 1D surface size */
    int    maxSurface2D[2];            /**< Maximum 2D surface dimensions */
    int    maxSurface3D[3];            /**< Maximum 3D surface dimensions */
    int    maxSurface1DLayered[2];     /**< Maximum 1D layered surface dimensions */
    int    maxSurface2DLayered[3];     /**< Maximum 2D layered surface dimensions */
    int    maxSurfaceCubemap;          /**< Maximum Cubemap surface dimensions */
    int    maxSurfaceCubemapLayered[2];/**< Maximum Cubemap layered surface dimensions */
    size_t surfaceAlignment;           /**< Alignment requirements for surfaces */
    int    concurrentKernels;          /**< Device can possibly execute multiple kernels concurrently */
    int    ECCEnabled;                 /**< Device has ECC support enabled */
    int    pciBusID;                   /**< PCI bus ID of the device */
    int    pciDeviceID;                /**< PCI device ID of the device */
    int    pciDomainID;                /**< PCI domain ID of the device */
    int    tccDriver;                  /**< 1 if device is a Tesla device using TCC driver, 0 otherwise */
    int    asyncEngineCount;           /**< Number of asynchronous engines */
    int    unifiedAddressing;          /**< Device shares a unified address space with the host */
    int    memoryClockRate;            /**< Peak memory clock frequency in kilohertz */
    int    memoryBusWidth;             /**< Global memory bus width in bits */
    int    l2CacheSize;                /**< Size of L2 cache in bytes */
    int    maxThreadsPerMultiProcessor;/**< Maximum resident threads per multiprocessor */
    int    streamPrioritiesSupported;  /**< Device supports stream priorities */
    int    globalL1CacheSupported;     /**< Device supports caching globals in L1 */
    int    localL1CacheSupported;      /**< Device supports caching locals in L1 */
    size_t sharedMemPerMultiprocessor; /**< Shared memory available per multiprocessor in bytes */
    int    regsPerMultiprocessor;      /**< 32-bit registers available per multiprocessor */
    int    managedMemory;              /**< Device supports allocating managed memory on this system */
    int    isMultiGpuBoard;            /**< Device is on a multi-GPU board */
    int    multiGpuBoardGroupID;       /**< Unique identifier for a group of devices on the same multi-GPU board */
};
#endif /* !__DRIVER_TYPES_H__ */


struct smcsi_pfparams;
struct gpu5s_problem;
struct gpu5s_devicedata;
struct gpu5s_hostdata;

int setgridsize(gpu5s_problem *g, int np);
void gpu5s_initialstates(gpu5s_problem *g);
int allocate_gpupfdata(gpu5s_problem *g, int maxT, int np, int maxnsteps, int maxK, int maxKp1substeps,
		int maxtotalsubsteps, int maxpresubsteps, unsigned long long seedval = 1006);
gpu5s_problem* init_gpu5sproblem_fromraw(FILE *pFile);
int save_gpu5sresults_toraw(FILE *pFile, gpu5s_problem * g);
void free_gpupfdata(gpu5s_problem *g); //only frees GPU arrays and sets the pointers to NULL
void destroy_gpu5s_problem(gpu5s_problem *g);
int cudainitrng(int nblocks, void *&rngstates, unsigned long long seedval);
void cudaseedrng(int nblocks, void *rngstates, int nrng, unsigned long long seedval = 0);
void free_gpupfdata(gpu5s_problem *g);
int runpfmainloop(gpu5s_problem *g);
void pushparameterstodevice(gpu5s_problem *g);
float gpu5s_marglik(gpu5s_problem *g);
void cudaMemcpy_h2d_wrapper(void *d, void *h, int nbytes);
void cudaMemcpy_d2h_wrapper(void *h, void *d, int nbytes);
void fillarray(float *x, float v, int np);
void cudaDeviceReset_wrapper();
int cudaDeviceList_wrapper(cudaDeviceProp *pDeviceList, int MaxDevices);

typedef struct smcsi_options {
	int nparticles, nsteps, ntimepoints_pre, nblocks, nt, vt, nblocks_pf, nt_pf, vt_pf, nblocks_rs, nt_rs, vt_rs, K, n_newton_iterations;
    float resamplethreshold;
    bool computenmean;
} smcsi_options;

typedef struct smcsi_pfparams {
    float db1, db2, db3, db4;
    float vF, sigma_r, lambda, S, FBGp1, gain, maxex, kd_ex, c0, fdc, Btot0, Btot1;
    float kon0, kon1, kon2, kon3;
    float koff0, koff1, koff2, koff3;
    float kon_B0, koff_B0, kon_B1, koff_B1;
} smcsi_pfparams;

typedef struct gpu5s_devicedata {
	int maxK, maxKp1substeps, maxtotalsubsteps, maxnsteps, maxpresubsteps, maxT, maxnparticles, maxnblocks;
    void *rngstates;

    //temporary storage arrays for CUB reductions, and their sizes
    void *d_temp_storage_float_max_reduction, *d_temp_storage_float_add_reduction, *d_temp_storage_float_add_scan, *d_temp_storage_uint_max_scan;
    size_t temp_storage_bytes_float_max_reduction, temp_storage_bytes_float_add_reduction, temp_storage_bytes_float_add_scan, temp_storage_bytes_uint_max_scan;

    void *states; // states 0-4 are stored
    void *cbr; //calcium, buffers and log(baseline)

    bool *ns, *ns_pre; //spiking. buffered nsteps * K times and nsteps * ntimepoints_pre times, respectively

    float *gp; //gain times photon flux. expected value of F - F_{DC}. buffered (K + 1) times

	float *log_w, *w, *w_corrected, *log_w_corrected, *W; //weights

	//arrays to store per-block reduction values
    float *block_lw_max, *block_sw_over_max, *block_swsq_over_maxsq, *block_lw_c_max, *block_sw_c_over_max, *block_sw_spike, *block_sw_spike_pre, *block_wgp, *block_wgpsq;

    float *log_sum_raw_w, *neff; //arrays with T elements, filled as we run the algorithm

    float *gpmean, *gpsqmean; //arrays with ntimepoints_pre + T elements

    float *q_spike, *log_pq_spike, *log_pq_nospike, *nmean; //arrays with (ntimepoints_pre + T) * nsteps elements

    int *O, *Oi; //cumulative offspring before and after max-int-scan to force it to be increasing

    int *ancestor; //nparticles * (maxK + 1) elements. important that it's buffered maxK and not K times so we can do fast modular arithmetic using & on the GPU (maxK + 1 is a power of 2)

} gpu5s_devicearraystruct;

typedef struct gpu5s_hostdata {
    float *fobs, *u; //arrays with T elements. u is for systematic resampling
} gpu5s_hostarraystruct;

typedef struct gpu5s_problem { //particle filtering problem
    smcsi_pfparams params;
	smcsi_options options;
    int T;
    float dt;
    gpu5s_hostdata h;
    gpu5s_devicedata d;
} smcsi_pfproblem;

#endif /* CUDA5s_H_ */
