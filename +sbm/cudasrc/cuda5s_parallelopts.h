# ifndef CUDA5S_PARALLELOPTS_H_
# define CUDA5S_PARALLELOPTS_H_


#define WARPSIZE 32
#define BLOCKSIZE_PF 128
#define BLOCKSIZE 256
#define BLOCKSIZE_RS 256

#define NT BLOCKSIZE //threads per block (CTA)
#define VT 2 //values processed per thread (in most cases)
#define VT_PF 2
#define VT_RS 1
#define NV (NT * VT) //value processed per block
#define NV_RS (BLOCKSIZE_RS * VT_RS)
#define NV_PF (BLOCKSIZE_PF * VT_PF)

# endif /* CUDA5S_PARALLELOPTS_H_ */
