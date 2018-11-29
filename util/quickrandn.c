/* all ideas and most code from Marsaglia & Tsang 2002 */
#include <math.h>
#include <mex.h>
#include <stdint.h>

float nfix(void);
void zigset();

#define SHR3 (jz=jsr, jsr^=(jsr<<13), jsr^=(jsr>>17), jsr^=(jsr<<5),jz+jsr)
#define UNI (.5 + (signed) SHR3*.2328306e-9)
#define IUNI SHR3

uint32_t jz, jsr;
static int32_t hz;
static uint32_t iz, kn[128];
static float wn[128],fn[128];

#define RNOR (hz=SHR3, iz=hz&127, (fabs(hz)<kn[iz])? hz*wn[iz] : nfix())

void mexFunction(
        int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    
    uint32_t n, j;
    float *w;
    
    if (nrhs != 2 || nlhs != 1) {
        mexErrMsgTxt("must have two inputs and one output");
    }
    if (mxGetClassID(prhs[0]) != mxSINGLE_CLASS) {
        mexErrMsgTxt("first input must be single class");
    }
    if ((mxGetNumberOfElements(prhs[1]) != 1) | (mxGetClassID(prhs[1]) != mxUINT32_CLASS)) {
        mexErrMsgTxt("second input must be a uint32 scalar");
    }
    plhs[0] = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
    n = (uint32_t) mxGetNumberOfElements(prhs[0]);
    jsr = *((uint32_t *) mxGetData(prhs[1]));
    
    zigset();
    w = (float *) mxGetData(prhs[0]);
    for (j = 0; j < n; j++) {
        w[j] = RNOR;
    }
    *((uint32_t *) mxGetData(plhs[0])) = jsr; /* set the second input as the new seed */    
}

float nfix(void)
{
    const float r = 3.442620f;     /* The start of the right tail */
    static float x, y;
    for(;;)
    {  x=hz*wn[iz];      /* iz==0, handles the base strip */
       if(iz==0)
       { do{ x=-log(UNI)*0.2904764; y=-log(UNI);}	/* .2904764 is 1/r */
         while(y+y<x*x);
         return (hz>0)? r+x : -r-x;
       }
       /* iz>0, handle the wedges of other strips */
       if( fn[iz]+UNI*(fn[iz-1]-fn[iz]) < exp(-.5*x*x) ) return x;
       
       /* initiate, try to exit for(;;) for loop*/
       hz=SHR3;
       iz=hz&127;
       if(fabs(hz)<kn[iz]) return (hz*wn[iz]);
    }
}

void zigset()
{  const double m1 = 2147483648.0;
   double dn=3.442619855899,tn=dn,vn=9.91256303526217e-3, q;
   int32_t i;
   
   /* Set up tables for RNOR */
   q=vn/exp(-.5*dn*dn);
   kn[0]=(dn/q)*m1;
   kn[1]=0;
   
   wn[0]=q/m1;
   wn[127]=dn/m1;
   
   fn[0]=1.;
   fn[127]=exp(-.5*dn*dn);
   
   for(i=126;i>=1;i--)
   {dn=sqrt(-2.*log(vn/dn+exp(-.5*dn*dn)));
    kn[i+1]=(dn/tn)*m1;
    tn=dn;
    fn[i]=exp(-.5*dn*dn);
    wn[i]=dn/m1;
   }
}