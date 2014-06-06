#include <NTL/lzz_p.h>
#include <NTL/FFT.h>
#include "CTFT.h"


NTL_CLIENT

/*---------------------------------------------------------*/
/* applies the premultiplier before negacyclic convolution */
/* input length = 2^k                                      */
/* input is word-size, output is in [0,p)                  */
/*---------------------------------------------------------*/
void CTFT_rescale(long *y, const long* x, const long k, const CTFT_multipliers& mult){
  long K = 1L << k;
  
  const long p = mult.p;
  const long * powers = mult.z[k].elts();
  const unsigned long * powers_precomp = mult.z_precomp[k].elts();

  for (long i = 0; i < K; i++){
    y[i] = LazyReduce(LazyMulModPrecon(x[i], powers[i], p, powers_precomp[i]), p);
  }
}

/*---------------------------------------------------------*/
/* applies postmultiplier after negacyclic convolution     */
/* input length = 2^k                                      */
/* input is word-size, output is in [0,p)                  */
/* incorporates the division by 2^k for FFT^(-1)           */
/*---------------------------------------------------------*/
void CTFT_rescale_inverse(long *y, const long* x, const long k, const CTFT_multipliers& mult){
  long K = 1L << k;
  
  const long p = mult.p;
  const long * powers = mult.invz[k].elts();
  const unsigned long * powers_precomp = mult.invz_precomp[k].elts();

  for (long i = 0; i < K; i++)
    y[i] = LazyReduce(LazyMulModPrecon(x[i], powers[i], p, powers_precomp[i]), p);
}

