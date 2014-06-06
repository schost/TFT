#include <NTL/lzz_pX.h>
#include <NTL/FFT.h>
#include "CTFT.h"


NTL_CLIENT


/*---------------------------------------------------------*/
/* transpose negacyclic convolution in length 2^k w.r.t. b */
/* wk and wk2 are workspaces of length 2^k                 */
/* done by replacing e.g. [a0 a1 a2 a3] by [a0 -a3 -a2 -a1]*/
/*---------------------------------------------------------*/
void CTFT_negacyclic_convolution_t(long *b, const long *a, const long *c, const long k, const CTFT_multipliers& mult, long *wk, long *wk2){
  
  const long index = mult.index;
  const long q = mult.p;
  const double qinv = mult.pinv;

  long K = 1L << k;

  switch(k){
  case 0:
    b[0] = MulMod(a[0], c[0], q);
    return;
  case 1:
    b[0] = AddMod(MulMod(a[0], c[0], q), MulMod(a[1], c[1], q), q);
    b[1] = SubMod(MulMod(a[0], c[1], q), MulMod(a[1], c[0], q), q);
  default:

    const long * powers = mult.z[k].elts();
    const unsigned long * powers_precomp = mult.z_precomp[k].elts();

    wk[0] = LazyReduce(LazyMulModPrecon(a[0], powers[0], q, powers_precomp[0]), q);
    for (long i = 1; i < K; i++)
      wk[i] = LazyReduce(LazyMulModPrecon(NegateMod(LazyReduce(a[K-i], q), q), powers[i], q, powers_precomp[i]), q);
    FFTFwd(wk, wk, k, index);

    CTFT_rescale(wk2, c, k, mult);
    FFTFwd(wk2, wk2, k, index);
    
    for (long j = 0; j < K; j++)
      b[j] = MulMod(wk[j], wk2[j], q, qinv);
    
    FFTRev(b, b, k, index);
    CTFT_rescale_inverse(b, b, k, mult);
  }
}

