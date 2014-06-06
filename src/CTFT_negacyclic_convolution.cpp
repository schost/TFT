#include <NTL/lzz_pX.h>
#include <NTL/FFT.h>
#include "CTFT.h"


NTL_CLIENT


/*---------------------------------------------------------*/
/* negacyclic convolution in length 2^k                    */
/* wk and wk2 are workspaces of length 2^k                 */
/*---------------------------------------------------------*/
void CTFT_negacyclic_convolution(long *c, const long *a, const long *b, const long k, const CTFT_multipliers& mult, long *wk, long *wk2){
  
  const long index = mult.index;
  const long q = mult.p;
  const double qinv = mult.pinv;

  long K = 1L << k;

  switch(k){
  case 0:
    c[0] = MulMod(a[0], b[0], q);
    return;
  case 1:
    c[0] = SubMod(MulMod(a[0], b[0], q), MulMod(a[1], b[1], q), q);
    c[1] = AddMod(MulMod(a[0], b[1], q), MulMod(a[1], b[0], q), q);
    return;
  default:
    CTFT_rescale(wk, a, k, mult);
    FFTFwd(wk, wk, k, index);
    CTFT_rescale(wk2, b, k, mult);
    FFTFwd(wk2, wk2, k, index);
    
    for (long j = 0; j < K; j++)
      c[j] = MulMod(wk[j], wk2[j], q, qinv);
    
    FFTRev(c, c, k, index);
    CTFT_rescale_inverse(c, c, k, mult);
  }
}

