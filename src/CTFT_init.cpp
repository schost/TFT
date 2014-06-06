#include <NTL/lzz_p.h>
#include <NTL/FFT.h>
#include "CTFT.h"


NTL_CLIENT

/*--------------------------------------------------------------*/
/* initializes the k-th row of pre-multipliers for negacyclic   */
/* convolution (== in size 2^k)                                 */
/*--------------------------------------------------------------*/
void CTFT_init_multipliers(CTFT_multipliers& mult, const long k){

   if (k <= mult.MaxK) 
     return;

   mult.z.SetLength(k+1);
   mult.z_precomp.SetLength(k+1);
   mult.invz.SetLength(k+1);
   mult.invz_precomp.SetLength(k+1);

   const long p = mult.p;
   const double pinv = mult.pinv;
   const long index = mult.index;
   const FFTPrimeInfo& info = FFTTables[index];
   const long * root = info.RootTable.elts();

   // root[i] is a root of order 2^i
   // so root[0] = 1, root[1] = -1
   //
   // root[i+1]^j, j = 0..2^i-1
   for (long i = mult.MaxK+1; i <= k; i++){
     mult.z[i].SetLength(1 << i);
     mult.z_precomp[i].SetLength(1 << i);
     mult.invz[i].SetLength(1 << i);
     mult.invz_precomp[i].SetLength(1 << i);

     zz_p tmpz = to_zz_p(root[i+1]);
     zz_p tmp_invz = 1/tmpz;

     mult.z[i][0] = 1;
     mult.z_precomp[i][0] = LazyPrepMulModPrecon(mult.z[i][0], p, pinv);
     mult.invz[i][0] = info.TwoInvTable[k];
     mult.invz_precomp[i][0] = LazyPrepMulModPrecon(mult.invz[i][0], p, pinv);

     for (long j = 1; j < (1L << i); j++){
       mult.z[i][j] = (to_zz_p(mult.z[i][j-1]) * tmpz).LoopHole();
       mult.z_precomp[i][j] = LazyPrepMulModPrecon(mult.z[i][j], p, pinv);
       mult.invz[i][j] = (to_zz_p(mult.invz[i][j-1]) * tmp_invz).LoopHole();
       mult.invz_precomp[i][j] = LazyPrepMulModPrecon(mult.invz[i][j], p, pinv);
     }
   }

   mult.MaxK = k;
}


