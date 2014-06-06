#include <NTL/lzz_p.h>
#include <NTL/lzz_pX.h>
#include <assert.h>

#include "CTFT.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/* if opt = 1, runs a check                                   */
/* else, runs timings                                         */
/*------------------------------------------------------------*/
void check(int opt){

  CTFT_multipliers mult = CTFT_multipliers(0);
  
  for (long i = 0; i < 15; i++){
    long len = 1L << i;
    CTFT_init_multipliers(mult, i);
    long *src = new long[len];
    long *dest = new long[len];

    for (long j = 0; j < len; j++)
      src[j] = random_zz_p().LoopHole();

    long p = zz_p::modulus();

    if (opt == 1){
      CTFT_rescale(dest, src, i, mult);
      for (long j = 0; j < len; j++)
	assert (dest[j] == (to_zz_p(src[j])*to_zz_p(mult.z[i][j])).LoopHole());
      CTFT_rescale_inverse(dest, dest, i, mult);

      long pow2 = power(to_zz_p(2), i).LoopHole();
      for (long j = 0; j < len; j++)
	assert (MulMod(pow2, dest[j], p) == src[j]);

      cout << i << endl;
    }
    else {
      cout << i;
	    
      double t = GetTime();
      for (long j = 0; j < 100000; j++)
	CTFT_rescale(dest, src, i, mult);
      t = GetTime()-t;
      cout << " " << t;

      double u = GetTime();
      for (long j = 0; j < 100000; j++)
	CTFT_rescale_inverse(dest, src, i, mult);
      u = GetTime()-u;
      cout << " " << u;

      cout << endl;
    }

    delete[] src;
    delete[] dest;
  }

}

/*------------------------------------------------------------*/
/* main just calls check()                                    */
/* if not argument is given, runs timings                     */
/* if the argument 1 is given, runs check                     */
/*------------------------------------------------------------*/
int main(int argc, char **argv){
  zz_p::FFTInit(0);

  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);

  return 0;
}

