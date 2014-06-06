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

  for (long i = 0; i < 1000; i+=1){

    long *src = new long[2*i];
    long *tmp = new long[2*i];
    long *dest = new long[2*i];

    for (long j = 0; j < i; j++){
      src[j] = random_zz_p().LoopHole();
      dest[j] = src[j];
    }
    for (long j = i; j < 2*i; j++){
      src[j] = 0;
      dest[j] = 0;
    }

    CTFT_reduce_t(src, i, zz_p::modulus());

    if (opt == 1){
      CTFT_CRT_t(src, i, zz_p::modulus(), zz_p::ModulusInverse(), tmp);

      for (long k = 0; k < i; k++)
	assert (dest[k] ==  src[k]);
      cout << i << endl;
    }
    else{
      long p = zz_p::modulus();
      double pinv = zz_p::ModulusInverse();

      cout << i << " ";

      double v = GetTime();
      for (long j = 0; j < 100000; j++)
	CTFT_CRT(dest, i, p, pinv, tmp);
      v = GetTime()-v;
      cout << v << " ";

      double u = GetTime();
      for (long j = 0; j < 100000; j++)
	CTFT_CRT_t(dest, i, p, pinv, tmp);
      u = GetTime()-u;
      cout << u << " ";

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
