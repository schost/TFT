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
    long *dest = new long[i];

    for (long j = 0; j < i; j++)
      src[j] = random_zz_p().LoopHole();
    for (long j = i; j < 2*i; j++)
      src[j] = 0;

    if (opt == 1){

      zz_pX res;
      clear(res);
      long *degrees = CTFT_exponents(i);
      long j = 0;
      long idx = 0;
      while (degrees[j] != -1){
      	long ell = 1 << degrees[j];
      	zz_pX b;
	for (long k = 0; k < ell; k++)
	  SetCoeff(b, k, src[idx+k]);
      	idx += ell;
	long pow = 0;
	while (pow < i){
	  res += LeftShift(b, pow);
	  b = -b;
	  pow += ell;
	}
	res = trunc(res, i);
      	j++;
      }

      CTFT_reduce_t(src, i, zz_p::modulus());
      
      for (long k = 0; k < i; k++)
	assert (src[k] == coeff(res, k));

      free(degrees);

      cout << i << endl;
    }
    else{
      cout << i;

      double t = GetTime();
      for (long j = 0; j < 100000; j++)
      	CTFT_reduce(src, i, zz_p::modulus());
      t = GetTime()-t;
      cout << " " << t;

      double u = GetTime();
      for (long j = 0; j < 100000; j++)
      	CTFT_reduce_t(src, i, zz_p::modulus());
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
