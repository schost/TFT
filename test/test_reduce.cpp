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

    for (long j = 0; j < i; j++)
      src[j] = random_zz_p().LoopHole();
    for (long j = i; j < 2*i; j++)
      src[j] = 0;

    if (opt == 1){
      zz_pX a;
      for (long k = 0; k < i; k++)
	SetCoeff(a, k, src[k]);
      CTFT_reduce(src, i, zz_p::modulus());

      long *degrees = CTFT_exponents(i);

      long j = 0;
      long idx = 0;
      while (degrees[j] != -1){
	long ell = 1 << degrees[j];
	zz_pX b;
	SetCoeff(b, 0, 1);
	SetCoeff(b, ell, 1);
	zz_pX c = a % b;

	for (long k = 0; k < ell; k++)
	  assert (coeff(c, k) == (src[idx+k]) % zz_p::modulus());

	idx += ell;
	j++;
	
      }
      free(degrees);

      cout << i << endl;
    }
    else{
      cout << i;

      double t = GetTime();
      for (long j = 0; j < 1000000; j++)
      	CTFT_reduce(src, i, zz_p::modulus());
      t = GetTime()-t;
      cout << " " << t;

      cout << endl;
    }

    delete[] src;
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
