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

  for (long i = 1; i < 15; i++){

    long len = 1L << i;
    long *a = new long[len];
    long *b = new long[len];
    long *a2 = new long[len];
    long *b2 = new long[len];
    long *c = new long[len];
    long *c2 = new long[len];
    long *wk = new long[2*len];

    CTFT_init_multipliers(mult, i);

    for (long j = 0; j < len; j++){
      a[j] = random_zz_p().LoopHole();
      b[j] = random_zz_p().LoopHole();

      a2[j] = a[j];
      b2[j] = b[j];
    }

    if (opt == 1){

      zz_pX A, B, C, M;
      for (long j = 0; j < len; j++){
	SetCoeff(A, j, a[j]);
	SetCoeff(B, j, b[j]);
      }
      SetCoeff(M, 0, 1);
      SetCoeff(M, len, 1);

      CTFT_negacyclic_convolution(c, a, b, i, mult, wk, wk+len);

      C = (A*B) % M;

      for (long j = 0; j < len; j++)
	assert  (c[j] == coeff(C, j));

      cout << i << endl;
    }
    else{
      cout << i;
	    
      double t = GetTime();
      for (long j = 0; j < 100000; j++)
	CTFT_negacyclic_convolution(c, a, b, i, mult, wk, wk+len);
      t = GetTime()-t;
      cout << " " << t;

      zz_pX A, B, C;
      for (long j = 0; j < len/2; j++){
	SetCoeff(A, j, a[j]);
	SetCoeff(B, j, b[j]);
      }
      double v = GetTime();
      for (long j = 0; j < 100000; j++)
	FFTMul(C, A, B);
	
      v = GetTime()-v;
      cout << " " << v;

      cout << endl;
    }

    delete[] a;
    delete[] b;
    delete[] a2;
    delete[] b2;
    delete[] c;
    delete[] c2;
    delete[] wk;
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

