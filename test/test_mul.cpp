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

  for (long i = 1; i < 770; i+=1){

    long j = 0;
    while ((1 << j) < (2*i+1)){
      CTFT_init_multipliers(mult, j);
      j++;
    }

    zz_pX A, B, C1, C2;
    for (long j = 0; j < i; j++){
      SetCoeff(A, j, random_zz_p());
      SetCoeff(B, j, random_zz_p());
    }
    for (long j = 0; j < 2*i+1; j++){
      SetCoeff(C1, j, random_zz_p());
      SetCoeff(C2, j, random_zz_p());
    }
    
    if (opt == 1){
      CTFT_mul(C1, A, B, i, mult);
      FFTMul(C2, A, B);
      assert (C1 == C2);
      cout << i << endl;
    }
    else{
      cout << i;

      double t = GetTime();
      for (long j = 0; j < 10000; j++)
	CTFT_mul(C1, A, B, i, mult);
      t = GetTime()-t;
      cout << " " << t;

      FFTMul(C2, A, B);
      double v = GetTime();
      for (long j = 0; j < 10000; j++)
	FFTMul(C2, A, B);
      v = GetTime()-v;
      cout << " " << v;

      cout << endl;
    }
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

