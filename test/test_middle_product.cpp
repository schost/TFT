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

    zz_pX A, B1, B2, C;
    for (long j = 0; j < i; j++){
      SetCoeff(A, j, random_zz_p());
      SetCoeff(B1, j, random_zz_p());
      SetCoeff(B2, j, random_zz_p());
    }
    for (long j = 0; j < 2*i-1; j++)
      SetCoeff(C, j, random_zz_p());

    if (opt == 1){
      CTFT_middle_product(B1, A, C, i, mult);
      
      B2 = A*C;
      
      for (long j = 0; j < i; j++)
	assert (coeff(B1, j) == coeff(B2, j+i-1));

      cout << i << endl;
		
    }
    else{
      cout << i;

      double u = GetTime();
      for (long j = 0; j < 10000; j++)
     	CTFT_middle_product(B1, A, C, i, mult);
      u = GetTime()-u;
      cout << " " << u;

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

