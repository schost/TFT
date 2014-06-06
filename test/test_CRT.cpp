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

  for (long i = 1; i < 1000; i+=1){

    long *src = new long[2*i];
    long *tmp = new long[2*i];
    long *dest = new long[2*i];

    for (long j = 0; j < i; j++){
      src[j] = random_zz_p().LoopHole();
      dest[j] = src[j];
    }
    for (long j = i; j < 2*i; j++)
      src[j] = 0;

    if (opt == 1){
      CTFT_reduce(src, i, zz_p::modulus());
      for (long j = 0; j < i; j++)   // warning, src is in [0,2p)
	src[j] = src[j] % zz_p::modulus();

      CTFT_CRT(src, i, zz_p::modulus(), zz_p::ModulusInverse(), tmp);

      for (long j = 0; j < i; j++)
	assert (dest[j] == src[j]);

      cout << i << endl;
    }
    else{
      cout << i << " ";

      double t;
      t = GetTime();
      for (long j = 0; j < 100000; j++)
	CTFT_CRT(dest, i, zz_p::modulus(), zz_p::ModulusInverse(), tmp);
      t = GetTime()-t;
      cout << t << " ";

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
