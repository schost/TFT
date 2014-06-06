#include <NTL/lzz_p.h>
#include "CTFT.h"


/*-----------------------------------------------------------*/
/* transpose multi-reduction (= extension of sequences)      */
/* in place                                                  */
/* input has length N, but needs space 2N                    */
/*-----------------------------------------------------------*/
void CTFT_reduce_t(long* x, const long N, const long p){
  if (N <= 1)
    return;

  long a = 1;
  long k = 0;

  while (! (a & N)){
    k++;
    a <<= 1;
  }

  if (a == N)
    return;

  for (long i = 0; i < a; i++)
    x[N+i] = NegateMod(x[N-a+i], p);

  long n = a;
  while (n != N){
    long b = a;
    long ell = k;
    a <<= 1;
    k++;
    while (!(a & N)){
      a <<= 1;
      k++;
    }

    const long b2 = b << 1;
    long lambda = 1L << (k-ell-1);

    n = n+a;
    
    long t = b2;
    // we have to deal with i=0 separately
    // (since otherwise we would erase the source)
    //
    // we also write special cases for b=1 (and should do b=2)
    if (b == 1){
      for (long i = 1; i < lambda; i++){
	long u, v;
	u = x[N-n+t];
	v = x[N-n+a];
	x[N-n+t] = AddMod(v, u, p);
	x[N-n+a+t] = SubMod(v, u, p);
	t++;
	u = x[N-n+t];
	v = x[N-n+a+1];
	x[N-n+t] = AddMod(v, u, p);
	x[N-n+a+t] = SubMod(v, u, p);
	t++;
      }
      
      long u, v;
      u = x[N-n];
      v = x[N-n+a];
      x[N-n] = AddMod(v, u, p);
      x[N-n+a] = SubMod(v, u, p);
      u = x[N-n+1];
      v = x[N-n+a+1];
      x[N-n+1] = AddMod(v, u, p);
      x[N-n+a+1] = SubMod(v, u, p);

    }
    else{
      for (long i = 1; i < lambda; i++)
	for (long j = 0; j < b2; j++){
	  long u = x[N-n+t];
	  long v = x[N-n+a+j];
	  x[N-n+t] = AddMod(v, u, p);
	  x[N-n+a+t] = SubMod(v, u, p);
	  t++;
	}
      
      for (long j = 0; j < b2; j++){
	long u = x[N-n+j];
	long v = x[N-n+a+j];
	x[N-n+j] = AddMod(v, u, p);
	x[N-n+a+j] = SubMod(v, u, p);
      }
    }
  }

}
