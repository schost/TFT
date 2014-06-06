#include <NTL/lzz_p.h>
#include "CTFT.h"

NTL_CLIENT

/*-----------------------------------------------------------*/
/* applies the multi-reduction map to x (mod p), in place    */
/* input is in [0,p), output is in [0,2p)                    */
/* x must have size at least 2N, and upper half must be zero */
/*-----------------------------------------------------------*/
void CTFT_reduce(long* x, const long N, const long p){
  if (N <= 1)
    return;

  long n = N;
  long k = 0;
  long a = 1L;
  while (a <= n/2){
    k++;
    a <<= 1;
  }

  // if N is a power of 2, nothing to do
  if (a == n)
    return;

  // else, loop through the reduction process
  while (a != n){

    long b = a >> 1;
    long ell = k-1;
    while (! (n & b)){
      b >>= 1;
      ell--;
    }

    long t = 0;
    const long b2 = b << 1;

    // b = 1: unroll the loop
    if (b == 1){
      long u, v;
      u = x[0];
      v = x[a+0];
      x[0] = u - v + p;
      x[a+0] = AddMod(u, v, p);
      u = x[1];
      v = x[a+1];
      x[1] = u - v + p;
      x[a+1] = AddMod(u, v, p);

      long t = 2;
      long lambda = 1L << (k-ell-1);
      for (long i = 1; i < lambda; i++){
	long u, v;
	u = x[t];
	v = x[a+t];
	x[t] = u - v + p;
	x[a] = AddMod(AddMod(u, v, p), x[a], p);
	t++;
	u = x[t];
	v = x[a+t];
	x[t] = u - v + p;
	x[a+1] = AddMod(AddMod(u, v, p), x[a+1], p);
	t++;
      }
    }
    else{
      // i = 0 to 2^(k-ell-1) = a/(2b), j = 0 to 2b
      // i = 0 : special case
      for (; t < b2; t++){
	long u, v;
	u = x[t];
	v = x[a+t];
	x[t] = u - v + p;
	x[a+t] = AddMod(u, v, p);
      }
      
      long lambda = 1L << (k-ell-1);
      for (long i = 1; i < lambda; i++){
	for (long j = 0; j < b2; j++){
	  long u, v;
	  u = x[t];
	  v = x[a+t];
	  x[t] = u - v + p;
	  x[a+j] = AddMod(AddMod(u, v, p), x[a+j], p);
	  t++;
	}
      }
    }

    x += a;
    n -= a;
    a = b;
    k = ell;
  }

  for (long t = 0; t < a; t++){
    long u, v;
    u = x[t];
    v = x[a+t];
    x[t] = u - v + p;
  }
}



