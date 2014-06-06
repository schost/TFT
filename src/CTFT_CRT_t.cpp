#include <NTL/lzz_p.h>
#include "CTFT.h"

NTL_CLIENT

/*-----------------------------------------------------------*/
/* transpose CRT (= decomposition of sequences)              */
/* in place                                                  */
/* input has length N, but needs workspace tmp of size N     */
/*-----------------------------------------------------------*/
void CTFT_CRT_t(long* x, const long N, const long p, const double pinv, long *tmp){

  // p must be odd!
  // half = 1/2 mod p
  const long half = NegateMod(p >> 1, p);
  const mulmod_precon_t half_pinv = PrepMulModPrecon(half, p, pinv);

  if (N <= 1)
    return;

  long n = N;
  long a, b, b2, n2;
  a = 1L;
  n2 = (n >> 1);
  while (a <= n2)
    a <<= 1;

  while (n != a){
    n = n-a;
    b = 1L;
    n2 = (n >> 1);
    while (b <= n2)
      b <<= 1;
    b2 = b << 1;
    
    for (long i = 0; i < b; i++){
      long u = MulModPrecon(AddMod(x[i], x[a+i], p), half, p, half_pinv);
      x[a+i] = u;
      tmp[i] = u;
    }

    for (long i = b; i < n; i++)
      x[a+i] = MulModPrecon(AddMod(x[i], x[a+i], p), half, p, half_pinv);

    tmp += b2;
    x += a;
    a = b;
  }

  while (n != N){
    b = a;
    a <<= 1;
    b2 = a;
    while (!(a & N))
      a <<= 1;
    tmp -= b2;
    x -= a;
    
    for (long i = 0; i < b; i++)
      tmp[i+b] = SubMod(tmp[i], AddMod(x[a+i], x[a+i], p), p);

    long t = 0;
    while (t < a)
      for (long i = 0; i < b2; i++, t++)
    	x[t] = SubMod(x[t], tmp[i], p);

    n = n+a;

  }

}

