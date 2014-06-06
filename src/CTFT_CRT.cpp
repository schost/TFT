#include <NTL/lzz_p.h>
#include "CTFT.h"

NTL_CLIENT

/*---------------------------------------------------------*/
/* reduces s mod (X^sz-1), assuming s has length N         */
/* assumes sz divides N                                    */
/* all calculations are mod p                              */
/*---------------------------------------------------------*/
static inline 
void fold_minus(long *x, const long *s, const long sz, const long N, const long p){
  long i;

  for (i = 0; i < sz; i++)
    x[i] = s[i];

  if (sz == 2){
    while (i < N){
      x[0] = AddMod(x[0], s[i], p);
      i++;
      x[1] = AddMod(x[1], s[i], p);
      i++;
    }
  }
  else{
    while (i < N)
      for (long j = 0; j < sz; j++, i++)
	x[j] = AddMod(x[j], s[i], p);
  }
}

/*-----------------------------------------------------------*/
/* applies the CRT map to x (mod p), in place                */
/* x has length N                                            */
/* result is in y                                            */
/* input is in [0,p), output is in [0,p)                     */
/* tmp is a temporary workspace, of length at least 2N       */
/*-----------------------------------------------------------*/
void CTFT_CRT(long* x, const long N, const long p, const double pinv, long *tmp){

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
    
    fold_minus(tmp, x, b2, a, p);

    for(long i = 0; i < b; i++)
      x[a+i] = AddMod(AddMod(tmp[b+i], tmp[b+i], p), x[a+i], p);

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
    
    for (long i = 0; i < b; i++){
      long u = AddMod(tmp[i], tmp[b+i], p);
      u = SubMod(x[a+i], u, p);
      x[a+i] = MulModPrecon(u, half, p, half_pinv);
      x[i] = AddMod(x[i], x[a+i], p);
    }
      
    for (long i = b; i < n; i++){
      x[a+i] = MulModPrecon(x[a+i], half, p, half_pinv);
      x[i] = AddMod(x[i], x[a+i], p);
    }

    n = n+a;

  }
}
