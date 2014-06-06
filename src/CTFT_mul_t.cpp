#include <NTL/lzz_p.h>
#include "CTFT.h"


NTL_CLIENT

/*---------------------------------------------------------*/
/* transpose multiple negacyclic convolutions in length N  */
/* wk and wk2 are workspaces of length 2^k                 */
/*---------------------------------------------------------*/
void CTFT_mul_remainders_t(long *y, const long *x, const long *z, const long N, const CTFT_multipliers& mult, long *wk, long *wk2){

  if (N == 0)
    return;

  if (N == 1){
    y[0] = MulMod(x[0], z[0], mult.p);
    return;
  }
  
  long n = N;
  long k = 0;
  long a = 1L;

  while (a <= n/2){
    k++;
    a <<= 1;
  }

  do {
    CTFT_negacyclic_convolution_t(y, x, z, k, mult, wk, wk2);

    z += a;
    x += a;
    y += a;

    n = n-a;

    a = 1L;
    k = 0;
    while (a <= n/2){
      k++;
      a <<= 1;
    }
  }
  while (n != 0);
}


/*---------------------------------------------------------*/
/* middle product of A and C                               */
/* inputs have length <= n (A) and <= 2n-1 (C)             */
/* output has length <= n                                  */
/*---------------------------------------------------------*/
void CTFT_middle_product(zz_pX& B, const zz_pX& A, const zz_pX& C, const long n, const CTFT_multipliers& mult){
  const long N = 2*n-1;
  const long N2 = 2*N;
  const long p = mult.p;
  const double pinv = mult.pinv;

  long *wk = new long[5*N];

  const zz_p *a = A.rep.elts();
  long n0 = A.rep.length();
  memset(wk, 0, (n-n0)*sizeof(long));
  for (long i = 0; i < n0; i++)
    wk[n-1-i] = a[i]._zz_p__rep;
  memset(wk+n0, 0, (N2-n0)*sizeof(long));
  CTFT_reduce(wk, N, p);

  long *wkN = wk+N;
  const zz_p *c = C.rep.elts();
  long n1 = C.rep.length();
  for (long i = 0; i < n1; i++)
    wkN[i] = c[i]._zz_p__rep;
  memset(wkN+n1, 0, (N2-n1)*sizeof(long));
  CTFT_CRT_t(wkN, N, p, pinv, wk+2*N);

  long *wkN2 = wk+N2;

  CTFT_mul_remainders_t(wkN2, wk, wkN, N, mult, wk+3*N, wk+4*N);
  CTFT_reduce_t(wkN2, N, p);

  clear(B);
  B.rep.SetLength(n);
  zz_p *b = B.rep.elts();
  for (long i = 0; i < n; i++)
    b[i] = zz_p(wkN2[i], INIT_LOOP_HOLE);

  B.normalize();

  delete[] wk;
}
