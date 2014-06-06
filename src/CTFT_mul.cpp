#include <NTL/lzz_p.h>
#include "CTFT.h"

NTL_CLIENT

/*---------------------------------------------------------*/
/* multiple negacyclic convolutions in length N            */
/* wk and wk2 are workspaces of length 2^k                 */
/*---------------------------------------------------------*/
void CTFT_mul_remainders(long *z, const long *x, const long *y, const long N, const CTFT_multipliers& mult, long *wk, long *wk2){
  if (N == 0)
    return;

  if (N == 1){
    z[0] = MulMod(x[0], y[0], mult.p);
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
    CTFT_negacyclic_convolution(z, x, y, k, mult, wk, wk2);

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
/* multiplies a and b                                      */
/* both inputs have length <= n, output has length 2n-1    */
/*---------------------------------------------------------*/
void CTFT_mul(zz_pX& C, const zz_pX& A, const zz_pX& B, const long n, const CTFT_multipliers& mult){
  const long N = 2*n-1;
  const long N2 = 2*N;
  const long p = mult.p;
  const double pinv = mult.pinv;

  long *wk = new long[5*N];

  const zz_p *a = A.rep.elts();
  long n0 = A.rep.length();
  for (long i = 0; i < n0; i++)
    wk[i] = a[i]._zz_p__rep;
  memset(wk+n0, 0, (2*N-n0)*sizeof(long));
  CTFT_reduce(wk, N, p);

  const zz_p *b = B.rep.elts();
  long n1 = B.rep.length();
  long *wkN2 = wk + N2;
  for (long i = 0; i < n1; i++)
    wkN2[i] = b[i]._zz_p__rep;
  memset(wkN2+n1, 0, (N2-n1)*sizeof(long));
  CTFT_reduce(wkN2, N, p);

  long *wkN4 = wk + 4*N;
  CTFT_mul_remainders(wkN4, wk, wkN2, N, mult, wk+N, wk+3*N);
  CTFT_CRT(wkN4, N, p, pinv, wk);

  clear(C);
  C.rep.SetLength(N);
  zz_p *c = C.rep.elts();
  for (long i = 0; i < N; i++)
    c[i] = zz_p(wkN4[i], INIT_LOOP_HOLE);
  C.normalize();

  delete[] wk;
}
