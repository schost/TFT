#ifndef __CTFT_H
#define __CTFT_H

#include <cstring>
#include <NTL/lzz_p.h>
#include <NTL/lzz_pX.h>

NTL_CLIENT

/*---------------------------------------------------------*/
/*---------------------------------------------------------*/
/* some routines for fast modular arithmetic from NTL      */
/* based on David Harvey's idea of lazy reductions         */
/*---------------------------------------------------------*/
/*---------------------------------------------------------*/

/*---------------------------------------------------------*/
/* prepares a multiplier-by-b mod n                        */
/*---------------------------------------------------------*/
static inline 
unsigned long LazyPrepMulModPrecon(long b, long n, double ninv)
{
   unsigned long q, r;

   q = (long) ( (((double) b) * NTL_SP_FBOUND) * ninv ); 
   r = (((unsigned long) b) << NTL_SP_NBITS ) - q * ((unsigned long) n);

   if (r >> (NTL_BITS_PER_LONG-1)) {
      q--;
      r += n;
   }
   else if (((long) r) >= n) {
      q++;
      r -=n;
   }

   unsigned long res = q << (NTL_BITS_PER_LONG - NTL_SP_NBITS);
   long qq;
   MulDivRem(qq, (long) r, 4, n, 4*ninv);
   res = res + (qq << (NTL_BITS_PER_LONG - NTL_SP_NBITS-2));
   return res;
}

/*---------------------------------------------------------*/
/* multplies a by b mod n, using the premultiplier         */
/* the output is in [0,2n)                                 */
/*---------------------------------------------------------*/
static inline 
unsigned long LazyMulModPrecon(unsigned long a, unsigned long b, 
                               unsigned long n, unsigned long bninv)
{
   unsigned long q = MulHiUL(a, bninv);
   unsigned long res = a*b - q*n;
   return res;
}

/*---------------------------------------------------------*/
/* reduction modulo q, assuming a is in [0,2q)             */
/*---------------------------------------------------------*/
static inline 
unsigned long LazyReduce(unsigned long a, unsigned long q)
{
  unsigned long res;
#if (NTL_ARITH_RIGHT_SHIFT && defined(NTL_AVOID_BRANCHING) && !defined(NTL_CLEAN_INT))
  res = a - q;
  res  += (((long) res) >> (NTL_BITS_PER_LONG-1)) & q; 
#elif (defined(NTL_AVOID_BRANCHING))
  res = a - q;
  res  += (-(res >> (NTL_BITS_PER_LONG-1))) & q; 
#else
  if (a >= q)
    res = a - q;
  else
    res = a;
#endif

  return res;
}

/*---------------------------------------------------------*/
/*---------------------------------------------------------*/
/* a class that contains precomputed information           */
/* for negacyclic convolution                              */
/* index is the index of the FFT prime, p is the prime,    */
/* pinv = 1/(double)p.                                     */
/* MaxK is the maximum K for which we know the multipliers */
/* z, z_precomp are multipliers for forward transform      */
/* invz, invz_precomp are for the inverse transform        */
/*---------------------------------------------------------*/
/*---------------------------------------------------------*/

class CTFT_multipliers {
 public:
  long index;
  long p;
  double pinv;
  long MaxK;
  Vec< Vec<long> > z;
  Vec< Vec<mulmod_precon_t> > z_precomp;
  Vec< Vec<long> > invz;
  Vec< Vec<mulmod_precon_t> > invz_precomp;
  
 CTFT_multipliers(long index_in) { 
   index = index_in;
   p = FFTPrime[index];
   pinv = FFTPrimeInv[index];
   MaxK = -1;
  }
};

/*--------------------------------------------------------------*/
/* initializes the k-th row of pre-multipliers for negacyclic   */
/* convolution (== in size 2^k)                                 */
/*--------------------------------------------------------------*/
void CTFT_init_multipliers(CTFT_multipliers& mult, const long k);

/*-----------------------------------------*/
/* returns the exponents ki such that      */
/*    n = sum_i 2^k_i                      */
/*-----------------------------------------*/
long * CTFT_exponents(const long nn);


/*-----------------------------------------------------------*/
/*-----------------------------------------------------------*/
/* direct algorithms                                         */
/*-----------------------------------------------------------*/
/*-----------------------------------------------------------*/

/*-----------------------------------------------------------*/
/* applies the multi-reduction map to x (mod p), in place    */
/* input is in [0,p), output is in [0,2p)                    */
/* x must have size at least 2N, and upper half must be zero */
/*-----------------------------------------------------------*/
void CTFT_reduce(long* f, const long N, const long p);

/*-----------------------------------------------------------*/
/* applies the CRT map to x (mod p), in place                */
/* x has length N                                            */
/* input is in [0,p), output is in [0,p)                     */
/* tmp is a temporary workspace, of length at least 2N       */
/*-----------------------------------------------------------*/
void CTFT_CRT(long* x, const long N, const long p, const double pinv, long *tmp);

/*-----------------------------------------------------------*/
/* applies the premultiplier before negacyclic convolution   */
/* input length = 2^k                                        */
/* input is word-size, output is in [0,p)                    */
/*-----------------------------------------------------------*/
void CTFT_rescale(long *y, const long* x, const long k, const CTFT_multipliers& mult);

/*-----------------------------------------------------------*/
/* applies postmultiplier after negacyclic convolution       */
/* input length = 2^k                                        */
/* input is word-size, output is in [0,p)                    */
/* incorporates the division by 2^k for FFT^(-1)             */
/*-----------------------------------------------------------*/
void CTFT_rescale_inverse(long *y, const long* x, const long k, const CTFT_multipliers& mult);

/*-----------------------------------------------------------*/
/* negacyclic convolution in length 2^k                      */
/* wk and wk2 are workspaces of length 2^k                   */
/*-----------------------------------------------------------*/
void CTFT_negacyclic_convolution(long *c, const long *a, const long *b, const long k, const CTFT_multipliers& mult, long *wk, long *wk2);

/*-----------------------------------------------------------*/
/* multiple negacyclic convolutions in length N              */
/* wk and wk2 are workspaces of length N                     */
/*-----------------------------------------------------------*/
void CTFT_mul_remainders(long *z, const long *x, const long *y, const long N, const CTFT_multipliers& mult, long *wk, long *wk2);

/*-----------------------------------------------------------*/
/* multiplies a and b                                        */
/* both inputs have length n, output has length 2n+1         */
/*-----------------------------------------------------------*/
void CTFT_mul(zz_pX& C, const zz_pX& A, const zz_pX& B, const long n, const CTFT_multipliers& mult);

/*-----------------------------------------------------------*/
/*-----------------------------------------------------------*/
/* tranposed algorithms                                      */
/*-----------------------------------------------------------*/
/*-----------------------------------------------------------*/

/*-----------------------------------------------------------*/
/* transpose multi-reduction (= extension of sequences)      */
/* in place                                                  */
/* input has length N, but needs space 2N                    */
/*-----------------------------------------------------------*/
void CTFT_reduce_t(long* x, const long N, const long p);

/*-----------------------------------------------------------*/
/* transpose CRT (= decomposition of sequences)              */
/* in place                                                  */
/* input has length N, but needs workspace tmp of size N     */
/*-----------------------------------------------------------*/
void CTFT_CRT_t(long* x, const long N, const long p, const double pinv, long *tmp);

/*---------------------------------------------------------*/
/* transpose negacyclic convolution in length 2^k w.r.t.  b*/
/* done by replacing e.g. [a0 a1 a2 a3] by [a0 -a3 -a2 -a1]*/
/* wk and wk2 are workspaces of length 2^k                 */
/*---------------------------------------------------------*/
void CTFT_negacyclic_convolution_t(long *b, const long *a, const long *c, const long k, const CTFT_multipliers& mult, long *wk, long *wk2);

/*---------------------------------------------------------*/
/* transpose multiple negacyclic convolutions in length N  */
/* wk and wk2 are workspaces of length 2^k                 */
/*---------------------------------------------------------*/
void CTFT_mul_remainders_t(long *y, const long *x, const long *z, const long N, const CTFT_multipliers& mult, long *wk, long *wk2);

/*---------------------------------------------------------*/
/* middle product of A and C                               */
/* inputs have length <= n (A) and <= 2n-1 (C)             */
/* output has length <= n                                  */
/*---------------------------------------------------------*/
void CTFT_middle_product(zz_pX& B, const zz_pX& a, const zz_pX& c, const long n, const CTFT_multipliers& mult);

#endif
