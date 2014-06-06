#include "CTFT.h"

/*-----------------------------------------*/
/* returns the exponents ki such that      */
/*    n = sum_i 2^k_i                      */
/*-----------------------------------------*/
long * CTFT_exponents(const long nn){
  long i = weight(nn);
  long * vec = new long[i+1];

  if (nn == 0){
    vec[0] = -1;
    return vec;
  }
  
  long j;
  for (j = 0; j < i; j++)
    vec[j] = 0;
  vec[i] = -1;

  long n = nn;
  long a = 1;
  i = 0;
  j = 0;
  while (a <= n/2){
    a <<= 1;
    j++;
  }
  vec[i++] = j;

  while (a != n){
    long k = j-1;
    long b = a >> 1;
    while (! (n & b)){
      k--;
      b >>= 1;
    }
    n -= a;
    a = b;
    vec[i++] = k;
    j = k;
  }

  return vec;
}
