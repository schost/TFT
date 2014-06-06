#include <NTL/lzz_p.h>
#include "CTFT.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/* if opt = 1, runs a check                                   */
/* else, runs timings                                         */
/*------------------------------------------------------------*/
void check(int opt){
  long i;

  for (i = 0; i < 20; i+=1){
    cout << i << ":   ";
    long * degrees = CTFT_exponents(i);
    long j = 0;
    while (degrees[j] != -1){
      cout << (1 << degrees[j]) << " ";
      j++;
    }
    cout << endl;
    delete[] degrees;
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

      
