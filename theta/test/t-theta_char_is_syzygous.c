#include <stdlib.h>

#include "theta.h"


int main()
{
  flint_printf("theta_char_is_syzygous....");
  fflush(stdout);

  {
    slong g = 2;
    ulong ch1, ch2, ch3;
    slong n = n_pow(2, 2*g);
    slong count = 0;
    
    for (ch1 = 0; ch1 < n; ch1++) {
      for (ch2 = ch1 + 1; ch2 < n; ch2++) {
	for (ch3 = ch2 + 1; ch3 < n; ch3++) {
	  if (theta_char_is_syzygous(ch1, ch2, ch3, g))
	    {
	      count += 1;
	    }
	}
      }
    }
    if (count != 60)
      {
	flint_printf("FAIL\n");
	fflush(stdout);
	flint_abort();
      }
    
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
  }
}
