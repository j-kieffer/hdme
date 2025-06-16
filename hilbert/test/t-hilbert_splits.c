#include <stdlib.h>

#include "hilbert.h"

int main()
{
  slong delta;
  slong ell;
  int splits;
  fmpz_poly_t beta;
  
  flint_printf("hilbert_splits....");
  fflush(stdout);

  fmpz_poly_init(beta);
    
  for (delta = 5; delta < 100; delta++)
    {
      if (hilbert_is_fundamental(delta))
	{
	  for (ell = 2; ell < 200; ell++)
	    {
	      if (n_is_prime(ell))
		{
		  splits = hilbert_splits(beta, ell, delta);
		  /*
		  flint_printf("%wd splits correctly in Q(sqrt(%wd)): %d\n",
		  ell, delta, splits); */
		  if (splits)
		    {
		      /* flint_printf("%wd splits correctly in Q(sqrt(%wd)), generator: ",
				   ell, delta);
		      fmpz_poly_print_pretty(beta, "x");
		      flint_printf("\n"); */
		    }
		  
		}
	    }
	}
    }

  fmpz_poly_clear(beta);
  
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

