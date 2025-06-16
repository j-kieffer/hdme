#include <stdlib.h>

#include "igusa.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("cov_rescale_fmpz_si....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
      slong nb = 1 + n_randint(state, 10);
      slong mag_bits = 10;
      slong* weights;
      fmpz* I;
      fmpz* S;
      slong scal = 1 + n_randint(state, 1000);
      slong k;

      I = _fmpz_vec_init(nb);
      S = _fmpz_vec_init(nb);
      weights = flint_malloc(nb * sizeof(slong));
      
      for (k = 0; k < nb; k++)
	{
	  weights[k] = 1 + n_randint(state, 10);
	  fmpz_randtest(&I[k], state, mag_bits);
	}
      
      cov_rescale_fmpz_si(S, I, scal, nb, weights);
      cov_divexact_fmpz_si(S, S, scal, nb, weights);
      
      if (!_fmpz_vec_equal(S, I, nb))
	{
	  flint_printf("FAIL\n"); 
	  for (k = 0; k < nb; k++)
	    {
	      flint_printf("Scal = %wd, weight %wd: ", scal, weights[k]);
	      fmpz_print(&S[k]); flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_abort();
	}
      
      _fmpz_vec_clear(I, nb);
      _fmpz_vec_clear(S, nb);
      flint_free(weights);
    }

  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
