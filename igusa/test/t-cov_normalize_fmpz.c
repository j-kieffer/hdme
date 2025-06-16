#include <stdlib.h>

#include "igusa.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("cov_normalize_fmpz....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
      slong nb = 1 + n_randint(state, 10);
      slong mag_bits = 10;
      slong* weights;
      fmpz* I;
      fmpz* S;
      fmpz_t scal;
      slong k;

      I = _fmpz_vec_init(nb);
      S = _fmpz_vec_init(nb);
      fmpz_init(scal);
      weights = flint_malloc(nb * sizeof(slong));
      
      for (k = 0; k < nb; k++)
	{
	  weights[k] = 1 + n_randint(state, 10);
	  fmpz_randtest_not_zero(&I[k], state, mag_bits);
	}
      
      fmpz_set_si(scal, 1 + n_randint(state, n_pow(2, COV_FACTOR_BITS-5)));
      cov_rescale_fmpz(S, I, scal, nb, weights);
      cov_normalize_fmpz(I, I, nb, weights);
      cov_normalize_fmpz(S, S, nb, weights);
     
      if (!_fmpz_vec_equal(S, I, nb))
	{
	  flint_printf("FAIL\n");
	  flint_printf("Scalar: "); fmpz_print(scal); flint_printf("\n");
	  for (k = 0; k < nb; k++)
	    {
	      flint_printf("Weight %wd:\n", weights[k]);
	      fmpz_print(&I[k]); flint_printf("\n");
	      fmpz_print(&S[k]); flint_printf("\n");	      
	    }
	  fflush(stdout);
	  flint_abort();
	}

      _fmpz_vec_clear(I, nb);
      _fmpz_vec_clear(S, nb);
      fmpz_clear(scal);
      flint_free(weights);
    }

  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

      
