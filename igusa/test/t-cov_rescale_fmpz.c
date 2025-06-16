#include <stdlib.h>

#include "igusa.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("cov_rescale_fmpz....");
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
	  fmpz_randtest(&I[k], state, mag_bits);
	}
      fmpz_randtest_not_zero(scal, state, mag_bits);
      cov_rescale_fmpz(S, I, scal, nb, weights);

      if (!cov_divisible_fmpz(S, scal, nb, weights))
	{	  
	  flint_printf("FAIL\n");
	  flint_printf("Scalar: "); fmpz_print(scal); flint_printf("\n");
	  for (k = 0; k < nb; k++)
	    {
	      fmpz_print(&S[k]); flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_abort();
	}

      cov_divexact_fmpz(S, S, scal, nb, weights);
      if (!_fmpz_vec_equal(S, I, nb))
	{
	  flint_printf("FAIL\n");
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
