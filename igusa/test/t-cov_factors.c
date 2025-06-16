#include <stdlib.h>

#include "igusa.h"

int main()
{  
  slong iter;
  flint_rand_t state;
  
  flint_printf("cov_factors....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 500 * flint_test_multiplier(); iter++)
    {
      slong nb = 1 + n_randint(state, 10);
      slong mag_bits = 10;
      slong* weights;
      fmpz* I;
      fmpz_factor_t fac;
      slong p1, p2, p3;
      slong k;
      int res1, res2;
      
      weights = flint_malloc(nb * sizeof(slong));
      I = _fmpz_vec_init(nb);
      fmpz_factor_init(fac);
      
      for (k = 0; k < nb; k++)
	{
	  weights[k] = 1 + n_randint(state, 5);
	  fmpz_randtest(&I[k], state, mag_bits);
	}
      k = n_randint(state, nb);
      fmpz_randtest_not_zero(&I[k], state, mag_bits);

      p1 = n_randprime(state, COV_FACTOR_BITS/2, 1);
      cov_rescale_fmpz_si(I, I, p1, nb, weights);
      p2 = n_randprime(state, COV_FACTOR_BITS + 5, 1);
      cov_rescale_fmpz_si(I, I, p2, nb, weights);
      p3 = n_randprime(state, COV_FACTOR_BITS + 5, 1);
      cov_rescale_fmpz_si(I, I, p3, nb, weights);

      cov_factors(fac, I, nb);
      res1 = 0;
      res2 = 0;
      for (k = 0; k < cov_factor_nb(fac); k++)
	{
	  if (fmpz_equal_si(cov_factor_p(fac, k), p1)) res1++;
	  if (fmpz_equal_si(cov_factor_p(fac, k), p2)) res2++;
	}
      if (res1 != 1 || (p2 != p3 && res2 != 1))
	{
	  flint_printf("FAIL\n");
	  flint_printf("p1 = %wd, p2 = %wd, res1 = %wd, res2 = %wd\n", p1, p2, res1, res2);
	  fflush(stdout);
	  flint_abort();
	}      

      flint_free(weights);
      _fmpz_vec_clear(I, nb);
      fmpz_factor_clear(fac);
    }

  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
