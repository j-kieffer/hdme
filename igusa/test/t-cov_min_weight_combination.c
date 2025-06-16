#include <stdlib.h>

#include "igusa.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("cov_min_weight_combination....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
      slong nb = 1 + n_randint(state, 5);
      slong* weights;
      fmpz* I;
      slong k;
      slong wt;
      slong test;
      slong* exps;
      int res;

      exps = flint_malloc(nb * sizeof(slong));
      weights = flint_malloc(nb * sizeof(slong));
      I = _fmpz_vec_init(nb);

      for (k = 0; k < nb; k++)
	{
	  weights[k] = 1 + n_randint(state, 10);
	  if (n_randint(state, 2) == 0) fmpz_one(&I[k]);
	  else fmpz_zero(&I[k]);
	}
      k = n_randint(state, nb);
      fmpz_one(&I[k]);

      cov_min_weight_combination(&wt, exps, I, nb, weights);

      res = 1;
      test = 0;
      for (k = 0; k < nb; k++)
	{
	  if (fmpz_is_zero(&I[k]) && (exps[k] != 0)) res = 0;
	  if (!fmpz_is_zero(&I[k]) && (weights[k] % wt != 0)) res = 0;
	  test += exps[k] * weights[k];
	}
      if (!res || (test != wt))
	{
	  flint_printf("FAIL\n");
	  flint_printf("wt = %wd\n", wt);
	  for (k = 0; k < nb; k++)
	    {
	      flint_printf("Weight %wd, is zero? %wd, exponent %wd\n",
			   weights[k], fmpz_is_zero(&I[k]), exps[k]);
	    }
	  fflush(stdout);
	  flint_abort();
	}

      flint_free(exps);
      flint_free(weights);
      _fmpz_vec_clear(I, nb);
    }

  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
