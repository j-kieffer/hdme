#include <stdlib.h>
#include "igusa.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("cov_rescale....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
      slong nb = 1 + n_randint(state, 10);
      slong prec = 500;
      slong mag_bits = 10;
      slong* weights;
      acb_ptr I, S;
      acb_t scal;
      slong k;

      I = _acb_vec_init(nb);
      S = _acb_vec_init(nb);
      acb_init(scal);
      weights = flint_malloc(nb * sizeof(slong));
      
      for (k = 0; k < nb; k++)
	{
	  weights[k] = 1 + n_randint(state, 10);
	  acb_randtest_precise(&I[k], state, prec, mag_bits);
	}
      acb_randtest_precise(scal, state, prec, mag_bits);
      cov_rescale(S, I, scal, nb, weights, prec);

      if (cov_distinct(S, I, nb, weights, prec))
	{
	  flint_printf("FAIL\n");
	  flint_printf("Scalar: "); acb_printd(scal, 10); flint_printf("\n");
	  for (k = 0; k < nb; k++)
	    {
	      flint_printf("Weight %wd: ", weights[k]);
	      acb_div(scal, &S[k], &I[k], prec);
	      acb_printd(scal, 10); flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_abort();
	}

      _acb_vec_clear(I, nb);
      _acb_vec_clear(S, nb);
      acb_clear(scal);
      flint_free(weights);
    }

  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
