#include <stdlib.h>
#include "igusa.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("cov_no_rescale_to_one....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
      slong nb = 2 + n_randint(state, 10);
      slong prec = 500;
      slong mag_bits = 10;
      slong* weights;
      acb_ptr I, test;
      acb_t root;
      acb_t scal;
      slong k;

      I = _acb_vec_init(nb);
      test = _acb_vec_init(nb);
      acb_init(root);
      acb_init(scal);
      weights = flint_malloc(nb * sizeof(slong));
      
      for (k = 0; k < nb; k++)
	{
	  weights[k] = 1 + n_randint(state, 10);
	  acb_one(&I[k]);
	}
      acb_randtest_precise(scal, state, prec, mag_bits);

      cov_rescale(I, I, scal, nb, weights, prec);      
      
      if (cov_no_rescale_to_one(I, nb, weights, prec))
	{
	  flint_printf("FAIL\n");
	  flint_printf("Scalar: ");
	  acb_printd(scal, 10); flint_printf("\n");
	  for (k = 0; k < nb; k++)
	    {
	      flint_printf("Weight %wd: ", weights[k]);
	      acb_printd(&I[k], 10); flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_abort();
	}
      
      for (k = 0; k < nb; k++)
	{
	  acb_one(&I[k]);
	}
      acb_zero(&I[0]);
      if (!cov_no_rescale_to_one(I, nb, weights, prec))
	{
	  flint_printf("FAIL\n");
	  flint_printf("Scalar: ");
	  acb_printd(scal, 10); flint_printf("\n");
	  for (k = 0; k < nb; k++)
	    {
	      flint_printf("Weight %wd: ", weights[k]);
	      acb_printd(&I[k], 10); flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_abort();
	}

      _acb_vec_clear(I, nb);
      _acb_vec_clear(test, nb);
      acb_clear(root);
      acb_clear(scal);
      flint_free(weights);
    }

  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
