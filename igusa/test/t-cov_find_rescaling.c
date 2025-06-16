#include <stdlib.h>

#include "igusa.h"

int main()
{  
  slong iter;
  flint_rand_t state;
  
  flint_printf("cov_find_rescaling....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {      
      slong nb = 1 + n_randint(state, 4);
      slong prec = 1000;
      slong mag_bits = 10;
      slong* weights;
      fmpz* I;
      acb_ptr S;
      acb_t scal;
      slong k;
      int res;
      
      I = _fmpz_vec_init(nb);
      S = _acb_vec_init(nb);
      acb_init(scal);
      weights = flint_malloc(nb * sizeof(slong));

      for (k = 0; k < nb; k++)
	{
	  weights[k] = 1 + n_randint(state, 10);
	  if (k == 0) fmpz_randtest_not_zero(&I[k], state, mag_bits);
	  else if (n_randint(state, 2) == 0) fmpz_randtest(&I[k], state, mag_bits);
	  else fmpz_zero(&I[k]);
	  acb_set_fmpz(&S[k], &I[k]);
	}
      acb_randtest_precise(scal, state, prec, mag_bits);
      cov_rescale(S, S, scal, nb, weights, prec);
      /*for (k = 0; k < nb; k++)
	{
	  flint_printf("Weight %wd: ", weights[k]);
	  acb_printd(&S[k], 10); flint_printf("\n");
	  fmpz_print(&I[k]); flint_printf("\n");
	  }*/
      
      cov_find_rescaling(scal, S, I, nb, weights, prec);
      acb_inv(scal, scal, prec);
      cov_rescale(S, S, scal, nb, weights, prec);

      res = 1;
      for (k = 0; k < nb; k++)
	{
	  if (!acb_contains_fmpz(&S[k], &I[k])) res = 0;
	}
      if (!res)
	{	  
	  flint_printf("FAIL\n");
	  flint_printf("Scalar: "); acb_printd(scal, 10); flint_printf("\n");
	  for (k = 0; k < nb; k++)
	    {
	      flint_printf("Weight %wd: ", weights[k]);
	      acb_printd(&S[k], 10); flint_printf("\n");
	      fmpz_print(&I[k]); flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_abort();
	}      
      
      _fmpz_vec_clear(I, nb);
      _acb_vec_clear(S, nb);
      acb_clear(scal);
      flint_free(weights);
    }

  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
