#include <stdlib.h>

#include "igusa.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("tau_from_igusa....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 5 * flint_test_multiplier(); iter++)
    {
      slong prec = 2000 + n_randint(state, 4000);
      slong weights[4] = IGUSA_WEIGHTS;
      acb_mat_t tau;
      acb_ptr I, I_test;
      slong mag_bits = 1 + n_randint(state, 10);
      slong k;
      int res;
      
      acb_mat_init(tau, 2, 2);
      I = _acb_vec_init(4);
      I_test = _acb_vec_init(4);

      /* Generate Igusa invariants */
      for (k = 0; k < 4; k++) acb_randtest_precise(&I_test[k], state, prec, mag_bits);
      res = tau_from_igusa(tau, I_test, prec);

      if (!res)
	{ 
	  flint_printf("FAIL (tau)\n");
	  for (k = 0; k < 4; k++)
	    {
	      acb_printd(&I_test[k], 30); flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_abort();
	}

      /* Compute Igusa invariants and check */
      res = igusa_from_tau(I, tau, prec);
      if (!res)
	{ 
	  flint_printf("FAIL (I)\n");
	  for (k = 0; k < 4; k++)
	    {
	      acb_printd(&I_test[k], 30); flint_printf("\n");
	    }
	  acb_mat_printd(tau, 10);
	  fflush(stdout);
	  flint_abort();
	}

      if (cov_distinct(I, I_test, 4, weights, prec))
	{
	  flint_printf("FAIL (overlap)\n");
	  for (k = 0; k < 4; k++)
	    {
	      acb_printd(&I_test[k], 30); flint_printf("\n");
	      acb_printd(&I[k], 30); flint_printf("\n");
	    }
	  acb_mat_printd(tau, 10);
	  fflush(stdout);
	  flint_abort();
	}

      acb_mat_clear(tau);
      _acb_vec_clear(I_test, 4);
      _acb_vec_clear(I, 4);
    }
  
  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
      
