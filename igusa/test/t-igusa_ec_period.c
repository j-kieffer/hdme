#include <stdlib.h>
#include "igusa.h"

int main()
{
  
  slong iter;
  flint_rand_t state;
  
  flint_printf("igusa_ec_period....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
      acb_t j;
      acb_t tau;
      acb_t j_test;
      arf_t tol;
      slong prec = 1000;
      slong mag_bits = 5;
      int r;
      int valid;
      
      acb_init(j);
      acb_init(tau);
      acb_init(j_test);
      arf_init(tol);

      arf_one(tol);
      arf_mul_2exp_si(tol, tol, -prec/2);

      r = n_randint(state, 10);
      if (r == 0) acb_zero(j);
      else if (r == 1) acb_set_si(j, 1728);
      else
	{
	  /* Not too close to 0, 1728 */
	  valid = 0;
	  while (!valid)
	    {
	      acb_randtest_precise(j, state, prec, mag_bits);
	      acb_set(j_test, j);
	      acb_add_error_arf(j_test, tol);
	      valid = !acb_contains_zero(j_test);
	      acb_sub_si(j_test, j, 1728, prec);
	      acb_add_error_arf(j_test, tol);
	      valid = valid && !acb_contains_zero(j_test);
	    }
	}
      
      r = igusa_ec_period(tau, j, prec);
      acb_modular_j(j_test, tau, prec);

      if (!r || !acb_overlaps(j_test, j))
	{
	  flint_printf("FAIL\n");
	  flint_printf("r = %wd, j, tau, j_test:\n", r);
	  acb_printd(j, 10); flint_printf("\n");
	  acb_printd(tau, 10); flint_printf("\n");
	  acb_printd(j_test, 10); flint_printf("\n");
	  fflush(stdout);
	  flint_abort();
	}
      
      acb_clear(j);
      acb_clear(tau);
      acb_clear(j_test);
      arf_clear(tol);
    }
     
  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}


      
		
