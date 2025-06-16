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
      acb_t tau;
      acb_t j, test;
      acb_ptr kp2;
      acb_ptr th;
      arf_t tol;
      int r;
      slong k;
      slong prec = 1000;
      slong mag_bits = 2;

      acb_init(tau);
      acb_init(j);
      acb_init(test);
      th = _acb_vec_init(4);
      kp2 = _acb_vec_init(6);
      arf_init(tol);

      arb_randtest_precise(acb_realref(tau), state, prec, mag_bits);
      arb_randtest_precise(acb_imagref(tau), state, prec, 1);
      arb_exp(acb_imagref(tau), acb_imagref(tau), prec);
      arf_one(tol);
      arf_mul_2exp_si(tol, tol, -prec/2);

      acb_modular_j(j, tau, prec);

      acb_set(test, j);
      acb_add_error_arf(test, tol);
      r = !acb_contains_zero(test);      
      acb_sub_si(test, j, 1728, prec);
      acb_add_error_arf(test, tol);
      r = r && !acb_contains_zero(test);

      if (r)
	{
	    
	  acb_zero(test);
	  acb_modular_theta(&th[0], &th[1], &th[2], &th[3], test, tau, prec);
	  acb_div(test, &th[3], &th[2], prec);
	  acb_sqr(test, test, prec);
	  acb_sqr(test, test, prec);
	  
	  r = igusa_ec_possible_kp2(kp2, j, prec);
	  if (!r)
	    {
	      flint_printf("FAIL (kp2)\n");
	      flint_printf("tau, j, kp2:\n");
	      acb_printd(tau, 10); flint_printf("\n");
	      acb_printd(j, 10); flint_printf("\n");
	      acb_printd(test, 10); flint_printf("\n");
	      fflush(stdout);
	      flint_abort();
	    }

	  r = 0;
	  for (k = 0; k < 6; k++)
	    {
	      if (acb_overlaps(test, &kp2[k])) r = 1;
	    }
	  if (!r)
	    {
	      flint_printf("FAIL (kp2 not found)\n");
	      flint_printf("tau, j, kp2:\n");
	      acb_printd(tau, 10); flint_printf("\n");
	      acb_printd(j, 10); flint_printf("\n");
	      acb_printd(test, 10); flint_printf("\n");
	      flint_printf("kp2 vector:\n");
	      for (k = 0; k < 6; k++)
		{
		  acb_printd(&kp2[k], 10); flint_printf("\n");
		}
	      fflush(stdout);
	      flint_abort();
	    }	  
	}
      
      acb_clear(tau);
      acb_clear(j);
      acb_clear(test);
      _acb_vec_clear(th, 4);
      _acb_vec_clear(kp2, 6);
      arf_clear(tol);
    }

  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

