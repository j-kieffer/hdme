#include <stdlib.h>
#include "theta.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("theta_newton_k1_k2....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
      slong g = 2;
      slong bits = 10 + n_randint(state, 1000);
      slong mult = 1 + n_randint(state, 10);
      slong prec = 5 * bits;
      int res;
      slong k1, k2;
      int flag_naive, flag_newton;
      
      arb_t tol;
      acb_mat_t tau;
      fmpz_mat_t m;
      
      arb_init(tol);
      acb_mat_init(tau, g, g);
      fmpz_mat_init(m, 2*g, 2*g);

      arb_set_si(tol, 1);
      arb_mul_2exp_si(tol, tol, -bits); /* tol is larger than 2^(-prec) */

      siegel_halfspace_randtest(tau, state, prec);
      /* Closer to the cusp for testing purposes */
      acb_mat_scalar_mul_si(tau, tau, mult, prec); 
      
      res = siegel_fundamental_domain(tau, m, tau, tol, prec);
      if (!res)
	{
	  flint_printf("FAIL (fundamental domain)\n");
	  flint_printf("res = %wd\n", res);
	  flint_printf("tau = "); acb_mat_printd(tau, 30); flint_printf("\n");
	  flint_abort();
	}



      k2 = theta_newton_k2(tau, tau, prec);
      siegel_reduce_real(tau, m, tau, tol, prec);
      res = siegel_is_weakly_reduced(tau, tol, prec);

      if (!res)
	{
	  flint_printf("FAIL (k2)\n");
	  flint_printf("res = %wd, k2 = %wd\n", res, k2);
	  flint_printf("tau = "); acb_mat_printd(tau, 30); flint_printf("\n");
	  flint_abort();
	}

      k1 = theta_newton_k1(tau, tau, prec);
      res = siegel_is_weakly_reduced(tau, tol, prec);

      if (!res)
	{
	  flint_printf("FAIL (k1)\n");
	  flint_printf("res = %wd, k1 = %wd,\n", res, k1);
	  flint_printf("tau = "); acb_mat_printd(tau, 30); flint_printf("\n");
	  flint_abort();
	}
      
      flag_naive = theta_use_naive(tau, prec);
      flag_newton = theta_use_newton(tau, prec);

      if (!flag_naive && !flag_newton)
	{
	  flint_printf("FAIL (k1+k2)\n");
	  flint_printf("res = %wd\n", res);
	  flint_printf("tau = "); acb_mat_printd(tau, 30); flint_printf("\n");
	  flint_abort();
	}
      
      /* flint_printf("prec = %wd\n", prec); */
      /* flint_printf("tau = "); acb_mat_printd(tau, 30); flint_printf("\n"); */
      /* flint_printf("k2 = %wd, k1 = %wd\n", k2, k1); */
      /* flint_printf("flag_naive = %i, flag_newton = %i\n\n", flag_naive, flag_newton); */
      
      arb_clear(tol);
      acb_mat_clear(tau);
      fmpz_mat_clear(m);
    }

  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

  
