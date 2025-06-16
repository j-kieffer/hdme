#include <stdlib.h>
#include "theta.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("theta_0123half_inverse....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 500 * flint_test_multiplier(); iter++)
    {
      slong g = 2;
      slong bits = 10 + n_randint(state, 500);
      slong prec = 5 * bits;
      int res;
      int i;

      arb_t tol;
      acb_mat_t tau;
      acb_mat_t tau_half;
      acb_mat_t tau_test;
      fmpz_mat_t m;
      acb_ptr th_half;

      arb_init(tol);
      acb_mat_init(tau, g, g);
      acb_mat_init(tau_half, g, g);
      acb_mat_init(tau_test, g, g);
      fmpz_mat_init(m, 2*g, 2*g);
      th_half = _acb_vec_init(4);
      
      arb_set_si(tol, 1);
      arb_mul_2exp_si(tol, tol, -bits); /* tol is larger than 2^(-prec) */
      
      siegel_halfspace_randtest(tau, state, prec);
      res = siegel_fundamental_domain(tau, m, tau, tol, prec);
      if (!res)
	{
	  flint_printf("FAIL (fundamental domain)\n");
	  flint_printf("res = %wd\n", res);
	  flint_printf("tau = "); acb_mat_printd(tau, 30); flint_printf("\n");
	  flint_abort();
	}
      
      acb_mat_scalar_mul_2exp_si(tau_half, tau, -1);

      res = theta_0123_naive(th_half, tau_half, prec);
      if (!res)
	{
	  flint_printf("FAIL (naive theta)\n");
	  flint_printf("res = %wd\n", res);
	  flint_printf("tau_half = "); acb_mat_printd(tau_half, 30); flint_printf("\n");
	  flint_abort();
	}

      res = theta_0123half_inverse(tau_test, th_half, prec);
      if (!res || !acb_mat_overlaps(tau_test, tau))
	{
	  flint_printf("FAIL (inverse theta)\n");
	  flint_printf("res = %wd\n", res);
	  flint_printf("tau = "); acb_mat_printd(tau, 30); flint_printf("\n");
	  flint_printf("tau_half = "); acb_mat_printd(tau_half, 30); flint_printf("\n");
	  for (i = 0; i < 4; i++)
	    {
	      flint_printf("th_half[%wd] = ", i); acb_printd(&th_half[i], 30); flint_printf("\n");
	    }
	  flint_printf("tau_test = "); acb_mat_printd(tau_test, 30); flint_printf("\n");
	  flint_abort();
	}
      
      arb_clear(tol);
      acb_mat_clear(tau);
      acb_mat_clear(tau_half);
      acb_mat_clear(tau_test);
      fmpz_mat_clear(m);
      _acb_vec_clear(th_half, 4);
    }

  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
