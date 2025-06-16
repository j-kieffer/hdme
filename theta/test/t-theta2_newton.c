#include <stdlib.h>
#include "theta.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("theta2_newton....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
      slong g = 2;
      slong bits = 10 + n_randint(state, 2500);
      slong prec = 5 * bits;
      int res;
      int i;

      arb_t tol;
      acb_mat_t tau;
      fmpz_mat_t m;
      acb_ptr th;
      acb_ptr th_test;
      acb_t th0;

      arb_init(tol);
      acb_mat_init(tau, g, g);
      fmpz_mat_init(m, 2*g, 2*g);
      th = _acb_vec_init(16);
      th_test = _acb_vec_init(16);
      acb_init(th0);
      
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

      res = theta2_naive(th, tau, prec);
      if (!res)
	{
	  flint_printf("FAIL (naive theta)\n");
	  flint_printf("res = %wd\n", res);
	  flint_printf("tau = "); acb_mat_printd(tau, 30); flint_printf("\n");
	  for (i = 0; i < 16; i++)
	    {
	      flint_printf("th[%wd] = ", i); acb_printd(&th[i], 30); flint_printf("\n");
	    }
	  flint_abort();
	}
      acb_set(th0, &th[0]);
      _acb_vec_scalar_div(th, th, 16, th0, prec);

      res = theta2_newton(th_test, tau, prec);
      if (!res)
	{
	  flint_printf("FAIL (newton)\n");
	  flint_printf("res = %wd\n", res);
	  flint_printf("tau = "); acb_mat_printd(tau, 30); flint_printf("\n");
	  for (i = 0; i < 16; i++)
	    {
	      flint_printf("th[%wd] = ", i); acb_printd(&th[i], 30); flint_printf("\n");
	    }
	  flint_abort();
	}
      acb_set(th0, &th_test[0]);
      _acb_vec_scalar_div(th_test, th_test, 16, th0, prec);
      
      /* flint_printf("prec = %wd\n", prec);
      flint_printf("tau = "); acb_mat_printd(tau, 30); flint_printf("\n");
      acb_sub(th0, &th_test[1], &th[1], prec);
      flint_printf("th_test[%wd] - th[%wd] = ", 1, 1);
      acb_printd(th0, 30); flint_printf("\n"); */
      
      for (i = 0; i < 16; i++)
	{
	  if (!acb_overlaps(&th[i], &th_test[i])) res = 0;
	}
      if (res == 0)
	{
	  flint_printf("FAIL\n");
	  flint_printf("tau = "); acb_mat_printd(tau, 30); flint_printf("\n");
	  for (i = 0; i < 16; i++)
	    {
	      flint_printf("th[%wd] = ", i);
	      acb_printd(&th[i], 30); flint_printf("\n");
	      flint_printf("th_test[%wd] = ", i);
	      acb_printd(&th_test[i], 30); flint_printf("\n");
	      acb_sub(th0, &th_test[i], &th[i], prec);
	      flint_printf("th_test[%wd] - th[%wd] = ", i, i);
	      acb_printd(th0, 30); flint_printf("\n");
	    }
	  flint_abort();
	}
      
      arb_clear(tol);
      acb_mat_clear(tau);
      fmpz_mat_clear(m);
      _acb_vec_clear(th, 16);
      _acb_vec_clear(th_test, 16);
      acb_clear(th0);
    }
  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
