#include <stdlib.h>
#include "theta.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("theta_0123half_inverse_diff....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
      slong g = 2;
      slong bits = 20 + n_randint(state, 500);
      slong prec = 5 * bits;
      int res;
      slong i, j;
      
      arb_t tol;
      acb_mat_t tau;
      acb_mat_t tau_half;
      acb_mat_t tau_test;
      fmpz_mat_t m;
      acb_ptr th_half;
      acb_mat_t dtau;
      acb_mat_t dtau_inv;
      acb_mat_t dth_naive;
      acb_t th0;
      mag_t abs;

      arb_init(tol);
      acb_mat_init(tau, g, g);
      acb_mat_init(tau_half, g, g);
      acb_mat_init(tau_test, g, g);
      fmpz_mat_init(m, 2*g, 2*g);
      th_half = _acb_vec_init(4);
      acb_mat_init(dtau, 3, 3);
      acb_mat_init(dtau_inv, 3, 3);
      acb_mat_init(dth_naive, 3, 3);
      acb_init(th0);
      mag_init(abs);
      
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
      
      /* Custom value */
      /* acb_set_si_si(acb_mat_entry(tau, 0, 0), 0, 2);
      acb_set_si_si(acb_mat_entry(tau, 1, 1), 0, 3);
      acb_one(acb_mat_entry(tau, 0, 1));
      acb_mul_2exp_si(acb_mat_entry(tau, 0, 1), acb_mat_entry(tau, 0, 1), -200);
      acb_set(acb_mat_entry(tau, 1, 0), acb_mat_entry(tau, 0, 1)); */
      
      acb_mat_scalar_mul_2exp_si(tau_half, tau, -1);

      res = theta_0123_naive(th_half, tau_half, prec);
      if (!res)
	{
	  flint_printf("FAIL (naive theta)\n");
	  flint_printf("res = %wd\n", res);
	  flint_printf("tau_half = "); acb_mat_printd(tau_half, 30); flint_printf("\n");
	  flint_abort();
	}
      acb_set(th0, &th_half[0]);
      _acb_vec_scalar_div(th_half, th_half, 4, th0, prec);

      theta_0123half_inverse_diff(dtau, tau, th_half, prec);
      acb_mat_inv(dtau_inv, dtau, prec);
      theta_0123half_diff_naive(dth_naive, tau, prec);

      /* dtau_inv should be small. */
      for (i = 0; i < 3; i++)
	{
	  for (j = 0; j < 3; j++)
	    {
	      acb_get_mag(abs, acb_mat_entry(dtau_inv, i, j));
	      if (mag_cmp_2exp_si(abs, 0) > 0) res = 0;
	    }
	}
      
      if (res == 0)
	{
	  flint_printf("FAIL (dtau_inv too large)\n");
	  flint_printf("prec = %wd\n", prec);
	  flint_printf("tau = "); acb_mat_printd(tau, 30); flint_printf("\n");
	  for (i = 0; i < 4; i++)
	    {
	      flint_printf("th_half[%wd] = ", i); acb_printd(&th_half[i], 30);
	      flint_printf("\n");
	    }
	  flint_printf("dtau = "); acb_mat_printd(dtau, 30); flint_printf("\n");
	  flint_printf("dtau_inv = "); acb_mat_printd(dtau_inv, 30); flint_printf("\n");
	  flint_printf("dth_naive = "); acb_mat_printd(dth_naive, 30); flint_printf("\n");
	  acb_mat_sub(dtau_inv, dtau_inv, dth_naive, prec);
	  flint_abort();
	}
	      
      /* dtau_inv and dth_naive should coincide. */
      if (!acb_mat_overlaps(dtau_inv, dth_naive))
	{
	  flint_printf("FAIL (dth_naive)\n");
	  flint_printf("prec = %wd\n", prec);
	  flint_printf("tau = "); acb_mat_printd(tau, 30); flint_printf("\n");
	  for (i = 0; i < 4; i++)
	    {
	      flint_printf("th_half[%wd] = ", i); acb_printd(&th_half[i], 30);
	      flint_printf("\n");
	    }
	  flint_printf("dtau = "); acb_mat_printd(dtau, 30); flint_printf("\n");
	  flint_printf("dtau_inv = "); acb_mat_printd(dtau_inv, 30); flint_printf("\n");
	  flint_printf("dth_naive = "); acb_mat_printd(dth_naive, 30); flint_printf("\n");
	  acb_mat_sub(dtau_inv, dtau_inv, dth_naive, prec);
	  flint_printf("error = "); acb_mat_printd(dtau_inv, 30); flint_printf("\n");
	  flint_abort();
	}
      
      arb_clear(tol);
      acb_mat_clear(tau);
      acb_mat_clear(tau_half);
      acb_mat_clear(tau_test);
      fmpz_mat_clear(m);
      _acb_vec_clear(th_half, 4);
      acb_mat_clear(dtau);
      acb_mat_clear(dtau_inv);
      acb_mat_clear(dth_naive);
      acb_clear(th0);
      mag_clear(abs);
    }

  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

