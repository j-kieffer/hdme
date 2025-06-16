#include <stdlib.h>
#include "theta.h"


int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("theta2_unif....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 5 * flint_test_multiplier(); iter++)
    {
      slong g = 2;
      slong bits = 1000 + n_randint(state, 5000);
      slong prec = 5 * bits;
      slong mult1 = 1 + n_randint(state, 5);
      slong mult2 = 1 + n_randint(state, 5);
      slong n = n_pow(2, 2*g);
      slong k1, k2;
      slong i;
      int res, skip;

      arb_t tol;
      acb_mat_t tau;
      acb_mat_t w;
      acb_ptr th2; /* Computed using theta2_unif */
      acb_ptr th2_test; /* Computed using the naive algorithm */
      fmpz_mat_t m;
      
      arb_init(tol);
      acb_mat_init(tau, g, g);
      acb_mat_init(w, g, g);
      th2 = _acb_vec_init(n);
      th2_test = _acb_vec_init(n);
      fmpz_mat_init(m, 2*g, 2*g);
      
      arb_set_si(tol, 1);
      arb_mul_2exp_si(tol, tol, -bits); /* tol is larger than 2^(-prec) */
      
      siegel_halfspace_randtest(tau, state, prec);
      /* Closer to the cusp for testing purposes */
      acb_mat_scalar_mul_si(tau, tau, mult1, prec);
      res = siegel_fundamental_domain(tau, m, tau, tol, prec);
      acb_mul_si(acb_mat_entry(tau, 1, 1), acb_mat_entry(tau, 1, 1), mult2, prec);
      res = res && siegel_fundamental_domain(tau, m, tau, tol, prec);
      
      if (!res)
	{
	  flint_printf("FAIL (fundamental domain)\n");
	  flint_printf("prec = %wd\n", prec);
	  flint_printf("tau = "); acb_mat_printd(tau, 30); flint_printf("\n");
	  flint_abort();
	}

      k2 = theta_newton_k2(w, tau, prec);
      k1 = theta_newton_k1(w, w, prec);
      /* flint_printf("prec = %wd, k1 = %wd, k2 = %wd\n", prec, k1, k2); */
      skip = (k1 == 0) && (k2 == 0);
      if (skip) iter -= 1;

      if (!skip) res = theta2_naive(th2_test, tau, prec);
      if (!res)
	{
	  flint_printf("FAIL (naive theta)\n");
	  flint_printf("prec = %wd\n", prec);
	  flint_printf("tau = "); acb_mat_printd(tau, 30); flint_printf("\n");
	  flint_abort();
	}
      
      if (!skip) res = theta2_unif(th2, tau, prec);
      if (!res)
	{
	  flint_printf("FAIL (theta2_unif)\n");
	  flint_printf("prec = %wd, k1 = %wd, k2 = %wd\n", prec, k1, k2);
	  flint_printf("tau = "); acb_mat_printd(tau, 30); flint_printf("\n");
	  flint_abort();
	}
      
      /* Rescale: theta0 cannot be exactly zero */
      if (!skip)
	{
	  for (i = 1; i < n; i++)
	    {
	      acb_div(&th2[i], &th2[i], &th2[0], prec);
	      acb_div(&th2_test[i], &th2_test[i], &th2_test[0], prec);
	    }
	  acb_one(&th2[0]);
	  acb_one(&th2_test[0]);
	}
      
      /* Check overlap */
      for (i = 0; i < n; i++)
	{
	  /* acb_printd(&th2[i], 30); flint_printf("\n");
	     acb_printd(&th2_test[i], 30); flint_printf("\n\n"); */
	  if (!acb_overlaps(&th2[i], &th2_test[i]))
	    {
	      flint_printf("FAIL (theta overlap)\n");
	      flint_printf("prec = %wd, k1 = %wd, k2 = %wd\n", prec, k1, k2);
	      flint_printf("tau = "); acb_mat_printd(tau, 30); flint_printf("\n");
	      flint_printf("i = %wd\n", i);
	      flint_printf("th2[i] = ");
	      acb_printd(&th2[i], 30); flint_printf("\n");
	      flint_printf("th2_test[i] = ");
	      acb_printd(&th2_test[i], 30); flint_printf("\n");
	      flint_abort();
	    }
	}
      
      arb_clear(tol);
      acb_mat_clear(tau);
      acb_mat_clear(w);
      _acb_vec_clear(th2, n);
      _acb_vec_clear(th2_test, n);
      fmpz_mat_clear(m);
    }
  
  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
