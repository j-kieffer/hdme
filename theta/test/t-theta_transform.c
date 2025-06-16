#include <stdlib.h>
#include "theta.h"


int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("theta_transform....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
      slong g = 2;
      slong bits = 10 + n_randint(state, 500);
      slong eta_bits = 1;
      slong prec = 5 * bits;
      slong n = n_pow(2, 2*g);
      slong i;
      int res;

      arb_t tol;
      acb_mat_t tau;
      fmpz_mat_t eta;
      acb_ptr th2; /* At eta*tau, computed using theta2_transform */
      acb_ptr th2_test; /* At eta*tau, computed using the naive algorithm */
      
      arb_init(tol);
      acb_mat_init(tau, g, g);
      fmpz_mat_init(eta, 2*g, 2*g);
      th2 = _acb_vec_init(n);
      th2_test = _acb_vec_init(n);
      
      arb_set_si(tol, 1);
      arb_mul_2exp_si(tol, tol, -bits); /* tol is larger than 2^(-prec) */
      
      siegel_halfspace_randtest(tau, state, prec);
      res = siegel_fundamental_domain(tau, eta, tau, tol, prec);
      if (!res)
	{
	  flint_printf("FAIL (fundamental domain)\n");
	  flint_printf("prec = %wd\n", prec);
	  flint_printf("tau = "); acb_mat_printd(tau, 30); flint_printf("\n");
	  flint_abort();
	}
      
      fmpz_mat_randtest_symplectic(eta, state, eta_bits);
      /* fmpz_mat_print(eta); flint_printf("\n\n"); */

      res = theta2_naive(th2, tau, prec);
      if (!res)
	{
	  flint_printf("FAIL (naive theta)\n");
	  flint_printf("prec = %wd\n", prec);
	  flint_printf("tau = "); acb_mat_printd(tau, 30); flint_printf("\n");
	  flint_abort();
	}
      theta2_transform(th2, eta, th2, prec);

      res = siegel_transform(tau, eta, tau, prec);
      if (!res)
	{
	  flint_printf("FAIL (Siegel transform)\n");
	  flint_printf("prec = %wd\n", prec);
	  flint_printf("tau = "); acb_mat_printd(tau, 30); flint_printf("\n");
	  flint_printf("eta = \n"); fmpz_mat_print(eta); flint_printf("\n");
	  flint_abort();
	}

      res = theta2_naive(th2_test, tau, prec);
      if (!res)
	{
	  flint_printf("FAIL (naive theta at eta*tau)\n");
	  flint_printf("prec = %wd\n", prec);
	  flint_printf("tau = "); acb_mat_printd(tau, 30); flint_printf("\n");
	  flint_printf("eta = \n"); fmpz_mat_print(eta); flint_printf("\n");
	  flint_abort();
	}

      /* Rescale: theta0 cannot be exactly zero */
      for (i = 1; i < n; i++)
	{
	  acb_div(&th2[i], &th2[i], &th2[0], prec);
	  acb_div(&th2_test[i], &th2_test[i], &th2_test[0], prec);
	}
      acb_one(&th2[0]);
      acb_one(&th2_test[0]);

      /* Check overlap */
      for (i = 0; i < n; i++)
	{
	  /* acb_printd(&th2[i], 30); flint_printf("\n");
	     acb_printd(&th2_test[i], 30); flint_printf("\n\n"); */
	  if (!acb_overlaps(&th2[i], &th2_test[i]))
	    {
	      flint_printf("FAIL (theta overlap)\n");
	      flint_printf("prec = %wd\n", prec);
	      flint_printf("tau = "); acb_mat_printd(tau, 30); flint_printf("\n");
	      flint_printf("eta = \n"); fmpz_mat_print(eta); flint_printf("\n");
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
      fmpz_mat_clear(eta);
      _acb_vec_clear(th2, n);
      _acb_vec_clear(th2_test, n);
    }
  
  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
