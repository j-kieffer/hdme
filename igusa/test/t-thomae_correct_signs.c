#include <stdlib.h>

#include "igusa.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("thomae_correct_signs....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 5 * flint_test_multiplier(); iter++)
    {
      slong prec = 1000 + n_randint(state, 5000);
      acb_mat_t tau;
      arb_t abs;
      acb_ptr th2;
      acb_ptr I, IC;
      acb_poly_t crv;
      acb_ptr roots;
      slong perm[1];
      slong signs[1];
      slong k;
      int res;

      acb_mat_init(tau, 2, 2);
      arb_init(abs);
      th2 = _acb_vec_init(16);
      I = _acb_vec_init(4);
      IC = _acb_vec_init(4);
      acb_poly_init(crv);
      roots = _acb_vec_init(6);

      /* Choose tau not too close to the diagonal for the root finding method */
      do
	{
	  siegel_fundamental_domain_randtest(tau, state, prec);
	  acb_abs(abs, acb_mat_entry(tau, 0, 1), prec);
	  arb_mul_2exp_si(abs, abs, 50);
	  arb_sub_si(abs, abs, 1, prec);
	} while (!arb_is_positive(abs));
      
      
      res = theta2_unif(th2, tau, prec);
      if (!res)
	{
	  flint_printf("FAIL (theta2_unif)\n");
	  fflush(stdout);
	  flint_abort();
	}
      igusa_from_theta2(I, th2, prec);
      igusa_IC(IC, I, prec);
      res = mestre(crv, IC, prec);
      if (!res)
	{
	  flint_printf("FAIL (mestre)\n");
	  fflush(stdout);
	  flint_abort();
	}
      res = thomae_roots(roots, crv, prec);
      if (!res)
	{
	  flint_printf("FAIL (roots)\n");
	  acb_mat_printd(tau, 30);
	  flint_printf("I: \n");
	  for (k = 0; k < 4; k++)
	    {
	      acb_printd(&I[k], 30); flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_printf("Curve: ");
	  acb_poly_printd(crv, 30);
	  flint_printf("\n");
	  flint_printf("Roots: \n");
	  for (k = 0; k < 6; k++)
	    {
	      acb_printd(&roots[k], 30); flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_abort();
	}
      res = thomae_correct_signs(perm, signs, roots, I, prec);
      if (!res)
	{
	  flint_printf("FAIL (no candidate found)\n");
	  acb_mat_printd(tau, 30);
	  flint_printf("Roots: \n");
	  for (k = 0; k < 6; k++)
	    {
	      acb_printd(&roots[k], 30); flint_printf("\n");
	    }
	  flint_printf("I: \n");
	  for (k = 0; k < 4; k++)
	    {
	      acb_printd(&I[k], 30); flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_abort();
	}

      acb_mat_clear(tau);
      arb_clear(abs);
      _acb_vec_clear(th2, 16);
      _acb_vec_clear(I, 4);
      _acb_vec_clear(IC, 4);
      acb_poly_clear(crv);
      _acb_vec_clear(roots, 6);
    }
  
  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
