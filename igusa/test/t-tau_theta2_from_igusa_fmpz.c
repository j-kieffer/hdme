#include <stdlib.h>

#include "igusa.h"

int main()
{
  
  slong iter;
  flint_rand_t state;
  
  flint_printf("tau_theta2_from_igusa_fmpz....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 15 * flint_test_multiplier(); iter++)
    {
      fmpz* I;
      acb_ptr I_acb;
      acb_mat_t tau;
      acb_ptr theta2;
      acb_ptr I_test;
      acb_poly_t crv;
      acb_t scal;
      slong k;
      int r;
      slong prec = 1000;
      slong mag_bits = 10;
      slong weights[4] = IGUSA_WEIGHTS;
      int print = 0;

      I = _fmpz_vec_init(4);
      I_acb = _acb_vec_init(4);
      acb_mat_init(tau, 2, 2);
      theta2 = _acb_vec_init(16);
      I_test = _acb_vec_init(4);
      acb_poly_init(crv);
      acb_init(scal);

      for (k = 0; k < 4; k++) fmpz_randtest_not_zero(&I[k], state, mag_bits);
      if (iter % 2 == 0) fmpz_zero(igusa_psi4(I));
      if (iter % 3 == 0) fmpz_zero(igusa_psi6(I));
      if (iter % 5 == 0) fmpz_zero(igusa_chi10(I));
      if (iter % 5 == 1) fmpz_zero(igusa_chi12(I)); /* chi10, chi12 have no common zeroes */
      for (k = 0; k < 4; k++) acb_set_fmpz(&I_acb[k], &I[k]);

      if (print)
	{
	  flint_printf("I:\n");
	  for (k = 0; k < 4; k++)
	    {
	      fmpz_print(&I[k]); flint_printf("\n");
	    }
	  fflush(stdout);
	}
      if (print && igusa_is_g2_curve_fmpz(I))
	{
	  igusa_IC_fmpz(I, I);
	  mestre_fmpz(crv, I, prec);
	  flint_printf("Curve:\n");
	  acb_poly_printd(crv, 10); flint_printf("\n");
	  igusa_from_IC_fmpz(I, I);
	}	

      r = tau_theta2_from_igusa_fmpz(tau, theta2, I, prec);
      igusa_from_theta2(I_test, theta2, prec);

      if (print)
	{
	  flint_printf("I_test:\n");
	  for (k = 0; k < 4; k++)
	    {
	      acb_printd(&I_test[k], 10); flint_printf("\n");
	    }
	  fflush(stdout);
	}

      if (!r || cov_distinct(I_test, I_acb, 4, weights, prec))
	{
	  flint_printf("FAIL\n");
	  fflush(stdout);
	  flint_abort();
	}
      cov_find_rescaling(scal, I_test, I, 4, weights, prec);

      _fmpz_vec_clear(I, 4);
      _acb_vec_clear(I_acb, 4);
      acb_mat_clear(tau);
      _acb_vec_clear(theta2, 16);
      _acb_vec_clear(I_test, 4);
      acb_poly_clear(crv);
      acb_clear(scal);
    }

  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

