#include <stdlib.h>
#include "igusa.h"

int main()
{
  
  slong iter;
  flint_rand_t state;
  
  flint_printf("igusa_ec_j1j2....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
      fmpz* I;
      acb_ptr j;
      acb_t tau1, tau2;
      acb_mat_t tau;
      acb_ptr I_acb;
      acb_ptr I_test;
      acb_t scal;
      slong prec = 1000;
      slong mag_bits = 20;
      slong weights[4] = IGUSA_WEIGHTS;
      int r;
      slong k;

      I = _fmpz_vec_init(4);
      j = _acb_vec_init(2);
      acb_init(tau1);
      acb_init(tau2);
      acb_mat_init(tau, 2, 2);
      I_acb = _acb_vec_init(4);
      I_test = _acb_vec_init(4);
      acb_init(scal);

      for (k = 0; k < 4; k++) fmpz_randtest_not_zero(&I[k], state, mag_bits);	
      fmpz_zero(igusa_chi10(I));
      if (n_randint(state, 3) == 0) fmpz_zero(igusa_psi6(I));
      if (n_randint(state, 3) == 0) fmpz_zero(igusa_psi4(I));

      for (k = 0; k < 4; k++) acb_set_fmpz(&I_acb[k], &I[k]);
      igusa_ec_j1j2(j, I, prec);
      igusa_ec_period(tau1, &j[0], prec);
      igusa_ec_period(tau2, &j[1], prec);
      if (arb_lt(acb_imagref(tau2), acb_imagref(tau1))) acb_swap(tau1, tau2);      
      acb_mat_zero(tau);
      acb_set(acb_mat_entry(tau, 0, 0), tau1);
      acb_set(acb_mat_entry(tau, 1, 1), tau2);
      r = igusa_from_tau(I_test, tau, prec);
      /*
      flint_printf("I_test:\n");
      for (k = 0; k < 4; k++)
	{
	  acb_printd(&I_test[k], 10); flint_printf("\n");
	}
      fflush(stdout); */

      if (!r || cov_distinct(I_test, I_acb, 4, weights, prec))
	{
	  flint_printf("FAIL\n");
	  flint_printf("r = %i, I:\n", r);
	  for (k = 0; k < 4; k++)
	    {
	      fmpz_print(&I[k]); flint_printf("\n");
	    }
	  flint_printf("j1, j2:\n");
	  for (k = 0; k < 2; k++)
	    {
	      acb_printd(&j[k], 10); flint_printf("\n");
	    }
	  flint_printf("tau:\n");
	  acb_mat_printd(tau, 10); flint_printf("\n");
	  flint_printf("I_test:\n");
	  for (k = 0; k < 4; k++)
	    {
	      acb_printd(&I_test[k], 10); flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_abort();
	}
      cov_find_rescaling(scal, I_test, I, 4, weights, prec);

      _fmpz_vec_clear(I, 4);
      _acb_vec_clear(j, 2);
      acb_clear(tau1);
      acb_clear(tau2);
      acb_mat_clear(tau);
      _acb_vec_clear(I_acb, 4);
      _acb_vec_clear(I_test, 4);
      acb_clear(scal);
    }

  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

