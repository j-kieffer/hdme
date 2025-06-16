#include <stdlib.h>

#include "igusa.h"

int main()
{
  
  slong iter;
  flint_rand_t state;
  
  flint_printf("mestre_fmpz....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
      fmpz_poly_t crv0;
      fmpz_t c;
      acb_t scal;
      acb_poly_t crv;      
      fmpz* I;
      acb_ptr I_acb;
      acb_ptr I_test;
      fmpz* IC;
      slong mag_bits = 10;
      slong prec = 1000;
      slong weights[4] = IGUSA_WEIGHTS;
      slong k;
      int valid;
      int print = 0;

      fmpz_poly_init(crv0);
      fmpz_init(c);
      acb_init(scal);
      acb_poly_init(crv);
      I = _fmpz_vec_init(4);
      I_acb = _acb_vec_init(4);
      I_test = _acb_vec_init(4);
      IC = _fmpz_vec_init(4);
      
      /* Generate curves of different Bolza types */
      fmpz_poly_zero(crv0);
      for (k = 0; k < 7; k++)
	{
	  fmpz_randtest_not_zero(c, state, mag_bits);
	  fmpz_poly_set_coeff_fmpz(crv0, k, c);
	}
      
      igusa_from_curve_fmpz(I, crv0);
      for (k = 0; k < 4; k++) acb_set_fmpz(&I_acb[k], &I[k]);
      valid = igusa_is_g2_curve_fmpz(I);
      if (valid)
	{
	  igusa_IC_fmpz(IC, I);
	  mestre_fmpz(crv, IC, prec);
	  igusa_from_curve(I_test, crv, prec);
	}
      if (valid && print)
	{
	  flint_printf("I_test (generic):\n");
	  for (k = 0; k < 4; k++)
	    {
	      acb_printd(&I_test[k], 10); flint_printf("\n");
	    }
	  fflush(stdout);
	}
      if (valid && cov_distinct(I_acb, I_test, 4, weights, prec))
	{
	  flint_printf("FAIL (generic)\n");
	  flint_printf("crv0:\n"); fmpz_poly_print(crv0); flint_printf("\n");
	  flint_printf("I:\n");
	  for (k = 0; k < 4; k++)
	    {
	      fmpz_print(&I[k]); flint_printf("\n");	      
	    }
	  flint_printf("I_test:\n");
	  for (k = 0; k < 4; k++)
	    {
	      acb_printd(&I_test[k], 10); flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_abort();
	}
      if (valid) cov_find_rescaling(scal, I_test, I, 4, weights, prec);
      
      fmpz_poly_zero(crv0);
      fmpz_randtest_not_zero(c, state, mag_bits);
      fmpz_poly_set_coeff_fmpz(crv0, 6, c);
      fmpz_poly_set_coeff_fmpz(crv0, 0, c);
      
      fmpz_randtest_not_zero(c, state, mag_bits);
      fmpz_poly_set_coeff_fmpz(crv0, 4, c);
      fmpz_randtest_not_zero(c, state, mag_bits);
      fmpz_poly_set_coeff_fmpz(crv0, 2, c);
      
      igusa_from_curve_fmpz(I, crv0);
      for (k = 0; k < 4; k++) acb_set_fmpz(&I_acb[k], &I[k]);
      valid = igusa_is_g2_curve_fmpz(I);
      if (valid)
	{
	  igusa_IC_fmpz(IC, I);
	  mestre_fmpz(crv, IC, prec);
	  igusa_from_curve(I_test, crv, prec);
	}
      if (valid && print)
	{
	  flint_printf("I_test (I):\n");
	  for (k = 0; k < 4; k++)
	    {
	      acb_printd(&I_test[k], 10); flint_printf("\n");
	    }
	  fflush(stdout);
	}
      if (valid && cov_distinct(I_acb, I_test, 4, weights, prec))
	{
	  flint_printf("FAIL (type I)\n");
	  fflush(stdout);
	  flint_abort();
	}
      if (valid) cov_find_rescaling(scal, I_test, I, 4, weights, prec);
      
      fmpz_poly_zero(crv0);
      fmpz_randtest_not_zero(c, state, mag_bits);
      fmpz_poly_set_coeff_fmpz(crv0, 6, c);
      fmpz_poly_set_coeff_fmpz(crv0, 1, c);
      igusa_from_curve_fmpz(I, crv0);
      for (k = 0; k < 4; k++) acb_set_fmpz(&I_acb[k], &I[k]);
      valid = igusa_is_g2_curve_fmpz(I);
      if (valid)
	{
	  igusa_IC_fmpz(IC, I);
	  mestre_fmpz(crv, IC, prec);
	  igusa_from_curve(I_test, crv, prec);
	}
      if (valid && print)
	{
	  flint_printf("I_test (II):\n");
	  for (k = 0; k < 4; k++)
	    {
	      acb_printd(&I_test[k], 10); flint_printf("\n");
	    }
	  fflush(stdout);
	}
      if (valid && cov_distinct(I_acb, I_test, 4, weights, prec))
	{
	  flint_printf("FAIL (type II)\n");
	  fflush(stdout);
	  flint_abort();
	}
      if (valid) cov_find_rescaling(scal, I_test, I, 4, weights, prec);
      
      fmpz_poly_zero(crv0);
      fmpz_randtest_not_zero(c, state, mag_bits);
      fmpz_poly_set_coeff_fmpz(crv0, 5, c);
      fmpz_poly_set_coeff_fmpz(crv0, 1, c);
      
      fmpz_randtest_not_zero(c, state, mag_bits);
      fmpz_poly_set_coeff_fmpz(crv0, 3, c);
      
      igusa_from_curve_fmpz(I, crv0);
      for (k = 0; k < 4; k++) acb_set_fmpz(&I_acb[k], &I[k]);
      valid = igusa_is_g2_curve_fmpz(I);
      if (valid)
	{
	  igusa_IC_fmpz(IC, I);
	  mestre_fmpz(crv, IC, prec);
	  igusa_from_curve(I_test, crv, prec);
	}
      if (valid && print)
	{
	  flint_printf("I_test (III):\n");
	  for (k = 0; k < 4; k++)
	    {
	      acb_printd(&I_test[k], 10); flint_printf("\n");
	    }
	  fflush(stdout);
	}
      if (valid && cov_distinct(I_acb, I_test, 4, weights, prec))
	{
	  flint_printf("FAIL (type III)\n");
	  fflush(stdout);
	  flint_abort();
	}
      if (valid) cov_find_rescaling(scal, I_test, I, 4, weights, prec);
      
      fmpz_poly_zero(crv0);
      fmpz_randtest_not_zero(c, state, mag_bits);
      fmpz_poly_set_coeff_fmpz(crv0, 6, c);
      fmpz_poly_set_coeff_fmpz(crv0, 0, c);
      
      fmpz_randtest_not_zero(c, state, mag_bits);
      fmpz_poly_set_coeff_fmpz(crv0, 3, c);
      
      igusa_from_curve_fmpz(I, crv0);
      for (k = 0; k < 4; k++) acb_set_fmpz(&I_acb[k], &I[k]);
      valid = igusa_is_g2_curve_fmpz(I);
      if (valid)
	{
	  igusa_IC_fmpz(IC, I);
	  mestre_fmpz(crv, IC, prec);
	  igusa_from_curve(I_test, crv, prec);
	}
      if (valid && print)
	{
	  flint_printf("I_test (IV):\n");
	  for (k = 0; k < 4; k++)
	    {
	      acb_printd(&I_test[k], 10); flint_printf("\n");
	    }
	  fflush(stdout);
	}
      if (valid && cov_distinct(I_acb, I_test, 4, weights, prec))
	{
	  flint_printf("FAIL (type IV)\n");
	  fflush(stdout);
	  flint_abort();
	}
      if (valid) cov_find_rescaling(scal, I_test, I, 4, weights, prec);
      
      fmpz_poly_zero(crv0);
      fmpz_randtest_not_zero(c, state, mag_bits);
      fmpz_poly_set_coeff_fmpz(crv0, 6, c);
      fmpz_poly_set_coeff_fmpz(crv0, 0, c);
      
      igusa_from_curve_fmpz(I, crv0);
      for (k = 0; k < 4; k++) acb_set_fmpz(&I_acb[k], &I[k]);
      valid = igusa_is_g2_curve_fmpz(I);
      if (valid)
	{      
	  igusa_IC_fmpz(IC, I);
	  mestre_fmpz(crv, IC, prec);
	  igusa_from_curve(I_test, crv, prec);
	}
      if (valid && print)
	{
	  flint_printf("I_test (V):\n");
	  for (k = 0; k < 4; k++)
	    {
	      acb_printd(&I_test[k], 10); flint_printf("\n");
	    }
	  fflush(stdout);
	}
      if (valid && cov_distinct(I_acb, I_test, 4, weights, prec))
	{
	  flint_printf("FAIL (type V)\n");
	  fflush(stdout);
	  flint_abort();
	}
      if (valid) cov_find_rescaling(scal, I_test, I, 4, weights, prec);

      fmpz_poly_zero(crv0);
      fmpz_randtest_not_zero(c, state, mag_bits);
      fmpz_poly_set_coeff_fmpz(crv0, 5, c);
      fmpz_poly_set_coeff_fmpz(crv0, 1, c);
      
      igusa_from_curve_fmpz(I, crv0);
      for (k = 0; k < 4; k++) acb_set_fmpz(&I_acb[k], &I[k]);
      valid = igusa_is_g2_curve_fmpz(I);
      if (valid)
	{
	  igusa_IC_fmpz(IC, I);
	  mestre_fmpz(crv, IC, prec);
	  igusa_from_curve(I_test, crv, prec);
	}
      if (valid && print)
	{
	  flint_printf("I_test (VI):\n");
	  for (k = 0; k < 4; k++)
	    {
	      acb_printd(&I_test[k], 10); flint_printf("\n");
	    }
	  fflush(stdout);
	}
      if (valid && cov_distinct(I_acb, I_test, 4, weights, prec))
	{
	  flint_printf("FAIL (type VI)\n");
	  fflush(stdout);
	  flint_abort();
	}      
      if (valid) cov_find_rescaling(scal, I_test, I, 4, weights, prec);
      
      fmpz_poly_clear(crv0);
      fmpz_clear(c);
      acb_clear(scal);
      acb_poly_clear(crv);
      _fmpz_vec_clear(I, 4);
      _acb_vec_clear(I_acb, 4);
      _acb_vec_clear(I_test, 4);
      _fmpz_vec_clear(IC, 4);
    }
     
  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
