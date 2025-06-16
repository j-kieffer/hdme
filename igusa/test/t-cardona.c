#include <stdlib.h>
#include "igusa.h"

int main()
{
  
  slong iter;
  flint_rand_t state;
  
  flint_printf("cardona....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
      acb_poly_t crv;
      acb_t c;
      acb_ptr I, I_test;
      acb_ptr IC;
      slong mag_bits = 5;
      slong prec = 1000;
      slong weights[4] = IGUSA_WEIGHTS;
      slong k;

      acb_poly_init(crv);
      I = _acb_vec_init(4);
      I_test = _acb_vec_init(4);
      IC = _acb_vec_init(4);
      acb_init(c);

      for (k = 0; k < 4; k++)
	{
	  acb_randtest_precise(c, state, prec, mag_bits);
	  acb_poly_set_coeff_acb(crv, 2*k, c);
	}

      /*
      acb_poly_set_coeff_si(crv, 6, 1);
      acb_poly_set_coeff_si(crv, 4, 1);
      acb_poly_set_coeff_si(crv, 2, 2);
      acb_poly_set_coeff_si(crv, 0, 1);*/
      
      igusa_from_curve(I, crv, prec);
      igusa_IC(IC, I, prec);
      igusa_R2_from_IC(c, IC, prec);

      if (!acb_contains_zero(c))
	{
	  flint_printf("FAIL (R2)\n");
	  fflush(stdout);
	  flint_abort();
	}

      cardona(crv, IC, prec);
      igusa_from_curve(I_test, crv, prec);
      igusa_ABCD_from_IC(IC, IC, prec);
	  
      if (cov_distinct(I_test, I, 4, weights, prec))
	{
	  flint_printf("FAIL\n");
	  flint_printf("I:\n");
	  for (k = 0; k < 4; k++)
	    {
	      acb_printd(&I[k], 10); flint_printf("\n");
	      acb_printd(&I_test[k], 10); flint_printf("\n");
	    }
	  flint_printf("ABCD:\n");
	  for (k = 0; k < 4; k++)
	    {
	      acb_printd(&IC[k], 10); flint_printf("\n");
	    }
	  flint_printf("Curve:\n");
	  acb_poly_printd(crv, 10); flint_printf("\n");

	  flint_printf("Normalized covariants:\n");
	  cov_normalize(I, I, 4, weights, prec);
	  for (k = 0; k < 4; k++)
	    {
	      acb_printd(&I[k], 10); flint_printf("\n");
	    }
	  cov_normalize(I_test, I_test, 4, weights, prec);
	  for (k = 0; k < 4; k++)
	    {
	      acb_printd(&I_test[k], 10); flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_abort();	  
	}
      
      acb_poly_clear(crv);
      _acb_vec_clear(I, 4);
      _acb_vec_clear(I_test, 4);
      _acb_vec_clear(IC, 4);
      acb_clear(c);
    }
     
  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
