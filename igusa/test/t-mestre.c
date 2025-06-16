#include <stdlib.h>

#include "igusa.h"

int main()
{
  
  slong iter;
  flint_rand_t state;
  
  flint_printf("mestre....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
      slong weights[4] = IC_WEIGHTS;
      slong prec = 500 + n_randint(state, 1000);
      acb_poly_t crv_test, crv;
      acb_ptr IC_test, I, IC;
      slong k;

      acb_poly_init(crv_test);
      acb_poly_init(crv);
      IC_test =  _acb_vec_init(4);
      IC = _acb_vec_init(4);
      I =  _acb_vec_init(4);

      igusa_generic_randtest(crv_test, IC_test, state, prec);
      if (acb_contains_zero(&IC_test[3]))
	{
	  flint_printf("FAIL (chi10 = 0)\n");
	  flint_printf("Curve: "); acb_poly_printd(crv_test, 30); flint_printf("\n");
	  fflush(stdout);
	  flint_abort();
	}
      
      mestre(crv, IC_test, prec);
      igusa_from_curve(I, crv, prec);      
      if (acb_contains_zero(igusa_chi10(I)))
	{
	  flint_printf("FAIL (chi10 = 0)\n");
	  flint_printf("Curve: "); acb_poly_printd(crv, 30); flint_printf("\n");
	  fflush(stdout);
	  flint_abort();
	}
      igusa_IC(IC, I, prec);
      
      if (cov_distinct(IC, IC_test, 4, weights, prec))
	{
	  flint_printf("FAIL (overlap)\n");
	  flint_printf("crv_test: "); acb_poly_printd(crv_test, 30); flint_printf("\n");
	  flint_printf("crv: "); acb_poly_printd(crv, 30); flint_printf("\n");
	  for (k = 0; k < 4; k++)
	    {
	      flint_printf("IC_test[%d]: ", k); acb_printd(&IC_test[k], 30); flint_printf("\n");
	      flint_printf("IC[%d]: ", k); acb_printd(&I[k], 30); flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_abort();
	}

      acb_poly_clear(crv_test);
      acb_poly_clear(crv);
      _acb_vec_clear(IC_test, 4);
      _acb_vec_clear(IC, 4);
      _acb_vec_clear(I, 4);
    }
  
  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
