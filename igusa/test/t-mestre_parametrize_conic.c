#include <stdlib.h>

#include "igusa.h"

int main()
{
  
  slong iter;
  flint_rand_t state;
  
  flint_printf("mestre_parametrize_conic....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
      slong prec = 50 + n_randint(state, 1000);
      acb_ptr conic;
      acb_ptr pt;
      acb_poly_t x1, x2, x3;
      acb_t param;
      acb_ptr test_pt;
      int res = 1;
      int k;
      slong mag_bits = 2;
      
      conic = _acb_vec_init(6);
      pt = _acb_vec_init(3);
      acb_poly_init(x1);
      acb_poly_init(x2);
      acb_poly_init(x3);
      acb_init(param);
      test_pt = _acb_vec_init(3);

      mestre_conic_randtest(conic, state, prec);
      mestre_point_on_conic(pt, conic, prec);
      mestre_parametrize_conic(x1, x2, x3, pt, conic, prec);

      /* A random point on the parametrized P^1 should lie on the conic */
      acb_randtest_precise(param, state, prec, mag_bits);
      acb_poly_evaluate(&test_pt[0], x1, param, prec);
      acb_poly_evaluate(&test_pt[1], x2, param, prec);
      acb_poly_evaluate(&test_pt[2], x3, param, prec);
      res = !mestre_point_is_outside_conic(test_pt, conic, prec);

      if (!res)
	{
	  flint_printf("FAIL (Point outside conic)\n");
	  flint_printf("Conic: \n");
	  for (k = 0; k < 6; k++)
	    {
	      flint_printf("    "); acb_printd(&conic[k], 30); flint_printf("\n");
	    }
	  flint_printf("Parametrization: \n");
	  flint_printf("    "); acb_poly_printd(x1, 30); flint_printf("\n");
	  flint_printf("    "); acb_poly_printd(x2, 30); flint_printf("\n");
	  flint_printf("    "); acb_poly_printd(x3, 30); flint_printf("\n");
	  flint_printf("Parameter: "); acb_printd(param, 30); flint_printf("\n");
	  flint_printf("Point: \n");
	  for (k = 0; k < 3; k++)
	    {
	      flint_printf("    "); acb_printd(&test_pt[k], 30); flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_abort();
	}
      
      _acb_vec_clear(conic, 6);
      _acb_vec_clear(pt, 3);
      acb_poly_clear(x1);
      acb_poly_clear(x2);
      acb_poly_clear(x3);
      acb_clear(param);
      _acb_vec_clear(test_pt, 3);
    }

  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
