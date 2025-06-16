#include <stdlib.h>
#include "polynomials.h"

int main(void)
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("acb_poly_product_tree_1....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
      acb_poly_t P, P_test;
      acb_ptr xi, yi;
      slong d = n_randint(state, 100);
      slong prec = 50 + n_randint(state, 200);
      slong mag_bits = 1 + n_randint(state, 4);
      slong k;

      acb_poly_init(P);
      acb_poly_init(P_test);
      xi = _acb_vec_init(d);
      yi = _acb_vec_init(d);

      for (k = 0; k < d; k++) acb_randtest_precise(&yi[k], state, prec, mag_bits);
      for (k = 0; k < d; k++) acb_one(&xi[k]);

      acb_poly_product_roots(P_test, yi, d, prec);
      for (k = 0; k < d; k++) acb_neg(&yi[k], &yi[k]);
      acb_poly_product_tree_1(P, xi, yi, d, prec);

      if (!acb_poly_overlaps(P, P_test))
	{
	  flint_printf("FAIL (overlap)\n");
	  acb_poly_printd(P_test, 10);
	  acb_poly_printd(P, 10);
	  /* for (k = 0; k < d; k++)
	    {
	      acb_printd(&yi[k], 30); flint_printf("\n");
	      } */
	  fflush(stdout);
	  flint_abort();
	}

      
      acb_poly_clear(P);
      acb_poly_clear(P_test);
      _acb_vec_clear(xi, d);
      _acb_vec_clear(yi, d);
    }
  
  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
      
