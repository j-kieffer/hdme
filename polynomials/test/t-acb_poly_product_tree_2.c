#include <stdlib.h>
#include "polynomials.h"
#include <flint/fmpq.h>

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("acb_poly_product_tree_2....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
      acb_poly_t Q, Q_test;
      acb_ptr xi, yi, zi;
      slong d = n_randint(state, 100);
      slong prec = 50 + n_randint(state, 200);
      slong mag_bits = 1 + n_randint(state, 4);
      slong k;

      acb_poly_init(Q);
      acb_poly_init(Q_test);
      xi = _acb_vec_init(d);
      yi = _acb_vec_init(d);
      zi = _acb_vec_init(d);

      for (k = 0; k < d; k++) acb_randtest_precise(&yi[k], state, prec, mag_bits);
      for (k = 0; k < d; k++) acb_randtest_precise(&xi[k], state, prec, mag_bits);
      for (k = 0; k < d; k++) acb_set(&zi[k], &xi[k]);

      acb_poly_product_tree_1(Q_test, xi, yi, d, prec);
      acb_poly_derivative(Q_test, Q_test, prec);
      acb_poly_product_tree_2(Q, xi, yi, zi, d, prec);
      
      if (!acb_poly_overlaps(Q, Q_test))
	{
	  flint_printf("FAIL (overlap)\n");
	  acb_poly_printd(Q_test, 10);
	  acb_poly_printd(Q, 10);
	  /* for (k = 0; k < d; k++)
	    {
	      acb_printd(&xi[k], 30); flint_printf("\n");
	      acb_printd(&yi[k], 30); flint_printf("\n\n");
	      } */
	  fflush(stdout);
	  flint_abort();
	}
      
      acb_poly_clear(Q);
      acb_poly_clear(Q_test);
      _acb_vec_clear(xi, d);
      _acb_vec_clear(yi, d);
      _acb_vec_clear(zi, d);
    }
  
  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
      
  
