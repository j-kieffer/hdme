#include <stdlib.h>

#include "siegel.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("siegel_halfspace_randtest....");
  fflush(stdout);
  
  flint_rand_init(state);
  
  for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
      acb_mat_t z;
      arb_mat_t re, im, ret, imt;
      arb_t det, tr;
      slong g, prec;
      
      g = 1 + n_randint(state, 10);

      acb_mat_init(z, g, g);
      arb_mat_init(re, g, g);
      arb_mat_init(ret, g, g);
      arb_mat_init(im, g, g);
      arb_mat_init(imt, g, g);
      arb_init(det);
      arb_init(tr);

      prec = 20 + n_randint(state, 200);

      /* Check that re and im are symmetric, and det, tr of im are >0 */
      siegel_halfspace_randtest(z, state, prec);
      acb_mat_get_real(re, z);
      acb_mat_get_imag(im, z);
      arb_mat_transpose(ret, re);
      arb_mat_transpose(imt, im);

      if (!arb_mat_overlaps(ret, re) || !arb_mat_overlaps(imt, im))
	{
	  flint_printf("FAIL\n");
	  flint_printf("g = %wd\n", g);
	  flint_printf("z = "); acb_mat_printd(z, 30); flint_printf("\n\n");
	  flint_printf("re = "); arb_mat_printd(re, 30); flint_printf("\n\n");
	  flint_printf("im = "); arb_mat_printd(im, 30); flint_printf("\n\n");
	  flint_abort();
	}

      arb_mat_det(det, im, prec);
      arb_mat_trace(tr, im, prec);
      
      if (!arb_is_positive(det) || !arb_is_positive(tr))
	{
	  flint_printf("FAIL\n");
	  flint_printf("g = %wd\n", g);
	  flint_printf("z = "); acb_mat_printd(z, 30); flint_printf("\n\n");
	  flint_printf("im = "); arb_mat_printd(im, 30); flint_printf("\n\n");
	  flint_printf("det = "); arb_printd(det, 30); flint_printf("\n\n");
	  flint_printf("tr = "); arb_printd(tr, 30); flint_printf("\n\n");
	  flint_abort();
	}

      acb_mat_clear(z);
      arb_mat_clear(re);
      arb_mat_clear(ret);
      arb_mat_clear(im);
      arb_mat_clear(imt);
      arb_clear(det);
      arb_clear(tr);
    }
  
  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
