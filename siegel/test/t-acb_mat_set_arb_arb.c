#include <stdlib.h>
#include "siegel.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("acb_mat_set_arb_arb....");
  fflush(stdout);
  
  flint_rand_init(state);
  
  for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
      acb_mat_t z;
      arb_mat_t re, im, re2, im2;
      slong nrows = 1 + n_randint(state, 10);
      slong ncols = 1 + n_randint(state, 10);

      slong prec = 2 + n_randint(state, 200);
      slong mag_bits = 2 + n_randint(state, 20);

      acb_mat_init(z, nrows, ncols);
      arb_mat_init(re, nrows, ncols);
      arb_mat_init(re2, nrows, ncols);
      arb_mat_init(im, nrows, ncols);
      arb_mat_init(im2, nrows, ncols);

      arb_mat_randtest(re, state, prec, mag_bits);
      arb_mat_randtest(im, state, prec, mag_bits);
      acb_mat_set_arb_arb(z, re, im);
      acb_mat_get_real(re2, z);
      acb_mat_get_imag(im2, z);

      if (!arb_mat_overlaps(re, re2) || !arb_mat_overlaps(im, im2))
	{
	  flint_printf("FAIL\n");
	  flint_printf("nrows = %wd\n", nrows);
	  flint_printf("ncols = %wd\n", ncols);
	  flint_printf("z = "); acb_mat_printd(z, 30); flint_printf("\n\n");
	  flint_printf("re = "); arb_mat_printd(re, 30); flint_printf("\n\n");
	  flint_printf("im = "); arb_mat_printd(im, 30); flint_printf("\n\n");
	  flint_printf("re2 = "); arb_mat_printd(re2, 30); flint_printf("\n\n");
	  flint_printf("im2 = "); arb_mat_printd(im2, 30); flint_printf("\n\n");
	  flint_abort();
	}

      acb_mat_clear(z);
      arb_mat_clear(re);
      arb_mat_clear(re2);
      arb_mat_clear(im);
      arb_mat_clear(im2);
    }
  
  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
