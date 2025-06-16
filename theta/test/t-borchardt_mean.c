#include <stdlib.h>
#include "theta.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("borchardt_mean....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
      slong mag_bits = n_randint(state, 10);
      slong prec = 10 + n_pow(2, mag_bits) + n_randint(state, 200);
      
      acb_ptr a;
      arb_t x, y, z;
      acb_t r;

      int res;
      int i;

      a = _acb_vec_init(4);
      arb_init(x);
      arb_init(y);
      arb_init(z);
      acb_init(r);

      /* Check that result coincides with real AGM */
      arb_randtest_precise(x, state, prec, mag_bits);
      arb_randtest_precise(y, state, prec, mag_bits);
      while (arb_contains_zero(x)) arb_randtest_precise(x, state, prec, mag_bits);
      while (arb_contains_zero(y)) arb_randtest_precise(y, state, prec, mag_bits);
      if (!arb_is_positive(x)) arb_neg(x, x);
      if (!arb_is_positive(y)) arb_neg(y, y);
      
      arb_agm(z, x, y, prec);

      acb_set_arb(&a[0], x);
      acb_set_arb(&a[1], x);
      acb_set_arb(&a[2], y);
      acb_set_arb(&a[3], y);
      res = borchardt_mean(r, a, prec);

      if (!res || !arb_overlaps(z, acb_realref(r)) || !arb_contains_zero(acb_imagref(r)))
	{
	  flint_printf("FAIL\n");
	  flint_printf("res = %wd\n", res);
	  flint_printf("x = "); arb_printd(x, 100); flint_printf("\n");
	  flint_printf("y = "); arb_printd(y, 100); flint_printf("\n");
	  flint_printf("z = "); arb_printd(z, 100); flint_printf("\n");
	  for (i = 0; i < 4; i++)
	    {
	      flint_printf("a[%wd] = ", i); acb_printd(&a[i], 30); flint_printf("\n");
	    }
	  flint_printf("r = "); acb_printd(r, 30); flint_printf("\n");
	  flint_abort();
	}

      _acb_vec_clear(a, 4);
      arb_clear(x);
      arb_clear(y);
      arb_clear(z);
      acb_clear(r);
    }
  
  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
