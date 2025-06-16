#include <stdlib.h>
#include "theta.h"


int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("borchardt_mean_m0_M0....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
      slong mag_bits = n_randint(state, 5);
      slong prec = 32 + n_randint(state, 200);
      
      acb_ptr a;
      acb_ptr b;
      acb_t scal;
      arb_t m0;
      arb_t M0;
      arb_t abs;
      int i;

      a = _acb_vec_init(4);
      b = _acb_vec_init(4);
      acb_init(scal);
      arb_init(m0);
      arb_init(M0);
      arb_init(abs);

      /* Generate four complex numbers with positive real part */
      for (i = 0; i < 4; i++)
	{
	  acb_randtest_precise(&b[i], state, prec, mag_bits);
	  while (arb_contains_zero(acb_realref(&b[i])))
	    {
	      acb_randtest_precise(&b[i], state, prec, mag_bits);
	    }
	  if (arb_is_negative(acb_realref(&b[i])))
	    {
	      acb_neg(&b[i], &b[i]);
	    } 
	}
      /* Scale so that first element becomes 1 */
      acb_inv(scal, &b[0], prec);
      _acb_vec_scalar_mul(a, b, 4, scal, prec);

      /* Compute m0, M0 */
      borchardt_mean_m0(m0, a, prec);
      borchardt_mean_M0(M0, a, prec);

      if (!arb_is_positive(m0) || !arb_is_positive(M0))
	{
	  flint_printf("FAIL (not positive)\n");
	  for (i = 0; i < 4; i++)
	    {
	      flint_printf("a[%wd] = ", i); acb_printd(&a[i], 30); flint_printf("\n");
	      flint_printf("b[%wd] = ", i); acb_printd(&b[i], 30); flint_printf("\n");
	    }
	  flint_printf("m0 = "); arb_printd(m0, 30); flint_printf("\n");
	  flint_printf("M0 = "); arb_printd(M0, 30); flint_printf("\n");
	  flint_abort();
	}

      for (i = 0; i < 4; i++)
	{
	  acb_abs(abs, &a[i], prec);
	  if (arb_lt(abs, m0) || arb_lt(M0, abs))
	    {
	      flint_printf("FAIL (comparison)\n");
	      flint_printf("a[%wd] = ", i); acb_printd(&a[i], 30); flint_printf("\n");
	      flint_printf("abs = "); arb_printd(abs, 30); flint_printf("\n");
	      flint_printf("m0 = "); arb_printd(m0, 30); flint_printf("\n");
	      flint_printf("M0 = "); arb_printd(M0, 30); flint_printf("\n");
	      flint_abort();
	    }
	}
      
      _acb_vec_clear(a, 4);
      _acb_vec_clear(b, 4);
      acb_clear(scal);
      arb_clear(m0);
      arb_clear(M0);
      arb_clear(abs);
    }
  
  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
