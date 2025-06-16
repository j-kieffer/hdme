#include <stdlib.h>

#include "igusa.h"

int main()
{
  
  slong iter;
  flint_rand_t state;
  
  flint_printf("cardona_conic_cubic....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
      acb_ptr ABCD;
      acb_ptr Aij;
      acb_ptr aijk;
      acb_t scal;
      acb_ptr test;
      slong wts[4] = IC_WEIGHTS;
      slong conic_wts[4] = {6,8,10,16};
      slong cubic_wts[6] = {10,12,14,20,16,22};
      slong prec = 1000;
      slong mag_bits = 5;
      slong k;

      ABCD = _acb_vec_init(4);
      Aij = _acb_vec_init(4);
      aijk = _acb_vec_init(6);
      acb_init(scal);
      test = _acb_vec_init(6);

      for (k = 0; k < 4; k++) acb_randtest_precise(&ABCD[k], state, prec, mag_bits);
      acb_randtest_precise(scal, state, prec, mag_bits);

      cardona_conic(Aij, ABCD, prec);
      cov_rescale(test, ABCD, scal, 4, wts, prec);
      cardona_conic(test, test, prec);
      if (cov_distinct(Aij, test, 4, conic_wts, prec))
	{
	  flint_printf("FAIL (conic)\n");
	  fflush(stdout);
	  flint_abort();
	}
      
      cardona_cubic(aijk, ABCD, prec);
      cov_rescale(test, ABCD, scal, 4, wts, prec);
      cardona_cubic(test, test, prec);
      if (cov_distinct(aijk, test, 6, cubic_wts, prec))
	{
	  flint_printf("FAIL (cubic)\n");
	  /* There is a term of degree 18 in a233 */
	  fflush(stdout);
	  flint_abort();
	}

      _acb_vec_clear(ABCD, 4);
      _acb_vec_clear(Aij, 4);
      _acb_vec_clear(aijk, 6);
      acb_clear(scal);
      _acb_vec_clear(test, 6);
    }
     
  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

