#include <stdlib.h>

#include "igusa.h"

int main()
{
  
  slong iter;
  flint_rand_t state;
  
  flint_printf("mestre_subst_in_cubic....");
  fflush(stdout);
  
  flint_rand_init(state);
  
  for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
      slong prec = 50 + n_randint(state, 1000);
      
      acb_ptr cubic;
      acb_poly_t x1, x2, x3, subst;
      acb_ptr pt;
      acb_t param, ev, ev_test;
      int k;
      slong mag_bits = 1 + n_randint(state, 20);
      
      cubic = _acb_vec_init(10);
      pt = _acb_vec_init(3);
      acb_poly_init(x1);
      acb_poly_init(x2);
      acb_poly_init(x3);
      acb_poly_init(subst);
      acb_init(param);
      acb_init(ev);
      acb_init(ev_test);
      
      acb_poly_randtest(x1, state, 3, prec, mag_bits);
      acb_poly_randtest(x2, state, 3, prec, mag_bits);
      acb_poly_randtest(x3, state, 3, prec, mag_bits);
      for (k = 0; k < 10; k++)
	{
	  acb_randtest_precise(&cubic[k], state, prec, mag_bits);
	}
      acb_randtest_precise(param, state, prec, mag_bits);
      
      acb_poly_evaluate(&pt[0], x1, param, prec);
      acb_poly_evaluate(&pt[1], x2, param, prec);
      acb_poly_evaluate(&pt[2], x3, param, prec);
      mestre_eval_cubic(ev, pt, cubic, prec);
      
      mestre_subst_in_cubic(subst, x1, x2, x3, cubic, prec);
      acb_poly_evaluate(ev_test, subst, param, prec);
      
      if (!acb_overlaps(ev, ev_test))
	{
	  flint_printf("FAIL (overlap)\n");
	  fflush(stdout);
	  flint_abort();
	}
      
      _acb_vec_clear(cubic, 10);
      _acb_vec_clear(pt, 3);
      acb_poly_clear(x1);
      acb_poly_clear(x2);
      acb_poly_clear(x3);
      acb_poly_clear(subst);
      acb_clear(param);
      acb_clear(ev);
      acb_clear(ev_test);
    }

  
  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

