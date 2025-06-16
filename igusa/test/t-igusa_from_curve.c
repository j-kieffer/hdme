#include <stdlib.h>

#include "igusa.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("igusa_from_curve....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
      slong prec = 50 + n_randint(state, 1000);
      acb_ptr roots;
      acb_ptr I;
      acb_ptr IC;
      acb_t R2;
      acb_poly_t crv;
      slong k;

      /* Generate crv with multiple roots; check that covariants are zero */
      roots = _acb_vec_init(6);
      I = _acb_vec_init(4);
      IC = _acb_vec_init(4);
      acb_poly_init(crv);
      acb_init(R2);

      acb_randtest_precise(&roots[0], state, prec, 1);
      acb_set(&roots[1], &roots[0]);
      acb_set(&roots[2], &roots[0]);
      acb_set(&roots[3], &roots[0]);
      acb_randtest_precise(&roots[4], state, prec, 1);
      acb_randtest_precise(&roots[5], state, prec, 1);

      acb_poly_product_roots(crv, roots, 6, prec);

      igusa_from_curve(I, crv, prec);

      if (!acb_contains_zero(&I[0])
	  || !acb_contains_zero(&I[1])
	  || !acb_contains_zero(&I[2])
	  || !acb_contains_zero(&I[3]))
	{
	  flint_printf("FAIL\n");
	  flint_printf("Curve: "), acb_poly_printd(crv, 30), flint_printf("\n");
	  for (k = 0; k < 4; k++)
	    {	      
	      acb_printd(&I[k], 30); flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_abort();
	}
      
      /* Generate crv with less multiple roots; check that I10, R2 is zero */
      acb_randtest_precise(&roots[0], state, prec, 1);
      acb_set(&roots[1], &roots[0]);
      acb_randtest_precise(&roots[2], state, prec, 1);
      acb_set(&roots[3], &roots[2]);
      acb_randtest_precise(&roots[4], state, prec, 1);
      acb_randtest_precise(&roots[5], state, prec, 1);

      acb_poly_product_roots(crv, roots, 6, prec);

      igusa_from_curve(I, crv, prec);
      igusa_IC(IC, I, prec);
      igusa_R2_from_IC(R2, IC, prec);
      
      if (!acb_contains_zero(&IC[3])
	  || !acb_contains_zero(R2))
	{
	  flint_printf("FAIL\n");
	  flint_printf("Curve: "), acb_poly_printd(crv, 30), flint_printf("\n");
	  for (k = 0; k < 4; k++)
	    {
	      acb_printd(&I[k], 30); flint_printf("\n");
	    }
	  flint_printf("R2 = "); acb_printd(R2, 30); flint_printf("\n");
	  fflush(stdout);
	  flint_abort();
	}
      
      _acb_vec_clear(roots, 6);
      _acb_vec_clear(I, 4);
      _acb_vec_clear(IC, 4);
      acb_poly_clear(crv);
      acb_clear(R2);
    }
  
  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}


      
