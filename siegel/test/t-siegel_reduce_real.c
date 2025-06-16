#include <stdlib.h>

#include "siegel.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("siegel_reduce_real....");
  fflush(stdout);
  
  flint_rand_init(state);
  
  for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
      fmpz_mat_t u;
      acb_mat_t z, w1, w2;
      arb_t tol;
      int res;
      slong g = 1 + n_randint(state, 10);
      slong bits = 1 + n_randint(state, 100);
      slong prec = 5 * bits;

      acb_mat_init(z, g, g);
      acb_mat_init(w1, g, g);
      acb_mat_init(w2, g, g);
      fmpz_mat_init(u, 2*g, 2*g);
      arb_init(tol);
      
      arb_set_si(tol, 1);
      arb_mul_2exp_si(tol, tol, -bits); /* tol is larger than 2^(-prec) */
      
      siegel_halfspace_randtest(z, state, prec);
      res = siegel_reduce_real(w1, u, z, tol, prec);
      siegel_transform(w2, u, z, prec);
      
      if (!res || !siegel_is_real_reduced(w1, tol, prec) || !acb_mat_overlaps(w1, w2))
	{
	  flint_printf("FAIL\n");
	  flint_printf("g = %wd\n", g);
	  flint_printf("z = "); acb_mat_printd(z, 30); flint_printf("\n\n");
	  flint_printf("u = "); fmpz_mat_print(u); flint_printf("\n\n");
	  flint_printf("w1 = "); acb_mat_printd(w1, 30); flint_printf("\n\n");
	  flint_printf("w2 = "); acb_mat_printd(w2, 30); flint_printf("\n\n");
	  flint_abort();
	}

      acb_mat_clear(z);
      acb_mat_clear(w1);
      acb_mat_clear(w2);
      fmpz_mat_clear(u);
      arb_clear(tol);
    }
  
  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

