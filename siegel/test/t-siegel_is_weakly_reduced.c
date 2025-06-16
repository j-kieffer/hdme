#include <stdlib.h>

#include "siegel.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("siegel_is_weakly_reduced....");
  fflush(stdout);
  
  flint_rand_init(state);

  for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
      slong g = 2;
      slong bits = 10 + n_randint(state, 500);
      slong prec = 10 * g * bits;
      int res;

      arb_t tol;
      acb_mat_t z;
      fmpz_mat_t m;

      arb_init(tol);
      acb_mat_init(z, g, g);
      fmpz_mat_init(m, 2*g, 2*g);
      
      arb_set_si(tol, 1);
      arb_mul_2exp_si(tol, tol, -bits); /* tol is larger than 2^(-prec) */

      siegel_halfspace_randtest(z, state, prec);
      siegel_fundamental_domain(z, m, z, tol, prec);
      res = siegel_is_weakly_reduced(z, tol, prec);

      if (!res)
	{ 
	  flint_printf("FAIL (not weakly reduced)\n");
	  flint_printf("res = %wd\n", res);
	  flint_printf("z = "); acb_mat_printd(z, 30); flint_printf("\n");
	  flint_abort();
	}

      arb_clear(tol);
      acb_mat_clear(z);
      fmpz_mat_clear(m);
    }

  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
