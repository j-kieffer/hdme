#include <stdlib.h>

#include "siegel.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("siegel_transform....");
  fflush(stdout);
  
  flint_rand_init(state);
  
  /* Check associativity */
  for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
      fmpz_mat_t m, n, mn;
      acb_mat_t z1, z2, z3;
      slong g, prec;
      int valid_transform;

      g = 1 + n_randint(state, 10);
      fmpz_mat_init(m, 2*g, 2*g);
      fmpz_mat_init(n, 2*g, 2*g);
      fmpz_mat_init(mn, 2*g, 2*g);
      acb_mat_init(z1, g, g);
      acb_mat_init(z2, g, g);
      acb_mat_init(z3, g, g);

      fmpz_mat_randtest_symplectic(m, state, n_randint(state, 20));
      fmpz_mat_randtest_symplectic(n, state, n_randint(state, 20));
      fmpz_mat_mul(mn, m, n);

      prec = 200 + n_randint(state, 200);

      siegel_halfspace_randtest(z1, state, prec);
      valid_transform = siegel_transform(z2, mn, z1, prec);

      if (!valid_transform)
	{
	  flint_printf("FAIL (insufficient precision to invert)\n");
	  flint_printf("g = %wd\n", g);
	  flint_printf("mn = "); fmpz_mat_print(mn); flint_printf("\n\n");
	  flint_printf("z1 = "); acb_mat_printd(z1, 30); flint_printf("\n\n");
	  flint_abort();
	}

      valid_transform = siegel_transform(z3, n, z1, prec);
      
      if (!valid_transform)
	{
	  flint_printf("FAIL (insufficient precision to invert)\n");
	  flint_printf("g = %wd\n", g);
	  flint_printf("n = "); fmpz_mat_print(n); flint_printf("\n\n");
	  flint_printf("z1 = "); acb_mat_printd(z1, 30); flint_printf("\n\n");
	  flint_abort();
	}

      valid_transform = siegel_transform(z3, m, z3, prec);

      if (!valid_transform)
	{
	  flint_printf("FAIL (insufficient precision to invert)\n");
	  flint_printf("g = %wd\n", g);
	  flint_printf("m = "); fmpz_mat_print(m); flint_printf("\n\n");
	  flint_printf("z3 = "); acb_mat_printd(z3, 30); flint_printf("\n\n");
	  flint_abort();
	}

      if (!acb_mat_overlaps(z2, z3))
	{
	  
	  flint_printf("FAIL\n");
	  flint_printf("g = %wd\n", g);
	  flint_printf("m = "); fmpz_mat_print(m); flint_printf("\n\n");
	  flint_printf("n = "); fmpz_mat_print(n); flint_printf("\n\n");
	  flint_printf("z1 = "); acb_mat_printd(z1, 30); flint_printf("\n\n");
	  flint_printf("z2 = "); acb_mat_printd(z2, 30); flint_printf("\n\n");
	  flint_printf("z3 = "); acb_mat_printd(z3, 30); flint_printf("\n\n");
	  flint_abort();
	}

      fmpz_mat_clear(m);
      fmpz_mat_clear(n);
      fmpz_mat_clear(mn);
      acb_mat_clear(z1);
      acb_mat_clear(z2);
      acb_mat_clear(z3);
    }
  
  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
