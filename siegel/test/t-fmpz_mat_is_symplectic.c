#include <stdlib.h>

#include "siegel.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("fmpz_mat_is_symplectic....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
      fmpz_mat_t m;
      slong g;

      g = 1 + n_randint(state, 10);
      fmpz_mat_init(m, 2*g, 2*g);

      fmpz_mat_randtest_triangular_symplectic(m, state, n_randint(state, 50));
      if (!fmpz_mat_is_symplectic(m))
        {
	  flint_printf("FAIL\n");
	  flint_printf("g = %wd\n", g);
	  flint_printf("m = "); fmpz_mat_print(m); flint_printf("\n");
	  flint_abort();
	}
      fmpz_mat_randtest_diagonal_symplectic(m, state, n_randint(state, 50));
      if (!fmpz_mat_is_symplectic(m))
        {
	  flint_printf("FAIL\n");
	  flint_printf("g = %wd\n", g);
	  flint_printf("m = "); fmpz_mat_print(m); flint_printf("\n");
	  flint_abort();
	}
      fmpz_mat_randtest_symplectic(m, state, n_randint(state, 50));
      if (!fmpz_mat_is_symplectic(m))
        {
	  flint_printf("FAIL\n");
	  flint_printf("g = %wd\n", g);
	  flint_printf("m = "); fmpz_mat_print(m); flint_printf("\n");
	  flint_abort();
	}

      fmpz_mat_clear(m);
    }
  
  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

