#include <stdlib.h>

#include "siegel.h"

int main()
{
  slong i, j, g;
  fmpz_mat_t ui, uj;
  
  flint_printf("siegel_test_matrix....");
  fflush(stdout);

  /* Test matrices are distinct and symplectic */
  for (g = 1; g <= 2; g++)
    {
      fmpz_mat_init(ui, 2*g, 2*g);
      fmpz_mat_init(uj, 2*g, 2*g);
      for (i = 0; i < siegel_nb_test_matrices(g); i++)
	{
	  siegel_test_matrix(ui, i);
	  /* fmpz_mat_print(ui); flint_printf("\n\n"); */
	  if (!fmpz_mat_is_symplectic(ui))
	    {
	      flint_printf("FAIL (not symplectic)\n");
	      flint_printf("i = %wd\n", i);
	      fmpz_mat_print(ui); flint_printf("\n\n");
	      flint_abort();
	    }
	  for (j = 0; j < siegel_nb_test_matrices(g); j++)
	    {
	      siegel_test_matrix(uj, j);
	      if (i != j && fmpz_mat_equal(ui, uj))
		{
		  flint_printf("FAIL (equal matrices)\n");
		  flint_printf("(i, j) = (%wd, %wd)\n", i, j);
		  fmpz_mat_print(ui); flint_printf("\n\n");
		  flint_abort();
		}
	    }
	}
      fmpz_mat_clear(ui);
      fmpz_mat_clear(uj);
    }
  
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
