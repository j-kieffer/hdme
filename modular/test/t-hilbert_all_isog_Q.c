#include <stdlib.h>

#include "modular.h"

int main()
{
  slong iter;
  
  flint_printf("hilbert_all_isog_Q....");
  fflush(stdout);

  for (iter = 0; iter < 2; iter++)
    {
      fmpz* I1;
      slong ell;
      slong delta;
      slong max_nb_roots = 3;
      slong nb_roots = 1;
      fmpz* all_I;
      int print = 1;
      slong k, j;
            
      I1 = _fmpz_vec_init(4);
      all_I = _fmpz_vec_init(4 * max_nb_roots);
            
      if (iter == 0)
	{
	  /* https://beta.lmfdb.org/Genus2Curve/Q/529/a/529/1 */
	  fmpz_set_si(&I1[0], 284);
	  fmpz_set_si(&I1[1], 2401);
	  fmpz_set_si(&I1[2], 246639);
	  fmpz_set_si(&I1[3], -67712);
	  ell = 11;
	  delta = 5;
	}
      else
	{
	  /* https://beta.lmfdb.org/Genus2Curve/Q/841/a/841/1 */
	  fmpz_set_si(&I1[0], 1420);
	  fmpz_set_si(&I1[1], 4201);
	  fmpz_set_si(&I1[2], 1973899);
	  fmpz_set_si(&I1[3], 107648);
	  ell = 7;
	  delta = 8;
	}
      
      igusa_from_IC_fmpz(I1, I1);
      hilbert_all_isog_Q(&nb_roots, all_I, I1, ell, delta);
      if (nb_roots == 0)
	{
	  flint_printf("FAIL (roots)\n");
	  flint_printf("nb_roots = %wd\n", nb_roots);
	  fflush(stdout);
	  flint_abort();
	}

      if (print)
	{
	  flint_printf("Number of roots: %wd\n", nb_roots);
	  for (k = 0; k < nb_roots; k++)
	    {
	      flint_printf("Igusa covariants at root number %wd:\n", k+1);
	      for (j = 0; j < 4; j++)
		{
		  fmpz_print(&all_I[4*k + j]);
		  flint_printf("\n");
		}
	    }
	}

      _fmpz_vec_clear(I1, 4);
      _fmpz_vec_clear(all_I, 4 * max_nb_roots);
    }

  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

