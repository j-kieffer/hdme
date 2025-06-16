#include <stdlib.h>
#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>

#include "hecke.h"

int main()
{
  slong p;
  
  flint_printf("siegel_T1_coset....");
  fflush(stdout);

  for (p = 2; p < 6; p++)
    {
      slong nb = siegel_nb_T1_cosets(p);
      slong i, j;
      fmpz_mat_t m1, m2;
      fmpq_mat_t n, test;
      fmpz_t det;

      if (!n_is_prime(p)) continue;

      fmpz_mat_init(m1, 4, 4);
      fmpz_mat_init(m2, 4, 4);
      fmpq_mat_init(n, 4, 4);
      fmpq_mat_init(test, 4, 4);
      fmpz_init(det);

      for (i = 0; i < nb; i++)
	{
	  siegel_T1_coset(m1, i, p);
	  if (!fmpz_mat_is_general_symplectic(m1))
	    {
	      flint_printf("FAIL (not symplectic)\n");
	      flint_printf("i = %wd, p = %wd\n", i, p);
	      fmpz_mat_print(m1); flint_printf("\n");
	      fflush(stdout);
	      flint_abort();
	    }
	  
	  fmpz_mat_det(det, m1);
	  if (!fmpz_equal_si(det, n_pow(p, 4)))
	    {
	      flint_printf("FAIL (determinant)\n");
	      flint_printf("p = %wd, i = %wd\n", p, i);
	      fmpz_mat_print(m1);
	      fflush(stdout);
	      flint_abort();
	    }
	  
	  fmpq_mat_set_fmpz_mat(n, m1);
	  fmpq_mat_inv(n, n);

	  for (j = i+1; j < nb; j++)
	    {
	      siegel_T1_coset(m2, j, p);
	      fmpq_mat_mul_r_fmpz_mat(test, m2, n);
	      if (fmpq_mat_is_integral(test))
		{
		  flint_printf("FAIL (same cosets)\n");
		  flint_printf("i = %wd, j = %wd, p = %wd\n", i, j, p);
		  fmpz_mat_print_pretty(m1); flint_printf("\n");
		  fmpz_mat_print_pretty(m2); flint_printf("\n");
		  fmpq_mat_print(test); flint_printf("\n");
		  fflush(stdout);
		  flint_abort();
		}
	    }
	}

      fmpz_mat_clear(m1);
      fmpz_mat_clear(m2);
      fmpq_mat_clear(n);
      fmpq_mat_clear(test);
      fmpz_clear(det);
    }

  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
      
