#include <stdlib.h>
#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>

#include "hecke.h"

int main()
{
  slong ell;
  
  flint_printf("siegel_coset....");
  fflush(stdout);


  for (ell = 2; ell < 4; ell++)
    {
      slong nb = siegel_nb_cosets(ell);
      fmpz_mat_t m1, m2;
      fmpq_mat_t n, test;
      
      fmpz_mat_t c;
      fmpz_t det;
      slong k, i;

      fmpz_mat_init(m1, 4, 4);
      fmpz_mat_init(m2, 4, 4);
      fmpq_mat_init(n, 4, 4);
      fmpq_mat_init(test, 4, 4);      
      fmpz_mat_init(c, 2, 2);
      fmpz_init(det);
      
      /* Check:
	 - Determinant ell^2
	 - General symplectic
	 - Lower left is zero
	 - All distinct cosets */
      
      for (k = 0; k < nb; k++)
	{
	  siegel_coset(m1, k, ell);
	  
	  fmpz_mat_get_c(c, m1);
	  if (!fmpz_mat_is_zero(c))
	    {
	      flint_printf("FAIL (nonzero c)\n");
	      flint_printf("ell = %wd, k = %wd\n", ell, k);
	      fmpz_mat_print(m1);
	      fflush(stdout);
	      flint_abort();
	    }
	  
	  fmpz_mat_det(det, m1);
	  if (!fmpz_equal_si(det, ell*ell))
	    {
	      flint_printf("FAIL (determinant)\n");
	      flint_printf("ell = %wd, k = %wd\n", ell, k);
	      fmpz_mat_print(m1);
	      fflush(stdout);
	      flint_abort();
	    }

	  if (!fmpz_mat_is_general_symplectic(m1))
	    {
	      flint_printf("FAIL (not symplectic)\n");
	      flint_printf("ell = %wd, k = %wd\n", ell, k);
	      fmpz_mat_print(m1);
	      fflush(stdout);
	      flint_abort();
	    }
	  
	  fmpq_mat_set_fmpz_mat(n, m1);
	  fmpq_mat_inv(n, n);
	  
	  for (i = k+1; i < nb; i++)
	    {
	      siegel_coset(m2, i, ell);
	      fmpq_mat_mul_r_fmpz_mat(test, m2, n);
	      	      
	      if (fmpq_mat_is_integral(test))
		{
		  flint_printf("FAIL (same cosets)\n");
		  flint_printf("ell = %wd, k = %wd, i = %wd\n", ell, k, i);
		  fmpz_mat_print(m1);
		  fmpz_mat_print(m2);
		  fmpq_mat_print(test);
		  fflush(stdout);
		  flint_abort();
		}
	    }
	}
	 
      fmpz_mat_clear(m1);
      fmpz_mat_clear(m2);
      fmpq_mat_clear(n);
      fmpq_mat_clear(test);
      fmpz_mat_clear(c);
      fmpz_clear(det);
    }
  
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
      
