#include <stdlib.h>
#include "hecke.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("hecke_operator....");
  fflush(stdout);
  flint_rand_init(state);

  for (iter = 0; iter < 1 * flint_test_multiplier(); iter++)
    {
      slong ell = 2;
      slong ellmax = 10;
      acb_mat_t tau;
      hecke_t H;
      acb_ptr val;
      acb_t r;
      fmpz_t eig;
      slong k;
      slong nb;
      slong prec = 500 + n_randint(state, 1000);
      int res;

      acb_mat_init(tau, 2, 2);
      acb_init(r);
      fmpz_init(eig);
      
      siegel_fundamental_domain_randtest(tau, state, prec);
            
      while (ell < ellmax)
	{
	  nb = siegel_nb_cosets(ell);
	  hecke_init(H, nb);
	  val = _acb_vec_init(nb);

	  res = hecke_set_tau(H, tau, prec);
	  if (res) res = hecke_collect_siegel(H, ell, prec);
	  if (!res)
	    {
	      flint_printf("FAIL (Hecke)\n");
	      hecke_print(H, 10);
	      fflush(stdout);
	      flint_abort();
	    }

	  /* Test expected eigenvalues for I4, I6' */
	  for (k = 0; k < nb; k++) acb_set(&val[k], igusa_psi4(hecke_I(H, k)));
	  hecke_operator(r, H, val, ell, 4, 0, prec);
	  acb_div(r, r, &hecke_I_tau(H)[0], prec);
	  
	  hecke_eigenvalue_eisenstein(eig, 4, ell);
	  if (!acb_contains_fmpz(r, eig))
	    {
	      flint_printf("FAIL (I4 eigenvalue)\n");
	      flint_printf("ell = %wd\n", ell);
	      acb_printd(r, 30); flint_printf("\n");
	      fmpz_print(eig); flint_printf("\n");
	      fflush(stdout);
	      flint_abort();
	    }
	  
	  for (k = 0; k < nb; k++) acb_set(&val[k], igusa_psi6(hecke_I(H, k))); /* I6' */
	  hecke_operator(r, H, val, ell, 6, 0, prec);
	  acb_div(r, r, &hecke_I_tau(H)[1], prec);
	  
	  hecke_eigenvalue_eisenstein(eig, 6, ell);
	  if (!acb_contains_fmpz(r, eig))
	    {
	      flint_printf("FAIL (I6' eigenvalue)\n");
	      flint_printf("ell = %wd\n", ell);
	      acb_printd(r, 30); flint_printf("\n");
	      fmpz_print(eig); flint_printf("\n");
	      fflush(stdout);
	      flint_abort();
	    }

	  hecke_clear(H);
	  _acb_vec_clear(val, nb);
	  ell = n_nextprime(ell, 1);
	}

      acb_mat_clear(tau);
      acb_clear(r);
      fmpz_clear(eig);
    }
  
  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
