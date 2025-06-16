#include <stdlib.h>

#include "igusa.h"

int main()
{  
  slong iter;
  flint_rand_t state;
  
  flint_printf("igusa_from_monomials....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 500 * flint_test_multiplier(); iter++)
    {
      fmpz* I;
      fmpz* M;
      fmpz* test;
      slong nb;
      fmpz_mpoly_ctx_t ctx;
      fmpz_mpoly_t mon;
      slong k;
      slong wt;
      slong mag_bits = 10;
      slong wts[3] = {20, 30, 60};
      slong weights[4] = IGUSA_HALFWEIGHTS;
      int print = 0;

      I = _fmpz_vec_init(4);
      test = _fmpz_vec_init(4);
      /* Init M later when wt is known */
      fmpz_mpoly_ctx_init(ctx, 4, ORD_LEX);
      fmpz_mpoly_init(mon, ctx);

      for (k = 0; k < 4; k++) fmpz_randtest_not_zero(&I[k], state, mag_bits);
      if (iter % 2 == 0) fmpz_zero(igusa_psi4(I));
      if (iter % 3 == 0) fmpz_zero(igusa_psi6(I));
      if (iter % 5 == 0) fmpz_zero(igusa_chi10(I));
      if (iter % 5 == 1) fmpz_zero(igusa_chi12(I));

      cov_adjust_weights(weights, weights, I, 4);
      cov_normalize_fmpz(I, I, 4, weights);

      k = n_randint(state, 3);
      wt = wts[k];
      
      if (fmpz_is_zero(igusa_psi4(I)))
	{
	  wt = FLINT_MAX(wt, 30);
	}
      if (wt == 30 && fmpz_is_zero(igusa_psi6(I)) &&
	  (fmpz_is_zero(igusa_psi4(I)) || fmpz_is_zero(igusa_chi10(I))))
	{
	  wt = 60;
	}

      nb = igusa_nb_base_monomials(wt);
      M = _fmpz_vec_init(nb);
      for (k = 0; k < nb; k++)
	{
	  igusa_base_monomial(mon, wt, k, ctx);
	  cov_mpoly_eval_fmpz(&M[k], mon, I, ctx);
	}

      if (print)
	{
	  flint_printf("I:\n");
	  for (k = 0; k < 4; k++)
	    {
	      fmpz_print(&I[k]); flint_printf("\n");
	    }
	  flint_printf("Adjusted weights: ");
	  for (k = 0; k < 4; k++) flint_printf("%wd ", weights[k]);
	  flint_printf("\n");
	  flint_printf("Weight %wd, monomials:\n", wt);
	  for (k = 0; k < nb; k++)
	    {
	      fmpz_print(&M[k]); flint_printf("\n");
	    }
	}

      igusa_from_monomials(test, M, wt);
      if (!_fmpz_vec_equal(test, I, 4))
	{
	  flint_printf("FAIL\n");
	  flint_printf("Weight %wd, I:\n", wt);
	  for (k = 0; k < 4; k++)
	    {
	      fmpz_print(&I[k]); flint_printf("\n");
	    } 
	  flint_printf("test:\n");
	  for (k = 0; k < 4; k++)
	    {
	      fmpz_print(&test[k]); flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_abort();
	}
      
      _fmpz_vec_clear(I, 4);
      _fmpz_vec_clear(test, 4);
      _fmpz_vec_clear(M, nb);
      fmpz_mpoly_clear(mon, ctx);
      fmpz_mpoly_ctx_clear(ctx);      
    }
  
  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

