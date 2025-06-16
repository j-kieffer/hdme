#include <stdlib.h>

#include "igusa.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("cov_eval_all_monomials....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
      slong nb = 1 + n_randint(state, 4);
      slong mag_bits = 10;
      slong prec = 200;
      slong* weights;
      slong* exps;
      slong wt;
      acb_ptr I;
      acb_ptr ev;
      fmpz_mpoly_ctx_t ctx;
      fmpz_mpoly_t mon;
      acb_t test;
      slong m;
      slong k, j;
      int print = 0;
      int res;

      weights = flint_malloc(nb * sizeof(slong));
      I = _acb_vec_init(nb);      
      for (k = 0; k < nb; k++)
	{
	  weights[k] = 1 + n_randint(state, 10);
	  acb_randtest_precise(&I[k], state, prec, mag_bits);
	}
      wt = n_randint(state, 20);
      m = cov_nb_monomials(wt, nb, weights);
      ev = _acb_vec_init(m);
      exps = flint_malloc(m*nb * sizeof(slong));
      fmpz_mpoly_ctx_init(ctx, nb, ORD_LEX);
      fmpz_mpoly_init(mon, ctx);
      acb_init(test);

      cov_all_exps(exps, wt, nb, weights);

      if (print)
	{
	  flint_printf("Weights:\n");
	  for (k = 0; k < nb; k++) flint_printf("%wd ", weights[k]);
	  flint_printf("\nTarget weight %wd, number of monomials: %wd\n", wt, m);
	  flint_printf("Exponents:\n");
	  for (k = 0; k < m; k++)
	    {
	      for (j = 0; j < nb; j++) flint_printf("%wd ", exps[k*nb + j]);
	      flint_printf("\n");
	    }
	}

      cov_eval_all_monomials(ev, I, wt, nb, weights, prec);

      if (m >= 1)
	{
	  k = n_randint(state, m);
	  cov_monomial(mon, &exps[k*nb], ctx);
	  cov_mpoly_eval(test, mon, I, ctx, prec);
	  
	  res = 0;
	  for (k = 0; k < m; k++)
	    {
	      if (acb_overlaps(&ev[k], test)) res = 1;
	    }
	  
	  if (!res)
	    {
	      flint_printf("FAIL\n");
	      fflush(stdout);
	      flint_abort();
	    }
	}

      _acb_vec_clear(ev, m);
      _acb_vec_clear(I, nb);
      flint_free(weights);
      flint_free(exps);
      fmpz_mpoly_clear(mon, ctx);
      fmpz_mpoly_ctx_clear(ctx);
      acb_clear(test);
    }

  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

      
