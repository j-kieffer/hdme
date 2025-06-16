#include <stdlib.h>

#include "igusa.h"

int main()
{  
  slong iter;
  flint_rand_t state;
  slong wts[3] = {20, 30, 60};
  
  flint_printf("igusa_base_monomial....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 3; iter++)
    {
      slong wt = wts[iter];
      slong weights[4] = IGUSA_WEIGHTS;
      slong k;
      fmpz_mpoly_ctx_t ctx;
      fmpz_mpoly_t mon;
      slong exps[4];
      slong test;
      slong j;

      fmpz_mpoly_ctx_init(ctx, 4, ORD_LEX);
      fmpz_mpoly_init(mon, ctx);

      for (k = 0; k < igusa_nb_base_monomials(wt); k++)
	{
	  igusa_base_monomial(mon, wt, k, ctx);
	  cov_monomial_degrees(exps, mon, ctx);
	  test = 0;
	  for (j = 0; j < 4; j++) test += weights[j] * exps[j];
	  if (test != wt)
	    {
	      flint_printf("FAIL\n");
	      flint_printf("Weight %wd, k = %wd, exponents:\n", wt, k);
	      for (j = 0; j < 4; j++) flint_printf("%wd ", exps[j]);
	      flint_printf("\n");
	      igusa_base_exps(exps, wt, k);
	      for (j = 0; j < 4; j++) flint_printf("%wd ", exps[j]);
	      flint_printf("\n");
	      fflush(stdout);
	      flint_abort();
	    }
	}

      fmpz_mpoly_clear(mon, ctx);
      fmpz_mpoly_ctx_clear(ctx);
    }

  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
