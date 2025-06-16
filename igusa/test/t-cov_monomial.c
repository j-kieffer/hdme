#include <stdlib.h>

#include "igusa.h"

int main()
{  
  slong iter;
  flint_rand_t state;
  
  flint_printf("cov_monomial....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
      slong nb = 1 + n_randint(state, 5);
      slong* exps;
      slong* test;
      fmpz_mpoly_ctx_t ctx;
      fmpz_mpoly_t mon;
      acb_ptr val;
      acb_t ev1, ev2;
      slong mag_bits = 5;
      slong prec = 500;
      slong j, k;
      int res;

      exps = flint_malloc(nb * sizeof(slong));
      test = flint_malloc(nb * sizeof(slong));
      fmpz_mpoly_ctx_init(ctx, nb, ORD_LEX);
      fmpz_mpoly_init(mon, ctx);
      val = _acb_vec_init(nb);
      acb_init(ev1);
      acb_init(ev2);

      for (k = 0; k < nb; k++) exps[k] = n_randint(state, 10);
      cov_monomial(mon, exps, ctx);
      cov_monomial_degrees(test, mon, ctx);
      res = 1;
      for (k = 0; k < nb; k++)
	{
	  if (exps[k] != test[k]) res = 0;
	}
      if (!res)
	{
	  flint_printf("FAIL (exponents)\n");
	  for (k = 0; k < nb; k++) flint_printf("%wd, %wd\n", exps[k], test[k]);
	  fflush(stdout);
	  flint_abort();
	}

      for (k = 0; k < nb; k++) acb_randtest_precise(&val[k], state, prec, mag_bits);
      for (k = 0; k < nb; k++)
	{
	  for (j = 0; j < nb; j++) exps[j] = 0;
	  exps[k] = 1 + n_randint(state, 10);
	  cov_monomial(mon, exps, ctx);
	  cov_mpoly_eval(ev1, mon, val, ctx, prec);
	  acb_pow_si(ev2, &val[k], exps[k], prec);
	  if (!acb_overlaps(ev1, ev2))
	    {
	      flint_printf("FAIL (eval)\n");
	      flint_printf("exps[%wd] = %wd\n", k, exps[k]);
	      acb_printd(&val[k], 10); flint_printf("\n");
	      acb_printd(ev1, 10); flint_printf("\n");
	      acb_printd(ev2, 10); flint_printf("\n");
	    }
	}      

      fmpz_mpoly_clear(mon, ctx);
      fmpz_mpoly_ctx_clear(ctx);
      flint_free(exps);
      flint_free(test);
      _acb_vec_clear(val, nb);
      acb_clear(ev1);
      acb_clear(ev2);
    }

  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
