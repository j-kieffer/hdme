#include <stdlib.h>
#include "modular.h"

int main()
{
  
  slong iter;
  flint_rand_t state;
  
  flint_printf("modeq_ctx_choose....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
      slong nb = 1 + n_randint(state, 100);
      slong i1, i2;
      acb_ptr I;
      modeq_ctx_t ctx;
      acb_t ev;
      slong prec = 1000;
      slong mag_bits = 10;
      slong k;
      int res;

      I = _acb_vec_init(4*nb);
      modeq_ctx_init(ctx);
      acb_init(ev);

      for (k = 0; k < 4*nb; k++) acb_randtest_precise(&I[k], state, prec, mag_bits);
      for (k = 0; k < nb; k++) acb_one(&I[4*k + n_randint(state, 3)]);
      if (iter % 2 == 0) acb_zero(&I[0]);
      if (iter % 3 == 0) acb_zero(&I[1]);
      if (iter % 3 == 1) acb_zero(&I[2]);

      i1 = n_randint(state, nb);
      i2 = n_randint(state, nb);
      for (k = 0; k < 4; k++) acb_set(&I[4*i1 + k], &I[4*i2 + k]);
      
      modeq_ctx_choose(ctx, I, nb, prec);

      if (i1 != i2 && !modeq_ctx_is_pair(i1, i2, ctx)
	  && !modeq_ctx_is_pair(i2, i1, ctx))
	{
	  flint_printf("FAIL (pair)\n");
	  fflush(stdout);
	  flint_abort();
	}

      res = 1;
      for (k = 0; k < nb; k++)
	{
	  cov_mpoly_eval(ev, modeq_ctx_den(ctx), &I[4*k], modeq_ctx_ctx(ctx), prec);
	  if (acb_contains_zero(ev)) res = 0;	    
	}

      if (!res)
	{
	  flint_printf("FAIL (denominator)\n");
	  fflush(stdout);
	  flint_abort();
	}

      _acb_vec_clear(I, 4*nb);
      modeq_ctx_clear(ctx);
      acb_clear(ev);
    }
  
  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}


