#include <stdlib.h>
#include "modular.h"

int main()
{
  
  slong iter;
  flint_rand_t state;
  
  flint_printf("modeq_ctx_add_pair....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
      modeq_ctx_t ctx;
      slong nb = 1+n_randint(state, 200);
      slong k;
      slong* i1;
      slong* i2;
      int res;

      i1 = flint_malloc(nb * sizeof(slong));
      i2 = flint_malloc(nb * sizeof(slong));
      modeq_ctx_init(ctx);

      for (k = 0; k < nb; k++)
	{
	  i1[k] = n_randint(state, 1000);
	  i2[k] = n_randint(state, 1000);
	  modeq_ctx_add_pair(ctx, i1[k], i2[k]);
	}
      res = 1;
      for (k = 0; k < nb; k++)
	{
	  if (!modeq_ctx_is_pair(i1[k], i2[k], ctx))
	    {
	      res = 0;
	      break;
	    }
	}
      if (modeq_ctx_is_pair(-1, -1, ctx)) res = 0;
      
      if (!res)
	{
	  flint_printf("FAIL\n");
	  fflush(stdout);
	  flint_abort();
	}

      flint_free(i1);
      flint_free(i2);
      modeq_ctx_clear(ctx);
    }
  
  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

