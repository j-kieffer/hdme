#include <stdlib.h>
#include "modular.h"

int main()
{
  
  slong iter;
  flint_rand_t state;
  
  flint_printf("siegel_modeq_eval....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 2 * flint_test_multiplier(); iter++)
    {
      fmpz* I;
      modeq_t R;
      slong ell = 2;
      modeq_ctx_t ctx;
      slong k;
      slong mag_bits = 20;
      int print = 1;
      int zeroentries = 0;

      I = _fmpz_vec_init(4);
      modeq_init(R);
      modeq_ctx_init(ctx);

      for (k = 0; k < 4; k++) fmpz_randtest_not_zero(&I[k], state, mag_bits);
      if (zeroentries)
	{
	  if (iter % 2 == 0) fmpz_zero(igusa_psi4(I));
	  if (iter % 3 == 0) fmpz_zero(igusa_psi6(I));
	  if (iter % 3 == 1) fmpz_zero(igusa_chi10(I));
	}

      if (print)
	{
	  flint_printf("I:\n");
	  for (k = 0; k < 4; k++)
	    {
	      fmpz_print(&I[k]);
	      flint_printf("\n");
	    }
	}

      siegel_modeq_eval(R, ctx, I, ell);

      if (fmpz_poly_is_zero(modeq_equation(R)))
	{
	  flint_printf("FAIL\n");
	  fflush(stdout);
	  flint_abort();
	}

      _fmpz_vec_clear(I, 4);
      modeq_clear(R);
      modeq_ctx_clear(ctx);      
    }
  
  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}



