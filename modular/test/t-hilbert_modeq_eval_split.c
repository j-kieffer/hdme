#include <stdlib.h>
#include "modular.h"

int main()
{
  
  slong iter;
  flint_rand_t state;
  
  flint_printf("hilbert_modeq_eval_split....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 1 * flint_test_multiplier(); iter++)
    {
      fmpq* mn;
      fmpz* G;
      fmpz* I;
      modeq_t R1, R2;
      modeq_ctx_t ctx1, ctx2;
      slong mag_bits = 3;
      int print = 1;
      slong delta = 5;
      slong primes[3] = {11,19,29};
      slong ell = primes[iter%3];
      slong k;
      int res;

      mn = _fmpq_vec_init(2);
      G = _fmpz_vec_init(3);
      I = _fmpz_vec_init(4);
      modeq_init(R1);
      modeq_init(R2);
      modeq_ctx_init(ctx1);
      modeq_ctx_init(ctx2);

      res = 0;
      while (!res)
	{
	  for (k = 0; k < 2; k++) fmpq_randtest(&mn[k], state, mag_bits);
	  gundlach_from_hilbert_param(G, mn, delta);
	  igusa_from_gundlach_fmpz(I, G, delta);
	  res = !fmpz_is_zero(igusa_chi10(I)) || !fmpz_is_zero(igusa_chi12(I));
	}

      if (print)
	{
	  flint_printf("Igusa covariants:\n");
	  for (k = 0; k < 4; k++)
	    {
	      fmpz_print(&I[k]); flint_printf("\n");
	    }
	}
	  
      hilbert_modeq_eval_split(R1, R2, ctx1, ctx2, I, ell, delta);

      _fmpq_vec_clear(mn, 2);
      _fmpz_vec_clear(G, 3);
      _fmpz_vec_clear(I, 4);      
      modeq_clear(R1);
      modeq_clear(R2);
      modeq_ctx_clear(ctx1);
      modeq_ctx_clear(ctx2);
    }
  
  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

