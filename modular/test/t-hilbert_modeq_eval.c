#include <stdlib.h>

#include "modular.h"

int main()
{
  slong iter;
  flint_rand_t state;

  flint_rand_init(state);
  
  flint_printf("hilbert_modeq_eval....");
  fflush(stdout);

  for (iter = 0; iter < 10; iter++)
    {
      fmpz* G;
      fmpz* I;
      slong primes[3] = {11, 19, 29};
      slong ell = primes[iter%3];
      slong delta = 5;
      slong k;
      slong mag_bits = 10;
      modeq_t E;
      modeq_ctx_t ctx;
      int print = 1;
      int valid = 0;

      G = _fmpz_vec_init(3);
      I = _fmpz_vec_init(4);
      modeq_ctx_init(ctx);
      modeq_init(E);

      while (!valid)
	{
	  for (k = 0; k < 3; k++) fmpz_randtest_not_zero(&G[k], state, mag_bits);
	  if (iter % 2 == 0) fmpz_zero(&G[0]);
	  if (iter % 3 == 0) fmpz_zero(&G[1]);
	  if (iter % 5 == 0) fmpz_zero(&G[2]);
	  k = n_randint(state, 3);
	  fmpz_randtest_not_zero(&G[k], state, mag_bits);
	  
	  igusa_from_gundlach_fmpz(I, G, delta);
	  valid = !fmpz_is_zero(igusa_chi12(I)) || !fmpz_is_zero(igusa_chi10(I));
	}
	
      if (print)
	{
	  flint_printf("Igusa covariants:\n");
	  for (k = 0; k < 4; k++)
	    {
	      fmpz_print(&I[k]); flint_printf("\n");
	    }
	}
      hilbert_modeq_eval(E, ctx, I, ell, delta);

      _fmpz_vec_clear(G, 3);
      _fmpz_vec_clear(I, 4);
      modeq_ctx_clear(ctx);
      modeq_clear(E);
    }
  
  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}



