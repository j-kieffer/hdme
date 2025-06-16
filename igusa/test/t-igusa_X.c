#include <stdlib.h>

#include "igusa.h"

int main()
{  
  slong iter;
  flint_rand_t state;
  
  flint_printf("igusa_X....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
      fmpz* I;
      fmpq* X;
      fmpz* XZ;
      fmpz* test;
      fmpz_t scal;
      slong mag_bits = 10;
      int res;
      int print = 0;
      slong k;
      slong Xweights[X_NB] = X_WEIGHTS;
      slong Iweights[4] = IGUSA_WEIGHTS;
      slong half[4] = IGUSA_HALFWEIGHTS;

      I = _fmpz_vec_init(4);
      X = _fmpq_vec_init(X_NB);
      XZ = _fmpz_vec_init(X_NB);
      test = _fmpz_vec_init(X_NB);
      fmpz_init(scal);

      for (k = 0; k < 4; k++) fmpz_randtest(&I[k], state, mag_bits);
      fmpz_randtest_not_zero(scal, state, mag_bits);
      cov_rescale_fmpz_si(I, I, 2*3, 4, half);
      igusa_X(X, I);

      if (print)
	{
	  flint_printf("I:\n");
	  for (k = 0; k < 4; k++)
	    {
	      fmpz_print(&I[k]); flint_printf("\n");
	    }
	  flint_printf("X:\n");
	  for (k = 0; k < X_NB; k++)
	    {
	      fmpq_print(&X[k]); flint_printf("\n");
	    }
	}

      res = 1;
      for (k = 0; k < X_NB; k++)
	{
	  if (!fmpz_is_one(fmpq_denref(&X[k]))) res = 0;
	}
      if (!res)
	{
	  flint_printf("FAIL (not integral)\n");
	  fflush(stdout);
	  flint_abort();
	}
      for (k = 0; k < X_NB; k++) fmpz_set(&XZ[k], fmpq_numref(&X[k]));
      cov_rescale_fmpz(XZ, XZ, scal, X_NB, Xweights);

      cov_rescale_fmpz(I, I, scal, 4, Iweights);
      igusa_X(X, I);
      for (k = 0; k < X_NB; k++) fmpz_set(&test[k], fmpq_numref(&X[k]));

      res = 1;
      for (k = 0; k < X_NB; k++)
	{
	  if (!fmpz_equal(&XZ[k], &test[k])) res = 0;
	}
      if (!res)
	{
	  flint_printf("FAIL (rescaling)\n");
	  for (k = 0; k < X_NB; k++)
	    {
	      flint_printf("k = %wd:\n", k);
	      fmpz_print(&XZ[k]); flint_printf("\n");
	      fmpz_print(&test[k]); flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_abort();
	}      

      _fmpz_vec_clear(I, 4);
      _fmpq_vec_clear(X, X_NB);
      _fmpz_vec_clear(XZ, X_NB);
      _fmpz_vec_clear(test, X_NB);
      fmpz_clear(scal);
    }

  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

