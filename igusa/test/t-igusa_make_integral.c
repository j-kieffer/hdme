#include <stdlib.h>

#include "igusa.h"

int main()
{  
  slong iter;
  flint_rand_t state;
  
  flint_printf("igusa_make_integral....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 500 * flint_test_multiplier(); iter++)
    {
      slong mag_bits = 20;
      slong weights[X_NB] = X_WEIGHTS;
      fmpz* I;
      fmpq* X;
      fmpq_t scal;
      fmpq_t test;
      fmpz_t gcd;
      slong k;
      int res;
      int print = 0;

      I = _fmpz_vec_init(4);
      X = _fmpq_vec_init(X_NB);
      fmpq_init(scal);
      fmpq_init(test);
      fmpz_init(gcd);

      for (k = 0; k < 4; k++) fmpz_randtest_not_zero(&I[k], state, mag_bits);
      if (iter % 2 == 0) fmpz_zero(&I[0]);
      if (iter % 3 == 0) fmpz_zero(&I[1]);
      if (iter % 5 == 0) fmpz_zero(&I[2]);
      if (iter % 5 == 1) fmpz_zero(&I[3]);

      igusa_X(X, I);
      res = 1;
      fmpz_zero(gcd);

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
	  fmpz_set_si(gcd, 2);
	  cov_valuation_fmpq(test, X, gcd, X_NB, weights);
	  flint_printf("Weighted valuation at 2: ");
	  fmpq_print(test); flint_printf("\n");
	}
      
      for (k = 0; k < X_NB; k++)
	{
	  igusa_make_integral(scal, I, weights[k], 0);
	  fmpq_mul(test, &X[k], scal);
	  if (print)
	    {
	      flint_printf("k = %wd: ", k); fmpq_print(test); flint_printf("\n");
	    }
	  if (!fmpz_is_one(fmpq_denref(test)))
	    {
	      res = 0;
	      break;
	    }	  
	  fmpz_gcd(gcd, gcd, fmpq_numref(test));
	}

      if (!res || !fmpz_is_one(gcd))
	{
	  flint_printf("FAIL\n");
	  fflush(stdout);
	  flint_abort();
	}
      
      _fmpz_vec_clear(I, 4);
      _fmpq_vec_clear(X, X_NB);
      fmpq_clear(scal);
      fmpq_clear(test);
      fmpz_clear(gcd);      
    }

  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

      
