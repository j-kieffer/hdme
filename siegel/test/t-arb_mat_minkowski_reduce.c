#include <stdlib.h>
#include <flint/fmpq.h>

#include "siegel.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("arb_mat_minkowski_reduce....");
  fflush(stdout);
  
  flint_rand_init(state);

  for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
      arb_t y1;
      arb_t y2;
      arb_t y3;
      arb_t tol;
      fmpq_t rand;
      fmpz_t den;
      
      arb_mat_t r;
      arb_mat_t m;
      arb_mat_t r2;
      fmpz_mat_t u;
      fmpz_mat_t uinv;
      fmpz_mat_t u2;
      
      int res;

      slong g = 2;
      slong bits = 1 + n_randint(state, 500);
      slong prec = 10 + 5 * bits;
      slong chance = n_randint(state, 10);

      arb_init(y1);
      arb_init(y2);
      arb_init(y3);
      arb_init(tol);
      fmpq_init(rand);
      fmpz_init(den);
      
      arb_mat_init(r, g, g);
      arb_mat_init(m, g, g);
      arb_mat_init(r2, g, g);
      fmpz_mat_init(u, g, g);
      fmpz_mat_init(uinv, g, g);
      fmpz_mat_init(u2, g, g);
      
      arb_set_si(tol, 1);
      arb_mul_2exp_si(tol, tol, -bits); /* tol is larger than 2^(-prec) */
      
      /* Generate random Minkowski-reduced matrix */
      fmpq_randtest_not_zero(rand, state, bits);
      fmpq_abs(rand, rand);
      arb_set_fmpq(y1, rand, prec);
      fmpq_randtest_not_zero(rand, state, bits);
      fmpq_abs(rand, rand);
      arb_set_fmpq(y2, rand, prec);
      if (chance == 0) arb_set(y2, y1); /* sometimes on the boundary */
      if (arb_lt(y2, y1)) arb_swap(y1, y2);
      
      fmpq_randtest(rand, state, bits);
      fmpq_abs(rand, rand);
      while (fmpq_cmp_si(rand, 1) > 0) fmpq_div_2exp(rand, rand, 1);
      if (chance == 1) fmpq_set_si(rand, 1, 1); /* sometimes on the boundary */
      fmpq_div_2exp(rand, rand, 1);
      arb_set_fmpq(y3, rand, prec);
      arb_mul(y3, y3, y1, prec);

      arb_set(arb_mat_entry(r, 0, 0), y1);
      arb_set(arb_mat_entry(r, 0, 1), y3);
      arb_set(arb_mat_entry(r, 1, 0), y3);
      arb_set(arb_mat_entry(r, 1, 1), y2);

      if (!arb_mat_is_minkowski_reduced(r, tol, prec))
	{	  
	  flint_printf("FAIL\n");
	  flint_printf("r = "); arb_mat_printd(r, 30); flint_printf("\n\n");
	  flint_abort();
	}

      /* Generate random transformation in SL_2(ZZ) */
      fmpz_mat_one(uinv);
      fmpz_mat_randops(uinv, state, bits);
      arb_mat_congr_fmpz_mat(m, uinv, r, prec);
      fmpz_mat_inv(u, den, uinv);
      if (!fmpz_is_one(den)) fmpz_mat_neg(u, u); /* in case den == -1 */

      /* Reduce */
      res = arb_mat_minkowski_reduce(r2, u2, m, tol, prec);
      if (!res || !arb_mat_overlaps(r, r2))
	{
	  flint_printf("FAIL (res=%d)\n", res);
	  flint_printf("r = "); arb_mat_printd(r, 30); flint_printf("\n\n");
	  flint_printf("u = "); fmpz_mat_print_pretty(u); flint_printf("\n\n");
	  flint_printf("uinv = "); fmpz_mat_print_pretty(uinv); flint_printf("\n\n");
	  flint_printf("m = "); arb_mat_printd(m, 30); flint_printf("\n\n");
	  flint_printf("u2 = "); fmpz_mat_print_pretty(u2); flint_printf("\n\n");
	  flint_printf("r2 = "); arb_mat_printd(r2, 30); flint_printf("\n\n");
	  flint_abort();
	}
      /* u and u2 may be different when r lies on the boundary. */
      
      arb_clear(y1);
      arb_clear(y2);
      arb_clear(y3);
      arb_clear(tol);
      fmpq_clear(rand);
      fmpz_clear(den);
      
      arb_mat_clear(r);
      arb_mat_clear(m);
      arb_mat_clear(r2);
      fmpz_mat_clear(u);
      fmpz_mat_clear(uinv);
      fmpz_mat_clear(u2);
    }
  
  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
