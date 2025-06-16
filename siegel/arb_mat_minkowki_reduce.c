
#include "siegel.h"

static int
arb_round_fmpz(fmpz_t z, const arb_t x, const arb_t tol)
{
  arb_t rad;
  int res;
  arb_init(rad);
  
  res = 1;
  arb_get_rad_arb(rad, x);
  
  if (!arb_lt(rad, tol)) res = 0;
  if (res)
    {
      arf_get_fmpz(z, arb_midref(x), ARF_RND_NEAR);
    }
  
  arb_clear(rad);
  return res;
}

static int
minkowski_g2_red_y3(fmpz_mat_t u, const arb_mat_t m,
		    const arb_t tol, slong prec)
{
  arb_t quo;
  int res;
  fmpz_t round;

  arb_init(quo);
  fmpz_init(round);

  arb_inv(quo, arb_mat_entry(m, 0, 0), prec);
  arb_mul(quo, quo, arb_mat_entry(m, 0, 1), prec);
  res = arb_round_fmpz(round, quo, tol);

  fmpz_mat_one(u);
  fmpz_neg(fmpz_mat_entry(u, 1, 0), round);

  arb_clear(quo);
  fmpz_clear(round);
  return res;
}

static void
minkowski_g2_sgn_y3(fmpz_mat_t u, const arb_mat_t m,
		    const arb_t tol, slong prec)
{
  fmpz_mat_one(u);
  if (arf_sgn(arb_midref(arb_mat_entry(m, 0, 1))) < 0)
    {
      fmpz_neg(fmpz_mat_entry(u, 1, 1), fmpz_mat_entry(u, 1, 1));
    }
}

static void
minkowski_g2_swap(fmpz_mat_t u, slong prec)
{
  fmpz_mat_zero(u);
  fmpz_one(fmpz_mat_entry(u, 0, 1));
  fmpz_one(fmpz_mat_entry(u, 1, 0));
}

static int
minkowski_reduce_g2(arb_mat_t r, fmpz_mat_t u, const arb_mat_t m,
		    const arb_t tol, slong prec)
{
  int res = 1;
  int iter = 0;
  int stop = 0;
  arb_mat_t r_cur;
  fmpz_mat_t u_step;
  slong g = arb_mat_nrows(m);

  arb_mat_init(r_cur, g, g);
  fmpz_mat_init(u_step, g, g);
  
  arb_mat_set(r_cur, m);
  fmpz_mat_one(u);
  
  while (!stop)
    {
      iter++;
      (void)iter; /* Used in debug code */
      
      /* Recompute r */
      arb_mat_congr_fmpz_mat(r_cur, u, m, prec);
      
      /*flint_printf("Step %d\nr = ", iter); arb_mat_printd(r_cur, 30); flint_printf("\n\n");
	flint_printf("u = "); fmpz_mat_print_pretty(u); flint_printf("\n\n");
	if (iter > 20) flint_abort();*/
      
      /* Reduce y3 */
      res = minkowski_g2_red_y3(u_step, r_cur, tol, prec);
      if (!res) break;
      fmpz_mat_mul(u, u_step, u);
      arb_mat_congr_fmpz_mat(r_cur, u_step, r_cur, prec);

      /*flint_printf("Reduce y3\nr = "); arb_mat_printd(r_cur, 30); flint_printf("\n\n");
	flint_printf("u = "); fmpz_mat_print_pretty(u); flint_printf("\n\n");*/

      /* Sign of y3 */
      minkowski_g2_sgn_y3(u_step, r_cur, tol, prec);
      fmpz_mat_mul(u, u_step, u);
      arb_mat_congr_fmpz_mat(r_cur, u_step, r_cur, prec);
      
      /*flint_printf("Sign of y3\nr = "); arb_mat_printd(r_cur, 30); flint_printf("\n\n");
	flint_printf("u = "); fmpz_mat_print_pretty(u); flint_printf("\n\n");*/

      /* Swap? */
      if (arb_lt(arb_mat_entry(r_cur, 1, 1), arb_mat_entry(r_cur, 0, 0)))
	{
	  minkowski_g2_swap(u_step, prec);
	  fmpz_mat_mul(u, u_step, u);
	  arb_mat_congr_fmpz_mat(r_cur, u_step, r_cur, prec);
	}
      else stop = 1;
      
      /*flint_printf("Swap\nr = "); arb_mat_printd(r_cur, 30); flint_printf("\n\n");
	flint_printf("u = "); fmpz_mat_print_pretty(u); flint_printf("\n\n");*/
      
    }
  /* Control result */
  arb_mat_congr_fmpz_mat(r_cur, u, m, prec);
  if (res)
    {
      res = arb_mat_is_minkowski_reduced(r_cur, tol, prec);
    }
  arb_mat_set(r, r_cur);

  arb_mat_clear(r_cur);
  fmpz_mat_clear(u_step);
  return res;
}


int
arb_mat_minkowski_reduce(arb_mat_t r, fmpz_mat_t u, const arb_mat_t m,
			 const arb_t tol, slong prec)
{
  slong g = arb_mat_nrows(m);
  int res;
  
  switch(g)
    {
    case 1:
      res = 1;
      break;
    case 2:
      res = minkowski_reduce_g2(r, u, m, tol, prec);
      break;
    default:
      flint_fprintf(stderr, "Error: Minkowski reduction not implemented for g=%wd\n", g);
      flint_abort();
    }
  return res;
}
