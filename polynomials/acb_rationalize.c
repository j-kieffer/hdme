#include <gmp.h>
#include <flint/fmpz_vec.h>
#include "polynomials.h"

static int
cont_frac_step(fmpz_t r, arf_t next, const arf_t current, slong prec, slong tol_exp)
{
  int res = 0;    
  arf_get_fmpz(r, current, ARF_RND_FLOOR);
  /* arf_printd(current, 30); flint_printf("\n"); */
  arf_sub_fmpz(next, current, r, prec, ARF_RND_NEAR);
  /* arf_printd(next, 30); flint_printf("\n"); */
  if (arf_cmp_2exp_si(next, tol_exp) < 0)
    {
      /* flint_printf("Stop cont_frac with value ");
	 arf_printd(next, 30); flint_printf("\n"); */
      res = 1;
    }
  else
    {
      arf_ui_div(next, 1, next, prec, ARF_RND_NEAR);
    }
  return res;
}

static slong
cont_frac_max_nb_steps(const acb_t x, slong prec)
{
  return prec/2;
}

static void
cont_frac_get_fmpq(fmpq_t c, fmpz* r_vec, slong nb_steps)
{
  slong k;
  fmpq_zero(c);
  fmpq_add_fmpz(c, c, &r_vec[nb_steps-1]);
  for (k = nb_steps-2; k >= 0; k--)
    {
      fmpq_inv(c, c);
      fmpq_add_fmpz(c, c, &r_vec[k]);
    }
}

int acb_rationalize(fmpq_t c, fmpz_t den, const acb_t x,
		    const fmpz_t probable_den, slong prec)
{
  acb_t z;
  arf_t current;
  mpz_t n, d;
  fmpz_t d_fmpz;
  int res = 1;
  int stop = 0;
  slong max_steps = cont_frac_max_nb_steps(x, prec);
  fmpz* r_vec;
  slong k;
  slong tol_exp;

  acb_init(z);
  arf_init(current);
  mpz_init(n);
  mpz_init(d);
  fmpz_init(d_fmpz);
  r_vec = _fmpz_vec_init(max_steps);
  
  acb_mul_fmpz(z, x, probable_den, prec);
  
  if (!arb_contains_zero(acb_imagref(z)))
    {      
      flint_printf("(acb_rationalize) Error: contains no real number\n");
      acb_printd(x, 10); flint_printf("\n");
      fflush(stdout);
      flint_abort();
    }
  
  if (mag_cmp_2exp_si(arb_radref(acb_realref(z)), 0) > 0)
    {
      /* Too imprecise */
      res = 0;
    }
  
  if (res)
    {
      arf_set(current, arb_midref(acb_realref(z)));
      k = 0;
      while (!stop && (k < max_steps))
	{
	  if (k < 2)
	    {
	      tol_exp = -prec/2;
	    }
	  else
	    {
	      tol_exp = -prec/8;
	    }
	  stop = cont_frac_step(&r_vec[k], current, current, prec, tol_exp);
	  k++;
	}
      if (k == max_steps)
	{
	  res = 0;
	}
      else
	{
	  res = 1;
	  cont_frac_get_fmpq(c, r_vec, k);
	  fmpq_div_fmpz(c, c, probable_den);
	  fmpq_get_mpz_frac(n, d, c);
	  fmpz_set_mpz(d_fmpz, d);
	  fmpz_lcm(den, d_fmpz, probable_den);
	}
    }

  acb_clear(z);
  arf_clear(current);
  mpz_clear(n);
  mpz_clear(d);
  fmpz_clear(d_fmpz);
  _fmpz_vec_clear(r_vec, max_steps);
  return res;
}
