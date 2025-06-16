#include "polynomials.h"
#include <flint/fmpz_mod_poly_factor.h>
#include <flint/fmpz_mod.h>

void pol_roots_Fp(slong* nb_roots, fmpz* roots, slong* mults,
		  const fmpz_mod_poly_t pol, const fmpz_mod_ctx_t ctx)
{
  fmpz_mod_poly_factor_t fac;
  fmpz_mod_poly_t lin;
  fmpz_t num, den;
  slong k;
  int v = get_pol_verbose();

  fmpz_mod_poly_factor_init(fac, ctx);
  fmpz_mod_poly_init(lin, ctx);
  fmpz_init(num);
  fmpz_init(den);
  
  if (v) flint_printf("(modeq_roots_Fp) Finding roots...\n");
  fmpz_mod_poly_roots(fac, pol, 1, ctx);
  *nb_roots = fac->num;
  if (v) flint_printf("(siegel_modeq_roots_Fp) Found %wd.\n", *nb_roots);

  for (k = 0; k < *nb_roots; k++)
    {
      fmpz_mod_poly_set(lin, &((fac->poly)[k]), ctx);
      fmpz_mod_poly_get_coeff_fmpz(num, lin, 0, ctx);
      fmpz_mod_neg(num, num, ctx);
      fmpz_mod_poly_get_coeff_fmpz(den, lin, 1, ctx);
      fmpz_mod_inv(den, den, ctx);
      fmpz_mod_mul(num, num, den, ctx);
      
      fmpz_set(&roots[k], num);
      mults[k] = (fac->exp)[k];
    }
  
  fmpz_mod_poly_factor_clear(fac, ctx);
  fmpz_mod_poly_clear(lin, ctx);
  fmpz_clear(num);
  fmpz_clear(den);
}
