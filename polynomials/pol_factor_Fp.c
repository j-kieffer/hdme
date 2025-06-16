#include "polynomials.h"
#include <flint/fmpz_mod_poly_factor.h>

void pol_factor_Fp(slong* nb_factors, fmpz_mod_poly_struct* factors, slong* exps,
		   const fmpz_mod_poly_t pol, const fmpz_mod_ctx_t ctx)
{
  fmpz_mod_poly_factor_t fac;
  slong k;
  int v = get_pol_verbose();

  fmpz_mod_poly_factor_init(fac, ctx);
  if (v) flint_printf("(pol_factor_Fp) Factoring...\n");
  fmpz_mod_poly_factor(fac, pol, ctx);
  *nb_factors = fac->num;

  for (k = 0; k < *nb_factors; k++)
    {
      fmpz_mod_poly_set(&factors[k], (&(fac->poly))[k], ctx);
      exps[k] = (fac->exp)[k];
    }
  
  if (v)
    {
      flint_printf("(pol_factor_Fp) Factorization pattern: %wd",
		   fmpz_mod_poly_degree(&factors[0], ctx));
      for (k = 1; k < *nb_factors; k++)
	{
	  flint_printf(" + %wd", fmpz_mod_poly_degree(&factors[k], ctx));
	}
      flint_printf("\n");
    }
      
  fmpz_mod_poly_factor_clear(fac, ctx);
}
