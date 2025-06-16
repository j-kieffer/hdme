#include "polynomials.h"
#include <flint/fmpz_poly_factor.h>

void pol_factor_Q(slong* nb_factors, fmpz_poly_struct* factors, slong* exps,
		  const fmpz_poly_t pol)
{
  fmpz_poly_factor_t fac;
  slong k;
  int v = get_pol_verbose();

  fmpz_poly_factor_init(fac);
  
  if (v) flint_printf("(pol_factor_Q) Factoring...\n");
  
  fmpz_poly_factor(fac, pol);
  *nb_factors = fac->num;

  for (k = 0; k < *nb_factors; k++)
    {
      fmpz_poly_set(&factors[k], &((fac->p)[k]));
      exps[k] = (fac->exp)[k];
    }
  
  if (v)
    {
      flint_printf("(pol_factor_Q) Factorization pattern: %wd",
		   fmpz_poly_degree(&factors[0]));
      for (k = 1; k < *nb_factors; k++)
	{
	  flint_printf(" + %wd", fmpz_poly_degree(&factors[k]));
	}
      flint_printf("\n");
    }
  
  fmpz_poly_factor_clear(fac);
}
