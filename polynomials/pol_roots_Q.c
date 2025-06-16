#include "polynomials.h"
#include <flint/fmpq.h>

void pol_roots_Q(slong* nb_roots, fmpq* roots, slong* mults,
		 const fmpz_poly_t pol)
{
  slong deg = fmpz_poly_degree(pol);
  slong nb_factors = 0;
  slong* exps;
  fmpz_poly_struct* factors;
  fmpz_t num, den;
  slong k;

  factors = flint_malloc(deg * sizeof(fmpz_poly_struct));
  for (k = 0; k < deg; k++) fmpz_poly_init(&factors[k]);
  
  fmpz_init(num);
  fmpz_init(den);
  exps = flint_malloc(deg * sizeof(slong));

  pol_factor_Q(&nb_factors, factors, exps, pol);
  *nb_roots = 0;
  for (k = 0; k < nb_factors; k++)
    {
      if (fmpz_poly_degree(&factors[k]) == 1)
	{
	  fmpz_poly_get_coeff_fmpz(num, &factors[k], 0);
	  fmpz_poly_get_coeff_fmpz(den, &factors[k], 1);
	  fmpz_neg(num, num);
	  fmpq_set_fmpz_frac(&roots[*nb_roots], num, den);
	  mults[*nb_roots] = exps[*nb_roots];
	  (*nb_roots)++;
	}
    }

  for (k = 0; k < deg; k++) fmpz_poly_clear(&factors[k]);
  flint_free(factors);
  fmpz_clear(num);
  fmpz_clear(den);
  flint_free(exps);
}
