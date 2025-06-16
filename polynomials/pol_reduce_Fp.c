#include "polynomials.h"
#include <flint/fmpz_mod.h>

int pol_reduce_Fp(fmpz_mod_poly_t red, const fmpz_poly_t num,
		  const fmpz_t den, const fmpz_mod_ctx_t ctx)
{
  int res;
  fmpz_t one;
  fmpz_t invden;
  
  fmpz_init(one);
  fmpz_init(invden);
  
  fmpz_one(one);
  res = fmpz_mod_divides(invden, one, den, ctx);

  if (res)
    {
      fmpz_mod_poly_set_fmpz_poly(red, num, ctx);
      fmpz_mod_poly_scalar_mul_fmpz(red, red, invden, ctx);
    }

  fmpz_clear(one);
  fmpz_clear(invden);
  return res;  
}
