#include "polynomials.h"
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>

int acb_poly_rationalize(fmpq_poly_t pol, const acb_poly_t pol_acb,
			 slong degree, slong prec)
{  
  fmpq_t c;
  fmpz_t c_num;
  fmpz_t current_den;
  
  acb_t x;
  int success = 1;
  slong k;
  
  fmpq_init(c);
  fmpz_init(c_num);
  fmpz_init(current_den);
  acb_init(x);
  
  fmpz_one(current_den);

  for (k = 0; (k <= degree) && success; k++)
    {      
      acb_poly_get_coeff_acb(x, pol_acb, k);
      success = acb_rationalize(c, current_den, x, current_den, prec);
      if (success)
	{
	  fmpq_poly_set_coeff_fmpq(pol, k, c);
	}
    }
  
  fmpq_clear(c);
  fmpz_clear(c_num);
  fmpz_clear(current_den);
  acb_clear(x);

  return success;
}
