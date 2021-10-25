
#include "modular.h"

int hilbert_modeq_gundlach_round(fmpz_poly_t num1, fmpz_poly_t num2,
				 fmpz_t den, const acb_poly_t num1_acb,
				 const acb_poly_t num2_acb,
				 const acb_t den_acb, slong ell, slong delta)
{  
  acb_t coeff;
  arf_t radius, max_radius;
  fmpz_t rd;
  slong radius_prec = 100;
  slong radius_bits;
  slong d = 2*hilbert_nb_cosets(ell, delta);
  int res = 1;

  acb_init(coeff);
  fmpz_init(rd);
  arf_init(radius);
  arf_init(max_radius);
	  
  fmpz_poly_zero(num1);
  fmpz_poly_zero(num2);
  
  arf_zero(max_radius);
  acb_get_rad_ubound_arf(radius, den_acb, radius_prec);
  arf_max(max_radius, max_radius, radius);
  
  res = siegel_modeq_round_coeff(den, den_acb);

  if (res)
    {
      res = siegel_modeq_round_poly(num1, radius, num1_acb, d);
      arf_max(max_radius, max_radius, radius);      
    }
  if (res)
    {
      res = siegel_modeq_round_poly(num2, radius, num2_acb, d-1);
      arf_max(max_radius, max_radius, radius);      
    }
  
  radius_bits = arf_abs_bound_lt_2exp_si(max_radius);
  if (res) flint_printf("(hilbert_modeq_gundlach_round) Excess precision: %wd\n", -radius_bits);

  acb_clear(coeff);
  fmpz_clear(rd);
  arf_clear(radius);
  arf_clear(max_radius);
  return res;
}
