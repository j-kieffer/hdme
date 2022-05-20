
#include "igusa.h"

int cov_find_rescaling(acb_t scal, acb_srcptr I, fmpz* S,
		       slong nb, slong* weights, slong prec)
{
  slong wt;
  slong* exponents;  
  acb_t r, c;
  acb_ptr rescale;
  slong j;
  int res = 1;

  exponents = flint_malloc(nb * sizeof(slong));
  acb_init(r);
  acb_init(c);
  rescale = _acb_vec_init(nb);

  /* Compute combination of I/S with minimum weight */
  cov_min_weight_combination(&wt, exponents, S, nb, weights);
  acb_one(c);
  for (j = 0; j < nb; j++)
    {
      acb_pow_si(r, &I[j], exponents[j], prec);
      acb_mul(c, c, r, prec);
      acb_set_fmpz(r, &S[j]);
      acb_pow_si(r, r, exponents[j], prec);
      acb_div(c, c, r, prec);
    }
  borchardt_root_ui(c, c, wt, prec);

  /* Check that rescaling works */
  for (j = 0; j < nb; j++) acb_set_fmpz(&rescale[j], &S[j]);
  cov_rescale(rescale, rescale, c, nb, weights, prec);
  
  for (j = 0; j < 4; j++)
    {
      if (!acb_overlaps(&rescale[j], &I[j]))
	{
	  res = 0;
	  break;
	}
    }

  flint_free(exponents);
  acb_clear(r);
  acb_clear(c);
  _acb_vec_clear(rescale, nb);
  return res;
}
