
#include <string.h>
#include "hilbert.h"

void hilbert_parametrize(acb_ptr I, acb_srcptr rs, slong delta, slong prec)
{
  acb_t x, y, z;
  acb_ptr gh;
  fmpq_mpoly_t pol;
  fmpq_mpoly_ctx_t ctx;
  char** vars;
  acb_ptr vals;

  fmpq_mpoly_ctx_init(ctx, 2, ORD_LEX);
  fmpq_mpoly_init(pol, ctx);
  acb_init(x);
  acb_init(y);
  acb_init(z);
  gh = _acb_vec_init(2);
  vars = hdme_data_vars_init(2);
  vals = _acb_vec_init(2);

  _acb_vec_set(vals, rs, 2);

  if (delta == 5)
    {
      strcpy(vars[0], "m");
      strcpy(vars[1], "n");
      hilbert_get_mpoly(pol, (const char**) vars, "g", delta, ctx);
      hdme_data_evaluate_acb(&gh[0], pol, vals, ctx, prec);

      strcpy(vars[1], "g");
      acb_set(&vals[1], &gh[0]);
      
      hilbert_get_mpoly(pol, (const char**) vars, "h", delta, ctx);
      hdme_data_evaluate_acb(&gh[1], pol, vals, ctx, prec);      
    }

  else if (delta == 8)
    {      
      strcpy(vars[0], "m");
      strcpy(vars[1], "n");
      hilbert_get_mpoly(pol, (const char**) vars, "rnum", delta, ctx);
      hdme_data_evaluate_acb(&gh[0], pol, vals, ctx, prec);      
      hilbert_get_mpoly(pol, (const char**) vars, "rden", delta, ctx);
      hdme_data_evaluate_acb(x, pol, vals, ctx, prec);
      acb_div(&gh[0], &gh[0], x, prec);

      strcpy(vars[1], "r");
      acb_set(&vals[1], &gh[0]);
      
      hilbert_get_mpoly(pol, (const char**) vars, "snum", delta, ctx);
      hdme_data_evaluate_acb(&gh[1], pol, vals, ctx, prec);      
      hilbert_get_mpoly(pol, (const char**) vars, "sden", delta, ctx);
      hdme_data_evaluate_acb(x, pol, vals, ctx, prec);
      acb_div(&gh[1], &gh[1], x, prec);
    }

  else if (delta == 12)
    {
      strcpy(vars[0], "f");
      strcpy(vars[1], "g");
      hilbert_get_mpoly(pol, (const char**) vars, "enum", delta, ctx);
      hdme_data_evaluate_acb(&gh[0], pol, vals, ctx, prec);      
      hilbert_get_mpoly(pol, (const char**) vars, "eden", delta, ctx);
      hdme_data_evaluate_acb(x, pol, vals, ctx, prec);
      acb_div(&gh[0], &gh[0], x, prec); /* &gh[0] is e */

      acb_set(&gh[1], &vals[0]);
    }

  else if (delta == 13)
    {    
      strcpy(vars[0], "m");
      strcpy(vars[1], "n");
      hilbert_get_mpoly(pol, (const char**) vars, "g2num", delta, ctx);
      hdme_data_evaluate_acb(x, pol, vals, ctx, prec);      
      hilbert_get_mpoly(pol, (const char**) vars, "g2den", delta, ctx);
      hdme_data_evaluate_acb(y, pol, vals, ctx, prec);
      acb_div(x, x, y, prec); /* x is g2 */
      acb_div_si(y, x, 54, prec); /* y is g1 */
      acb_mul(z, &vals[0], y, prec); /* z is h1 */
      acb_set_si(&gh[1], 64);
      acb_div_si(&gh[1], &gh[1], 27, prec);
      acb_add(&gh[1], &gh[1], z, prec);
      acb_set_si(&gh[0], -1);
      acb_div_si(&gh[0], &gh[0], 54, prec);
      acb_add(&gh[0], &gh[0], y, prec);
    }

  else if (delta == 17)
    {
      strcpy(vars[0], "m");
      strcpy(vars[1], "n");
      hilbert_get_mpoly(pol, (const char**) vars, "gnum", delta, ctx);
      hdme_data_evaluate_acb(&gh[0], pol, vals, ctx, prec);      
      hilbert_get_mpoly(pol, (const char**) vars, "gden", delta, ctx);
      hdme_data_evaluate_acb(x, pol, vals, ctx, prec);
      acb_div(&gh[0], &gh[0], x, prec);

      strcpy(vars[1], "g");
      acb_set(&vals[1], &gh[0]);

      hilbert_get_mpoly(pol, (const char**) vars, "h", delta, ctx);
      hdme_data_evaluate_acb(&gh[1], pol, vals, ctx, prec); 
    }

  else
    {
      flint_printf("(hilbert_parametrize) Error: no Hilbert parametrization for discriminant %wd\n", delta);
      fflush(stdout);
      flint_abort();
    }

  humbert_parametrize(I, gh, delta, prec);

  fmpq_mpoly_clear(pol, ctx);
  fmpq_mpoly_ctx_clear(ctx);
  acb_clear(x);
  acb_clear(y);
  acb_clear(z);
  _acb_vec_clear(gh, 2);
  hdme_data_vars_clear(vars, 2);
  _acb_vec_clear(vals, 2);
}
