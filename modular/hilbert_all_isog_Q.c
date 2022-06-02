
#include "modular.h"

int hilbert_all_isog_Q(slong* nb_roots, fmpz* all_I, fmpz* I,
		       slong ell, slong delta)
{
  modeq_t E;
  modeq_ctx_t ctx;
  int res;

  modeq_init(E);
  modeq_ctx_init(ctx);
  
  res = hilbert_modeq_eval(E, ctx, I, ell, delta);
  if (res) modeq_all_isog_Q(nb_roots, all_I, E, ctx);
  
  modeq_clear(E);
  modeq_ctx_clear(ctx);
  return res;
}
