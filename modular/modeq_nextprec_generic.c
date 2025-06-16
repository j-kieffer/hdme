
#include <math.h>
#include "modular.h"

slong modeq_nextprec_generic(slong current_prec)
{
  slong res = ceil(MODEQ_MUL_PREC * current_prec);
  return 100 * (res/100 + 1);
}
