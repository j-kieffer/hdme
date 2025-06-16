
#include <math.h>
#include "hecke.h"

slong hecke_charpoly_nextprec(slong prec)
{
  slong res = ceil(HECKE_CHARPOLY_PREC_MUL * prec);
  return 100*(res/100 + 1);  
}
