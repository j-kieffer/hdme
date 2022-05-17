
#include "igusa.h"

void cov_from_IC(acbptr I, acb_srcptr IC, slong prec)
{
  acb_t I12, I6prime;

  acb_init(I12);
  acb_init(I6prime);

  acb_mul(I12, &IC[0], &IC[3], prec);
  igusa_I6prime(I6prime, IC, prec);

  acb_set(cov_I4(I), &IC[1]);
  acb_set(cov_I6prime(I), I6prime);
  acb_set(cov_I10(I), &IC[3]);
  acb_set(cov_I12(I), I12);

  acb_clear(I12);
  acb_clear(I6prime);
}