#include "siegel.h"
#include <flint/acb.h>

void
acb_mat_get_real(arb_mat_t re, const acb_mat_t z)
{
  slong nrows = acb_mat_nrows(z);
  slong ncols = acb_mat_ncols(z);
  slong i, j;

  for (i = 0; i < nrows; i++)
    {
      for (j = 0; j < ncols; j++)
	{
	  acb_get_real(arb_mat_entry(re, i, j), arb_mat_entry(z, i, j));
	}
    }
}
