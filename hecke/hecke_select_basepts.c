
#include "hecke.h"

int hecke_select_basepts(acb_mat_struct* pts, acb_mat_t basis_inv,
			 slong wt, slong prec)
{
  flint_rand_t state;
  slong k, j;
  slong jmax = HECKE_SELECT_BASEPTS_TRIALS;
  slong weights[4] = IGUSA_WEIGHTS;
  slong nb = cov_nb_monomials(wt, 4, weights);
  int res = 0;
  
  flint_rand_init(state);
  
  for (j = 0; j < jmax; j++)
    {
      for (k = 0; k < nb; k++)
	{
	  siegel_fundamental_domain_randtest(&pts[k], state, prec);
	}
      res = hecke_basis_matrix(basis_inv, nb, pts, wt, prec);
      if (!res) continue;
      res = acb_mat_inv(basis_inv, basis_inv, prec);
      if (res) break;
    }

  flint_rand_clear(state);
  return res;  
}
