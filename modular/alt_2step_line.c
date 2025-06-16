
#include <flint/nmod_mat.h>
#include "modular.h"

int alt_2step_line(fmpz_mat_t L, slong* indices, slong nb, const hecke_t H)
{
  slong ell = hecke_ell(H);
  nmod_mat_t current;
  nmod_mat_t temp;
  nmod_mat_t rectangle;
  
  slong k, j;
  slong dim = 4;

  nmod_mat_init(current, 4, 4, ell);
  nmod_mat_init(temp, 4, 4, ell);
  nmod_mat_init(rectangle, 8, 4, ell);
  for (k = 0; k < 4; k++) nmod_mat_entry(current, k, k) = 1;
  

  for (k = 0; k < nb; k++)
    {
      j = indices[k];
      fmpz_mat_get_nmod_mat(temp, hecke_coset(H, j));
      nmod_mat_nullspace(temp, temp);
      /* Get duals */
      nmod_mat_transpose(current, current);
      nmod_mat_transpose(temp, temp);
      nmod_mat_nullspace(current, current);
      nmod_mat_nullspace(temp, temp);      
      nmod_mat_transpose(current, current);
      nmod_mat_transpose(temp, temp);
      /* Stack them */
      nmod_mat_concat_vertical(rectangle, current, temp);
      dim = nmod_mat_nullspace(current, rectangle);
    }

  fmpz_mat_set_nmod_mat(L, current);

  nmod_mat_clear(current);
  nmod_mat_clear(temp);
  nmod_mat_clear(rectangle);  
  return (dim == 1);
}
