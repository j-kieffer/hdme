
#include <flint/nmod_mat.h>
#include "hecke.h"

int siegel_T1_coset_contains_line_dual(const fmpz_mat_t m,
				       const fmpz_mat_t L, slong ell)
{
  fmpz_mat_t J;
  nmod_mat_t dual;
  nmod_mat_t red;
  int res;

  fmpz_mat_init(J, 4, 4);
  nmod_mat_init(red, 4, 4, ell);
  nmod_mat_init(dual, 4, 4, ell);

  fmpz_mat_J(J);

  fmpz_mat_get_nmod_mat(dual, L);
  nmod_mat_transpose(dual, dual);
  fmpz_mat_get_nmod_mat(red, J);
  nmod_mat_mul(dual, dual, red);  
  nmod_mat_nullspace(dual, dual);

  fmpz_mat_get_nmod_mat(red, m);
  nmod_mat_mul(red, red, dual);
  res = nmod_mat_is_zero(red);

  fmpz_mat_clear(J);
  nmod_mat_clear(red);
  nmod_mat_clear(dual);
  return res;
}
