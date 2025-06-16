#include <stdlib.h>

#include "hecke.h"

int main()
{
  slong wt;

  flint_printf("hecke_charpoly....");
  fflush(stdout);

  for (wt = 4; wt < 10; wt += 2) {
    slong ell = 2;
    fmpz_poly_t pol;
    fmpz_t eig;
    int print = 1;

    fmpz_poly_init(pol);
    hecke_charpoly(pol, ell, wt);
    fmpz_init(eig);

    if (print) {
      flint_printf("Charpoly for weight %wd:\n", wt);
      fmpz_poly_print_pretty(pol, "x"); flint_printf("\n");
    }

    hecke_eigenvalue_eisenstein(eig, wt, ell);
    fmpz_poly_evaluate_fmpz(eig, pol, eig);

    if (!fmpz_is_zero(eig)) {
      flint_printf("FAIL\n");
      fflush(stdout);
      flint_abort();
    }

    fmpz_poly_clear(pol);
  }

  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}


