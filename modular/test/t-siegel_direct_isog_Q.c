#include <stdlib.h>

#include "modular.h"
#include "hecke.h"

int main()
{
  slong iter;

  flint_printf("siegel_direct_isog_Q....");
  fflush(stdout);

  /* Up to iter < 3 */
  for (iter = 0; iter < 3; iter++) {
    fmpz* I1;
    fmpz* all_I;
    slong ell;
    slong max_nb_roots = 8;
    slong nb_roots = 1;

    I1 = _fmpz_vec_init(4);
    all_I = _fmpz_vec_init(4 * max_nb_roots);

    switch(iter) {
      case 0:
        /* https://beta.lmfdb.org/Genus2Curve/Q/990/a/8910/1 */
        fmpz_set_si(&I1[0], 3268);
        fmpz_set_si(&I1[1], 252577);
        fmpz_set_si(&I1[2], 318023313);
        fmpz_set_si(&I1[3], 1140480);
        ell = 2;
        break;

      case 1:
        /* https://beta.lmfdb.org/Genus2Curve/Q/676/b/17576/1 */
        fmpz_set_si(&I1[0], 1244);
        fmpz_set_si(&I1[1], 1249);
        fmpz_set_si(&I1[2], 129167);
        fmpz_set_si(&I1[3], 2249728);
        ell = 3;
        break;

      case 2:
        fmpz_set_str(&I1[0], "117322407584104", 10);
        fmpz_set_str(&I1[1], "166443185276308771469224348", 10);
        fmpz_set_str(&I1[2], "7813005098059830157559238453206675869220", 10);
        fmpz_set_str(&I1[3], "36493703564806262820376912749375714135215499089814947142944384", 10);
        ell = 3;
        break;

      default:
        flint_printf("iteration = %wd has not been initialized\n", iter);
        flint_abort();
  }

    igusa_from_IC_fmpz(I1, I1);
    siegel_direct_isog_Q(&nb_roots, all_I, I1, ell);

    if (nb_roots == 0) {
      flint_printf("FAIL (roots)\n");
      flint_printf("nb_roots = %wd\n", nb_roots);
      fflush(stdout);
      flint_abort();
    }
      _fmpz_vec_clear(I1, 4);
      _fmpz_vec_clear(all_I, 4 * max_nb_roots);
    }

  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
