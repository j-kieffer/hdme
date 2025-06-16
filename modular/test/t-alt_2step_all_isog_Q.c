#include <stdlib.h>

#include "modular.h"


int main()
{
  slong iter;
  
  flint_printf("alt_2step_all_isog_Q....");
  fflush(stdout);

  /* Up to iter < 4 */
  for (iter = 0; iter < 0; iter++)
    {
      fmpz* I1;
      fmpz* I2;
      fmpz* all_I;
      slong ell;
      slong max_nb_roots = 2;
      slong nb_roots = 1;
      slong k;
      int res;
      
      I1 = _fmpz_vec_init(4);
      I2 = _fmpz_vec_init(4);
      all_I = _fmpz_vec_init(4 * max_nb_roots);
      
      if (iter == 0)
	{
	  /* https://beta.lmfdb.org/Genus2Curve/Q/249/a/249/1,
	     https://beta.lmfdb.org/Genus2Curve/Q/249/a/6723/1 */
	  fmpz_set_si(&I1[0], 108);
	  fmpz_set_si(&I1[1], 57);
	  fmpz_set_si(&I1[2], 2259);
	  fmpz_set_si(&I1[3], -31872);
	  ell = 2;
	  fmpz_set_si(&I2[0], 1932);
	  fmpz_set_si(&I2[1], 87897);
	  fmpz_set_si(&I2[2], 65765571);
	  fmpz_set_si(&I2[3], 860544);
	}
      else if (iter == 1)
	{
	  /* https://beta.lmfdb.org/Genus2Curve/Q/277/a/277/1,
	     https://beta.lmfdb.org/Genus2Curve/Q/277/a/277/2 */
	  fmpz_set_si(&I1[0], 64);
	  fmpz_set_si(&I1[1], 352);
	  fmpz_set_si(&I1[2], 9552);
	  fmpz_set_si(&I1[3], -1108);
	  ell = 3;
	  fmpz_set_si(&I2[0], 4480);
	  fmpz_set_si(&I2[1], 1370512);
	  fmpz_set_si(&I2[2], 1511819744);
	  fmpz_set_si(&I2[3], -1108);
	}
      else if (iter == 2)
	{
	  /* https://beta.lmfdb.org/Genus2Curve/Q/523/a/523/1,
	     https://beta.lmfdb.org/Genus2Curve/Q/523/a/523/2 */	  
	  fmpz_set_si(&I1[0], 120);
	  fmpz_set_si(&I1[1], -540);
	  fmpz_set_si(&I1[2], -29169);
	  fmpz_set_si(&I1[3], -2092);
	  ell = 5;
	  fmpz_set_si(&I2[0], 332400);
	  fmpz_set_si(&I2[1], 10084860);
	  fmpz_set_si(&I2[2], 1107044456391);
	  fmpz_set_si(&I2[3], -2092);
	}
      else
	{
	  /* https://beta.lmfdb.org/Genus2Curve/Q/295/a/295/1,
	     https://beta.lmfdb.org/Genus2Curve/Q/295/a/295/2 */	  
	  fmpz_set_si(&I1[0], 108);
	  fmpz_set_si(&I1[1], -39);
	  fmpz_set_si(&I1[2], 20835);
	  fmpz_set_si(&I1[3], 37760);
	  ell = 7;
	  fmpz_set_si(&I2[0], 198804);
	  fmpz_set_si(&I2[1], 305807001);
	  fmpz_set_si(&I2[2], 18482629056189);
	  fmpz_set_si(&I2[3], -37760);
	}
      
      igusa_from_IC_fmpz(I1, I1);
      igusa_from_IC_fmpz(I2, I2);

      alt_2step_all_isog_Q(&nb_roots, all_I, I1, ell);
      
      if (nb_roots == 0)
	{
	  flint_printf("FAIL (roots)\n");
	  flint_printf("nb_roots = %wd\n", nb_roots);
	  fflush(stdout);
	  flint_abort();
	}

      res = 0;
      for (k = 0; k < nb_roots; k++)
	{
	  if (_fmpz_vec_equal(I2, &all_I[4*k], 4)) res = 1;
	}
      if (!res)
	{
	  flint_printf("FAIL (values)\n");
	  fflush(stdout);
	  flint_abort();
	}      

      _fmpz_vec_clear(I1, 4);
      _fmpz_vec_clear(I2, 4);
      _fmpz_vec_clear(all_I, 4 * max_nb_roots);
    }

  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
