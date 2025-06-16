#include <stdlib.h>

#include "hecke.h"

int main()
{
  slong iter;
  
  flint_printf("hecke_make_integral....");
  fflush(stdout);

  for (iter = 0; iter < 1; iter++)
    {
      fmpz* I;
      slong ell = 2;
      hecke_t H;
      slong nb = siegel_nb_T1_cosets(ell);
      slong k, j;
      int print = 0;
      slong prec = 1000;
      int res = 0;
      int is_int;

      I = _fmpz_vec_init(4);
      ell = 2;
      hecke_init(H, nb);
            
      fmpz_set_si(&I[0], 108);
      fmpz_set_si(&I[1], 57);
      fmpz_set_si(&I[2], 2259);
      fmpz_set_si(&I[3], -31872);
      igusa_from_IC_fmpz(I, I);
      
      hecke_set_I_fmpz(H, I, prec);
      hecke_collect_T1(H, ell, prec);
      hecke_make_integral(H, I, prec);

      for (k = 0; k < nb; k++)
	{
	  is_int = 1;
	  if (print)
	    {
	      flint_printf("Is this integral?\n");
	      acb_printd(&hecke_I_norm(H, k)[0], 20); flint_printf("\n");
	      acb_printd(&hecke_I_norm(H, k)[1], 20); flint_printf("\n");
	      acb_printd(&hecke_I_norm(H, k)[2], 20); flint_printf("\n");
	      acb_printd(&hecke_I_norm(H, k)[3], 20); flint_printf("\n");
	    }
	  for (j = 0; j < 4; j++)
	    {
	      is_int = is_int && acb_contains_int(&hecke_I_norm(H, k)[j]);
	    }
	  if (is_int)
	    {
	      res = 1;
	      break;
	    }
	}

      if (!res)
	{
	  flint_printf("FAIL\n");
	  fflush(stdout);
	  flint_abort();
	}

      _fmpz_vec_clear(I, 4);
      hecke_clear(H);
    }
  
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

      
