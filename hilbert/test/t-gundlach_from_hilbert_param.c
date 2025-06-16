#include <stdlib.h>

#include "hilbert.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("gundlach_from_hilbert_param....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
      slong delta = 5;
      slong prec = 1000 + n_randint(state, 2000);
      slong weights[3] = GUNDLACH_WEIGHTS_5;
      fmpq* mn;
      slong mn_bits = 1 + n_randint(state, 10);
      acb_ptr rs;
      acb_ptr I;
      fmpz* G;
      acb_ptr G_acb;
      acb_ptr G_test;
      slong k;

      mn = _fmpq_vec_init(2);
      rs = _acb_vec_init(2);
      I = _acb_vec_init(4);
      G = _fmpz_vec_init(3);
      G_acb = _acb_vec_init(3);
      G_test = _acb_vec_init(3);

      for (k = 0; k < 2; k++) fmpq_randbits(&mn[k], state, mn_bits);
      gundlach_from_hilbert_param(G, mn, delta);
      for (k = 0; k < 3; k++) acb_set_fmpz(&G_acb[k], &G[k]);
      
      acb_set_fmpq(&rs[0], &mn[0], prec);
      acb_set_fmpq(&rs[1], &mn[1], prec);
      
      hilbert_parametrize(I, rs, delta, prec);
      gundlach_from_igusa(G_test, I, delta, prec);

      if (cov_distinct(G_acb, G_test, 3, weights, prec))
	{
	  flint_printf("FAIL\n");
	  for (k = 0; k < 3; k++)
	    {
	      fmpz_print(&G[k]); flint_printf("\n");
	      acb_printd(&G_test[k], 300); flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_abort();
	}
      
      _fmpq_vec_clear(mn, 2);
      _acb_vec_clear(rs, 2);
      _acb_vec_clear(I, 4);
      _fmpz_vec_clear(G, 3);
      _acb_vec_clear(G_acb, 3);
      _acb_vec_clear(G_test, 3);
    }  

  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}



      
      
