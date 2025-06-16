#include <stdlib.h>

#include "hilbert.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("humbert_parametrize....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 1 * flint_test_multiplier(); iter++)
    {
      slong delta;
      fmpq* rs;
      acb_ptr rs_acb;
      acb_ptr I, I_test;
      acb_ptr t;
      acb_mat_t tau;
      fmpz_mat_t eta;
      slong rs_bits = 5 + n_randint(state, 10);
      slong prec = 1000;
      slong weights[4] = IGUSA_WEIGHTS;
      slong k;
      int res;
      slong delta_max = 100;
      int v = 0;
      
      rs = _fmpq_vec_init(2);
      rs_acb = _acb_vec_init(2);
      I = _acb_vec_init(4);
      I_test = _acb_vec_init(4);
      t = _acb_vec_init(2);
      acb_mat_init(tau, 2, 2);
      fmpz_mat_init(eta, 4, 4);
      
      for (delta = 5; delta < delta_max; delta++)
	{
	  if (hilbert_is_fundamental(delta))
	    {
	      for (k = 0; k < 2; k++)
		{
		  fmpq_randbits(&rs[k], state, rs_bits);
		}
	      if (v)
		{
		  flint_printf("delta = %wd; parameters are\n", delta);
		  fmpq_print(&rs[0]); flint_printf("\n");
		  fmpq_print(&rs[1]); flint_printf("\n");
		}
	      acb_set_fmpq(&rs_acb[0], &rs[0], prec);
	      acb_set_fmpq(&rs_acb[1], &rs[1], prec);
	      humbert_parametrize(I, rs_acb, delta, prec);
	      res = tau_from_igusa(tau, I, prec);

	      if (!res)
		{
		  flint_printf("FAIL (period matrix)\n");
		  for (k = 0; k < 2; k++)
		    {
		      fmpq_print(&rs[k]); flint_printf("\n");
		    }
		  for (k = 0; k < 4; k++)
		    {
		      acb_printd(&I[k], 30); flint_printf("\n");
		    }
		  fflush(stdout);
		  flint_abort();
		}

	      res = hilbert_inverse(t, eta, tau, delta, prec);
	      
	      if (!res)
		{
		  flint_printf("FAIL (Hilbert inversion)\n");
		  flint_printf("delta = %wd\n", delta);
		  for (k = 0; k < 2; k++)
		    {
		      fmpq_print(&rs[k]); flint_printf("\n");
		    }
		  fflush(stdout);
		  flint_abort();
		}

	      hilbert_map(tau, t, delta, prec);
	      igusa_from_tau(I_test, tau, prec);
	      
	      if (cov_distinct(I, I_test, 4, weights, prec))
		{
		  flint_printf("FAIL (overlap)\n");
		  for (k = 0; k < 2; k++)
		    {
		      fmpq_print(&rs[k]); flint_printf("\n");
		    }
		  fflush(stdout);
		  flint_abort();
		}
	    }
	}

      _fmpq_vec_clear(rs, 2);
      _acb_vec_clear(rs_acb, 2);
      _acb_vec_clear(I, 4);
      _acb_vec_clear(I_test, 4);
      _acb_vec_clear(t, 2);
      fmpz_mat_clear(eta);
    }

  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
