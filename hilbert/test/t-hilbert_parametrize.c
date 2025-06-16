#include <stdlib.h>

#include "hilbert.h"


int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("hilbert_parametrize....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 2 * flint_test_multiplier(); iter++)
    {
      slong delta;
      fmpq* rs;
      acb_ptr rs_acb;
      acb_ptr I;
      
      slong rs_bits = 5 + n_randint(state, 10);
      slong prec = 1000 + n_randint(state, 1000);
      slong delta_max = 18;
      slong k;

      acb_mat_t tau;
      acb_ptr t;
      fmpz_mat_t m;
      int res;

      rs = _fmpq_vec_init(2);
      rs_acb = _acb_vec_init(2);
      I = _acb_vec_init(4);
      acb_mat_init(tau, 2, 2);
      fmpz_mat_init(m, 4, 4);
      t = _acb_vec_init(3);
	
      for (delta = 5; delta < delta_max; delta++)
	{
	  if (hilbert_is_fundamental(delta))
	    {
	      for (k = 0; k < 2; k++)
		{
		  fmpq_randbits(&rs[k], state, rs_bits);
		}
	      acb_set_fmpq(&rs_acb[0], &rs[0], prec);
	      acb_set_fmpq(&rs_acb[1], &rs[1], prec);
	      hilbert_parametrize(I, rs_acb, delta, prec);

	      res = tau_from_igusa(tau, I, prec);
	      if (!res)
		{
		  flint_printf("FAIL (tau)\n");
		  fflush(stdout);
		  flint_abort();
		}
	      res = hilbert_inverse(t, m, tau, delta, prec);
	      if (!res)
		{
		  flint_printf("FAIL (inverse)\n");
		  flint_printf("delta = %wd, prec = %wd\n", delta, prec);
		  fflush(stdout);
		  flint_abort();
		}
	    }
	}

      _fmpq_vec_clear(rs, 2);
      _acb_vec_clear(rs_acb, 2);
      _acb_vec_clear(I, 4);
      acb_mat_clear(tau);
      fmpz_mat_clear(m);
      _acb_vec_clear(t, 2);
    }

  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
