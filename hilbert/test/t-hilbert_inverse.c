#include <stdlib.h>

#include "hilbert.h"

int main()
{
  slong iter;
  flint_rand_t state;
  slong delta;
  
  flint_printf("hilbert_inverse....");
  fflush(stdout);
  
  flint_rand_init(state);
  
  for (delta = 5; delta < 100; delta++)
    {
      if (hilbert_is_fundamental(delta))
	{
	  for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
	    {
	      acb_ptr t, t_test;
	      acb_mat_t tau, tau_test;
	      slong prec = 1000 + n_randint(state, 1000);
	      slong m_bits = 4;
	      fmpz_mat_t m;
	      fmpz_mat_t minv;
	      int res;

	      t = _acb_vec_init(2);
	      t_test = _acb_vec_init(2);
	      acb_mat_init(tau, 2, 2);
	      acb_mat_init(tau_test, 2, 2);
	      fmpz_mat_init(m, 4, 4);
	      fmpz_mat_init(minv, 4, 4);
	      
	      hilbert_halfspace_randtest(t, state, prec);
	      hilbert_map(tau, t, delta, prec);
	      res = hilbert_inverse(t_test, minv, tau, delta, prec);

	      fmpz_mat_direct_inv(m, minv);
	      hilbert_map(tau_test, t_test, delta, prec);
	      siegel_transform(tau_test, m, tau_test, prec);

	      if (!res || !acb_mat_overlaps(tau, tau_test))
		{
		  flint_printf("FAIL (wrong preimage)\n");
		  flint_printf("delta = %wd\n", delta);
		  flint_printf("(t1, t2):\n");
		  acb_printd(&t[0], 30); flint_printf("\n");
		  acb_printd(&t[1], 30); flint_printf("\n");
		  acb_mat_printd(tau, 10); flint_printf("\n");
		  flint_printf("Preimage:\n");
		  acb_printd(&t_test[0], 30); flint_printf("\n");
		  acb_printd(&t_test[1], 30); flint_printf("\n");
		  fmpz_mat_print(m);
		  acb_mat_printd(tau_test, 10); flint_printf("\n");
		  fflush(stdout);
		  flint_abort();
		}

	      fmpz_mat_randtest_symplectic(m, state, m_bits);
	      siegel_transform(tau, m, tau, prec);
	      res = hilbert_inverse(t_test, minv, tau, delta, prec);
	      
	      fmpz_mat_direct_inv(minv, minv);
	      hilbert_map(tau_test, t_test, delta, prec);
	      siegel_transform(tau_test, minv, tau_test, prec);
	      
	      if (!res || !acb_mat_overlaps(tau, tau_test))
		{
		  flint_printf("FAIL (wrong preimage)\n");
		  flint_printf("delta = %wd\n", delta);
		  flint_printf("(t1, t2):\n");
		  acb_printd(&t[0], 30); flint_printf("\n");
		  acb_printd(&t[1], 30); flint_printf("\n");
		  fmpz_mat_print(m);
		  acb_mat_printd(tau, 10); flint_printf("\n");
		  flint_printf("Preimage:\n");
		  acb_printd(&t_test[0], 30); flint_printf("\n");
		  acb_printd(&t_test[1], 30); flint_printf("\n");
		  fmpz_mat_print(minv);
		  acb_mat_printd(tau_test, 10); flint_printf("\n");
		  fflush(stdout);
		  flint_abort();
		}

	      _acb_vec_clear(t, 2);
	      _acb_vec_clear(t_test, 2);
	      acb_mat_clear(tau);
	      acb_mat_clear(tau_test);
	      fmpz_mat_clear(m);
	      fmpz_mat_clear(minv);
	    }
	}
    }

  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

