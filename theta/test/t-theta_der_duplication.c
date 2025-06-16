#include <stdlib.h>
#include "theta.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("theta2_der_duplication....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
      slong prec = 200 + n_randint(state, 1000);
      
      acb_mat_t tau;
      acb_mat_t dth_half;
      acb_mat_t dth2;
      acb_mat_t dth2_test;
      acb_ptr th2_tau;
      acb_ptr th_half;
      int res;
      
      acb_mat_init(tau, 2, 2);
      acb_mat_init(dth_half, 4, 3);
      acb_mat_init(dth2, 16, 3);
      acb_mat_init(dth2_test, 16, 3);
      th2_tau = _acb_vec_init(16);
      th_half = _acb_vec_init(4);

      siegel_fundamental_domain_randtest(tau, state, prec);
      res = theta2_der_naive(th2_tau, dth2, tau, prec);
      
      acb_mat_scalar_div_si(tau, tau, 2, prec);
      res = res && theta_0123_der_naive(th_half, dth_half, tau, prec);
      
      theta_der_duplication(th2_tau, dth2_test, th_half, dth_half, prec);

      if (!res || !acb_mat_overlaps(dth2_test, dth2))
	{
	  flint_printf("FAIL\n");
	  acb_mat_printd(tau, 5); flint_printf("\n");
	  acb_mat_printd(dth_half, 5); flint_printf("\n");
	  acb_mat_printd(dth2, 5); flint_printf("\n");
	  acb_mat_printd(dth2_test, 5); flint_printf("\n");
	  fflush(stdout);
	  flint_abort();
	}
      
      acb_mat_clear(tau);
      acb_mat_clear(dth_half);
      acb_mat_clear(dth2);
      acb_mat_clear(dth2_test);
      _acb_vec_clear(th2_tau, 16);
      _acb_vec_clear(th_half, 4);

    }
  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

