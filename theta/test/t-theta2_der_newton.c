#include <stdlib.h>
#include "theta.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("theta2_der_newton....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 5 * flint_test_multiplier(); iter++)
    {
      slong prec = 5000 + n_randint(state, 5000);
      acb_mat_t tau;
      acb_ptr th2, th2_test;
      acb_mat_t dth2, dth2_test;
      acb_t df, df_test;
      slong j, k;
      int res = 1;

      acb_mat_init(tau, 2, 2);
      th2 = _acb_vec_init(16);
      th2_test = _acb_vec_init(16);
      acb_mat_init(dth2, 16, 3);
      acb_mat_init(dth2_test, 16, 3);
      acb_init(df);
      acb_init(df_test);

      siegel_fundamental_domain_randtest(tau, state, prec);
      theta2_der_naive(th2_test, dth2_test, tau, prec);
      theta2_der_newton(th2, dth2, tau, prec);

      for (k = 1; k < 16; k++)
	{
	  /* Compute derivatives of theta_k/theta_0 */
	  for (j = 0; j < 3; j++)
	    {
	      acb_mul(df, acb_mat_entry(dth2, k, j), &th2[0], prec);
	      acb_submul(df, &th2[k], acb_mat_entry(dth2, 0, j), prec);
	      acb_div(df, df, &th2[0], prec);
	      acb_div(df, df, &th2[0], prec);
	      
	      acb_mul(df_test, acb_mat_entry(dth2_test, k, j), &th2_test[0], prec);
	      acb_submul(df_test, &th2_test[k], acb_mat_entry(dth2_test, 0, j), prec);
	      acb_div(df_test, df_test, &th2_test[0], prec);
	      acb_div(df_test, df_test, &th2_test[0], prec);

	      if (!acb_overlaps(df, df_test)) res = 0;
	    }			
	}

      if (!res)
	{
	  flint_printf("FAIL\n");
	  acb_mat_printd(dth2, 5); flint_printf("\n");
	  acb_mat_printd(dth2_test, 5); flint_printf("\n");
	  fflush(stdout);
	  flint_abort();
	}

      acb_mat_clear(tau);
      _acb_vec_clear(th2, 16);
      _acb_vec_clear(th2_test, 16);
      acb_mat_clear(dth2);
      acb_mat_clear(dth2_test);
      acb_clear(df);
      acb_clear(df_test);
    }
  
  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

