#include <stdlib.h>
#include "theta.h"


int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("theta0123_naive....");
  fflush(stdout);

  flint_rand_init(state);

  /* Check duplication formula */
  for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
      slong g = 2;
      slong prec = 10 + n_randint(state, 200);
      acb_mat_t tau;
      acb_ptr th_0123;
      acb_ptr th_0123_2tau;
      acb_ptr th2_0123_2tau;
      acb_ptr th2_2tau;
      fmpz_t B;
      
      int res;
      slong i;
      
      acb_mat_init(tau, g, g);
      th_0123 = _acb_vec_init(4);
      th_0123_2tau = _acb_vec_init(4);
      th2_0123_2tau = _acb_vec_init(4);
      th2_2tau = _acb_vec_init(16);
      fmpz_init(B);
      
      siegel_halfspace_randtest(tau, state, prec);
      res = theta_0123_naive_B(B, tau, prec);
      res = theta_0123_naive(th_0123, tau, prec);
      acb_mat_scalar_mul_2exp_si(tau, tau, 1);
      res = res && theta_0123_naive(th_0123_2tau, tau, prec);

      if (!res)
	{
	  flint_printf("FAIL\n");
	  flint_printf("res = %wd\n", res);
	  flint_printf("tau = "); acb_mat_printd(tau, 30); flint_printf("\n");
	  for (i = 0; i < 4; i++)
	    {
	      flint_printf("th_0123[%wd] = ", i);
	      acb_printd(&th_0123[i], 30); flint_printf("\n");
	    }
	  for (i = 0; i < 4; i++)
	    {
	      flint_printf("th_0123_2tau[%wd] = ", i);
	      acb_printd(&th_0123_2tau[i], 30); flint_printf("\n");
	    }
	  flint_abort();
	}
      
      theta_duplication(th2_2tau, th_0123, prec);
      for (i = 0; i < 4; i++)
	{
	  acb_sqr(&th2_0123_2tau[i], &th_0123_2tau[i], prec);
	}
      for (i = 0; i < 4; i++)
	{
	  if (!acb_overlaps(&th2_0123_2tau[i], &th2_2tau[i]))
	    {
	      flint_printf("FAIL\n");
	      flint_printf("tau = "); acb_mat_printd(tau, 30); flint_printf("\n");

	      flint_printf("th_0123[%wd] = ", i);
	      acb_printd(&th_0123[i], 30); flint_printf("\n");
	      
	      flint_printf("th2_0123_2tau[%wd] = ", i);
	      acb_printd(&th2_0123_2tau[i], 30); flint_printf("\n");
	      
	      flint_abort();
	    }
	}
      
      acb_mat_clear(tau);
      _acb_vec_clear(th_0123, 4);
      _acb_vec_clear(th_0123_2tau, 4);
      _acb_vec_clear(th2_0123_2tau, 4);
      _acb_vec_clear(th2_2tau, 16);
      fmpz_clear(B);
    }
  
  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
