#include <stdlib.h>
#include "theta.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("theta_unbalanced....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 4 * flint_test_multiplier(); iter++)
    {
      slong g = 2;
      slong bits = 500 + n_randint(state, 1000);
      slong jump = 4 + n_randint(state, 100);
      slong prec = 10 * bits;
      acb_t c;
      arb_t tol;
      acb_mat_t tau;
      acb_mat_t w;
      fmpz_mat_t real_red;
      fmpz_mat_t eta, eta_inv, M1;
      fmpz_t nb_real_red;
      acb_ptr th2;
      acb_ptr current_th2;
      int res;
      slong k;
      /* slong k1, k2;
	 slong i, j; */

      acb_mat_init(tau, g, g);
      acb_mat_init(w, g, g);
      fmpz_mat_init(real_red, 2*g, 2*g);
      fmpz_mat_init(eta, 2*g, 2*g);
      fmpz_mat_init(eta_inv, 2*g, 2*g);
      fmpz_mat_init(M1, 2*g, 2*g);
      fmpz_init(nb_real_red);
      acb_init(c);
      arb_init(tol);
      th2 = _acb_vec_init(16);
      current_th2 = _acb_vec_init(16);

      /* Set tau to reduced, unbalanced value */
      siegel_fundamental_domain_randtest(tau, state, prec);
      arb_add_si(acb_imagref(acb_mat_entry(tau, 1, 1)),
		 acb_imagref(acb_mat_entry(tau, 1, 1)), jump, prec);
      
      /* arb_set_d(acb_realref(c), -0.415); */
      /* arb_set_d(acb_imagref(c), 1.05); */
      /* acb_set(acb_mat_entry(tau, 0, 0), c); */

      /* arb_set_d(acb_realref(c), 0.318); */
      /* arb_set_d(acb_imagref(c), 0.384); */
      /* acb_set(acb_mat_entry(tau, 0, 1), c); */
      /* acb_set(acb_mat_entry(tau, 1, 0), c); */

      /* arb_set_d(acb_realref(c), -0.23); */
      /* arb_set_d(acb_imagref(c), 5.79); */
      /* acb_set(acb_mat_entry(tau, 1, 1), c); */
      
      /* acb_mat_printd(tau, 10); flint_printf("\n"); */

      arb_one(tol);
      arb_mul_2exp_si(tol, tol, -bits);      
      res = siegel_is_in_fundamental_domain(tau, tol, prec);
      if (!res)
	{	  
	  flint_printf("FAIL (fundamental domain)\n");
	  acb_mat_printd(tau, 10);
	  flint_abort();
	}

      theta2_unif(th2, tau, prec);
      theta2_renormalize(th2, th2, prec);
      theta2_naive(current_th2, tau, prec);

      for (k = 0; k < 16; k++)
	{
	  if (!acb_overlaps(&th2[k], &current_th2[k])) res = 0;
	}
      if (!res)
	{
	  flint_printf("FAIL (theta)\n");
	  acb_mat_printd(tau, 10);
	  for (k = 0; k < 16; k++)
	    {
	      acb_printd(&th2[k], 30); flint_printf("\n");
	      acb_printd(&current_th2[k], 30); flint_printf("\n\n");
	    }
	  flint_abort();
	}
      
      /* for (k = 0; k < 16; k++) */
      /* 	{ */
      /* 	  acb_printd(&th2[k], 30); flint_printf("\n"); */
      /* 	} */

      /* flint_printf("Renormalized values:\n"); */
      /* theta2_renormalize(th2, th2, prec); */
      /* for (k = 0; k < 16; k++) */
      /* 	{ */
      /* 	  acb_printd(&th2[k], 30); flint_printf("\n"); */
      /* 	} */
      
      /* k2 = theta_newton_k2(w, tau, prec); */
      /* siegel_reduce_real(w, real_red, w, tol, prec); */
      /* k1 = theta_newton_k1(w, w, prec); */
      /* flint_printf("k1, k2: %wd, %wd\n", k1, k2); */

      /* flint_printf("Reduced matrix:\n"); */
      /* acb_mat_printd(w, 10); */

      /* theta2_newton(current_th2, w, prec); */
      /* for (k = 0; k < 16; k++) */
      /* 	{ */
      /* 	  acb_printd(&current_th2[k], 30); flint_printf("\n"); */
      /* 	} */
      
      /* /\* Propagate along k1 sequence *\/ */
      /* for (i = 0; (i < k1) && res; i++) */
      /* 	{ */
      /* 	  /\* Each step: rescale; compute square roots in four first */
      /* 	     entries of current_th2; apply theta_duplication *\/ */
      /* 	  for (j = 1; j < 4; j++) */
      /* 	    { */
      /* 	      acb_div(&current_th2[j], &current_th2[j], &current_th2[0], prec); */
      /* 	      res = res && acb_sqrt_goodpos(&current_th2[j], &current_th2[j], prec); */
      /* 	    } */
      /* 	  acb_one(&current_th2[0]); */
      /* 	  /\* theta_duplication uses only the first 4 entries and allows */
      /* 	     aliasing *\/ */
      /* 	  theta_duplication(current_th2, current_th2, prec); */
      /* 	} */

      /* flint_printf("After k1 steps:\n"); */
      /* for (k = 0; k < 16; k++) */
      /* 	{ */
      /* 	  acb_printd(&current_th2[k], 30); flint_printf("\n"); */
      /* 	} */
      
      /* /\* We're back to the end of the k2 sequence. Propagate again *\/ */
      /* fmpz_mat_zero(eta); */
      /* fmpz_set_si(fmpz_mat_entry(eta, 0, 2), 1); */
      /* fmpz_set_si(fmpz_mat_entry(eta, 1, 1), 1); */
      /* fmpz_set_si(fmpz_mat_entry(eta, 2, 0), -1); */
      /* fmpz_set_si(fmpz_mat_entry(eta, 3, 3), 1); */
      /* fmpz_mat_direct_inv(eta_inv, eta); */
      /* fmpz_mat_one(M1); */
      
      /* /\* Real reduction happens only on the y1 coordinate *\/ */
      /* fmpz_mat_print_pretty(real_red); */
      
      /* fmpz_set(nb_real_red, fmpz_mat_entry(real_red, 0, 2)); */
      /* if (arb_is_positive(acb_realref(acb_mat_entry(tau, 0, 0)))) */
      /* 	{ */
      /* 	  fmpz_set_si(fmpz_mat_entry(M1, 0, 2), 1); */
      /* 	} */
      /* else */
      /* 	{ */
      /* 	  fmpz_set_si(fmpz_mat_entry(M1, 0, 2), -1); */
      /* 	} */
      
      /* for (i = 0; (i < k2) && res; i++) */
      /* 	{ */
      /* 	  /\* Each step: apply translation if nb_real_red is currently odd; */
      /* 	     apply theta_transform with eta; square roots and duplication; */
      /* 	     finally apply theta2_transform with eta_inv. *\/ */
      /* 	  if (fmpz_is_odd(nb_real_red)) */
      /* 	    { */
      /* 	      theta2_transform(current_th2, M1, current_th2, prec); */
      /* 	    } */
      /* 	  /\* Divide by 2, round towards zero *\/ */
      /* 	  fmpz_tdiv_q_2exp(nb_real_red, nb_real_red, 1); */
      /* 	  /\* Get theta constants at eta * current *\/ */
      /* 	  theta2_transform(current_th2, eta, current_th2, prec); */
      /* 	  /\* Rescale and square roots *\/ */
      /* 	  for (j = 1; j < 4; j++) */
      /* 	    { */
      /* 	      acb_div(&current_th2[j], &current_th2[j], &current_th2[0], prec); */
      /* 	      res = res && acb_sqrt_goodpos(&current_th2[j], &current_th2[j], prec); */
      /* 	    } */
      /* 	  acb_one(&current_th2[0]); */
      /* 	  /\* Duplication *\/ */
      /* 	  theta_duplication(current_th2, current_th2, prec); */
      /* 	  /\* Apply eta_inv *\/ */
      /* 	  theta2_transform(current_th2, eta_inv, current_th2, prec); */
      /* 	} */

      
      /* flint_printf("After k2 steps:\n"); */
      /* for (k = 0; k < 16; k++) */
      /* 	{ */
      /* 	  acb_printd(&current_th2[k], 30); flint_printf("\n"); */
      /* 	} */
      /* flint_printf("Renormalized:\n"); */
      /* theta2_renormalize(current_th2, current_th2, prec); */
      /* for (k = 0; k < 16; k++) */
      /* 	{ */
      /* 	  acb_printd(&current_th2[k], 30); flint_printf("\n"); */
      /* 	} */
      
      acb_mat_clear(tau);
      acb_mat_clear(w);
      fmpz_mat_clear(real_red);
      fmpz_mat_clear(eta);
      fmpz_mat_clear(eta_inv);
      fmpz_mat_clear(M1);
      fmpz_clear(nb_real_red);
      acb_clear(c);
      arb_clear(tol);
      _acb_vec_clear(th2, 16);
      _acb_vec_clear(current_th2, 16);

    }  
      
  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
