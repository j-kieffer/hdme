#include <stdlib.h>

#include "igusa.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("thomae_theta4....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
      slong prec = 500 + n_randint(state, 1000);

      acb_ptr th2;
      acb_ptr th4_test;
      acb_ptr I, IC;
      acb_poly_t crv;
      acb_ptr roots;
      acb_ptr new_roots;
      acb_ptr ros;
      acb_ptr th4;
      slong k;
      slong perm;
      int res, cnt;

      th2 = _acb_vec_init(16);
      th4_test = _acb_vec_init(16);
      I = _acb_vec_init(4);
      IC = _acb_vec_init(4);
      acb_poly_init(crv);
      roots = _acb_vec_init(6);
      new_roots = _acb_vec_init(6);
      ros = _acb_vec_init(3);
      th4 = _acb_vec_init(16);

      /* Generate theta constants and covariants */
      theta2_randtest(th2, state, prec);
      for (k = 0; k < 16; k++) acb_sqr(&th4_test[k], &th2[k], prec);
      for (k = 1; k < 16; k++) acb_div(&th4_test[k], &th4_test[k], &th4_test[0], prec);
      acb_one(&th4_test[0]);
      igusa_from_theta2(I, th2, prec);
      igusa_IC(IC, I, prec);

      /* Mestre */
      res = mestre(crv, IC, prec);
      if (!res)
	{
	  flint_printf("FAIL (mestre)\n");
	  flint_printf("I: \n");
	  for (k = 0; k < 4; k++)
	    {
	      acb_printd(&I[k], 30); flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_abort();
	}
      res = thomae_roots(roots, crv, prec);
      if (!res)
	{
	  flint_printf("FAIL (roots)\n");
	  flint_printf("Curve: ");
	  acb_poly_printd(crv, 30);
	  flint_printf("\n");
	  flint_printf("Roots: \n");
	  for (k = 0; k < 6; k++)
	    {
	      acb_printd(&roots[k], 30); flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_abort();
	}

      res = 0;
      for (perm = 0; perm < 720; perm++)
	{
	  thomae_reorder(new_roots, roots, perm);
	  thomae_rosenhain(ros, new_roots, prec);
	  thomae_theta4(th4, ros, prec);
	  cnt = 0;
	  for (k = 0; k < 16; k++)
	    {
	      if (acb_overlaps(&th4[k], &th4_test[k])) cnt++;
	      /*if (acb_overlaps(&th4[k], &th4_test[k])
		  && (k!=0) && (k!=5) && (k!=7) && (k!=10)
		  && (k!=11) && (k!=13) && (k!=14))
		{
		  flint_printf("perm no.%wd: Index %wd correct\n", perm, k);
		  }*/
	    }
	  /*if (cnt > 10)
	    {
	      flint_printf("Perm no.%wd, th4: \n", perm);
	      for (k = 0; k < 16; k++)
		{
		  acb_printd(&th4[k], 30); flint_printf("\n");
		}
		} */
	  if (cnt == 16)
	    {
	      res = 1;
	      /*flint_printf("Found correct permutation: %wd\n", perm);*/
	    }
	}
      if (!res)
	{
	  flint_printf("FAIL (no permutation found)\n");
	  flint_printf("roots: \n");
	  for (k = 0; k < 6; k++)
	    {
	      acb_printd(&roots[k], 30); flint_printf("\n");
	    }
	  flint_printf("th4_test: \n");
	  for (k = 0; k < 16; k++)
	    {
	      acb_printd(&th4_test[k], 30); flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_abort();
	}

      _acb_vec_clear(th2, 16);
      _acb_vec_clear(th4_test, 16);
      _acb_vec_clear(I, 4);
      _acb_vec_clear(IC, 4);
      acb_poly_clear(crv);
      _acb_vec_clear(roots, 6);
      _acb_vec_clear(new_roots, 6);
      _acb_vec_clear(ros, 3);
      _acb_vec_clear(th4, 16);
    }
  
  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
