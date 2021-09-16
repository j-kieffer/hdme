
#include "siegel.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("siegel_star....");
  fflush(stdout);
  
  flint_randinit(state);
  
  /* Check the cocycle relation: I_g = m^*(z) m^(-1)*(mz) */
  for (iter = 0; iter < 1000 * arb_test_multiplier(); iter++)
    {
      sp2gz_t m, minv;
      acb_mat_t z1, z2, w1, w2, i, t;
      slong g, prec;
      int valid_transform;
      
      g = 1 + n_randint(state, 10);
      sp2gz_init(m, g);
      sp2gz_init(minv, g);
      acb_mat_init(z1, g, g);
      acb_mat_init(z2, g, g);
      acb_mat_init(w1, g, g);
      acb_mat_init(w2, g, g);
      acb_mat_init(i, g, g);
      acb_mat_init(t, g, g);
      
      acb_mat_one(i);
      sp2gz_randtest(m, state, n_randint(state, 20));
      sp2gz_inv(minv, m);
      
      prec = 200 + n_randint(state, 200);
      
      siegel_halfspace_randtest(z1, state, prec);
      valid_transform = siegel_transform(z2, m, z1, prec);

      if (!valid_transform)
	{
	  flint_printf("FAIL (insufficient precision to invert)\n");
	  flint_printf("g = %wd\n", g);
	  flint_printf("m = "); sp2gz_print(m); flint_printf("\n\n");
	  flint_printf("z1 = "); acb_mat_printd(z1, 30); flint_printf("\n\n");
	  flint_abort();
	}

      siegel_star(w1, m, z1, prec);
      siegel_star(w2, minv, z2, prec);
      acb_mat_mul(t, w1, w2, prec);

      if (!acb_mat_overlaps(t, i))
        {
	  flint_printf("FAIL\n");
	  flint_printf("g = %wd\n", g);
	  flint_printf("m = "); sp2gz_print(m); flint_printf("\n\n");
	  flint_printf("minv = "); sp2gz_print(minv); flint_printf("\n\n");
	  flint_printf("z1 = "); acb_mat_printd(z1, 30); flint_printf("\n\n");
	  flint_printf("z2 = "); acb_mat_printd(z2, 30); flint_printf("\n\n");
	  flint_printf("w1 = "); acb_mat_printd(w1, 30); flint_printf("\n\n");
	  flint_printf("w2 = "); acb_mat_printd(w2, 30); flint_printf("\n\n");
	  flint_printf("t1 = "); acb_mat_printd(t, 30); flint_printf("\n\n");
	  flint_abort();
        }

      sp2gz_clear(m);
      sp2gz_clear(minv);
      acb_mat_clear(z1);
      acb_mat_clear(z2);
      acb_mat_clear(w1);
      acb_mat_clear(w2);
      acb_mat_clear(t);
      acb_mat_clear(i);
    }

  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
