#include <stdlib.h>

#include "igusa.h"

int main()
{
  
  slong iter;
  flint_rand_t state;
  
  flint_printf("mestre_point_on_conic....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
      slong prec = 50 + n_randint(state, 1000);
      acb_ptr conic;
      acb_ptr pt;
      int res = 1;
      int k;

      conic = _acb_vec_init(6);
      pt = _acb_vec_init(3);

      mestre_conic_randtest(conic, state, prec);
      res = mestre_point_on_conic(pt, conic, prec);

      if (!res)
	{
	  flint_printf("FAIL (Compute pt)\n");
	  flint_printf("Conic: \n");
	  for (k = 0; k < 6; k++)
	    {
	      flint_printf("    "); acb_printd(&conic[k], 30); flint_printf("\n");
	    }
	  flint_printf("Point: \n");
	  for (k = 0; k < 3; k++)
	    {
	      flint_printf("    "); acb_printd(&pt[k], 30); flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_abort();
	}

      res = !mestre_point_is_outside_conic(pt, conic, prec);

      if (!res)
	{
	  flint_printf("FAIL (Point outside conic)\n");
	  flint_printf("Conic: \n");
	  for (k = 0; k < 6; k++)
	    {
	      flint_printf("    "); acb_printd(&conic[k], 30); flint_printf("\n");
	    }
	  flint_printf("Point: \n");
	  for (k = 0; k < 3; k++)
	    {
	      flint_printf("    "); acb_printd(&pt[k], 30); flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_abort();
	}
      
      _acb_vec_clear(conic, 6);
      _acb_vec_clear(pt, 3);
    }

  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
