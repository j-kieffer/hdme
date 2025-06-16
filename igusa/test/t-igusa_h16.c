#include <stdlib.h>

#include "igusa.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("igusa_h....");
  fflush(stdout);

  flint_rand_init(state);

  /* Test formula I6 = h16/h10 */
  for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
      acb_ptr theta2;
      acb_ptr I, IC;
      acb_t h10, h16;
      slong k;

      slong prec = 50 + n_randint(state, 2000);

      theta2 = _acb_vec_init(16);
      I = _acb_vec_init(4);
      IC = _acb_vec_init(4);
      acb_init(h10);
      acb_init(h16);

      theta2_randtest(theta2, state, prec);
      igusa_from_theta2(I, theta2, prec);
      igusa_IC(IC, I, prec);

      igusa_h10(h10, theta2, prec);
      igusa_h16(h16, theta2, prec);
      acb_div(h16, h16, h10, prec);

      if (!acb_overlaps(h16, &IC[2]))
	{
	  flint_printf("FAIL (h16)");
	  flint_printf("Theta constants:\n");
	  for (k = 0; k < 16; k++)
	    {
	      acb_printd(&theta2[k], 30); flint_printf("\n");
	    }
	  flint_printf("\nh16/h10: "); acb_printd(h16, 30);
	  flint_printf("\nI6: "); acb_printd(&IC[1], 30);
	  fflush(stdout);
	  flint_abort();
	}

      _acb_vec_clear(theta2, 16);
      _acb_vec_clear(IC, 4);
      _acb_vec_clear(I, 4);
      acb_clear(h16);
      acb_clear(h10);
    }
 
  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

      

