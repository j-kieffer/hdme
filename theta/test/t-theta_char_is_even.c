#include <stdlib.h>

#include "theta.h"


int main()
{
  flint_printf("theta_char_is_even....");
  fflush(stdout);

  {
    slong g = 2;
    ulong ch, a, b;
    slong i, label;
    
    slong even[10] = {0, 1, 2, 3, 4, 6, 8, 9, 12, 15};
    slong odd[6] = {5, 7, 10, 11, 13, 14};
    
    for (i = 0; i < 10; i++)
      {
	ch = theta_char_set_label_g2(even[i]);
	a = theta_char_get_a(ch, g);
	b = theta_char_get_b(ch, g);
	label = theta_char_get_label_g2(ch);
	if (!theta_char_is_even(ch, g) || label != even[i])
	  {
	    flint_printf("FAIL\n");
	    flint_printf("i = %wd\n", i);
	    flint_printf("a = %wx\n", a);
	    flint_printf("b = %wx\n", b);
	    flint_printf("ch = %wx\n", ch);
	    flint_printf("even[i] = %wd\n", even[i]);
	    flint_printf("label = %wd\n", label);
	    flint_abort();
	  }
      }
    
    for (i = 0; i < 6; i++)
      {
	ch = theta_char_set_label_g2(odd[i]);
	a = theta_char_get_a(ch, g);
	b = theta_char_get_b(ch, g);
	label = theta_char_get_label_g2(ch);
	if (theta_char_is_even(ch, g) || label != odd[i])
	  {
	    flint_printf("FAIL\n");
	    flint_printf("i = %wd\n", i);
	    flint_printf("a = %wx\n", a);
	    flint_printf("b = %wx\n", b);
	    flint_printf("ch = %wx\n", ch);
	    flint_printf("odd[i] = %wd\n", odd[i]);
	    flint_printf("label = %wd\n", label);
	    flint_abort();
	  }
      }
  }
  
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
