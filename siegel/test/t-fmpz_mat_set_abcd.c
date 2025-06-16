#include <stdlib.h>

#include "siegel.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("fmpz_mat_set_abcd....");
  fflush(stdout);
  
  flint_rand_init(state);
  
  for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
      slong g = 1 + n_randint(state, 5);
      fmpz_mat_t m, a, b, c, d, m_test;
      slong bits = 1 + n_randint(state, 10);

      fmpz_mat_init(m, 2*g, 2*g);
      fmpz_mat_init(a, g, g);
      fmpz_mat_init(b, g, g);
      fmpz_mat_init(c, g, g);
      fmpz_mat_init(d, g, g);
      fmpz_mat_init(m_test, 2*g, 2*g);

      fmpz_mat_randtest(m, state, bits);
      fmpz_mat_get_a(a, m);
      fmpz_mat_get_b(b, m);
      fmpz_mat_get_c(c, m);
      fmpz_mat_get_d(d, m);
      fmpz_mat_set_abcd(m_test, a, b, c, d);

      if (!fmpz_mat_equal(m, m_test))
	{
	  flint_printf("FAIL\n");
	  fmpz_mat_print_pretty(m); flint_printf("\n");
	  fmpz_mat_print_pretty(a); flint_printf("\n");
	  fmpz_mat_print_pretty(b); flint_printf("\n");
	  fmpz_mat_print_pretty(c); flint_printf("\n");
	  fmpz_mat_print_pretty(d); flint_printf("\n");
	  fmpz_mat_print_pretty(m_test); flint_printf("\n");
	  fflush(stdout);
	  flint_abort();
	}
      
      fmpz_mat_clear(m);
      fmpz_mat_clear(a);
      fmpz_mat_clear(b);
      fmpz_mat_clear(c);
      fmpz_mat_clear(d);
      fmpz_mat_clear(m_test);
    }
  
  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}


