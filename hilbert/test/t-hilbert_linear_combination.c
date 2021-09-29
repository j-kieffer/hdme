
#include "hilbert.h"

int main()
{
  slong iter;
  flint_rand_t state;
  slong delta;
  int print=0;
  
  flint_printf("hilbert_linear_combination....");
  fflush(stdout);
  
  flint_randinit(state);
  
  for (delta = 5; delta < 100; delta++)
    {
      if (hilbert_is_fundamental(delta))
	{
	  for (iter = 0; iter < 5 * arb_test_multiplier(); iter++)
	    {
	      acb_t t1, t2;
	      acb_mat_t R;
	      acb_mat_t tau;
	      fmpz* abcde;
	      slong prec = 1000 + n_randint(state, 1000);
	      slong m_bits = 4;
	      sp2gz_t m;
	      int res;
	      
	      acb_init(t1);
	      acb_init(t2);
	      acb_mat_init(R, 2, 2);
	      acb_mat_init(tau, 2, 2);
	      abcde = _fmpz_vec_init(5);
	      sp2gz_init(m, 2);
	      
	      hilbert_halfspace_randtest(t1, t2, state, prec);
	      hilbert_map(tau, t1, t2, delta, prec);
	      res = hilbert_linear_combination(abcde, tau, delta, prec);

	      if (!res)
		{
		  flint_printf("FAIL (wrong discriminant)\n");
		  flint_printf("delta = %wd, LLL result:\n", delta);
		  flint_printf("a = "); fmpz_print(&abcde[0]); flint_printf("\n");
		  flint_printf("b = "); fmpz_print(&abcde[1]); flint_printf("\n");
		  flint_printf("c = "); fmpz_print(&abcde[2]); flint_printf("\n");
		  flint_printf("d = "); fmpz_print(&abcde[3]); flint_printf("\n");
		  flint_printf("e = "); fmpz_print(&abcde[4]); flint_printf("\n");
		  fflush(stdout);
		  flint_abort();
		}
	      
	      if (print)
		{
		  flint_printf("delta = %wd, LLL result:\n", delta);
		  flint_printf("a = "); fmpz_print(&abcde[0]); flint_printf("\n");
		  flint_printf("b = "); fmpz_print(&abcde[1]); flint_printf("\n");
		  flint_printf("c = "); fmpz_print(&abcde[2]); flint_printf("\n");
		  flint_printf("d = "); fmpz_print(&abcde[3]); flint_printf("\n");
		  flint_printf("e = "); fmpz_print(&abcde[4]); flint_printf("\n");
		}

	      sp2gz_randtest(m, state, m_bits);
	      siegel_transform(tau, m, tau, prec);
	      res = hilbert_linear_combination(abcde, tau, delta, prec);

	      if (!res)
		{
		  flint_printf("FAIL (wrong discriminant)\n");
		  sp2gz_print(m); flint_printf("\n");
		  acb_mat_printd(tau, 10); flint_printf("\n");
		  flint_printf("delta = %wd, LLL result:\n", delta);
		  flint_printf("a = "); fmpz_print(&abcde[0]); flint_printf("\n");
		  flint_printf("b = "); fmpz_print(&abcde[1]); flint_printf("\n");
		  flint_printf("c = "); fmpz_print(&abcde[2]); flint_printf("\n");
		  flint_printf("d = "); fmpz_print(&abcde[3]); flint_printf("\n");
		  flint_printf("e = "); fmpz_print(&abcde[4]); flint_printf("\n");
		  fflush(stdout);
		  flint_abort();
		}
	      
	      if (print)
		{
		  flint_printf("delta = %wd, LLL result:\n", delta);
		  flint_printf("a = "); fmpz_print(&abcde[0]); flint_printf("\n");
		  flint_printf("b = "); fmpz_print(&abcde[1]); flint_printf("\n");
		  flint_printf("c = "); fmpz_print(&abcde[2]); flint_printf("\n");
		  flint_printf("d = "); fmpz_print(&abcde[3]); flint_printf("\n");
		  flint_printf("e = "); fmpz_print(&abcde[4]); flint_printf("\n");
		}
	      
	      acb_clear(t1);
	      acb_clear(t2);
	      acb_mat_clear(R);
	      acb_mat_clear(tau);
	      _fmpz_vec_clear(abcde, 5);
	      sp2gz_clear(m);
	    }
	}
    }

  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

