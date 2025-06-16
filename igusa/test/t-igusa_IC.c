#include <stdlib.h>

#include "igusa.h"

int main()
{  
  slong iter;
  flint_rand_t state;
  
  flint_printf("igusa_IC....");
  fflush(stdout);

  flint_rand_init(state);

  for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
      acb_ptr I_acb;
      acb_ptr other_acb;
      acb_ptr resc;
      acb_ptr test;
      acb_t scal;
      fmpz* I;
      fmpz* other;
      acb_t R2_acb;
      fmpq* ABCD;
      fmpq_t R2;

      slong half[4] = IGUSA_HALFWEIGHTS;
      slong weights[4] = IGUSA_WEIGHTS;
      slong weights2[4] = IC_WEIGHTS;
      slong mag_bits = 10;
      slong prec = 1000;
      slong k;
      int res = 1;

      I_acb = _acb_vec_init(4);
      other_acb = _acb_vec_init(4);
      resc = _acb_vec_init(4);
      test = _acb_vec_init(4);
      acb_init(scal);
      I = _fmpz_vec_init(4);
      other = _fmpz_vec_init(4);
      ABCD = _fmpq_vec_init(4);
      acb_init(R2_acb);
      fmpq_init(R2);

      for (k = 0; k < 4; k++) fmpz_randtest_not_zero(&I[k], state, mag_bits);
      cov_normalize_fmpz(I, I, 4, half);
      for (k = 0; k < 4; k++) acb_set_fmpz(&I_acb[k], &I[k]);
      acb_randtest_precise(scal, state, prec, mag_bits);

      igusa_streng(other_acb, I_acb, prec);
      cov_rescale(resc, I_acb, scal, 4, weights, prec);
      igusa_streng(resc, resc, prec);
      cov_rescale(test, other_acb, scal, 4, weights, prec);
      for (k = 0; k < 4; k++)
	{ if (!acb_overlaps(&test[k], &resc[k])) res = 0; }
      if (!res)
	{
	  flint_printf("FAIL (Streng rescaling)\n");
	  fflush(stdout);
	  flint_abort();
	}
      
      igusa_streng(other_acb, I_acb, prec);
      igusa_streng_fmpz(other, I);
      cov_find_rescaling(scal, other_acb, other, 4, weights, prec);

      igusa_from_streng(test, other_acb, prec);
      igusa_from_streng_fmpz(other, other);
      for (k = 0; k < 4; k++)
	{ if (!acb_overlaps(&test[k], &I_acb[k])) res = 0; }
      if (!res)
	{
	  flint_printf("FAIL (Streng inversion)\n");
	  fflush(stdout);
	  flint_abort();
	}      
      for (k = 0; k < 4; k++)
	{ if (!fmpz_equal(&other[k], &I[k])) res = 0; }
      if (!res)
	{
	  flint_printf("FAIL (Streng_fmpz inversion)\n");
	  for (k = 0; k < 4; k++)
	    {
	      fmpz_print(&other[k]); flint_printf("\n");
	      fmpz_print(&I[k]); flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_abort();
	}
      
      igusa_IC(other_acb, I_acb, prec);
      cov_rescale(resc, I_acb, scal, 4, weights, prec);
      igusa_IC(resc, resc, prec);
      cov_rescale(test, other_acb, scal, 4, weights2, prec);
      for (k = 0; k < 4; k++)
	{ if (!acb_overlaps(&test[k], &resc[k])) res = 0; }
      if (!res)
	{
	  flint_printf("FAIL (IC rescaling)\n");
	  fflush(stdout);
	  flint_abort();
	}

      
      /* igusa_streng_fmpz(other, I);
      for (k = 0; k < 4; k++)
	{	  	  
	  fmpz_print(&other[k]); flint_printf("\n");
	}
      for (k = 0; k < 4; k++)
	{	  	  
	  fmpz_print(&I[k]); flint_printf("\n");
	  } */
      
      igusa_IC(other_acb, I_acb, prec);
      igusa_IC_fmpz(other, I);
      /*
      for (k = 0; k < 4; k++)
	{
	  fmpz_print(&other[k]); flint_printf("\n");
	  acb_printd(&other_acb[k], 10); flint_printf("\n");	  
	  fmpz_print(&I[k]); flint_printf("\n\n");
	  } */
      cov_find_rescaling(scal, other_acb, other, 4, weights2, prec);
      
      igusa_from_IC(test, other_acb, prec);
      igusa_from_IC_fmpz(other, other);
      for (k = 0; k < 4; k++)
	{ if (!acb_overlaps(&test[k], &I_acb[k])) res = 0; }
      if (!res)
	{
	  flint_printf("FAIL (IC inversion)\n");
	  fflush(stdout);
	  flint_abort();
	}      
      for (k = 0; k < 4; k++)
	{ if (!fmpz_equal(&other[k], &I[k])) res = 0; }
      if (!res)
	{
	  flint_printf("FAIL (IC_fmpz inversion)\n");
	  for (k = 0; k < 4; k++)
	    {
	      fmpz_print(&other[k]); flint_printf("\n");
	      fmpz_print(&I[k]); flint_printf("\n");
	    }
	  igusa_IC_fmpz(other, I);
	  for (k = 0; k < 4; k++)
	    {
	      fmpz_print(&other[k]); flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_abort();
	}      
            
      igusa_ABCD_from_IC(other_acb, I_acb, prec);
      cov_rescale(resc, I_acb, scal, 4, weights2, prec);
      igusa_ABCD_from_IC(resc, resc, prec);
      cov_rescale(test, other_acb, scal, 4, weights2, prec);
      for (k = 0; k < 4; k++)
	{ if (!acb_overlaps(&test[k], &resc[k])) res = 0; }
      if (!res)
	{
	  flint_printf("FAIL (ABCD rescaling)\n");
	  fflush(stdout);
	  flint_abort();
	}
      
      igusa_ABCD_from_IC(other_acb, I_acb, prec);
      igusa_ABCD_from_IC_fmpz(ABCD, I);
      for (k = 0; k < 4; k++)
	{ if (!acb_contains_fmpq(&other_acb[k], &ABCD[k])) res = 0; }
      if (!res)
	{
	  flint_printf("FAIL (ABCD_fmpz)\n");
	  fflush(stdout);
	  flint_abort();
	}
      

      igusa_R2_from_IC(R2_acb, I_acb, prec);
      cov_rescale(resc, I_acb, scal, 4, weights2, prec);
      igusa_R2_from_IC(&resc[0], resc, prec);
      acb_pow_si(scal, scal, 30, prec);
      acb_mul(&test[0], R2_acb, scal, prec);
      if (!acb_overlaps(&test[0], &resc[0]))
	{	  
	  flint_printf("FAIL (R2 rescaling)\n");
	  fflush(stdout);
	  flint_abort();
	}

      igusa_R2_from_IC(R2_acb, I_acb, prec);
      igusa_R2_from_IC_fmpz(R2, I);
      if (!acb_contains_fmpq(R2_acb, R2))
	{	  
	  flint_printf("FAIL (R2_fmpz)\n");
	  fflush(stdout);
	  flint_abort();
	}      
      
      _acb_vec_clear(I_acb, 4);
      _acb_vec_clear(other_acb, 4);
      _acb_vec_clear(resc, 4);
      _acb_vec_clear(test, 4);
      acb_clear(scal);
      _fmpz_vec_clear(I, 4);
      _fmpz_vec_clear(other, 4);
      _fmpq_vec_clear(ABCD, 4);
      acb_clear(R2_acb);
      fmpq_clear(R2);
    }

  flint_rand_clear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
