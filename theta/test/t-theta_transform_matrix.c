#include <stdlib.h>

#include "theta.h"

/* Check using Dupont's tables p. 146-147 (these tables only give
   epsilon mod 4, not 8) */

static void
check_dupont_init(fmpz_mat_t t)
{
  fmpz_mat_zero(t);
  
  fmpz_set_si(fmpz_mat_entry(t, 1, 0), 1);
  fmpz_set_si(fmpz_mat_entry(t, 2, 0), 2);
  fmpz_set_si(fmpz_mat_entry(t, 3, 0), 3);
  fmpz_set_si(fmpz_mat_entry(t, 4, 0), 4);
  fmpz_set_si(fmpz_mat_entry(t, 6, 0), 6);
  fmpz_set_si(fmpz_mat_entry(t, 8, 0), 8);
  fmpz_set_si(fmpz_mat_entry(t, 9, 0), 9);
  fmpz_set_si(fmpz_mat_entry(t, 12, 0), 12);
  fmpz_set_si(fmpz_mat_entry(t, 15, 0), 15);
}

static void
check_dupont_id(fmpz_mat_t t)
{
  check_dupont_init(t);

  fmpz_set_si(fmpz_mat_entry(t, 0, 2), 0);
  fmpz_set_si(fmpz_mat_entry(t, 1, 2), 1);
  fmpz_set_si(fmpz_mat_entry(t, 2, 2), 2);
  fmpz_set_si(fmpz_mat_entry(t, 3, 2), 3);
  fmpz_set_si(fmpz_mat_entry(t, 4, 2), 4);
  fmpz_set_si(fmpz_mat_entry(t, 6, 2), 6);
  fmpz_set_si(fmpz_mat_entry(t, 8, 2), 8);
  fmpz_set_si(fmpz_mat_entry(t, 9, 2), 9);
  fmpz_set_si(fmpz_mat_entry(t, 12, 2), 12);
  fmpz_set_si(fmpz_mat_entry(t, 15, 2), 15);
}

static void
check_dupont_J(fmpz_mat_t t)
{
  check_dupont_init(t);

  fmpz_set_si(fmpz_mat_entry(t, 15, 1), 4);
  
  fmpz_set_si(fmpz_mat_entry(t, 1, 2), 4);
  fmpz_set_si(fmpz_mat_entry(t, 2, 2), 8);
  fmpz_set_si(fmpz_mat_entry(t, 3, 2), 12);
  fmpz_set_si(fmpz_mat_entry(t, 4, 2), 1);
  fmpz_set_si(fmpz_mat_entry(t, 6, 2), 9);
  fmpz_set_si(fmpz_mat_entry(t, 8, 2), 2);
  fmpz_set_si(fmpz_mat_entry(t, 9, 2), 6);
  fmpz_set_si(fmpz_mat_entry(t, 12, 2), 3);
  fmpz_set_si(fmpz_mat_entry(t, 15, 2), 15);
}

static void
check_dupont_M1(fmpz_mat_t t)
{
  check_dupont_init(t);
  
  fmpz_set_si(fmpz_mat_entry(t, 4, 1), 1);
  fmpz_set_si(fmpz_mat_entry(t, 6, 1), 1);
  fmpz_set_si(fmpz_mat_entry(t, 12, 1), 1);
  fmpz_set_si(fmpz_mat_entry(t, 15, 1), 1);

  fmpz_set_si(fmpz_mat_entry(t, 0, 2), 1);
  fmpz_set_si(fmpz_mat_entry(t, 1, 2), 0);
  fmpz_set_si(fmpz_mat_entry(t, 2, 2), 3);
  fmpz_set_si(fmpz_mat_entry(t, 3, 2), 2);
  fmpz_set_si(fmpz_mat_entry(t, 4, 2), 4);
  fmpz_set_si(fmpz_mat_entry(t, 6, 2), 6);
  fmpz_set_si(fmpz_mat_entry(t, 8, 2), 9);
  fmpz_set_si(fmpz_mat_entry(t, 9, 2), 8);
  fmpz_set_si(fmpz_mat_entry(t, 12, 2), 12);
  fmpz_set_si(fmpz_mat_entry(t, 15, 2), 15);
}

static void
check_dupont_M2(fmpz_mat_t t)
{
  check_dupont_init(t);
  
  fmpz_set_si(fmpz_mat_entry(t, 8, 1), 1);
  fmpz_set_si(fmpz_mat_entry(t, 9, 1), 1);
  fmpz_set_si(fmpz_mat_entry(t, 12, 1), 1);
  fmpz_set_si(fmpz_mat_entry(t, 15, 1), 1);

  fmpz_set_si(fmpz_mat_entry(t, 0, 2), 2);
  fmpz_set_si(fmpz_mat_entry(t, 1, 2), 3);
  fmpz_set_si(fmpz_mat_entry(t, 2, 2), 0);
  fmpz_set_si(fmpz_mat_entry(t, 3, 2), 1);
  fmpz_set_si(fmpz_mat_entry(t, 4, 2), 6);
  fmpz_set_si(fmpz_mat_entry(t, 6, 2), 4);
  fmpz_set_si(fmpz_mat_entry(t, 8, 2), 8);
  fmpz_set_si(fmpz_mat_entry(t, 9, 2), 9);
  fmpz_set_si(fmpz_mat_entry(t, 12, 2), 12);
  fmpz_set_si(fmpz_mat_entry(t, 15, 2), 15);
}

static void
check_dupont_M3(fmpz_mat_t t)
{
  check_dupont_init(t);
  
  fmpz_set_si(fmpz_mat_entry(t, 12, 1), 6);
  fmpz_set_si(fmpz_mat_entry(t, 15, 1), 6);

  fmpz_set_si(fmpz_mat_entry(t, 0, 2), 0);
  fmpz_set_si(fmpz_mat_entry(t, 1, 2), 1);
  fmpz_set_si(fmpz_mat_entry(t, 2, 2), 2);
  fmpz_set_si(fmpz_mat_entry(t, 3, 2), 3);
  fmpz_set_si(fmpz_mat_entry(t, 4, 2), 6);
  fmpz_set_si(fmpz_mat_entry(t, 6, 2), 4);
  fmpz_set_si(fmpz_mat_entry(t, 8, 2), 9);
  fmpz_set_si(fmpz_mat_entry(t, 9, 2), 8);
  fmpz_set_si(fmpz_mat_entry(t, 12, 2), 15);
  fmpz_set_si(fmpz_mat_entry(t, 15, 2), 12);
}

static void
check_dupont_I(fmpz_mat_t t)
{
  check_dupont_init(t);

  fmpz_set_si(fmpz_mat_entry(t, 0, 2), 0);
  fmpz_set_si(fmpz_mat_entry(t, 1, 2), 1);
  fmpz_set_si(fmpz_mat_entry(t, 2, 2), 3);
  fmpz_set_si(fmpz_mat_entry(t, 3, 2), 2);
  fmpz_set_si(fmpz_mat_entry(t, 4, 2), 12);
  fmpz_set_si(fmpz_mat_entry(t, 6, 2), 15);
  fmpz_set_si(fmpz_mat_entry(t, 8, 2), 8);
  fmpz_set_si(fmpz_mat_entry(t, 9, 2), 9);
  fmpz_set_si(fmpz_mat_entry(t, 12, 2), 4);
  fmpz_set_si(fmpz_mat_entry(t, 15, 2), 6);
}

static void
check_dupont_Sigma(fmpz_mat_t t)
{
  check_dupont_init(t);

  fmpz_set_si(fmpz_mat_entry(t, 0, 2), 0);
  fmpz_set_si(fmpz_mat_entry(t, 1, 2), 2);
  fmpz_set_si(fmpz_mat_entry(t, 2, 2), 1);
  fmpz_set_si(fmpz_mat_entry(t, 3, 2), 3);
  fmpz_set_si(fmpz_mat_entry(t, 4, 2), 8);
  fmpz_set_si(fmpz_mat_entry(t, 6, 2), 9);
  fmpz_set_si(fmpz_mat_entry(t, 8, 2), 4);
  fmpz_set_si(fmpz_mat_entry(t, 9, 2), 6);
  fmpz_set_si(fmpz_mat_entry(t, 12, 2), 12);
  fmpz_set_si(fmpz_mat_entry(t, 15, 2), 15);
}

int main()
{
  
  slong g = 2;
  fmpz_mat_t transf_mat;
  fmpz_mat_t test_mat;
  fmpz_mat_t aux;
  fmpz_mat_t T;
  fmpz_mat_t eta;
  int res = 0;

  
  flint_printf("theta_transform_matrix....");
  fflush(stdout);

  fmpz_mat_init(transf_mat, n_pow(2, 2*g), 3);
  fmpz_mat_init(test_mat, n_pow(2, 2*g), 3);
  fmpz_mat_init(aux, 2*g, 2*g);
  fmpz_mat_init(T, g, g);
  fmpz_mat_init(eta, 2*g, 2*g);

  fmpz_mat_one(eta);
  theta_transform_matrix(transf_mat, eta);
  check_dupont_id(test_mat);
  res = fmpz_mat_equal(transf_mat, test_mat);
  if (!res)
    {
      flint_printf("FAIL: eta =\n");
      fmpz_mat_print(eta);
      flint_printf("\ntransf_mat =\n");
      fmpz_mat_print_pretty(transf_mat);
      flint_printf("\ntest_mat =\n");
      fmpz_mat_print_pretty(test_mat);
      flint_abort();
    }

  fmpz_mat_J(eta);
  theta_transform_matrix(transf_mat, eta);
  check_dupont_J(test_mat);
  res = fmpz_mat_equal(transf_mat, test_mat);
  if (!res)
    {
      flint_printf("FAIL: eta =\n");
      fmpz_mat_print(eta);
      flint_printf("\ntransf_mat =\n");
      fmpz_mat_print_pretty(transf_mat);
      flint_printf("\ntest_mat =\n");
      fmpz_mat_print_pretty(test_mat);
      flint_abort();
    }
  
  fmpz_mat_one(eta);
  fmpz_one(fmpz_mat_entry(eta, 0, 2));
  theta_transform_matrix(transf_mat, eta);
  check_dupont_M1(test_mat);
  res = fmpz_mat_equal(transf_mat, test_mat);
  if (!res)
    {
      flint_printf("FAIL: eta =\n");
      fmpz_mat_print(eta);
      flint_printf("\ntransf_mat =\n");
      fmpz_mat_print_pretty(transf_mat);
      flint_printf("\ntest_mat =\n");
      fmpz_mat_print_pretty(test_mat);
      flint_abort();
    }
  
  fmpz_mat_one(eta);
  fmpz_one(fmpz_mat_entry(eta, 1, 3));
  theta_transform_matrix(transf_mat, eta);
  check_dupont_M2(test_mat);
  res = fmpz_mat_equal(transf_mat, test_mat);
  if (!res)
    {
      flint_printf("FAIL: eta =\n");
      fmpz_mat_print(eta);
      flint_printf("\ntransf_mat =\n");
      fmpz_mat_print_pretty(transf_mat);
      flint_printf("\ntest_mat =\n");
      fmpz_mat_print_pretty(test_mat);
      flint_abort();
    }
  
  fmpz_mat_one(eta);
  fmpz_one(fmpz_mat_entry(eta, 0, 3));
  fmpz_one(fmpz_mat_entry(eta, 1, 2));
  theta_transform_matrix(transf_mat, eta);
  check_dupont_M3(test_mat);
  res = fmpz_mat_equal(transf_mat, test_mat);
  if (!res)
    {
      flint_printf("FAIL: eta =\n");
      fmpz_mat_print(eta);
      flint_printf("\ntransf_mat =\n");
      fmpz_mat_print_pretty(transf_mat);
      flint_printf("\ntest_mat =\n");
      fmpz_mat_print_pretty(test_mat);
      flint_abort();
    }

  fmpz_mat_one(T);
  fmpz_one(fmpz_mat_entry(T, 0, 1));
  fmpz_mat_diagonal_symplectic(eta, T);
  theta_transform_matrix(transf_mat, eta);
  check_dupont_I(test_mat);
  res = fmpz_mat_equal(transf_mat, test_mat);
  if (!res)
    {
      flint_printf("FAIL: eta =\n");
      fmpz_mat_print(eta);
      flint_printf("\ntransf_mat =\n");
      fmpz_mat_print_pretty(transf_mat);
      flint_printf("\ntest_mat =\n");
      fmpz_mat_print_pretty(test_mat);
      flint_abort();
    }
  
  fmpz_mat_zero(aux);
  fmpz_set_si(fmpz_mat_entry(aux, 0, 1), -1);
  fmpz_set_si(fmpz_mat_entry(aux, 1, 0), 1);
  fmpz_set_si(fmpz_mat_entry(aux, 2, 3), -1);
  fmpz_set_si(fmpz_mat_entry(aux, 3, 2), 1);
  fmpz_mat_set(eta, aux);
  res = fmpz_mat_is_symplectic(eta);
  if (!res)
    {
      flint_printf("FAIL (Not symplectic):\n");
      fmpz_mat_print_pretty(aux);
      flint_abort();
    }
  theta_transform_matrix(transf_mat, eta);
  check_dupont_Sigma(test_mat);
  res = fmpz_mat_equal(transf_mat, test_mat);
  if (!res)
    {
      flint_printf("FAIL: eta =\n");
      fmpz_mat_print(eta);
      flint_printf("\ntransf_mat =\n");
      fmpz_mat_print_pretty(transf_mat);
      flint_printf("\ntest_mat =\n");
      fmpz_mat_print_pretty(test_mat);
      flint_abort();
    }
 
  fmpz_mat_clear(transf_mat);
  fmpz_mat_clear(test_mat);
  fmpz_mat_clear(aux);
  fmpz_mat_clear(T);
  fmpz_mat_clear(eta);

  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
  
