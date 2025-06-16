/*
    Copyright (C) 2021 Jean Kieffer

    This file is part of the hdme library.

    hdme is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License (GPL v3). See
    LICENCE or <http://www.gnu.org/licenses/> for more details.
*/

#ifndef SIEGEL_H
#define SIEGEL_H

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/arb.h>
#include <flint/arb_mat.h>
#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/ulong_extras.h>

#include "flint_compat.h"


/* Additional functions for real and complex matrices */

void acb_mat_get_real(arb_mat_t re, const acb_mat_t z);

void acb_mat_get_imag(arb_mat_t im, const acb_mat_t z);

void acb_mat_set_arb_arb(acb_mat_t z, const arb_mat_t re, const arb_mat_t im);

void acb_mat_set_window(acb_mat_t z, slong j, slong k, const acb_mat_t w);

void arb_mat_randtest_precise(arb_mat_t r, flint_rand_t state, slong prec,
			      slong mag_bits);

void arb_mat_randtest_sym_precise(arb_mat_t r, flint_rand_t state, slong prec,
				  slong mag_bits);

void arb_mat_randtest_sym_pos(arb_mat_t r, flint_rand_t state, slong prec);

int arb_mat_is_nonsymmetric(const arb_mat_t m);

void arb_mat_congr_fmpz_mat(arb_mat_t r, const fmpz_mat_t u, const arb_mat_t m,
			    long prec);

int arb_mat_is_minkowski_reduced(const arb_mat_t r, const arb_t tol, slong prec);

int arb_mat_not_minkowski_reduced(const arb_mat_t r, slong prec);

int arb_mat_minkowski_reduce(arb_mat_t r, fmpz_mat_t u, const arb_mat_t m,
			 const arb_t tol, slong prec);

void arb_mat_lambda(arb_t lambda, const arb_mat_t m, slong prec);


/* Symplectic matrices */

slong fmpz_mat_half_dim(const fmpz_mat_t m);

void fmpz_mat_direct_inv(fmpz_mat_t minv, const fmpz_mat_t m);

void fmpz_mat_get_a(fmpz_mat_t a, const fmpz_mat_t m);

void fmpz_mat_get_b(fmpz_mat_t a, const fmpz_mat_t m);

void fmpz_mat_get_c(fmpz_mat_t a, const fmpz_mat_t m);

void fmpz_mat_get_d(fmpz_mat_t a, const fmpz_mat_t m);

void fmpz_mat_set_abcd(fmpz_mat_t m,
		       const fmpz_mat_t a, const fmpz_mat_t b,
		       const fmpz_mat_t c, const fmpz_mat_t d);

void fmpz_mat_J(fmpz_mat_t m);

int fmpz_mat_is_J(const fmpz_mat_t m);

int fmpz_mat_is_scalar(const fmpz_mat_t m);

int fmpz_mat_is_symplectic(const fmpz_mat_t m);

int fmpz_mat_is_general_symplectic(const fmpz_mat_t m);

void fmpz_mat_diagonal_symplectic(fmpz_mat_t m, const fmpz_mat_t u);

void fmpz_mat_randtest_triangular_symplectic(fmpz_mat_t m, flint_rand_t state, slong bits);

void fmpz_mat_randtest_diagonal_symplectic(fmpz_mat_t m, flint_rand_t state, slong bits);

void fmpz_mat_randtest_symplectic(fmpz_mat_t m, flint_rand_t state, slong bits);


/* The Siegel half space */

void siegel_halfspace_randtest(acb_mat_t z, flint_rand_t state, slong prec);

void siegel_star(acb_mat_t w, const fmpz_mat_t m, const acb_mat_t z, slong prec);

int siegel_transform(acb_mat_t w, const fmpz_mat_t m, const acb_mat_t z, slong prec);

int siegel_is_real_reduced(const acb_mat_t z, const arb_t tol, slong prec);

int siegel_not_real_reduced(const acb_mat_t z, slong prec);

int siegel_reduce_real(acb_mat_t w, fmpz_mat_t u, const acb_mat_t z,
		       const arb_t tol, slong prec);

slong siegel_nb_test_matrices(slong g);

void siegel_test_matrix(fmpz_mat_t u, slong j);

int siegel_fundamental_domain(acb_mat_t w, fmpz_mat_t m,
			      const acb_mat_t z, const arb_t tol, slong prec);

int siegel_is_in_fundamental_domain(const acb_mat_t z, const arb_t tol, slong prec);

int siegel_not_in_fundamental_domain(const acb_mat_t z, slong prec);

int siegel_is_weakly_reduced(const acb_mat_t z, const arb_t tol, slong prec);

void siegel_fundamental_domain_randtest(acb_mat_t z, flint_rand_t state, slong prec);


#endif 
