/*
    Copyright (C) 2021 Jean Kieffer

    This file is part of the hdme library.

    hdme is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License (GPL v3). See
    LICENCE or <http://www.gnu.org/licenses/> for more details.
*/

#ifndef IGUSA_H
#define IGUSA_H

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include <flint/fmpq.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpq_poly.h>
#include <flint/acb_poly.h>
#include <flint/acb_modular.h>
#include <flint/fmpq_mpoly.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpq_vec.h>
#include <flint/fmpz_factor.h>


#include "flint_compat.h"
#include "hdme_data.h"
#include "siegel.h"
#include "theta.h"
#include "verbose.h"

#define IGUSA_WEIGHTS {4,6,10,12}
#define IGUSA_HALFWEIGHTS {2,3,5,6}
#define IC_WEIGHTS {2,4,6,10}
#define IC_HALFWEIGHTS {1,2,3,5}
/* #define IGUSA_NB 4 would clutter the code. */
#define X_WEIGHTS {4,6,10,12,12,16,18,24,28,30,36,40,42,48}
#define X_NB 14

#define COV_FACTOR_BITS 25

#define THOMAE_LOWPREC 50
#define THOMAE_MULPREC 8


/* We call I the following vector of covariants:
   psi4 = I4/4
   psi6 = I6prime/4 (Streng's notation)
   chi10 = -I10/2^12
   chi12 = I12/2^15

   They correspond to the Siegel modular forms with the following
   normalized q-expansions:
   psi4 = 1 + 240(q1+q2) + ...
   psi6 = 1 - 504(q1+q1) + ...
   chi10 = (q3 - 2 + q3^-1) + ...
   chi12 = (q3 + 10 + q3^-1) + ...

   This normalization differs slightly from Igusa's, who divides
   further chi10 and chi12 by 4 and 12, respectively. */

/* Igusa covariants from theta constants */

void igusa_h4(acb_t h4, acb_srcptr theta2, slong prec);

void igusa_h6(acb_t h6, acb_srcptr theta2, slong prec);

void igusa_h10(acb_t h10, acb_srcptr theta2, slong prec);

void igusa_h12(acb_t h12, acb_srcptr theta2, slong prec);

void igusa_h16(acb_t h16, acb_srcptr theta2, slong prec);

#define igusa_psi4(I) &(I)[0]
#define igusa_psi6(I) &(I)[1]
#define igusa_chi10(I) &(I)[2]
#define igusa_chi12(I) &(I)[3]

void igusa_from_theta2(acb_ptr I, acb_srcptr theta2, slong prec);

int igusa_from_tau(acb_ptr I, const acb_mat_t tau, slong prec);


/* Handle general projective coordinates */

void cov_rescale(acb_ptr I, acb_srcptr S, const acb_t scal,
		 slong nb, slong* weights, slong prec);

void cov_rescale_fmpz(fmpz* I, fmpz* S, const fmpz_t scal,
		      slong nb, slong* weights);

int cov_divisible_fmpz(fmpz* I, const fmpz_t scal, slong nb, slong* weights);

void cov_divexact_fmpz(fmpz* I, fmpz* S, const fmpz_t scal,
		       slong nb, slong* weights);

void cov_rescale_fmpz_si(fmpz* I, fmpz* S, slong scal, slong nb, slong* weights);

void cov_divexact_fmpz_si(fmpz* I, fmpz* S, slong scal, slong nb, slong* weights);

void cov_normalize(acb_ptr I, acb_srcptr S, slong nb, slong* weights, slong prec);

#define cov_factor_nb(fac) ((fac)->num)
#define cov_factor_p(fac, k) (&(fac)->p[k])

void cov_factors(fmpz_factor_t fac, fmpz* I, slong nb);

void cov_valuation_fmpq(fmpq_t val, fmpq* I, const fmpz_t p, slong nb, slong* weights);

void cov_adjust_weights(slong* adj, slong* weights, fmpz* I, slong nb);

void cov_normalize_fmpz(fmpz* I, fmpz* S, slong nb, slong* weights);

void cov_normalize_fmpz_wt1(fmpz* I, fmpz* S, slong nb);

void cov_normalize_fmpq_wt1(fmpz* I, fmpq* S, slong nb);

void cov_min_weight_combination(slong* wt, slong* exponents, fmpz* I,
				slong nb, slong* weights);

void cov_find_rescaling(acb_t scal, acb_srcptr I, fmpz* S,
		       slong nb, slong* weights, slong prec);

int cov_no_rescale_to_one(acb_srcptr I, slong nb, slong* weights,
			  slong prec);

int cov_distinct(acb_srcptr I1, acb_srcptr I2,
		 slong nb, slong* weights, slong prec);

slong cov_height(fmpz* I, slong nb, slong* weights);


/* Weighted polynomials in general covariants */

void cov_mpoly_eval(acb_t ev, const fmpz_mpoly_t pol, acb_srcptr I,
		    const fmpz_mpoly_ctx_t ctx, slong prec);

void cov_mpoly_eval_fmpz(fmpz_t ev, const fmpz_mpoly_t pol, fmpz* I,
			 const fmpz_mpoly_ctx_t ctx);

void cov_monomial(fmpz_mpoly_t mon, slong* exps, const fmpz_mpoly_ctx_t ctx);

void cov_monomial_degrees(slong* exps, const fmpz_mpoly_t mon,
			  const fmpz_mpoly_ctx_t ctx);

slong cov_nb_monomials(slong wt, slong nb, slong* weights);

void cov_all_exps(slong* exps, slong wt, slong nb, slong* weights);

void cov_eval_all_monomials(acb_ptr ev, acb_srcptr I, slong wt,
			    slong nb, slong* weights, slong prec);


/* Weighted polynomials in Igusa covariants */

slong igusa_nb_base_monomials(slong wt);

void igusa_base_exps(slong* exps, slong wt, slong k);

void igusa_base_monomial(fmpz_mpoly_t mon, slong wt, slong k,
			 const fmpz_mpoly_ctx_t ctx);

void igusa_from_monomials_zeroes(int* z4, int* z6, int* z10, int* z12,
				 fmpz* M, slong wt);

void igusa_from_monomials_exps(slong* e4, slong* e6, slong* e10, slong* e12,
			       int z4, int z6, int z10, int z12, slong wt);

void igusa_from_monomials(fmpz* I, fmpz* M, slong wt);


void igusa_print_coordinate(const fmpz_mpoly_t pol, const fmpz_mpoly_ctx_t ctx);

void igusa_try_coordinate(fmpz_mpoly_t pol, slong wt, slong j,
			  const fmpz_mpoly_ctx_t ctx);

void igusa_make_integral(fmpq_t scal, fmpz* I, slong wt, slong corr_2_3);


/* Different covariants:
   - Streng I4, I6prime, I10, I12,
   - classical Igusa--Clebsch I2, I4, I6, I10,
   - Clebsch A, B, C, D
   - Igusa Y12, X16, ..., X48 */

void igusa_streng(acb_ptr S, acb_srcptr I, slong prec);

void igusa_streng_fmpz(fmpz* S, fmpz* I);

void igusa_from_streng(acb_ptr I, acb_srcptr S, slong prec);

void igusa_from_streng_fmpz(fmpz* I, fmpz* S);

void igusa_IC(acb_ptr IC, acb_srcptr I, slong prec);

void igusa_IC_fmpz(fmpz* IC, fmpz* I);

void igusa_from_IC(acb_ptr I, acb_srcptr IC, slong prec);

void igusa_from_IC_fmpz(fmpz* I, fmpz* IC);

void igusa_ABCD_from_IC(acb_ptr ABCD, acb_srcptr IC, slong prec);

void igusa_R2_from_IC(acb_t res, acb_srcptr IC, slong prec);

void igusa_ABCD_from_IC_fmpz(fmpq* ABCD, fmpz* IC);

void igusa_R2_from_IC_fmpz(fmpq_t R2, fmpz* IC);

void igusa_X(fmpq* X, fmpz* I);


/* Links with genus 2 curves */

void curve_coeffs(acb_ptr ai, const acb_poly_t crv);

void curve_coeffs_fmpz(fmpz* ai, const fmpz_poly_t crv);

void igusa_from_curve(acb_ptr I, const acb_poly_t crv, slong prec);

void igusa_from_curve_fmpz(fmpz* I, const fmpz_poly_t crv);

int igusa_has_generic_automorphisms(acb_srcptr IC, slong prec);

int igusa_is_g2_curve(acb_srcptr I);

int igusa_is_g2_curve_fmpz(fmpz* I);


/* Mestre's algorithm: use I2, I4, I6, I10 */

void igusa_generic_randtest(acb_poly_t crv, acb_ptr IC, flint_rand_t state, slong prec);

void mestre_conic(acb_ptr conic, acb_srcptr ABCD, const acb_t U, const acb_t I10, slong prec);

void mestre_conic_randtest(acb_ptr conic, flint_rand_t state, slong prec);

void mestre_cubic(acb_ptr cubic, acb_srcptr ABCD, const acb_t U, const acb_t I10, slong prec);

void mestre_subst_in_conic(acb_poly_t subst, const acb_poly_t x1, const acb_poly_t x2,
			   const acb_poly_t x3, acb_srcptr conic, slong prec);

void mestre_subst_in_cubic(acb_poly_t subst, const acb_poly_t x1, const acb_poly_t x2,
			   const acb_poly_t x3, acb_srcptr cubic, slong prec);

int mestre_point_on_conic(acb_ptr pt, acb_srcptr conic, slong prec);

int mestre_point_is_outside_conic(acb_srcptr pt, acb_srcptr conic, slong prec);

void mestre_eval_cubic(acb_t res, acb_srcptr pt, acb_srcptr cubic, slong prec);

void mestre_parametrize_conic(acb_poly_t x1, acb_poly_t x2, acb_poly_t x3,
			     acb_srcptr pt, acb_srcptr conic, slong prec);

int mestre(acb_poly_t crv, acb_srcptr IC, slong prec);


/* Cardona's algorithm for curves with extra involution, Bolza classification */

void cardona_conic(acb_ptr Aij, acb_srcptr ABCD, slong prec);

void cardona_cubic(acb_ptr aijk, acb_srcptr ABCD, slong prec);

void cardona(acb_poly_t crv, acb_srcptr IC, slong prec);

int mestre_fmpz(acb_poly_t crv, fmpz* IC, slong prec);


/* Thomae's formulae: back to I4, I6prime, I10, I12 */

int thomae_roots(acb_ptr roots, const acb_poly_t crv, slong prec);

void thomae_reorder(acb_ptr new_roots, acb_srcptr roots, slong perm);

void thomae_rosenhain(acb_ptr ros, acb_srcptr roots, slong prec);

void thomae_theta4(acb_ptr th4, acb_srcptr ros, slong prec);

void thomae_theta2(acb_ptr th2, acb_srcptr th4, acb_srcptr ros, slong signs, slong prec);

int thomae_discard(acb_srcptr th2, slong prec);

int thomae_keep_candidate(const acb_mat_t tau, acb_srcptr I, slong prec);

slong thomae_startprec(slong prec);

int thomae_correct_signs(slong* perm, slong* signs, acb_srcptr roots,
			 acb_srcptr I, slong prec);

int thomae_correct_signs_with_lowprec(slong* perm, slong* signs, acb_srcptr roots,
				      acb_srcptr I, acb_srcptr th2_lp, slong prec);

int tau_theta2_from_curve(acb_mat_t tau, acb_ptr th2,
			  const acb_poly_t crv, slong prec);

int tau_theta2_from_curve_with_lowprec(acb_mat_t tau, acb_ptr th2, const acb_poly_t crv,
				       acb_srcptr th2_lp, slong prec);


/* Period computations from Igusa covariants */

int tau_from_igusa(acb_mat_t tau, acb_srcptr I, slong prec);

int tau_theta2_from_igusa(acb_mat_t tau, acb_ptr th2, acb_srcptr I, slong prec);

void igusa_ec_j1j2(acb_ptr j, fmpz* I, slong prec);

int igusa_ec_possible_kp2(acb_ptr kp2, const acb_t j, slong prec);

int igusa_ec_period(acb_t tau, const acb_t j, slong prec);

int tau_theta2_from_igusa_ec(acb_mat_t tau, acb_ptr th2, fmpz* I, slong prec);

int tau_theta2_from_igusa_fmpz(acb_mat_t tau, acb_ptr th2, fmpz* I, slong prec);

int tau_theta2_from_igusa_fmpz_with_lowprec(acb_mat_t tau, acb_ptr th2,
					    fmpz* I, acb_srcptr th2_lp, slong prec);


#endif
