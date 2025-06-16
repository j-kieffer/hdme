/*
    Copyright (C) 2021 Jean Kieffer

    This file is part of the hdme library.

    hdme is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License (GPL v3). See
    LICENCE or <http://www.gnu.org/licenses/> for more details.
*/

#ifndef THETA_H
#define THETA_H

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>
#include <flint/ulong_extras.h>

#include "flint_compat.h"
#include "siegel.h"
#include "verbose.h"

#define THETA_NEWTON_MINPREC 100
#define THETA_NEWTON_LOSS 25
#define THETA_NEWTON_DERIVATIVE_OFFSET 10
#define THETA_NEWTON_BASEPREC 4000
#define THETA_NEWTON_Y1 5000
#define THETA_NEWTON_Y2MAX 10
#define THETA_NEWTON_TOL_EXP -5
#define THETA_DER_LOSS 25



/* Borchardt means */

int acb_sqrt_goodpos(acb_t r, const acb_t z, slong prec);

void borchardt_sqrt(acb_t r, const acb_t z, slong prec);

void borchardt_root_ui(acb_t r, const acb_t z, ulong e, slong prec);

int borchardt_step(acb_ptr b, acb_srcptr a, slong prec);

void borchardt_mean_m0(arb_t m0, acb_srcptr a, slong prec);

void borchardt_mean_M0(arb_t M0, acb_srcptr a, slong prec);

void borchardt_mean_Delta0(arb_t Delta0, acb_srcptr a, slong prec);

int borchardt_mean_nb_steps_before_quad_conv(fmpz_t nb, acb_srcptr a, slong prec);

void borchardt_mean_nb_steps_after_quad_conv(fmpz_t nb, acb_srcptr a, slong prec);

int borchardt_mean_quad_conv_is_reached(acb_srcptr a, slong prec);

int borchardt_mean(acb_t r, acb_srcptr a, slong prec);

void borchardt_excl_half_planes(arf_struct* b, const acb_t z, slong prec);

int borchardt_mean_invalid(acb_srcptr a, slong prec);


/* Theta characteristics */

ulong theta_char_get_a(ulong ch, slong g);

ulong theta_char_get_b(ulong ch, slong g);

ulong theta_char_set_ab(ulong a, ulong b, slong g);

int theta_char_dot_product(ulong a, ulong b, slong g);

int theta_char_is_even(ulong ch, slong g);

slong theta_char_get_label_g2(ulong ch);

ulong theta_char_set_label_g2(slong label);

int theta_char_is_goepel(ulong ch1, ulong ch2, ulong ch3, ulong ch4, slong g);

int theta_char_is_syzygous(ulong ch1, ulong ch2, ulong ch3, slong g);


/* Theta constants */

void theta_duplication(acb_ptr th2_2tau, acb_srcptr th_tau, slong prec);

int theta_0123_naive_B(fmpz_t B, const acb_mat_t tau, slong prec);

int theta_0123_naive(acb_ptr th, const acb_mat_t tau, slong prec);

int theta2_naive(acb_ptr th, const acb_mat_t tau, slong prec);

int theta2_inverse(acb_mat_t tau, acb_srcptr th, slong prec);

int theta2_inverse_no_sqrt(acb_mat_t tau, acb_srcptr th, slong prec);

int theta2_invalid(acb_srcptr th2, slong prec);

int theta_0123half_diff_naive(acb_mat_t dth, const acb_mat_t tau, slong prec);

int theta_0123half_inverse(acb_mat_t tau, acb_srcptr th_half, slong prec);

int theta_0123half_inverse_no_sqrt(acb_mat_t tau, acb_srcptr th_half, slong prec);

int theta_0123half_inverse_diff(acb_mat_t dtau, const acb_mat_t tau, acb_srcptr th_half,
				slong prec);

int theta2_newton_step(acb_ptr th_half, const acb_mat_t tau, acb_srcptr th_half_approx,
		       slong prec);

slong theta2_newton_start_prec(slong prec);

int theta2_newton(acb_ptr th2, const acb_mat_t tau, slong prec);

slong theta_newton_k2(acb_mat_t w, const acb_mat_t z, slong prec);

slong theta_newton_k1(acb_mat_t w, const acb_mat_t z, slong prec);

int theta_use_naive(const acb_mat_t tau, slong prec);

int theta_use_newton(const acb_mat_t tau, slong prec);

ulong theta_transform_image_char(fmpz_t epsilon, ulong ch, const fmpz_mat_t eta);

void theta_transform_matrix(fmpz_mat_t res, const fmpz_mat_t eta);

void theta_transform(acb_ptr th_eta, const fmpz_mat_t eta, acb_srcptr th, slong prec);

void theta2_transform(acb_ptr th2_eta, const fmpz_mat_t eta, acb_srcptr th2, slong prec);

int theta2_unif(acb_ptr th2, const acb_mat_t tau, slong prec);

int theta2_renormalize(acb_ptr th2, acb_srcptr th2_proj, slong prec);

void theta2_randtest(acb_ptr theta2, flint_rand_t state, slong prec);


/* Derivatives of theta constants */

void theta_der_set_pert(arb_t eps, slong prec);

int theta_der_set_error(mag_t error, const acb_mat_t tau, slong prec);

int theta2_der_naive(acb_ptr th2_tau, acb_mat_t dth2_tau,
		      const acb_mat_t tau, slong prec);

int theta_0123_der_naive(acb_ptr th, acb_mat_t dth,
			 const acb_mat_t tau, slong prec);

void theta_der_duplication(acb_ptr th2_2tau, acb_mat_t dth2_2tau,
			   acb_srcptr th_tau, const acb_mat_t dth_tau,
			   slong prec);

int theta2_der_newton_step(acb_ptr th_half, acb_mat_t dth_approx,
			   const acb_mat_t tau, acb_srcptr th_half_approx,
			   slong prec);

int theta2_der_newton(acb_ptr th2, acb_mat_t dth2, const acb_mat_t tau, slong prec);


#endif
