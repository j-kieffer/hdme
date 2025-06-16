/*
    Copyright (C) 2022 Jean Kieffer

    This file is part of the hdme library.

    hdme is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License (GPL v3). See
    LICENCE or <http://www.gnu.org/licenses/> for more details.
*/

#ifndef HECKE_H
#define HECKE_H

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_poly_mat.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/acb_mat.h>


#include "flint_compat.h"
#include "siegel.h"
#include "theta.h"
#include "igusa.h"
#include "hilbert.h"
#include "verbose.h"

#define HECKE_RED_TOL_BITS 50
#define HECKE_CHARPOLY_PREC_MUL 1.8
#define HECKE_SELECT_BASEPTS_TRIALS 3


/* Define a data structure hecke_t that contains all the necessary
   data related to Hecke correspondences over the complex numbers,
   i.e. data about all period matrices isogenous to a given one. */

typedef struct
{
  /* Context information */
  slong prec;
  acb_mat_t tau; /* Original period matrix */
  acb_ptr theta2_tau;
  acb_ptr I_tau;
  slong ell; /* Type of isogenies */

  /* Hilbert case */
  acb_ptr t1t2; /* Periods */
  fmpz_mat_t eta; /* Phi_R(t1, t2) is eta*tau */
  fmpz_poly_t beta; /* Type of isogenies */

  /* Data to be computed */
  slong nb; /* Number of isogenous matrices */
  fmpz_mat_struct* cosets; /* General symplectic matrices acting */
  acb_mat_struct* isog; /* Isogenous period matrices in fundamental domain */
  acb_mat_struct* stars; /* Cocycles C*tau+D */
  acb_ptr stardets; /* Determinants of stars */
  acb_ptr theta2; /* Values of theta constants at isog */
  acb_ptr I; /* Values of Igusa-Clebsch covariants at isog */

  slong norm_ind; /* Individual normalization factor */
  fmpz_t norm_all; /* Global normalization factor */
  slong prod_ec; /* When tau is a product of elliptic curves, number
		    of isogenous such products */
  acb_ptr I_norm; /* Normalized isogenous covariants */

} hecke_struct;

typedef hecke_struct hecke_t[1];


/* Access macros */

#define hecke_prec(H) ((H)->prec)
#define hecke_tau(H) ((H)->tau)
#define hecke_theta2_tau(H) ((H)->theta2_tau)
#define hecke_I_tau(H) ((H)->I_tau)
#define hecke_ell(H) ((H)->ell)
#define hecke_t1t2(H) ((H)->t1t2)
#define hecke_eta(H) ((H)->eta)
#define hecke_beta(H) ((H)->beta)
#define hecke_nb(H) ((H)->nb)

#define hecke_coset(H, k) (&(H)->cosets[(k)])
#define hecke_isog(H, k) (&(H)->isog[(k)])
#define hecke_star(H, k) (&(H)->stars[(k)])
#define hecke_stardet(H, k) (&(H)->stardets[(k)])
#define hecke_theta2(H, k) (&(H)->theta2[16*(k)])
#define hecke_all_I(H) ((H)->I)
#define hecke_I(H, k) (&(H)->I[4*(k)])
#define hecke_norm_ind(H) ((H)->norm_ind)
#define hecke_norm_all(H) ((H)->norm_all)
#define hecke_prod_ec(H) ((H)->prod_ec)
#define hecke_I_norm(H, k) (&(H)->I_norm[4*(k)])


/* Memory management */

void hecke_init(hecke_t H, slong nb);

void hecke_clear(hecke_t H);

void hecke_print(const hecke_t H, slong digits);

void hecke_check_nb(const hecke_t H, slong nb);


/* Set t1t2, tau, I_tau, theta2_tau */

int hecke_set_tau(hecke_t H, const acb_mat_t tau, slong prec);

int hecke_set_I_fmpz(hecke_t H, fmpz* I, slong prec);

int hecke_set_t1t2(hecke_t H, acb_srcptr t, slong delta, slong prec);

int hecke_set_I_fmpz_hilbert(hecke_t H, fmpz* I, slong delta, slong prec);

int hecke_set_I_fmpz_with_lowprec(hecke_t H, fmpz* I, acb_srcptr th2_lp,
				  slong prec);


/* Set whole Hecke correspondence */

int hecke_set_entry(hecke_t H, slong k, const fmpz_mat_t gamma, slong prec);

void hecke_collect_verbose_start(slong nb);

void hecke_collect_print_status(int res, slong k, slong nb);


slong siegel_nb_cosets(slong ell);

void siegel_coset(fmpz_mat_t m, slong k, slong ell);

int hecke_collect_siegel(hecke_t H, slong ell, slong prec);


slong hilbert_nb_cosets(slong ell, slong delta);

void hilbert_coset(fmpz_poly_mat_t m, slong k, const fmpz_poly_t beta,
		   slong ell, slong delta);

int hecke_collect_hilbert(hecke_t H, const fmpz_poly_t beta,
			  slong ell, slong delta, slong prec);

int hecke_collect_hilbert_sym(hecke_t H, slong ell, slong delta, slong prec);


slong siegel_nb_T1_cosets(slong ell);

void siegel_T1_coset(fmpz_mat_t m, slong k, slong ell);

int hecke_collect_T1(hecke_t H, slong ell, slong prec);


slong siegel_nb_T1_cosets_with_line(slong ell);

int siegel_T1_coset_contains_line_dual(const fmpz_mat_t m, const fmpz_mat_t L, slong ell);

int hecke_collect_T1_with_line(hecke_t H, const fmpz_mat_t L,
			       slong ell, slong prec);

/* Normalization */

void hecke_normalize_entry(hecke_t H, slong k,
			   fmpz* I, slong norm_ind, slong prec);

void hecke_make_integral(hecke_t H, fmpz* I, slong prec);

int hecke_has_integral_precision(hecke_t H, slong prec);

slong hecke_integral_highprec(const hecke_t H, slong prec);

int hecke_all_isog_Q(slong* nb_roots, fmpz* all_I, const hecke_t H, fmpz* I, slong prec);


/* Hecke operators */

void hecke_slash(acb_ptr im, const acb_mat_t star, acb_srcptr val,
		 slong k, slong j, slong prec);

void hecke_slash_scalar(acb_t im, const acb_t stardet, const acb_t val,
			slong k, slong prec);

void hecke_operator(acb_ptr im, const hecke_t H, acb_srcptr val,
		    slong ell, slong k, slong j, slong prec);

void hecke_eigenvalue_eisenstein(fmpz_t eig, slong k, slong ell);

void hecke_eigenvalues_eisenstein_2(fmpz* eig, slong k, slong ell);


/* Hecke tables */

void hecke_action_all_monomials(acb_ptr r, const hecke_t H, slong wt,
				slong ell, slong prec);

int hecke_basis_matrix(acb_mat_t basis, slong nb, const acb_mat_struct* pts,
		       slong wt, slong prec);

int hecke_image_matrix(acb_mat_t image, slong nb, const acb_mat_struct* pts,
		       slong wt, slong ell, slong prec);

int hecke_select_basepts(acb_mat_struct* pts, acb_mat_t basis_inv,
			 slong wt, slong prec);

slong hecke_charpoly_startprec(slong wt);

slong hecke_charpoly_nextprec(slong prec);

int hecke_charpoly(fmpz_poly_t pol, slong ell, slong wt);

#endif

