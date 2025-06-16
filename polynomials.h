/*
    Copyright (C) 2022 Jean Kieffer

    This file is part of the hdme library.

    hdme is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License (GPL v3). See
    LICENCE or <http://www.gnu.org/licenses/> for more details.
*/

#ifndef HDME_POLY_H
#define HDME_POLY_H

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mod.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/acb_poly.h>

#include "flint_compat.h"
#include "verbose.h"


#define HDME_RD_RADIUS_PREC 100

int acb_vec_contains_int(acb_srcptr v, slong nb);

void acb_poly_product_tree_1(acb_poly_t P, acb_srcptr xi,
			     acb_srcptr yi, slong d, slong prec);

void acb_poly_product_tree_2(acb_poly_t Q, acb_srcptr xi, acb_srcptr yi,
			     acb_srcptr zi, slong d, slong prec);

int acb_round(fmpz_t c, arf_t radius, const acb_t x);

int acb_poly_round(fmpz_poly_t pol, arf_t max_radius,
		   const acb_poly_t pol_acb, slong degree);

int acb_rationalize(fmpq_t c, fmpz_t den, const acb_t x,
		    const fmpz_t probable_den, slong prec);

int acb_poly_rationalize(fmpq_poly_t pol, const acb_poly_t pol_acb,
			 slong degree, slong prec);

void pol_simplify(fmpz_poly_struct* num_vec, fmpz_t den, slong degree, slong nb);

void pol_factor_Q(slong* nb_factors, fmpz_poly_struct* factors, slong* exps,
		  const fmpz_poly_t pol);

void pol_roots_Q(slong* nb_roots, fmpq* roots, slong* mults,
		 const fmpz_poly_t pol);

void pol_remove_factor_Q(fmpz_poly_t r, const fmpz_poly_t pol,
			 const fmpz_poly_t fac, slong mult);

void pol_remove_root_Q(fmpz_poly_t r, const fmpz_poly_t pol,
		       const fmpq_t root, slong mult);

int pol_reduce_Fp(fmpz_mod_poly_t red, const fmpz_poly_t num,
		  const fmpz_t den, const fmpz_mod_ctx_t ctx);

void pol_factor_Fp(slong* nb_factors, fmpz_mod_poly_struct* factors, slong* exps,
		   const fmpz_mod_poly_t pol, const fmpz_mod_ctx_t ctx);

void pol_roots_Fp(slong* nb_roots, fmpz* roots, slong* mults,
		  const fmpz_mod_poly_t pol, const fmpz_mod_ctx_t ctx);

void pol_remove_root_Fp(fmpz_mod_poly_t r, const fmpz_mod_poly_t pol,
			const fmpz_t root, slong mult, const fmpz_mod_ctx_t ctx);

#endif
