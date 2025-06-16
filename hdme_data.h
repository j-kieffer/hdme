/*
    Copyright (C) 2021 Jean Kieffer

    This file is part of the hdme library.

    hdme is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License (GPL v3). See
    LICENCE or <http://www.gnu.org/licenses/> for more details.
*/

#ifndef HDME_DATA_H
#define HDME_DATA_H

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/acb.h>
#include <flint/arb.h>
#include <flint/fmpq_mpoly.h>
#include <flint/mpoly.h>
#include <stdio.h>

#ifndef HDME_PATH
#error "HDME_PATH must be defined"
#endif

#define HDME_DATA_STR_LEN 10000
#define HDME_DATA_FILE_LEN 1000
#define HDME_DATA_VAR_LEN 10
#define HDME_DATA_PATH HDME_PATH"/hdme_data"

char** hdme_data_vars_init(slong nb);

void hdme_data_vars_clear(char** vars, slong nb);

void hdme_data_vars_set(char** vars, const char* name, slong k);

void hdme_data_get_string(char* s, const char* key);

void hdme_data_read(fmpq_mpoly_t pol, const char** vars, const char* name,
		    const fmpq_mpoly_ctx_t ctx);

void hdme_data_evaluate_acb(acb_t ev, const fmpq_mpoly_t pol, acb_srcptr vals,
			    const fmpq_mpoly_ctx_t ctx, slong prec);

void hdme_data_evaluate_fmpq(fmpq_t ev, const fmpq_mpoly_t pol,
			     fmpq* vals, const fmpq_mpoly_ctx_t ctx);

#endif
