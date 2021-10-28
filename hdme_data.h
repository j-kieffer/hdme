
#ifndef HDME_DATA_H
#define HDME_DATA_H

#include <stdio.h>
#include "flint.h"
#include "mpoly.h"
#include "fmpq_mpoly.h"
#include "arb.h"
#include "acb.h"

#define HDME_DATA_STR_LEN 10000
#define HDME_DATA_FILE_LEN 1000
#define HDME_DATA_VAR_LEN 10
#define HDME_DATA_PATH HDME_PATH"/hdme_data"

char** hdme_data_vars_init(slong nb);

void hdme_data_vars_clear(char** vars, slong nb);

void hdme_data_vars_set(char** vars, const char* name, slong k);

void hdme_data_read(fmpq_mpoly_t pol, const char** vars, const char* name,
		    const fmpq_mpoly_ctx_t ctx);

void hdme_data_evaluate_acb(acb_t ev, const fmpq_mpoly_t pol, acb_srcptr vals,
			    const fmpq_mpoly_ctx_t ctx, slong prec);

#endif