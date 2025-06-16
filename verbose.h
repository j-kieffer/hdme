#ifndef VERBOSE_H
#define VERBOSE_H

#include <stdio.h>
#include "timing.h"

#define report(message, s) {\
  double utime, wtime;\
  timestamp_report(&utime, &wtime, &s);\
  printf("(%s) Done in user: %.3f wall: %.3f cpu: %.2f\n", message, utime, wtime, utime/wtime);\
}
#define report_end(s) report(__func__, s)


// there must be a smart way to do this with the precompiler...
int set_borchardt_verbose(int i);
int get_borchardt_verbose(void);

int set_hecke_verbose(int i);
int get_hecke_verbose(void);

int set_hilbert_lll_verbose(int i);
int get_hilbert_lll_verbose(void);

int set_modeq_verbose(int i);
int get_modeq_verbose(void);

int set_pol_verbose(int i);
int get_pol_verbose(void);

int set_siegel_verbose(int i);
int get_siegel_verbose(void);

int set_theta_verbose(int i);
int get_theta_verbose(void);

int set_thomae_verbose(int i);
int get_thomae_verbose(void);


void progress_bar(slong current, slong total, const char prefix[]);

#endif


