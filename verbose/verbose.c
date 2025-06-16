// there must be a smart way to do this with the precompiler...

int set_borchardt_verbose(int i) {
#ifndef BORCHARDT_VERBOSE
  static int verbose = 0;
#else
  static int verbose = BORCHARDT_VERBOSE;
#endif
  if (i >= 0) verbose = i;
  return verbose;
}

int get_borchardt_verbose(void) {
  return set_borchardt_verbose(-1);
}


int set_hecke_verbose(int i) {
#ifndef HECKE_VERBOSE
  static int verbose = 1;
#else
  static int verbose = HECKE_VERBOSE;
#endif
  if (i >= 0) verbose = i;
  return verbose;
}

int get_hecke_verbose(void) {
  return set_hecke_verbose(-1);
}


int set_hilbert_lll_verbose(int i) {
#ifndef HILBERT_LLL_VERBOSE
  static int verbose = 0;
#else
  static int verbose = HILBERT_LLL_VERBOSE;
#endif
  if (i >= 0) verbose = i;
  return verbose;
}

int get_hilbert_lll_verbose(void) {
  return set_hilbert_lll_verbose(-1);
}


int set_modeq_verbose(int i) {
#ifndef MODEQ_VERBOSE
  static int verbose = 1;
#else
  static int verbose = MODEQ_VERBOSE;
#endif
  if (i >= 0) verbose = i;
  return verbose;
}

int get_modeq_verbose(void) {
  return set_modeq_verbose(-1);
}


int set_pol_verbose(int i) {
#ifndef POL_VERBOSE
  static int verbose = 1;
#else
  static int verbose = POL_VERBOSE;
#endif
  if (i >= 0) verbose = i;
  return verbose;
}

int get_pol_verbose(void) {
  return set_pol_verbose(-1);
}


int set_siegel_verbose(int i) {
#ifndef SIEGEL_VERBOSE
  static int verbose = 1;
#else
  static int verbose = SIEGEL_VERBOSE;
#endif
  if (i >= 0) verbose = i;
  return verbose;
}

int get_siegel_verbose(void) {
  return set_siegel_verbose(-1);
}


int set_theta_verbose(int i) {
#ifndef THETA_VERBOSE
  static int verbose = 0;
#else
  static int verbose = THETA_VERBOSE;
#endif
  if (i >= 0) verbose = i;
  return verbose;
}

int get_theta_verbose(void) {
  return set_theta_verbose(-1);
}


int set_thomae_verbose(int i) {
#ifndef THOMAE_VERBOSE
  static int verbose = 0;
#else
  static int verbose = THOMAE_VERBOSE;
#endif
  if (i >= 0) verbose = i;
  return verbose;
}

int get_thomae_verbose(void) {
  return set_thomae_verbose(-1);
}

