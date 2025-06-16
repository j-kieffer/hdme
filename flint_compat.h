#ifndef FLINT_COMPAT_H
#define FLINT_COMPAT_H

#include <flint/flint.h>
// see https://github.com/flintlib/flint/commit/b6bd88951edf6409d225491859da9bbfce859c07
// FLINT < 3.2.0 uses flint_randinit/flint_randclear
// FLINT >= 3.2.0 uses flint_rand_init/flint_rand_clear
#if __FLINT_VERSION == 3 && __FLINT_VERSION_MINOR < 2
#define flint_rand_init flint_randinit
#define flint_rand_clear flint_randclear
#endif
// For FLINT 3.2.0+, no mapping needed as function names are already correct

#endif /* FLINT_COMPAT_H */ 