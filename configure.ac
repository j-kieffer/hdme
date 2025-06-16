m4_define([HDME_VERSION_STRING], m4_esyscmd_s([cat VERSION | xargs]))
AC_INIT([hdme], [HDME_VERSION_STRING], [https://github.com/edgarcosta/hdme/issues], [hdme])
AC_PREREQ([2.62])
AC_CONFIG_FILES([Makefile hdme.pc])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_SRCDIR([hdme_data.h])

# Required programs
AC_PROG_CC
AC_PROG_LIBTOOL
LT_INIT([disable-static])
AC_PROG_MAKE_SET
AC_CHECK_TOOL([AR], [ar])
AC_PATH_PROG([PKG_CONFIG], [pkg-config])

# Version numbers
HDME_VERSION="HDME_VERSION_STRING"
tail="$HDME_VERSION"
HDME_MAJOR="${tail%%.*}"
tail="${tail#*.}"
HDME_MINOR="${tail%%.*}"
tail="${tail#*.}"
HDME_PATCH="${tail%%-*}"
AC_SUBST([HDME_MAJOR])
AC_SUBST([HDME_MINOR])
AC_SUBST([HDME_PATCH])
AC_SUBST([HDME_VERSION])

# Set default AR flags
ARFLAGS="${ARFLAGS:-rcs}"
AC_SUBST([ARFLAGS])

# System detection
host_os=`uname -s 2>/dev/null | tr 'A-Z' 'a-z'` || host_os="unknown"

# Initialize variables
LIBS_EXTRA=""
INCS_EXTRA=""

# Check for required libraries
AC_CHECK_LIB([m], [sin], [], [AC_MSG_ERROR([Math library not found])])

# Helper function to check library
# Parameters:
#   $1: Library name (e.g., "mpfr", "flint")
#   $2: Directory path variable (e.g., $MPFR_DIR, $FLINT_DIR)
#   $3: Function to check (e.g., "mpfr_init", "flint_version")
#   $4: Additional libraries needed for linking (e.g., "-lgmp", "-lgmp -lmpfr")
#   $5: PKG_CONFIG module name (e.g., "MPFR", "FLINT")
#
# This macro handles two cases:
#   1. When a specific directory is provided (--with-lib=PATH)
#   2. When no directory is provided (fallback to direct library check)
#
# For each case, it:
#   - Sets up appropriate library and include paths
#   - Performs linking tests with dependencies
#   - Updates LIBS_EXTRA and INCS_EXTRA variables
AC_DEFUN([HDME_CHECK_LIB], [
    if test "x$2" != "x"; then
        LIBS_EXTRA="$LIBS_EXTRA -L\"$2/lib\" -l$1"
        INCS_EXTRA="$INCS_EXTRA -I\"$2/include\""
        # Set flags for library test
        saved_LDFLAGS="$LDFLAGS"
        saved_CPPFLAGS="$CPPFLAGS"
        LDFLAGS="$LDFLAGS -L$2/lib"
        CPPFLAGS="$CPPFLAGS -I$2/include"
        AC_CHECK_LIB([$1], [$3], [], [
            AC_MSG_ERROR([$1 library found in $2 but linking test failed])
        ], [$4])
        LDFLAGS="$saved_LDFLAGS"
        CPPFLAGS="$saved_CPPFLAGS"
    else
        AC_CHECK_LIB([$1], [$3], [
            LIBS_EXTRA="$LIBS_EXTRA -l$1"
        ], [
            AC_MSG_ERROR([$1 is required but not found. Please install $1 or specify --with-$1=PATH])
        ], [$4])
    fi
])

# GMP library (required)
AC_ARG_WITH([gmp],
    [AS_HELP_STRING([--with-gmp=PATH], [specify location of GMP (required)])],
    [GMP_DIR=$withval], [GMP_DIR=""])

HDME_CHECK_LIB([gmp], [$GMP_DIR], [__gmpz_init], [], [GMP])

# MPFR library (required)
AC_ARG_WITH([mpfr],
    [AS_HELP_STRING([--with-mpfr=PATH], [specify location of MPFR (required)])],
    [MPFR_DIR=$withval], [MPFR_DIR=""])

HDME_CHECK_LIB([mpfr], [$MPFR_DIR], [mpfr_init], [-lgmp], [MPFR])

# FLINT library (required)
AC_ARG_WITH([flint],
    [AS_HELP_STRING([--with-flint=PATH], [specify location of FLINT (required)])],
    [FLINT_DIR=$withval], [FLINT_DIR=""])

HDME_CHECK_LIB([flint], [$FLINT_DIR], [flint_version], [-lgmp -lmpfr], [FLINT])

# Enable/disable options
AC_ARG_ENABLE([shared],
    [AS_HELP_STRING([--disable-shared], [disable building of a shared library (default is enabled)])],
    [], [enable_shared=yes])

AC_ARG_ENABLE([static],
    [AS_HELP_STRING([--enable-static], [enable building of a static library (default is disabled)])],
    [], [enable_static=no])

if test "x$enable_shared" = "xyes"; then
    AC_SUBST(BUILD_SHARED, 1)
else
    AC_SUBST(BUILD_SHARED, 0)
fi

if test "x$enable_static" = "xyes"; then
    AC_SUBST(BUILD_STATIC, 1)
else
    AC_SUBST(BUILD_STATIC, 0)
fi

# Enable/disable options
AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--enable-debug], [enable debug symbols])],
    [DEBUG=$enableval], [DEBUG=no])

AC_ARG_ENABLE([openmp],
    [AS_HELP_STRING([--enable-openmp], [enable OpenMP support (default)])
AS_HELP_STRING([--disable-openmp], [disable OpenMP support])],
    [OPENMP=$enableval], [OPENMP=yes])

# Set up CFLAGS
BASE_CFLAGS="-Wall -Wextra -O2 -std=c99"

if test "x$DEBUG" = "xyes"; then
    BASE_CFLAGS="$BASE_CFLAGS -g"
fi

# Handle OpenMP
SHARED_CFLAGS=""
if test "x$enable_shared" = "xyes"; then
    SHARED_CFLAGS="-fPIC"
fi
AC_SUBST(SHARED_CFLAGS)

if test "x$OPENMP" = "xyes"; then
    case $host_os in
        darwin*)
            BASE_CFLAGS="$BASE_CFLAGS -Xpreprocessor -fopenmp"
            # Try common Homebrew locations
            if test -d "/opt/homebrew/opt/libomp/lib"; then
                LIBS_EXTRA="$LIBS_EXTRA -L\"/opt/homebrew/opt/libomp/lib\" -lomp"
                INCS_EXTRA="$INCS_EXTRA -I\"/opt/homebrew/opt/libomp/include\""
            elif test -d "/usr/local/opt/libomp/lib"; then
                LIBS_EXTRA="$LIBS_EXTRA -L\"/usr/local/opt/libomp/lib\" -lomp"
                INCS_EXTRA="$INCS_EXTRA -I\"/usr/local/opt/libomp/include\""
            else
                AC_CHECK_LIB([omp], [omp_get_num_threads], [LIBS_EXTRA="$LIBS_EXTRA -lomp"])
            fi
            ;;
        *)
            BASE_CFLAGS="$BASE_CFLAGS -fopenmp"
            ;;
    esac
fi

# Set up library names based on OS
HDME_LIB_STATIC="libhdme.a"
case $host_os in
    darwin*)
        HDME_LIBNAME="libhdme.$HDME_VERSION.dylib"
        HDME_LIB="libhdme.dylib"
        EXEEXT=""
        SHARED_LDFLAGS="-dynamiclib -install_name '@rpath/$HDME_LIBNAME'"
        ;;
    cygwin*|mingw*)
        HDME_LIBNAME="libhdme.dll"
        HDME_LIB="libhdme.so.$HDME_VERSION"
        EXEEXT=""
        SHARED_LDFLAGS="-shared -Wl,-soname,${HDME_LIBNAME}.${HDME_MAJOR}"
        ;;
    *)
        HDME_LIBNAME="libhdme.so"
        HDME_LIB="libhdme.so.$HDME_VERSION"
        EXEEXT=""
        SHARED_LDFLAGS="-shared -Wl,-soname,${HDME_LIBNAME}.${HDME_MAJOR}"
        ;;
esac

# Build final flags
INCS="-I\"\$(CURDIR)\" $INCS_EXTRA"
LIBS="-L\"\$(CURDIR)\" -lm $LIBS_EXTRA"

# Set final CFLAGS if not provided by user
if test "x$CFLAGS" = "x"; then
    CFLAGS="$BASE_CFLAGS"
fi

# Substitute variables for Makefile
AC_SUBST([HDME_LIB])
AC_SUBST([HDME_LIBNAME])
AC_SUBST([HDME_LIB_STATIC])
AC_SUBST([SHARED_LDFLAGS])
AC_SUBST([INCS])
AC_SUBST([LIBS])
AC_SUBST([EXEEXT])
AC_SUBST([host_os])

AC_OUTPUT

dnl Shorten the original help message
m4_cleardivert([HELP_BEGIN])dnl
m4_divert_push([HELP_BEGIN])dnl

if test "$ac_init_help" = "long";
then
    cat <<_ACEOF
\`configure' configures m4_ifset([AC_PACKAGE_STRING],
			[AC_PACKAGE_STRING],
			[this package]) to adapt to many kinds of systems.

Usage: $[0] [[OPTION]]... [[VAR=VALUE]]...

To assign environment variables (e.g., CC, CFLAGS...), specify them as
VAR=VALUE.  See below for descriptions of some of the useful variables.

Defaults for the options are specified in brackets.

Configuration:
  -h, --help              display this help and exit
      --help=short        display options specific to this package
      --help=recursive    display the short help of all the included packages
  -V, --version           display version information and exit
  -q, --quiet, --silent   do not print \`checking ...' messages
      --cache-file=FILE   cache test results in FILE [disabled]
  -C, --config-cache      alias for \`--cache-file=config.cache'
  -n, --no-create         do not create output files
      --srcdir=DIR        find the sources in DIR [configure dir or \`..']

Installation directories:
]AS_HELP_STRING([--prefix=PREFIX],
  [install architecture-independent files in PREFIX [$ac_default_prefix]])
AS_HELP_STRING([--exec-prefix=EPREFIX],
  [install architecture-dependent files in EPREFIX [PREFIX]])[

By default, \`make install' will install all the files in
\`$ac_default_prefix/include' and \`$ac_default_prefix/lib' etc.
You can specify an installation prefix other than \`$ac_default_prefix'
using \`--prefix', for instance \`--prefix=\$HOME'.

For better control, use the options below.

Fine tuning of the installation directories:
  --bindir=DIR            user executables [EPREFIX/bin]
  --libdir=DIR            object code libraries [EPREFIX/lib]
  --includedir=DIR        C header files [PREFIX/include]
_ACEOF

    cat <<\_ACEOF]
m4_divert_pop([HELP_BEGIN])dnl

dnl NOTE: There is supposed to be a missing `fi' here.


