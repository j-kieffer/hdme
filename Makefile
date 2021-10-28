#
#   Copyright (C) 2021 Jean Kieffer
#
#   This file is part of the hdme library.
#
#   hdme is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License (GPL v3). See
#   LICENCE or <http://www.gnu.org/licenses/> for more details.
#


# Directories

# MODIFY THIS:
# CURDIR is where the source files lie
# MYDIR is where the other libraries are installed. See also ARB_DIR, etc.

CURDIR = /home/jean/code/hdme
MYDIR = /home/jean/install

PARITWINE_DIR = $(MYDIR)
PARITWINE_INCDIR = $(PARITWINE_DIR)/include
PARITWINE_LIBDIR = $(PARITWINE_DIR)/lib

ARB_DIR = $(MYDIR)
ARB_INCDIR = $(ARB_DIR)/include
ARB_LIBDIR = $(ARB_DIR)/lib

FLINT_DIR = $(MYDIR)
FLINT_INCDIR = $(FLINT_DIR)/include/flint
FLINT_LIBDIR = $(FLINT_DIR)/lib

PARI_DIR = $(MYDIR)
PARI_INCDIR = $(PARI_DIR)/include/pari
PARI_LIBDIR = $(PARI_DIR)/lib

MPFR_DIR = $(MYDIR)
MPFR_INCDIR = $(MPFR_DIR)/include
MPFR_LIBDIR = $(MPFR_DIR)/lib

GMP_DIR = $(MYDIR)
GMP_INCDIR = $(GMP_DIR)/include
GMP_LIBDIR = $(GMP_DIR)/lib

# DO NOT MODIFY
# Compiling and linking

SHELL = /bin/sh

CC = gcc
CFLAGS = -ansi -Wall -pedantic -Werror -g -DHDME_PATH=\"$(CURDIR)\"

AR = ar
ARFLAGS = -rcs

VALGRIND = valgrind

INCS = -I$(CURDIR) -I$(PARITWINE_INCDIR) -I$(ARB_INCDIR) -I$(FLINT_INCDIR) \
-I$(PARI_INCDIR) -I$(MPFR_INCDIR) -I$(GMP_INCDIR)

LLIBS = -L$(PARITWINE_LIBDIR) -L$(ARB_LIBDIR) -L$(FLINT_LIBDIR) -L$(PARI_LIBDIR) \
-L$(MPFR_INCDIR) -L$(GMP_LIBDIR)

lDEPS = arb flint mpfr gmp pthread m
lLIBS = $(patsubst %, -l%, $(lDEPS))


# Modules to be compiled

MODULES = hdme_data siegel theta igusa hilbert modular
SOURCES = $(wildcard $(patsubst %, %/*.c, $(MODULES)))
HEADERS = $(patsubst %, %.h, $(MODULES))
OBJS = $(patsubst %.c, build/%.o, $(SOURCES))

EXMP_SOURCES = $(wildcard examples/*.c)
EXMP_NAMES = $(patsubst %.c, %, $(EXMP_SOURCES))

TEST_SOURCES = $(wildcard $(patsubst %, %/test/t-*.c, $(MODULES)))
TESTS = $(patsubst %.c, build/test/%, $(TEST_SOURCES))
RUNTESTS = $(patsubst build/test/%, run/%, $(TESTS))
RUNVALGRIND = $(patsubst build/test/%, valgrind/%, $(TESTS))

TIMING_SOURCES = $(wildcard $(patsubst %, %/time/time-*.c, $(MODULES)))
TIMINGS = $(patsubst %.c, build/time/%, $(TIMING_SOURCES))

# Timing data is written in TIMEDIR

TIMEDIR = $(CURDIR)/time

# Targets

all: library

check: library tests $(RUNTESTS)

valgrind: library tests $(RUNVALGRIND)

tests: library $(TESTS)

timings: library $(TIMINGS)

clean:
	rm -rf build/
	rm libhdme.a

library: libhdme.a

libhdme.a: $(OBJS) $(SOURCES) $(HEADERS) | build
	@echo AR libhdme.a
	@$(AR) $(ARFLAGS) libhdme.a $(OBJS)

examples: library

build:
	@mkdir -p build build/test build/time \
	$(foreach mod, $(MODULES), build/$(mod)) \
	$(foreach mod, $(MODULES), build/test/$(mod)/test) \
	$(foreach mod, $(MODULES), build/time/$(mod)/time) 

build/%.o: %.c $(HEADERS) | build
	@echo CC $<
	@$(CC) $(CFLAGS) $(INCS) -c $< -o $@

run/%: build/test/%
	@$<

valgrind/%: build/test/%
	$(VALGRIND) --leak-check=full --track-origins=yes $<

build/test/%: %.c $(SOURCES) $(HEADERS) | build
	@echo CC $<
	@$(CC) $(CFLAGS) $(INCS) $(LLIBS) $< $(OBJS) -o $@ $(lLIBS) -Wl,-rpath $(ARB_LIBDIR)

build/time/%: %.c $(SOURCES) $(HEADERS) | build
	@echo CC $<
	@$(CC) -DTIMEDIR=\"$(TIMEDIR)\" $(CFLAGS) $(INCS) $(LLIBS) $< $(OBJS) -o $@ $(lLIBS) -Wl,-rpath $(ARB_LIBDIR)

print-%:
	@echo '$*=$($*)'

plots:
	python $(TIMEDIR)/make_plots.py

gp2c:
	sage $(CURDIR)/gp2c/gp_poly_to_acb.sage

.PHONY: all check clean tests library examples build valgrind gp2c

.PRECIOUS: Makefile
