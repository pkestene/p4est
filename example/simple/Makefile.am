
# This file is part of p4est.
# Makefile.am in example/simple
# included non-recursively from toplevel directory

if P4EST_ENABLE_BUILD_2D
bin_PROGRAMS += example/simple/p4est_simple
example_simple_p4est_simple_SOURCES = example/simple/simple2.c

bin_PROGRAMS += example/simple/p4est_simple_iterate
example_simple_p4est_simple_iterate_SOURCES = example/simple/simple2_iterate.c example/simple/ring_connectivity.c

LINT_CSOURCES += \
	$(example_simple_p4est_simple_SOURCES) \
	$(example_simple_p4est_simple_iterate_SOURCES)
endif

if P4EST_ENABLE_BUILD_3D
bin_PROGRAMS += example/simple/p8est_simple
example_simple_p8est_simple_SOURCES = example/simple/simple3.c

LINT_CSOURCES += $(example_simple_p8est_simple_SOURCES)
endif
