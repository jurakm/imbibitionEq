# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.62])
DUNE_AC_INIT # gets module version from dune.module file
AC_CONFIG_SRCDIR([src/imbibitionEq.cc])
AC_CONFIG_HEADERS([config.h])

# we need no more than the standard DE-stuff
# this module depends on dune-common dune-istl dune-geometry dune-grid dune-localfunctions dune-typetree dune-pdelab
# this implies checking for [dune-common], [dune-istl], [dune-geometry], [dune-grid], [dune-localfunctions], [dune-typetree], [dune-pdelab]
DUNE_CHECK_ALL

# implicitly set the Dune-flags everywhere
AC_SUBST(AM_CPPFLAGS, $DUNE_CPPFLAGS)
AC_SUBST(AM_LDFLAGS, $DUNE_LDFLAGS)
LIBS="$DUNE_LIBS /usr/local/lib/libo2scl.a -lgsl -lblas"

AC_CONFIG_FILES([
  Makefile
  src/Makefile
  cmake/Makefile
  cmake/modules/Makefile
  doc/Makefile
  doc/doxygen/Makefile
  doc/doxygen/Doxyfile
  dune/Makefile
  dune/imbibitionEq/Makefile
  m4/Makefile
  imbibitionEq.pc
])
AC_OUTPUT
# finally print the summary information
DUNE_SUMMARY_ALL
