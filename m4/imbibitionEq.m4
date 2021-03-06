dnl -*- autoconf -*-
# Macros needed to find imbibitionEq and dependent libraries.  They are called by
# the macros in ${top_src_dir}/dependencies.m4, which is generated by
# "dunecontrol autogen"

# Additional checks needed to build imbibitionEq
# This macro should be invoked by every module which depends on imbibitionEq, as
# well as by imbibitionEq itself
AC_DEFUN([IMBIBITIONEQ_CHECKS])

# Additional checks needed to find imbibitionEq
# This macro should be invoked by every module which depends on imbibitionEq, but
# not by imbibitionEq itself
AC_DEFUN([IMBIBITIONEQ_CHECK_MODULE],
[
  DUNE_CHECK_MODULES([imbibitionEq],[imbibitionEq/imbibitionEq.hh])
])
