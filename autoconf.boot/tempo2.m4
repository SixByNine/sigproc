#
# SWIN_LIB_TEMPO2([ACTION-IF-FOUND [,ACTION-IF-NOT-FOUND]])
#
# This m4 macro checks availability of the TEMPO2 Predictor library by
# Russell Edwards, George Hobbs, and Willem van Straten
#
# TEMPO2_CFLAGS - autoconfig variable with flags required for compiling
# TEMPO2_LIBS   - autoconfig variable with flags required for linking
# HAVE_TEMPO2   - automake conditional
# HAVE_TEMPO2   - pre-processor macro in config.h
#
# This macro tries to link a test program, using 
#
#    -I$TEMPO2/include -L$TEMPO2/lib -ltempo2pred
#
# Notice that the environment variable TEMPO2 is required.
#
# ----------------------------------------------------------
AC_DEFUN([SWIN_LIB_TEMPO2],
[
  AC_PROVIDE([SWIN_LIB_TEMPO2])

  AC_MSG_CHECKING([for TEMPO2 Predictor library])

  TEMPO2_CFLAGS="-I$TEMPO2/include"
  TEMPO2_LIBS="-L$TEMPO2/lib -ltempo2pred -lm"

  ac_save_CXXFLAGS="$CXXFLAGS"
  ac_save_LIBS="$LIBS"

  AC_LANG_PUSH(C++)

  LIBS="$ac_save_LIBS $TEMPO2_LIBS"
  CXXFLAGS="$ac_save_CXXFLAGS $TEMPO2_CFLAGS"

  AC_TRY_LINK([#include "tempo2pred.h"],
              [T2Predictor p; T2Predictor_Read (&p, 0);],
              have_tempo2=yes, have_tempo2=no)

  AC_MSG_RESULT($have_tempo2)

  LIBS="$ac_save_LIBS"
  CXXFLAGS="$ac_save_CXXFLAGS"

  AC_LANG_POP(C++)

  if test x"$have_tempo2" = xyes; then
    AC_DEFINE([HAVE_TEMPO2], [1], [Define to 1 if you have the TEMPO2 library])
    [$1]
  else
    AC_MSG_WARN([TEMPO2 code will not be compiled])
    if test x"$TEMPO2" = x; then
      AC_MSG_WARN([Please set the TEMPO2 environment variable])
    fi
    TEMPO2_CFLAGS=""
    TEMPO2_LIBS=""
    [$2]
  fi

  AC_SUBST(TEMPO2_CFLAGS)
  AC_SUBST(TEMPO2_LIBS)
  AM_CONDITIONAL(HAVE_TEMPO2, [test x"$have_tempo2" = xyes])

])

