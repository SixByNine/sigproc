#
# SWIN_LIB_CPGPLOT([ACTION-IF-FOUND [,ACTION-IF-NOT-FOUND]])
#
# This m4 macro checks availability of the PGPLOT Graphics Subroutine 
# Library by T. J. Pearson, setting the following variables:
#
# PGPLOT_CFLAGS - autoconfig variable with flags required for compiling
# PGPLOT_LIBS   - autoconfig variable with flags required for linking
# HAVE_PGPLOT   - automake conditional
# HAVE_PGPLOT   - pre-processor macro in config.h
#
# This macro tries to link a test program, first using only 
#
#    -L$PGPLOT_DIR -lcpgplot -lpgplot $FLIBS
#
# and, if this fails, by adding 
#
#    $X_LIBS -lX11
#
# Notice that the environment variable PGPLOT_DIR is used.
#
# Extra libraries may be required if certain PGPLOT drivers are enabled.
# These may be specified using the
#
#    --with-pgplot-extra
#
# argument to the configure script.
#
# ----------------------------------------------------------
AC_DEFUN([SWIN_LIB_PGPLOT],
[
  AC_PROVIDE([SWIN_LIB_PGPLOT])
  AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
  AC_REQUIRE([SWIN_LIB_X11])

  AC_ARG_WITH([pgplot-extra],
              AC_HELP_STRING([--with-pgplot-extra=LIBS],
                             [Specify LIBS required by PGPLOT drivers]))

  AC_MSG_CHECKING([for PGPLOT installation])

  PGPLOT_CFLAGS=""
  PGPLOT_LIBS=""

  if test x"$PGPLOT_DIR" != x; then
    PGPLOT_CFLAGS="-I$PGPLOT_DIR"
    PGPLOT_LIBS="-L$PGPLOT_DIR"
  fi

  PGPLOT_LIBS="$PGPLOT_LIBS -lcpgplot -lpgplot $FLIBS"

  # "yes" is not a specification
  if test x"$with_pgplot_extra" != xyes; then
    PGPLOT_LIBS="$PGPLOT_LIBS $with_pgplot_extra"
  fi

  ac_save_CFLAGS="$CFLAGS"
  ac_save_LIBS="$LIBS"
  LIBS="$ac_save_LIBS $PGPLOT_LIBS"
  CFLAGS="$ac_save_CFLAGS $PGPLOT_CFLAGS"

  AC_TRY_LINK([#include <cpgplot.h>],[cpgopen(""); cpgend();],
              have_pgplot=yes, have_pgplot=no)

  if test $have_pgplot = no; then
    PGPLOT_LIBS="$PGPLOT_LIBS $X_LIBS -lX11"
    LIBS="$ac_save_LIBS $PGPLOT_LIBS"
    AC_TRY_LINK([#include <cpgplot.h>],[cpgopen(""); cpgend();],
                have_pgplot=yes, have_pgplot=no)
  fi

  if test $have_pgplot = no; then
    PGPLOT_LIBS="$PGPLOT_LIBS $X_LIBS -lX11 -lpng"
    LIBS="$ac_save_LIBS $PGPLOT_LIBS"
    AC_TRY_LINK([#include <cpgplot.h>],[cpgopen(""); cpgend();],
                have_pgplot=yes, have_pgplot=no)
  fi

  if test $have_pgplot = no; then
# Blade libs
    PGPLOT_LIBS="$PGPLOT_LIBS -L/usr/X11R6/lib -lX11 -lcpgplot -lpgplot -lpng"
    LIBS="$ac_save_LIBS $PGPLOT_LIBS"
    AC_TRY_LINK([#include <cpgplot.h>],[cpgopen(""); cpgend();],
                have_pgplot=yes, have_pgplot=no)
  fi

  AC_MSG_RESULT($have_pgplot)

  LIBS="$ac_save_LIBS"
  CFLAGS="$ac_save_CFLAGS"

  if test x"$have_pgplot" = xyes; then
    AC_DEFINE([HAVE_PGPLOT], [1], [Define to 1 if you have the PGPLOT library])
    [$1]
  else
    AC_MSG_WARN([PGPLOT code will not be compiled])
    if test x"$PGPLOT_DIR" = x; then
      AC_MSG_WARN([Please set the PGPLOT_DIR environment variable])
    fi
    PGPLOT_CFLAGS=""
    PGPLOT_LIBS=""
    [$2]
  fi

  AC_SUBST(PGPLOT_CFLAGS)
  AC_SUBST(PGPLOT_LIBS)
  AM_CONDITIONAL(HAVE_PGPLOT, [test x"$have_pgplot" = xyes])

])

