#
# This m4 macro checks availability of the PSRXML io library
# created by M.Keith.
#
# PSRXML_CFLAGS - autoconfig variable with flags required for compiling
# PSRXML_LIBS   - autoconfig variable with flags required for linking
# HAVE_PSRXML   - automake conditional
#
# ----------------------------------------------------------
AC_DEFUN([MJK_LIB_FFTW],
[
  AC_PROVIDE([MJK_LIB_FFTW])
  AC_REQUIRE([ACX_PTHREAD])

  AC_ARG_WITH([fftw-dir],
	  AC_HELP_STRING([--with-fftw-dir=DIR],
		  [Specify dir where the FFTW3 single precision library is installed]))



  AC_MSG_CHECKING([for FFTW3 single precision libraries])


  FFTW_CFLAGS=""
  FFTW_LIBS="-lfftw3f -lm"

  if test x"$with_fftw_dir" != xyes; then
  	FFTW_LIBS="-L$with_fftw_dir/lib $FFTW_LIBS"
	FFTW_CFLAGS="$FFTW_CFLAGS -I$with_fftw_dir/include"
  fi


  ac_save_LIBS=$LIBS
  ac_save_CFLAGS=$CFLAGS
  LIBS="$LIBS $FFTW_LIBS"
  CFLAGS="$CFLAGS $FFTW_CFLAGS"
  AC_TRY_LINK([#include <fftw3.h>],[fftwf_plan_dft_1d(0,0,0,FFTW_FORWARD,FFTW_ESTIMATE);],
              have_fftw=yes, have_fftw=no)

  AC_MSG_RESULT($have_fftw)

  LIBS="$ac_save_LIBS"
  CFLAGS="$ac_save_CFLAGS"

  if test x"$have_fftw" = xyes; then
    AC_DEFINE([HAVE_FFTW], [1], [Define to 1 if you have the PGPLOT library])
    [$1]
  else
    AC_MSG_WARN([FFTW code will not be compiled])
    FFTW_CFLAGS=""
    FFTW_LIBS=""
    [$2]
  fi



  FFTW_THREAD_LIBS="-lfftw3f_threads $FFTW_LIBS $PTHREAD_LIBS"
  FFTW_THREAD_CFLAGS="$FFTW_CFLAGS $PTHREAD_CFLAGS"

  echo $FFTW_THREAD_LIBS
    echo $FFTW_THREAD_CFLAGS


  LIBS="$LIBS $FFTW_THREAD_LIBS"
  CFLAGS="$CFLAGS $FFTW_THREAD_CFLAGS"

  AC_MSG_CHECKING([for FFTW3 multi-threaded library])


  AC_TRY_LINK([#include <fftw3.h>],[fftwf_init_threads();],
              have_fftw_threads=yes, have_fftw_threads=no)

  AC_MSG_RESULT($have_fftw_threads)
  LIBS="$ac_save_LIBS"
  CFLAGS="$ac_save_CFLAGS"

  if test x"$have_fftw_threads" = xyes; then
    AC_DEFINE([HAVE_FFTW_THREADS], [1], [Define to 1 if you have the FFTW_THREADS library])
    [$1]
    FFTW_CFLAGS="$FFTW_THREAD_CFLAGS -DFFTW_THREADS"
    FFTW_LIBS="$FFTW_THREAD_LIBS"
  else
    AC_MSG_WARN([FFTW multi-threaded code will not be compiled])
    [$2]
  fi


  AC_SUBST(FFTW_CFLAGS)
  AC_SUBST(FFTW_LIBS)
  AM_CONDITIONAL(HAVE_FFTW, [test x"$have_fftw" = xyes])

])

