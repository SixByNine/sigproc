dnl @synopsis SWIN_LIB_CFITSIO
dnl 
AC_DEFUN([SWIN_LIB_CFITSIO],
[
  AC_PROVIDE([SWIN_LIB_CFITSIO])
  AC_REQUIRE([ETR_SOCKET_NSL])

  SWIN_PACKAGE_OPTIONS([cfitsio])

  AC_MSG_CHECKING([for CFITSIO installation])

  if test "$have_cfitsio" != "user disabled"; then

    SWIN_PACKAGE_FIND([cfitsio],[fitsio.h])
    SWIN_PACKAGE_TRY_COMPILE([cfitsio],[#include <fitsio.h>])

    SWIN_PACKAGE_FIND([cfitsio],[libcfitsio.*])
    SWIN_PACKAGE_TRY_LINK([cfitsio],[#include <fitsio.h>],
                          [fits_movnam_hdu(0,0,0,0,0);],
                          [-lcfitsio -lm $SOCKET_LIBS])

  fi

  AC_MSG_RESULT([$have_cfitsio])

  if test $have_cfitsio = yes; then
    AC_DEFINE([HAVE_CFITSIO], [1], [Define if the CFITSIO library is present])
    [$1]
  else
    AC_MSG_WARN([The PSRFITS code will not be compiled])
    [$2]
  fi

  CFITSIO_LIBS="$cfitsio_LIBS"
  CFITSIO_CFLAGS="$cfitsio_CFLAGS"

  AC_SUBST(CFITSIO_LIBS)
  AC_SUBST(CFITSIO_CFLAGS)
  AM_CONDITIONAL(HAVE_CFITSIO,[test $have_cfitsio = yes])

])

