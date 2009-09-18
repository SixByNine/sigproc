#
# SWIN_LIB_CPGPLOT([ACTION-IF-FOUND [,ACTION-IF-NOT-FOUND]])
#
# This m4 macro checks availability of the PSRXML io library
# created by M.Keith.
#
# PSRXML_CFLAGS - autoconfig variable with flags required for compiling
# PSRXML_LIBS   - autoconfig variable with flags required for linking
# HAVE_PSRXML   - automake conditional
#
# ----------------------------------------------------------
AC_DEFUN([MJK_LIB_PSRXML],
[
  AC_PROVIDE([MJK_LIB_PSRXML])
  AC_REQUIRE([MJK_LIB_CRYPTO])
  AC_REQUIRE([AM_PATH_XML2])

  AC_ARG_WITH([psrxml-dir],
	  AC_HELP_STRING([--with-psrxml-dir=DIR],
		  [Specify dir where the PsrXML io library is installed]))



  AC_MSG_CHECKING([for PsrXML libraries])


  PSRXML_CFLAGS="$XML_CPP_FLAGS"
  PSRXML_LIBS="$XML_LIBS -lpsrxml $OPENSSL_LIBS"

  if test x"$with_psrxml_dir" != xyes; then
  	PSRXML_LIBS="-L$with_psrxml_dir/lib $PSRXML_LIBS"
	PSRXML_CFLAGS="$PSRXML_CFLAGS -I$with_psrxml_dir/include"
  fi


  ac_save_LIBS=$LIBS
  ac_save_CFLAGS=$CFLAGS
  LIBS="$LIBS $PSRXML_LIBS"
  CFLAGS="$CFLAGS $PSRXML_CFLAGS"
  AC_TRY_LINK([#include <psrxml.h>],[writePsrXml((psrxml*)NULL, "");],
              have_psrxml=yes, have_psrxml=no)

  AC_MSG_RESULT($have_psrxml)

  LIBS="$ac_save_LIBS"
  CFLAGS="$ac_save_CFLAGS"

  if test x"$have_psrxml" = xyes; then
    AC_DEFINE([HAVE_PGPLOT], [1], [Define to 1 if you have the PGPLOT library])
    [$1]
  else
    AC_MSG_WARN([PsrXML code will not be compiled])
    PGPLOT_CFLAGS=""
    PGPLOT_LIBS=""
    [$2]
  fi

  AC_SUBST(PSRXML_CFLAGS)
  AC_SUBST(PSRXML_LIBS)
  AM_CONDITIONAL(HAVE_PSRXML, [test x"$have_psrxml" = xyes])

])

