#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.57)
AC_INIT([sixproc],[5.0.4],[mkeith@pulsarastronomy.net])

if test -n "${CXXFLAGS}"; then
	user_set_cxxflags=yes
fi

SWIN_OPTIONS_SET

AC_CONFIG_AUX_DIR([config])
#AC_CONFIG_SRCDIR([Makefile.am])
# Create static libraries by default
AC_DISABLE_SHARED

# Enable large-file support
AC_SYS_LARGEFILE

# Determine the machine endian
AC_C_BIGENDIAN

# Checks for programs.
AC_PROG_CXX
AC_PROG_CXXCPP
AC_PROG_CC
AC_PROG_CPP
AM_PROG_CC_C_O
AC_PROG_F77
AC_PROG_F77_C_O
AC_F77_LIBRARY_LDFLAGS

if test -n "${CXXFLAGS}"; then
	user_set_cxxflags=yes
fi
AC_PROG_CXX
if test X$user_set_cxxflags != Xyes; then
	CXXFLAGS="-O2"
fi


AC_PROG_INSTALL
AC_PROG_LIBTOOL

AM_PATH_XML2

SWIN_LIB_CFITSIO
SWIN_LIB_PGPLOT
MJK_LIB_PSRXML

AX_OPENMP([have_openmp=yes],[have_openmp=no])
if test X$have_openmp == Xyes; then
	CFLAGS="$CFLAGS $OPENMP_CFLAGS -DHAVE_OPENMP"
	CXXFLAGS="$CXXFLAGS $OPENMP_CFLAGS -DHAVE_OPENMP"
fi


CHECK_ZLIB
MJK_LIB_FFTW

if test X$have_fftw != Xyes; then
	AC_WARN("No FFTW3 found... you will have to use the slower 'singleton' FFT. Specify location with --with-fftw-dir")
fi
if test X$have_openmp == Xno; then
	AC_MSG_WARN([**OpenMP not detected, will compile only single threaded**])
fi



AC_SUBST(XML_CPPFLAGS)
AC_SUBST(PSRXML_LIBS)
AC_SUBST(PSRXML_CFLAGS)

AM_INIT_AUTOMAKE([foreign subdir-objects])



AM_CONFIG_HEADER([config.h])
AC_CONFIG_FILES([Makefile src/Makefile])

AC_OUTPUT
