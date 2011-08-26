AC_DEFUN([MJK_LIB_CRYPTO],
[
	AC_PROVIDE([MJK_LIB_CRYPTO])
	AC_CHECK_HEADERS([openssl/sha.h])
	AC_CHECK_LIB([crypto],[SHA1],have_openssl=yes, have_openssl=no)
	AM_CONDITIONAL(HAVE_OPENSSL, [test x"$have_openssl" = xyes])
	if test x"$have_openssl" = xyes; then
		OPENSSL_LIBS="-lcrypto"
		AC_SUBST(OPENSSL_LIBS)
		echo "Using openssl for SHA checksums."
	fi
])
