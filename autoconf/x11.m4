# SWIN_LIB_X11([search paths, [ACTION-IF-FOUND [,ACTION-IF-NOT-FOUND]]])
# ----------------------------------------------------------
AC_DEFUN([SWIN_LIB_X11],
[
  AC_PROVIDE([SWIN_LIB_X11])
  AC_REQUIRE([AC_PATH_XTRA])

  AC_MSG_CHECKING([for libX11.so])

  swin_x11_found="no"

  if test "x$X_LIBS" != "x"; then
    swin_x11_found="yes"
  else
    swin_lib_x11_path_list="$1 /usr/X11R6/lib"
    for swin_dir in $swin_lib_x11_path_list; do
      if test -r "$swin_dir/libX11.so"; then
        X_LIBS="-L$swin_dir"
        swin_x11_found="yes"
      fi
    done
  fi

  AC_MSG_RESULT($swin_x11_found)
  AC_MSG_RESULT([X_LIBS=$X_LIBS])
  AC_MSG_RESULT([X_PRE_LIBS=$X_PRE_LIBS])
  AC_MSG_RESULT([X_EXTRA_LIBS=$X_EXTRA_LIBS])
  AC_MSG_RESULT([X_CFLAGS=$X_CFLAGS])


])
