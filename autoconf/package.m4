dnl @synopsis SWIN_PACKAGE_OPTIONS(
dnl 
AC_DEFUN([SWIN_PACKAGE_OPTIONS],
[
  AC_PROVIDE([SWIN_PACKAGE_OPTIONS])

  AC_ARG_WITH([[$1]-dir],
              AC_HELP_STRING([--with-[$1]-dir=DIR],
                             [[$1] is in DIR]))

  AC_ARG_WITH([[$1]-include-dir],
              AC_HELP_STRING([--with-[$1]-include-dir=DIR],
                             [[$1] header files are in DIR]))

  AC_ARG_WITH([[$1]-lib-dir],
              AC_HELP_STRING([--with-[$1]-lib-dir=DIR],
                             [[$1] library is in DIR]))

  if test x"$with_[$1]_dir" = x"no" ||
     test x"$with_[$1]_include-dir" = x"no" ||
     test x"$with_[$1]_lib_dir" = x"no"; then

    # user disabled [$1]. Leave cache alone.
    have_[$1]="user disabled"

  else

    # "yes" is not a specification
    if test x"$with_[$1]_dir" = xyes; then
      with_[$1]_dir=""
    fi

    if test x"$with_[$1]_include_dir" = xyes; then
      with_[$1]_include_dir=""
    fi

    if test x"$with_[$1]_lib_dir" = xyes; then
      with_[$1]_lib_dir=""
    fi

    # _dir overrides include_dir and lib_dir    
    if test x"$with_[$1]_dir" != x; then
      with_[$1]_include_dir=$with_[$1]_dir/include
      with_[$1]_lib_dir=$with_[$1]_dir/lib
    fi

  fi

])

dnl @synopsis SWIN_PACKAGE_LIB_OPTIONS(
dnl 
AC_DEFUN([SWIN_PACKAGE_LIB_OPTIONS],
[
  AC_PROVIDE([SWIN_PACKAGE_LIB_OPTIONS])

  AC_ARG_WITH([[$1]-dir],
              AC_HELP_STRING([--with-[$1]-dir=DIR],
                             [[$1] is in DIR]))

  if test x"$with_[$1]_dir" = x"no"; then

    # user disabled [$1]. Leave cache alone.
    have_[$1]="user disabled"

  else

    # "yes" is not a specification
    if test x"$with_[$1]_dir" = xyes; then
      with_[$1]_dir=""
    fi

    # _dir overrides include_dir and lib_dir    
    if test x"$with_[$1]_dir" != x; then
      with_[$1]_lib_dir=$with_[$1]_dir
    fi

  fi

])

dnl @synopsis SWIN_PACKAGE_FIND(name,file)
dnl 
AC_DEFUN([SWIN_PACKAGE_FIND],
[
  AC_PROVIDE([SWIN_PACKAGE_FIND])

  swin_[$1]_found=""

  if test x"$PACKAGES" != x; then
    for cf_file in `find $PACKAGES -name "[$2]" 2> /dev/null`; do
      cf_path=`dirname $cf_file`
      swin_[$1]_found="$swin_[$1]_found $cf_path"
    done
  fi
])

dnl @synopsis SWIN_PACKAGE_TRY_COMPILE(name,includes,function body)
dnl 
AC_DEFUN([SWIN_PACKAGE_TRY_COMPILE],
[
  AC_PROVIDE([SWIN_PACKAGE_TRY_COMPILE])

  cf_include_path_list="$with_[$1]_include_dir $swin_[$1]_found [$4] ."

  ac_save_CFLAGS="$CFLAGS"

  for cf_dir in $cf_include_path_list; do
    CFLAGS="-I$cf_dir $ac_save_CFLAGS"
    AC_TRY_COMPILE([$2], [$3], have_[$1]=yes, have_[$1]=no)
    if test $have_[$1] = yes; then
      if test x"$cf_dir" != x.; then
        [$1]_CFLAGS="-I$cf_dir"
        swin_[$1]_include_dir="$cf_dir"
      fi
      break
    fi
  done

  CFLAGS="$ac_save_CFLAGS"
])

dnl @synopsis SWIN_PACKAGE_TRY_LINK(name,includes,function body,ldflags)
dnl 
AC_DEFUN([SWIN_PACKAGE_TRY_LINK],
[
  AC_PROVIDE([SWIN_PACKAGE_TRY_LINK])

  swin_search_path="$swin_[$1]_found"
  if test x"$swin_[$1]_include_dir" != x; then
    swin_base=`dirname $swin_[$1]_include_dir`
    if test -d $swin_base/lib; then
      swin_search_path="$swin_base/lib $swin_search_path"
    fi
  fi

  cf_lib_path_list="$with_[$1]_lib_dir $swin_search_path [$5] ."

  ac_save_CFLAGS="$CFLAGS"
  CFLAGS="$[$1]_CFLAGS $ac_save_CFLAGS"

  ac_save_LIBS="$LIBS"

  for cf_dir in $cf_lib_path_list; do
    LIBS="-L$cf_dir [$4] $ac_save_LIBS"
    AC_TRY_LINK([$2], [$3], have_[$1]=yes, have_[$1]=no)
    if test $have_[$1] = yes; then
      if test x"$cf_dir" = x.; then
        [$1]_LIBS="[$4]"
      else
        [$1]_LIBS="-L$cf_dir [$4]"
      fi
      break
    fi
  done

  CFLAGS="$ac_save_CFLAGS"
  LIBS="$ac_save_LIBS"

])

