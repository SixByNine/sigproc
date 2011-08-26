dnl @synopsis SWIN_OPTIONS_SET
dnl 
AC_DEFUN([SWIN_OPTIONS_SET],
[
  AC_PROVIDE([SWIN_OPTIONS_SET])

  if test "$CXXFLAGS"; then
    swin_cxxflags_set="yes"
  fi

  if test "$CFLAGS"; then
    swin_cflags_set="yes"
  fi

  if test "$FFLAGS"; then
    swin_fflags_set="yes"
  fi

])

dnl @synopsis SWIN_DEBUG
dnl 
AC_DEFUN([SWIN_DEBUG],
[
  AC_PROVIDE([SWIN_DEBUG])

  AC_ARG_ENABLE([debug],
                AC_HELP_STRING([--enable-debug],
                               [Enable debugging information]),
                [swin_debug=yes])

  if test x"$swin_debug" != xyes; then

    AC_MSG_NOTICE([Disabling debugging information compiler option (-g)])

    #
    # 'man grep' states that \> matches the empty string at the end of a word
    #
    # however, sed on OS X does not support this (not even sed -E)
    #
    # therefore, strip "-g " and "-g$" ($ = empty string at the end of a line)
    #
    echo 's/-g //g' > .ac_strip_g
    echo 's/-g$//g' >> .ac_strip_g

    if test x"$swin_cxxflags_set" != xyes; then
      CXXFLAGS="`echo $CXXFLAGS | sed -f .ac_strip_g` -Wall"
    else
      AC_MSG_NOTICE([   CXXFLAGS set by user.])
    fi

    if test x"$swin_cflags_set" != xyes; then
      CFLAGS="`echo $CFLAGS | sed -f .ac_strip_g` -Wall"
    else
      AC_MSG_NOTICE([   CFLAGS set by user.])
    fi

    if test x"$swin_fflags_set" != xyes; then
      FFLAGS=`echo $FFLAGS | sed -f .ac_strip_g`
    else
      AC_MSG_NOTICE([   FFLAGS set by user.])
    fi

    rm -f .ac_strip_g

  else
    AC_MSG_WARN([Debugging information not disabled.  Binaries may be large.])
  fi

])

