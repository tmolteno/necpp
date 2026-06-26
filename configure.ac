AC_INIT([necpp],[2.0.0])
AC_CONFIG_AUX_DIR([m4])
AC_CONFIG_MACRO_DIR([config/m4])
AC_CONFIG_HEADERS(config.h)
AM_INIT_AUTOMAKE([foreign subdir-objects])
LT_INIT
AM_SILENT_RULES([yes])
LT_INIT

CXXFLAGS="-O2 -Wall -Wextra -Wshadow -Werror"

dnl ************************************************
dnl *** Check for Array Bounds
dnl ************************************************

AC_ARG_WITH([bounds],
  [AS_HELP_STRING([--with-bounds],
  [Enable support for bounds-checking])],
  [],
  [with_bounds=no])

AS_IF([test "x$with_bounds" = xyes],
[
  AC_SUBST([AM_CPPFLAGS], ["-DNEC_ERROR_CHECK=1"])
  AC_SUBST([CXXFLAGS], ["-g3 -O0 -DNEC_ERROR_CHECK=1"])
])

dnl ************************************************
dnl *** Trigger typesafe integer checking
dnl ************************************************

AC_ARG_WITH([typecheck],
  [AS_HELP_STRING([--with-typecheck],
  [Enable support for typesafe integer checking])],
  [],
  [with_typecheck=no])

AS_IF([test "x$with_typecheck" = xyes],
[
  AC_SUBST([AM_CPPFLAGS], ["-DTYPESAFE_PEDANTIC=1"])
  AC_SUBST([CXXFLAGS], ["-g3 -O0 -DTYPESAFE_PEDANTIC=1"])
])

dnl ************************************************
dnl *** Add an automatic build date
dnl ************************************************
AC_PATH_PROG(DATE,date)
if test -z "$ac_cv_path_DATE" ; then
    AC_MSG_ERROR(** "date" not found. Can't determine build-date.)
    AC_SUBST(BUILD_DATE, unknown)
    AC_DEFINE(BUILD_DATE, "unknown", [The build date])
else
    AC_SUBST(BUILD_DATE, `date +"%Y-%m-%d"`)
    AC_MSG_NOTICE(Build date: $BUILD_DATE)
    AC_DEFINE_UNQUOTED(BUILD_DATE, "$BUILD_DATE", [The build date])
fi

dnl ************************************************
dnl *** Check for Eigen (mandatory)
dnl ************************************************

AC_LANG([C++])
AC_PROG_CXX
LT_INIT

CPPFLAGS="$CPPFLAGS -I/usr/include/eigen3"
AC_CHECK_HEADER(
  [Eigen/Dense],
  [],
  [AC_MSG_FAILURE([Eigen library not found. Install libeigen3-dev.])],
  [])

dnl Eigen is now the default for all components.
CPPFLAGS="$CPPFLAGS -DUSING_EIGEN_ARRAY=1 -DUSING_EIGEN_3VECT=1 -DEIGEN_LU=1"

AC_CHECK_LIB([m],[pow])

EXPLICIT_LIBS=""
PRIVATE_LIBS=" $LDFLAGS -lstdc++ -lpthread -lm"

AC_SUBST(EXPLICIT_LIBS)
AC_SUBST(PRIVATE_LIBS)

AC_CONFIG_FILES([Makefile src/Makefile necpp.pc])
AC_OUTPUT
