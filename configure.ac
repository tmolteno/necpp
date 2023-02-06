AC_INIT([necpp],[1.7.4])
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
dnl *** Check for Lapack library
dnl ************************************************
AC_PROG_F77

AC_ARG_WITH([lapack],
  [AS_HELP_STRING([--without-lapack],
  [disable support for lapack])],
  [],
  [with_lapack=yes])

LIBLAPACK=
  AS_IF([test "x$with_lapack" != xno],
  [
    AX_BLAS([], [AC_MSG_ERROR([BLAS library not found])])
    AX_LAPACK([], [AC_MSG_ERROR([LAPACK library not found])])

    LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"
    AC_CHECK_LIB(
    [lapack], [clapack_zgetrf], 
    [AC_SUBST([LIBLAPACK], ["$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"]) AC_DEFINE([LAPACK], [1], [Define if you have liblapack])],
    [AC_MSG_FAILURE([lapack library test failed (--without-lapack to disable)])],
    [])
  ])

dnl ************************************************
dnl *** Check for Eigen library
dnl ************************************************

AC_LANG([C++])
AC_PROG_CXX
LT_INIT


AC_ARG_WITH([eigen],
  [AS_HELP_STRING([--without-eigen],
  [disable support for eigen])],
  [],
  [with_eigen=no])

AC_ARG_WITH([eigenv],
  [AS_HELP_STRING([--without-eigenv],
  [disable support for eigen])],
  [],
  [with_eigenv=no])

LIBEIGEN=
  AS_IF([test "x$with_eigen" != xno],
  [
    CPPFLAGS="$CPPFLAGS -I/usr/include/eigen3"
    AC_CHECK_HEADER(
    [Eigen/Core],
    [AC_DEFINE([USING_EIGEN_ARRAY], [1], [Define if you have eigen])],
    [AC_MSG_FAILURE([Eigen library test failed (--without-eigen to disable)])],
    [])
  ])

LIBEIGENV=
  AS_IF([test "x$with_eigenv" != xno],
  [
    CPPFLAGS="$CPPFLAGS -I/usr/include/eigen3"
    AC_CHECK_HEADER(
    [Eigen/Core],
    [AC_DEFINE([USING_EIGEN_3VECT], [1], [Define if you have eigen])],
    [AC_MSG_FAILURE([Eigen library test failed (--without-eigen to disable)])],
    [])
  ])


AC_CHECK_LIB([m],[pow])

EXPLICIT_LIBS=""
PRIVATE_LIBS=" $LDFLAGS -lstdc++ -llapack -lblas -lpthread -lm"

AC_SUBST(EXPLICIT_LIBS)
AC_SUBST(PRIVATE_LIBS)

AC_CONFIG_FILES([Makefile src/Makefile necpp.pc])
AC_OUTPUT

