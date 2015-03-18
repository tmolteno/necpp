AC_INIT([necpp], [1.6.1])
AC_CONFIG_MACRO_DIR([m4])
AM_CONFIG_HEADER(config.h)
AM_INIT_AUTOMAKE

LT_INIT

CXXFLAGS="-O2 -Wall -Wextra"

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
dnl *** Add an automatic build date
dnl ************************************************
AC_PATH_PROG(DATE,date)
if test -z "$ac_cv_path_DATE" ; then
    AC_MSG_ERROR(** \"date\" not found. Can't determine build-date.)
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


AC_ARG_WITH([lapack],
	[AS_HELP_STRING([--without-lapack],
	[disable support for lapack])],
	[],
	[with_lapack=yes])

LIBLAPACK=
	AS_IF([test "x$with_lapack" != xno],
	[
		LDFLAGS="$LDFLAGS -L/usr/lib/atlas-base/atlas"
		AC_CHECK_LIB(
		[lapack], [clapack_zgetrf], 
		[AC_SUBST([LIBLAPACK], ["-L/usr/lib/atlas-base/atlas -llapack -lblas -lpthread"]) AC_DEFINE([LAPACK], [1], [Define if you have liblapack])],
		[AC_MSG_FAILURE([lapack library test failed (--without-lapack to disable)])],
		[])
	])

dnl ************************************************
dnl *** Check for Eigen library
dnl ************************************************

AC_LANG_CPLUSPLUS
AC_PROG_CXX
AM_PROG_LIBTOOL


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

AC_OUTPUT(Makefile src/Makefile necpp.pc)

