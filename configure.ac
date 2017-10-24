# prelude
AC_INIT([sid], [0.0], [lennart.weckeck@mailbox.org])
AM_INIT_AUTOMAKE([foreign subdir-objects -Wall -Werror])

# make sure we are in source directory
AC_CONFIG_SRCDIR(sid.cpp)

# check for c++ compiler
AC_PROG_CXX([g++ clang++])
AC_PROG_CXX
AC_OPENMP

# recommended check for GNU scientific library, from https://www.gnu.org/software/gsl/manual/html_node/Autoconf-Macros.html (6.6.2017)
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

AC_CONFIG_FILES([Makefile test/Makefile])
AC_OUTPUT
