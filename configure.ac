# Process this file with autoconf to produce a configure script.

AC_INIT(stamp, 4.6.01, stamp@russelllab.org, stamp, https://code.google.com/p/bio-stamp)
AC_PREREQ(2.59)
AC_CANONICAL_SYSTEM

# Needs to include the package name? 
AC_DEFINE_UNQUOTED([STAMPDIR], [$datadir]/stamp,
                   [Architecture-independent data directory.])

# 'foreign' doesn't require GNU template files like AUTHORS, etc
AM_INIT_AUTOMAKE(foreign)
AC_PROG_CC
# Create a header with config options
AC_CONFIG_HEADERS(config.h)
AC_CONFIG_FILES(Makefile src/Makefile)
AC_OUTPUT
