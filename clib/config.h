/* fftw/config.h.  Generated automatically by configure.  */
/* -*- C -*- */
/*
 * Copyright (c) 1997-1999 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/* fftw.h -- system-wide definitions */
/* $Id: config.h,v 1.1.1.1 2003-01-19 21:58:43 giannozz Exp $ */

/* configuration options (guessed by configure) */

/* Define to empty if the keyword does not work.  */
/* #undef const */

/* Define if you have the gettimeofday function.  */
#define HAVE_GETTIMEOFDAY 1

/* Define if you have the BSDgettimeofday function.  */
/* #undef HAVE_BSDGETTIMEOFDAY */

/* Define if you have the <sys/time.h> header file.  */
#define HAVE_SYS_TIME_H 1

/* Define if you have the <unistd.h> header file.  */
#define HAVE_UNISTD_H 1

/* Define if you have the <getopt.h> header file.  */
#define HAVE_GETOPT_H 1

/* Define if you have the <malloc.h> header file */
#define HAVE_MALLOC_H 1

/* Define if you have gethrtime() a la Solaris 2 */
/* #undef HAVE_GETHRTIME */
/* #undef HAVE_HRTIME_T */

/* Define to sizeof int and long long, if available: */
#define SIZEOF_INT 4
#define SIZEOF_LONG_LONG 8

#if (SIZEOF_INT != 0) && (SIZEOF_LONG_LONG >= 2 * SIZEOF_INT)
#  define LONGLONG_IS_TWOINTS
#endif

/* Define to use "unsafe" modular multiply (can cause integer overflow
   and errors for transforms of large prime sizes using Rader). */
/* #undef FFTW_ENABLE_UNSAFE_MULMOD */

/* Define if you have getopt() */
#define HAVE_GETOPT 1

/* Define if you have getopt_long() */
#define HAVE_GETOPT_LONG 1

/* Define if you have isnan() */
#define HAVE_ISNAN 1

/* Define for enabling the high resolution Pentium timer */
/* #undef FFTW_ENABLE_PENTIUM_TIMER */

/*
 * When using FFTW_ENABLE_PENTIUM_TIMER, set FFTW_CYCLES_PER_SEC 
 * to your real CPU clock speed! 
 */
/* This is for 200 MHz */
/* #define FFTW_CYCLES_PER_SEC 200000000L */

/*
 * Define to enable a gcc/x86 specific hack that aligns
 * the stack to an 8-byte boundary 
 */
/* #undef FFTW_ENABLE_I386_HACKS */

/* Define when using a version of gcc that aligns the stack properly */
#define FFTW_GCC_ALIGNS_STACK 1

/* Define to enable extra runtime checks for debugging. */
/* #undef FFTW_DEBUG */

/* Define to enable vector-recurse feature. */
/* #undef FFTW_ENABLE_VECTOR_RECURSE */

/*
 * Define to enable extra runtime checks for the alignment of variables
 * in the codelets (causes coredump for misaligned double on x86). 
 */
/* #undef FFTW_DEBUG_ALIGNMENT */

#define FFTW_VERSION "2.1.3" 

/* Use Win32 high-resolution timer */
#if defined(__WIN32__) || defined(WIN32) || defined(_WINDOWS)
#  define HAVE_WIN32_TIMER
#  define HAVE_WIN32
#endif

/* Use MacOS Time Manager timer */
#if defined(MAC) || defined(macintosh)
#  define HAVE_MAC_TIMER
#  define HAVE_MACOS

/* Define to use nanosecond timer on PCI PowerMacs: */
/* (WARNING: experimental, use at your own risk.) */
/* #undef HAVE_MAC_PCI_TIMER */
#endif

/* define if you have alloca.h: */
#define HAVE_ALLOCA_H 1

/* define if you have the alloca function: */
#define HAVE_ALLOCA 1

/************************** threads configuration ************************/

/* The following preprocessor symbols select which threads library
   to use when compiling the FFTW threads parallel libraries: */

/* #undef FFTW_USING_SOLARIS_THREADS */
/* #undef FFTW_USING_POSIX_THREADS */
/* #undef FFTW_USING_BEOS_THREADS */
/* #undef FFTW_USING_MACH_THREADS */

/* #undef HAVE_PTHREAD_CREATE_UNDETACHED */

/* #undef HAVE_MACH_CTHREADS_H */
/* #undef HAVE_CTHREADS_H */
/* #undef HAVE_CTHREAD_H */

#ifdef HAVE_WIN32
#define FFTW_USING_WIN32_THREADS
#endif

#ifdef HAVE_MACOS
#define FFTW_USING_MACOS_THREADS
#endif

/*********************** fortran wrapper configuration *********************/

/* These symbols select how to mangle function names so that they will
   be recognized by the linker.  If none of them are defined, then
   Fortran wrappers will not be compiled. */

/* #undef FFTW_FORTRANIZE_LOWERCASE */
#define FFTW_FORTRANIZE_LOWERCASE_UNDERSCORE 1
/* #undef FFTW_FORTRANIZE_UPPERCASE */
/* #undef FFTW_FORTRANIZE_UPPERCASE_UNDERSCORE */

/* define the following if names with an underscore get an extra one: */
/* #undef FFTW_FORTRANIZE_EXTRA_UNDERSCORE */

/* The following symbols control how MPI_Comm data structures are
   translated between Fortran and C for the fftw_mpi wrappers.  See
   the file mpi/fftw_f77_mpi.h for more information. */
/* #undef HAVE_MPI_COMM_F2C */
/* #undef FFTW_USE_F77_MPI_COMM */
/* #undef FFTW_USE_F77_MPI_COMM_P */
