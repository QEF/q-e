! Input/Output Tool Kit (IOTK)
! Copyright (C) 2004 Giovanni Bussi
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!------------------------------------------------------------------------------!
! CONFIGURATION FILE FOR IOTK 0.3.1
!------------------------------------------------------------------------------!
! The following lines map some commonly defined system macro to the internal
! iotk macros.
! Iotk macros which are not defined take their default values.
! See the manual for a list of iotk macros.

#define __IOTK_GENERIC_ARCHITECTURE

#ifdef __AIX
@PROCESS OPTIMIZE(0)
#   define __IOTK_BINARY_FORMAT "IBM-SP/XLF"
#   define __IOTK_LOGICAL1 1
#   define __IOTK_LOGICAL2 2
#   define __IOTK_LOGICAL3 4
#   define __IOTK_LOGICAL4 8
#   define __IOTK_INTEGER1 1
#   define __IOTK_INTEGER2 2
#   define __IOTK_INTEGER3 4
#   define __IOTK_INTEGER4 8
#   define __IOTK_REAL1    4
#   define __IOTK_REAL2    8
#   define __IOTK_REAL3    16
#   define __IOTK_AVOID_EMPTY_FILES
#   undef __IOTK_GENERIC_ARCHITECTURE
#endif


#ifdef __LINUX
#   ifdef __INTEL
#         define __IOTK_BINARY_FORMAT "PC-LINUX/IFC"
#         define __IOTK_LOGICAL1 1
#         define __IOTK_LOGICAL2 2
#         define __IOTK_LOGICAL3 4
#         define __IOTK_LOGICAL4 8
#         define __IOTK_INTEGER1 1
#         define __IOTK_INTEGER2 2
#         define __IOTK_INTEGER3 4
#         define __IOTK_INTEGER4 8
#         define __IOTK_REAL1    4
#         define __IOTK_REAL2    8
#         define __IOTK_REAL3    16
#         define __IOTK_AVOID_EMPTY_FILES
#         undef __IOTK_GENERIC_ARCHITECTURE
#   endif
#   ifdef __G95
#         define __IOTK_BINARY_FORMAT "PC-LINUX/G95"
#         define __IOTK_LOGICAL1 1
#         define __IOTK_LOGICAL2 2
#         define __IOTK_LOGICAL3 4
#         define __IOTK_LOGICAL4 8
#         define __IOTK_INTEGER1 1
#         define __IOTK_INTEGER2 2
#         define __IOTK_INTEGER3 4
#         define __IOTK_INTEGER4 8
#         define __IOTK_REAL1    4
#         define __IOTK_REAL2    8
#         define __IOTK_AVOID_EMPTY_FILES
#         undef __IOTK_GENERIC_ARCHITECTURE
#   endif
#endif


#if defined (__LAHEY)
#    define __IOTK_UNKNOWN_PROCESS
#    define __IOTK_BINARY_FORMAT "UNKNOWN"
#    define __IOTK_LOGICAL1 1
#    define __IOTK_LOGICAL2 2
#    define __IOTK_LOGICAL3 4
#    define __IOTK_LOGICAL4 8
#    define __IOTK_INTEGER1 4
#    define __IOTK_REAL1    8
#    define __IOTK_MAXRANK 5
#    undef __IOTK_GENERIC_ARCHITECTURE
#endif


#ifdef __SGI
#    define __IOTK_BINARY_FORMAT "SGI-ORIGIN"
#    define __IOTK_LOGICAL1 1
#    define __IOTK_LOGICAL2 2
#    define __IOTK_LOGICAL3 4
#    define __IOTK_LOGICAL4 8
#    define __IOTK_INTEGER1 1
#    define __IOTK_INTEGER2 2
#    define __IOTK_INTEGER3 4
#    define __IOTK_INTEGER4 8
#    define __IOTK_REAL1    4
#    define __IOTK_REAL2    8
#    define __IOTK_REAL3    16
#    define __IOTK_AVOID_EMPTY_FILES
#    undef __IOTK_GENERIC_ARCHITECTURE
#endif


#ifdef __IOTK_GENERIC_ARCHITECTURE
#    define __IOTK_UNKNOWN_PROCESS
#    define __IOTK_BINARY_FORMAT "UNKNOWN"
#    define __IOTK_LOGICAL1 1
#    define __IOTK_LOGICAL2 2
#    define __IOTK_LOGICAL3 4
#    define __IOTK_LOGICAL4 8
#    define __IOTK_INTEGER1 4
#    define __IOTK_INTEGER2 2
#    define __IOTK_INTEGER3 1
#    define __IOTK_INTEGER4 8
#    define __IOTK_REAL1    8
#    define __IOTK_REAL2    4
#    define __IOTK_MAXRANK 5
#endif


#ifdef __MPI
#  define __IOTK_MPI_ABORT
#endif

