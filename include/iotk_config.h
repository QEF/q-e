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
! CONFIGURATION FILE FOR IOTK 0.3.5
!------------------------------------------------------------------------------!
! The following lines map some commonly defined system macro to the internal
! iotk macros.
! Iotk macros which are not defined take their default values.
! See the manual for a list of iotk macros.

#ifndef __IOTK_CONFIG_H

#define __IOTK_CONFIG_H

#  define __IOTK_MAXRANK 4

#ifdef __AIX
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
#endif

#ifdef __LINUX
#   ifdef __INTEL
#         define __IOTK_BINARY_FORMAT "PC-LINUX/IFC"
#         define __IOTK_INTEGER1 1
#         define __IOTK_INTEGER2 2
#         define __IOTK_INTEGER3 4
#         define __IOTK_REAL1    4
#         define __IOTK_REAL2    8
#         define __IOTK_WORKAROUND1
#         define __IOTK_WORKAROUND3
#         define __IOTK_WORKAROUND5
#   endif
#   ifdef __G95
#         define __IOTK_BINARY_FORMAT "PC-LINUX/G95"
#         define __IOTK_INTEGER1 1
#         define __IOTK_INTEGER2 2
#         define __IOTK_INTEGER3 4
#         define __IOTK_REAL1    4
#         define __IOTK_REAL2    8
#   endif
#   ifdef __PGI
#         define __IOTK_BINARY_FORMAT "PC-LINUX/PGI"
#         undef __IOTK_LOGICAL1 1
#         undef __IOTK_LOGICAL2 2
#         define __IOTK_LOGICAL3 4
#         undef __IOTK_LOGICAL4 8
#         undef __IOTK_INTEGER1 1
#         undef __IOTK_INTEGER2 2
#         define __IOTK_INTEGER3 4
#         undef __IOTK_INTEGER4 8
#         undef __IOTK_REAL1    4
#         define __IOTK_REAL2    8
#         define __IOTK_WORKAROUND2
#         define __IOTK_WORKAROUND4
#   endif
#   ifdef __NAG
#         define __IOTK_INTEGER1 1
#         define __IOTK_INTEGER2 2
#         define __IOTK_INTEGER3 3
#         define __IOTK_INTEGER4 4
#         define __IOTK_LOGICAL1 1
#         define __IOTK_LOGICAL2 2
#         define __IOTK_LOGICAL3 3
#         define __IOTK_LOGICAL4 4
#         define __IOTK_REAL1 1
#         define __IOTK_REAL2 2
#         define __IOTK_WORKAROUND4
#   endif
#endif

#ifdef __LINUX64
#   ifdef __INTEL
#         define __IOTK_BINARY_FORMAT "PC-LINUX/IFC"
#         define __IOTK_INTEGER1 1
#         define __IOTK_INTEGER2 2
#         define __IOTK_INTEGER3 4
#         define __IOTK_REAL1    4
#         define __IOTK_REAL2    8
#         define __IOTK_WORKAROUND1
#         define __IOTK_WORKAROUND3
#         define __IOTK_WORKAROUND5
#   endif
#endif

#ifdef __SGI
#   define __IOTK_BINARY_FORMAT "SGI-ORIGIN"
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
#endif

#ifdef __PARA
#  define __IOTK_MPI_ABORT
#endif

#endif
