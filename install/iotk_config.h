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
! CONFIGURATION FILE FOR IOTK 1.1.0 for Quantum-Espresso
!------------------------------------------------------------------------------!
! The following lines map some commonly defined system macro to the internal
! iotk macros.
! Iotk macros which are not defined take their default values.
! See the manual for a list of iotk macros.

#ifndef __IOTK_CONFIG_H
#define __IOTK_CONFIG_H

! Generic options valid for quantum-espresso
! QE uses ranks up to four and default integer/logicals only

#define __IOTK_MAXRANK  4

! some compilers do not like the following
!    #define __IOTK_REAL1 selected_real_kind(6,30)
!    #define __IOTK_REAL2 selected_real_kind(14,200)
! so we use explicit kinds

! Something from an ancient past that probably is not valid anymore...
! #if defined(__NAG)
! #   define __IOTK_REAL1 1
! #   define __IOTK_REAL2 2
! #else
#if defined(__SX6)
#   define __IOTK_REAL2 8
#else
#   define __IOTK_REAL1 4
#   define __IOTK_REAL2 8
#endif

! Machine-dependent options
! Only for compilers that require some special tricks

#ifdef __IOTK_SAFEST
    !
    ! force to define all the workarounds
    !
#   define __IOTK_WORKAROUND1
#   define __IOTK_WORKAROUND2
#   define __IOTK_WORKAROUND3
#   define __IOTK_WORKAROUND4
#   define __IOTK_WORKAROUND5
#   define __IOTK_WORKAROUND6
#   define __IOTK_WORKAROUND7
#   define __IOTK_WORKAROUND9
#   define __IOTK_WORKAROUND10
#else
    !
    ! proceed with a machine dependent def where available
    !
#   if defined(__XLF)
#      define __IOTK_WORKAROUND5
#      define __IOTK_WORKAROUND9
#      define __IOTK_WORKAROUND10
#   elif defined(__INTEL)
#      define __IOTK_WORKAROUND1
#      define __IOTK_WORKAROUND3
#      define __IOTK_WORKAROUND5
#   elif defined(__PGI)
#      define __IOTK_WORKAROUND2
#      define __IOTK_WORKAROUND4
#   elif defined(__NAG)
#      define __IOTK_WORKAROUND4
#   elif defined(__ALPHA)
#      define __IOTK_WORKAROUND1
#      define __IOTK_WORKAROUND6
#      define __IOTK_WORKAROUND8
#   elif defined(__SX6)
#      define __IOTK_WORKAROUND5
#      define __IOTK_WORKAROUND7
#   else
#      define __IOTK_WORKAROUND1
#      define __IOTK_WORKAROUND2
#      define __IOTK_WORKAROUND3
#      define __IOTK_WORKAROUND4
#      define __IOTK_WORKAROUND5
#      define __IOTK_WORKAROUND6
#      define __IOTK_WORKAROUND7
#   endif
#endif


#if defined(__MPI)
#  define __IOTK_MPI_ABORT
#endif

#endif
