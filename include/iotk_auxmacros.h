# 5 "iotk_auxmacros.spp"

#ifndef __IOTK_AUXMACROS_H
#define __IOTK_AUXMACROS_H

! The macros are defined with -D option or inside iotk_config.h
! The default values are set here
! Maximum rank of an array
#ifndef __IOTK_MAXRANK
#  define __IOTK_MAXRANK 7
#endif
! Minimum value used in iotk_free_unit
#ifndef __IOTK_UNITMIN
#  define __IOTK_UNITMIN 90000
#endif
! Maximum value used in iotk_free_unit
#ifndef __IOTK_UNITMAX
#  define __IOTK_UNITMAX 99999
#endif
! Unit for errors
#ifndef __IOTK_ERROR_UNIT
#  define __IOTK_ERROR_UNIT 0
#endif
! Kind for header in binary files
#ifndef __IOTK_HEADER_KIND
#  define __IOTK_HEADER_KIND selected_int_kind(8)
#endif
! Character (or eventually string) for newline
! It may be adjusted for particular systems
! Unix    achar(10)
! Mac-OS  achar(13)
! Windows ? (now it should be a single byte)
#ifndef __IOTK_NEWLINE
#  define __IOTK_NEWLINE achar(10)
#endif
! Character for EOS
#ifndef __IOTK_EOS
#  define __IOTK_EOS achar(0)
#endif
! These are the default kinds, which depend on the options used
! during the library compilation
! Only default characters are implemented
#define __IOTK_CHARACTER1 iotk_defkind_character
! For logical, integer and real types, the c precompiler
! looks for defined kinds. If no kind is found, the default
! is used as __IOTK_type1
# 54 "iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL1
# 54 "iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL2
# 54 "iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL3
# 54 "iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL4
# 56 "iotk_auxmacros.spp"
#define __IOTK_LOGICAL1 iotk_defkind_LOGICAL
# 59 "iotk_auxmacros.spp"
#endif
# 59 "iotk_auxmacros.spp"
#endif
# 59 "iotk_auxmacros.spp"
#endif
# 59 "iotk_auxmacros.spp"
#endif
# 54 "iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER1
# 54 "iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER2
# 54 "iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER3
# 54 "iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER4
# 56 "iotk_auxmacros.spp"
#define __IOTK_INTEGER1 iotk_defkind_INTEGER
# 59 "iotk_auxmacros.spp"
#endif
# 59 "iotk_auxmacros.spp"
#endif
# 59 "iotk_auxmacros.spp"
#endif
# 59 "iotk_auxmacros.spp"
#endif
# 54 "iotk_auxmacros.spp"
#ifndef __IOTK_REAL1
# 54 "iotk_auxmacros.spp"
#ifndef __IOTK_REAL2
# 54 "iotk_auxmacros.spp"
#ifndef __IOTK_REAL3
# 54 "iotk_auxmacros.spp"
#ifndef __IOTK_REAL4
# 56 "iotk_auxmacros.spp"
#define __IOTK_REAL1 iotk_defkind_REAL
# 59 "iotk_auxmacros.spp"
#endif
# 59 "iotk_auxmacros.spp"
#endif
# 59 "iotk_auxmacros.spp"
#endif
# 59 "iotk_auxmacros.spp"
#endif
# 62 "iotk_auxmacros.spp"

! Some useful check follow
#if __IOTK_MAXRANK > 7
#  error
#endif
#if __IOTK_MAXRANK < 1
#  error
#endif
# 72 "iotk_auxmacros.spp"
#ifdef __IOTK_LOGICAL5
#  error
#endif
# 72 "iotk_auxmacros.spp"
#ifdef __IOTK_LOGICAL6
#  error
#endif
# 72 "iotk_auxmacros.spp"
#ifdef __IOTK_LOGICAL7
#  error
#endif
# 72 "iotk_auxmacros.spp"
#ifdef __IOTK_LOGICAL8
#  error
#endif
# 72 "iotk_auxmacros.spp"
#ifdef __IOTK_LOGICAL9
#  error
#endif
# 72 "iotk_auxmacros.spp"
#ifdef __IOTK_LOGICAL10
#  error
#endif
# 72 "iotk_auxmacros.spp"
#ifdef __IOTK_INTEGER5
#  error
#endif
# 72 "iotk_auxmacros.spp"
#ifdef __IOTK_INTEGER6
#  error
#endif
# 72 "iotk_auxmacros.spp"
#ifdef __IOTK_INTEGER7
#  error
#endif
# 72 "iotk_auxmacros.spp"
#ifdef __IOTK_INTEGER8
#  error
#endif
# 72 "iotk_auxmacros.spp"
#ifdef __IOTK_INTEGER9
#  error
#endif
# 72 "iotk_auxmacros.spp"
#ifdef __IOTK_INTEGER10
#  error
#endif
# 72 "iotk_auxmacros.spp"
#ifdef __IOTK_REAL5
#  error
#endif
# 72 "iotk_auxmacros.spp"
#ifdef __IOTK_REAL6
#  error
#endif
# 72 "iotk_auxmacros.spp"
#ifdef __IOTK_REAL7
#  error
#endif
# 72 "iotk_auxmacros.spp"
#ifdef __IOTK_REAL8
#  error
#endif
# 72 "iotk_auxmacros.spp"
#ifdef __IOTK_REAL9
#  error
#endif
# 72 "iotk_auxmacros.spp"
#ifdef __IOTK_REAL10
#  error
#endif
# 77 "iotk_auxmacros.spp"

! Complex are treated indentically to reals
! These lines map the definitions.
# 81 "iotk_auxmacros.spp"
#ifdef __IOTK_REAL1
#  define __IOTK_COMPLEX1 __IOTK_REAL1
#else
#  undef __IOTK_COMPLEX1
#endif
# 81 "iotk_auxmacros.spp"
#ifdef __IOTK_REAL2
#  define __IOTK_COMPLEX2 __IOTK_REAL2
#else
#  undef __IOTK_COMPLEX2
#endif
# 81 "iotk_auxmacros.spp"
#ifdef __IOTK_REAL3
#  define __IOTK_COMPLEX3 __IOTK_REAL3
#else
#  undef __IOTK_COMPLEX3
#endif
# 81 "iotk_auxmacros.spp"
#ifdef __IOTK_REAL4
#  define __IOTK_COMPLEX4 __IOTK_REAL4
#else
#  undef __IOTK_COMPLEX4
#endif
# 87 "iotk_auxmacros.spp"

#endif
