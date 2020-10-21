! Input/Output Tool Kit (IOTK)
! Copyright (C) 2004-2006 Giovanni Bussi
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

!------------------------------------------------------------------------------!
! Inclusion of configuration file
#include "iotk_config.h"
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
! Inclusion of the auxiliary macros
#include "iotk_auxmacros.h"
!------------------------------------------------------------------------------!

module iotk_base
implicit none
save

!------------------------------------------------------------------------------!
! In this module, all names are public
! For this reason, it should not be used directly by the end user.
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! Version strings and integer constants
character(5),      parameter :: iotk_version            = "1.2.0"
integer,                            parameter :: iotk_version_major      = 1
integer,                            parameter :: iotk_version_minor      = 2
integer,                            parameter :: iotk_version_patch      = 0
character(3), parameter :: iotk_file_version       = "1.0"
integer,                            parameter :: iotk_file_version_major = 1
integer,                            parameter :: iotk_file_version_minor = 0
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! Maximum allowed rank
integer, parameter :: iotk_maxrank_hard = 7       ! Controlled by sprep preprocessing
integer, parameter :: iotk_maxrank      = __IOTK_MAXRANK ! Controlled by cpp
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! Default kinds, depending on compilers and compilation options for the library source
integer, parameter :: iotk_character_defkind = kind("a")
integer, parameter :: iotk_logical_defkind   = kind(.true.)
integer, parameter :: iotk_integer_defkind   = kind(1)
integer, parameter :: iotk_real_defkind      = kind(1.0)
integer, parameter :: iotk_complex_defkind   = kind(1.0)
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! Kinds for the multiple interfaces
#ifdef __IOTK_LOGICAL1
integer, parameter :: iotk_LOGICAL1 = __IOTK_LOGICAL1
#endif
#ifdef __IOTK_LOGICAL2
integer, parameter :: iotk_LOGICAL2 = __IOTK_LOGICAL2
#endif
#ifdef __IOTK_LOGICAL3
integer, parameter :: iotk_LOGICAL3 = __IOTK_LOGICAL3
#endif
#ifdef __IOTK_LOGICAL4
integer, parameter :: iotk_LOGICAL4 = __IOTK_LOGICAL4
#endif
#ifdef __IOTK_INTEGER1
integer, parameter :: iotk_INTEGER1 = __IOTK_INTEGER1
#endif
#ifdef __IOTK_INTEGER2
integer, parameter :: iotk_INTEGER2 = __IOTK_INTEGER2
#endif
#ifdef __IOTK_INTEGER3
integer, parameter :: iotk_INTEGER3 = __IOTK_INTEGER3
#endif
#ifdef __IOTK_INTEGER4
integer, parameter :: iotk_INTEGER4 = __IOTK_INTEGER4
#endif
#ifdef __IOTK_REAL1
integer, parameter :: iotk_REAL1 = __IOTK_REAL1
#endif
#ifdef __IOTK_REAL2
integer, parameter :: iotk_REAL2 = __IOTK_REAL2
#endif
#ifdef __IOTK_REAL3
integer, parameter :: iotk_REAL3 = __IOTK_REAL3
#endif
#ifdef __IOTK_REAL4
integer, parameter :: iotk_REAL4 = __IOTK_REAL4
#endif
#ifdef __IOTK_COMPLEX1
integer, parameter :: iotk_COMPLEX1 = __IOTK_COMPLEX1
#endif
#ifdef __IOTK_COMPLEX2
integer, parameter :: iotk_COMPLEX2 = __IOTK_COMPLEX2
#endif
#ifdef __IOTK_COMPLEX3
integer, parameter :: iotk_COMPLEX3 = __IOTK_COMPLEX3
#endif
#ifdef __IOTK_COMPLEX4
integer, parameter :: iotk_COMPLEX4 = __IOTK_COMPLEX4
#endif
#ifdef __IOTK_CHARACTER1
integer, parameter :: iotk_CHARACTER1 = __IOTK_CHARACTER1
#endif
#ifdef __IOTK_CHARACTER2
integer, parameter :: iotk_CHARACTER2 = __IOTK_CHARACTER2
#endif
#ifdef __IOTK_CHARACTER3
integer, parameter :: iotk_CHARACTER3 = __IOTK_CHARACTER3
#endif
#ifdef __IOTK_CHARACTER4
integer, parameter :: iotk_CHARACTER4 = __IOTK_CHARACTER4
#endif
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! Kind for the header integer (number of digits in (iotk_ncontrol+1)*(iotk_taglenx+1))
integer, parameter :: iotk_header_kind = __IOTK_HEADER_KIND
!------------------------------------------------------------------------------!

#ifdef __IOTK_STREAMS
!------------------------------------------------------------------------------!
! Kinds for stream i/o compatibility
integer, parameter :: iotk_record_kind   = __IOTK_RECORD_KIND
integer, parameter :: iotk_record_length = __IOTK_RECORD_LENGTH
!------------------------------------------------------------------------------!
#endif

!------------------------------------------------------------------------------!
! Special characters
character, parameter :: iotk_newline = __IOTK_NEWLINE
character, parameter :: iotk_eos     = __IOTK_EOS
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! Max number of controls
integer, parameter :: iotk_ncontrol = 255 ! (2**8-1)
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! Map of controls into XML tags
! control = 1 <       >
! control = 2 </      >
! control = 3 <      />
! control = 4 <!--  -->
! control = 5 <?     ?>
! control = 128 is a special tag for binary files (continuation tag)
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! Max lengths for strings
integer, parameter :: iotk_taglenx =  65535 ! (2**16-1)
integer, parameter :: iotk_namlenx =  256
integer, parameter :: iotk_attlenx =  iotk_taglenx - iotk_namlenx - 1 ! for space
integer, parameter :: iotk_vallenx =  32768
integer, parameter :: iotk_linlenx =  4096
integer, parameter :: iotk_fillenx =  1024
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! Max number of arguments for the iotk tool
integer, parameter :: iotk_maxargs = __IOTK_MAXARGS
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! Alphabet
character(26), parameter :: lowalphabet = "abcdefghijklmnopqrstuvwxyz"
character(26), parameter :: upalphabet  = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
character(52), parameter :: alphabet    = lowalphabet//upalphabet
character(53), parameter :: alphabet_   = alphabet//"_"
character(10), parameter :: numbers     = "0123456789"
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! List of characters which are not separators in a dat or attribute array
character(66), parameter :: not_separator = alphabet_//numbers//"+-."
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! Rules for names
character(54), parameter :: iotk_namcharfirst = alphabet//"_:"
character(66), parameter :: iotk_namchar      = iotk_namcharfirst//numbers//".-"
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! Internal type dealing with io units
type iotk_unit
  integer                     :: unit  ! fortran unit
  character(iotk_namlenx)     :: root  ! name of the root tag
  logical                     :: skip_root ! if true, root tag is not written automatically
  logical                     :: raw   ! if true, the file is raw data
!-<
  logical                     :: qe_syntax ! if true, the file is read/written using the new syntax
!->
  integer                     :: level ! the hierarchical level inside the file
  logical                     :: close_at_end ! if true, the file has to be fortran-closed when iotk_close_* is called
  type (iotk_unit),   pointer :: son    ! a pointer to the son in the multi-file model
  type (iotk_unit),   pointer :: parent ! a pointer to the parent in the multi-file model
end type iotk_unit
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! Special type used to force optional argument labelling.
type iotk_dummytype
  integer :: dummy
end type iotk_dummytype
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! Internal type dealing with error messages
type iotk_error
  character, pointer :: str(:)
end type iotk_error
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! Max length of a line in the error message. Any longer line will be cut
integer, parameter :: iotk_error_linelength  = 120
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! Maximum number of errors which can be handled at the same time
integer, parameter :: iotk_error_pool_size   = 100
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! Static pool of errors
type(iotk_error) :: iotk_error_pool       (iotk_error_pool_size)
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! Flags concerning the error pool:
! If true, that element of the pool is in usage
logical          :: iotk_error_pool_used  (iotk_error_pool_size) = .false.
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! These integers are set in increasing order to trace the order errors are raised
! They are then used to eliminate old errors if the user forgets to do that
integer          :: iotk_error_pool_order (iotk_error_pool_size) = 0
!------------------------------------------------------------------------------!

! The following options can be modified runtime
! X_def is the default value for variable X

!------------------------------------------------------------------------------!
! Margins for unit search
integer, parameter :: iotk_unitmin_def     = __IOTK_UNITMIN
integer            :: iotk_unitmin         = iotk_unitmin_def
integer, parameter :: iotk_unitmax_def     = __IOTK_UNITMAX
integer            :: iotk_unitmax         = iotk_unitmax_def
integer, parameter :: iotk_error_unit_def  = __IOTK_ERROR_UNIT
integer            :: iotk_error_unit      = iotk_error_unit_def
integer, parameter :: iotk_output_unit_def = __IOTK_OUTPUT_UNIT
integer            :: iotk_output_unit     = iotk_output_unit_def
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! Size of the buffer for iotk_getline
! (it is intended for efficiency; the total length of a line should be <= iotk_linlenx)
integer, parameter :: iotk_getline_buffer_def = 1024
integer            :: iotk_getline_buffer     = iotk_getline_buffer_def
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! Parameters for text file beautyfier
! Length of a line:
integer, parameter :: iotk_linlen_def    = 128
integer            :: iotk_linlen        = iotk_linlen_def
! Number of spaces for each indentation
integer, parameter :: iotk_indent_def    = 2
integer            :: iotk_indent        = iotk_indent_def
! Maximum number of spaces in indentation
integer, parameter :: iotk_maxindent_def = 12
integer            :: iotk_maxindent     = iotk_maxindent_def
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! If true, exhausting the error pool causes an overflow warning
logical, parameter :: iotk_error_warn_overflow_def = .false.
logical            :: iotk_error_warn_overflow     = iotk_error_warn_overflow_def
!------------------------------------------------------------------------------!

end module iotk_base
