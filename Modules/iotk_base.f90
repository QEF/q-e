# 1 "iotk_base.spp"
! Input/Output Tool Kit (IOTK)
! Copyright (C) 2004,2005 Giovanni Bussi
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

# 28 "iotk_base.spp"
#include "iotk_auxmacros.h"
# 30 "iotk_base.spp"

module iotk_base
implicit none
save

! In this module, all names are public
! For this reason, it should not be used directly by the end user.

! This line set the version string
character(9),      parameter :: iotk_version            = "1.0.0test"
integer,                            parameter :: iotk_version_major      = 1
integer,                            parameter :: iotk_version_minor      = 0
integer,                            parameter :: iotk_version_patch      = 0
character(3), parameter :: iotk_file_version       = "1.0"
integer,                            parameter :: iotk_file_version_major = 1
integer,                            parameter :: iotk_file_version_minor = 0

character, parameter :: iotk_newline = __IOTK_NEWLINE
character, parameter :: iotk_eos     = __IOTK_EOS

! Max number of controls
integer, parameter :: iotk_ncontrol = 255 ! (2**8-1)

! Max length for strings
integer, parameter :: iotk_taglenx =  65535 ! (2**16-1)
integer, parameter :: iotk_namlenx =  256
integer, parameter :: iotk_attlenx =  iotk_taglenx - iotk_namlenx - 1 ! for space
integer, parameter :: iotk_vallenx =  32768
integer, parameter :: iotk_linlenx =  4096
integer, parameter :: iotk_fillenx =  1024
integer, parameter :: iotk_linlen  =  128
integer, parameter :: iotk_indent  =    2
integer, parameter :: iotk_maxindent = 12

! These options can be modified runtime
! Margins for unit search
integer :: iotk_unitmin = __IOTK_UNITMIN
integer :: iotk_unitmax = __IOTK_UNITMAX
integer :: iotk_error_unit = __IOTK_ERROR_UNIT
! Size of the buffer for iotk_getline
integer :: iotk_getline_buffer = 1024
! If true, exhausting the error pool causes an overflow warning
logical :: iotk_error_warn_overflow = .false.

! Kind for the header integer (number of digits in (iotk_ncontrol+1)*(iotk_taglenx+1))
integer, parameter :: iotk_header_kind = __IOTK_HEADER_KIND

! Map of controls into XML tags
! control = 1 <       >
! control = 2 </      >
! control = 3 <      />
! control = 4 <!--  -->
! control = 5 <?     ?>
! control = 128 is a special tag for binary files (continuation tag)

! Alphabet
character(26), parameter :: lowalphabet = "abcdefghijklmnopqrstuvwxyz"
character(26), parameter :: upalphabet  = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
character(52), parameter :: alphabet    = lowalphabet//upalphabet
character(53), parameter :: alphabet_   = alphabet//"_"
character(10), parameter :: numbers     = "0123456789"

! Rules for names
character(54), parameter :: iotk_namcharfirst = alphabet//"_:"
character(66), parameter :: iotk_namchar      = iotk_namcharfirst//numbers//".-"

! Default kinds, depending on compilers and compilation options for the library source
integer, parameter :: iotk_defkind_character = kind("a")
integer, parameter :: iotk_defkind_logical   = kind(.true.)
integer, parameter :: iotk_defkind_integer   = kind(1)
integer, parameter :: iotk_defkind_real      = kind(1.0)
integer, parameter :: iotk_defkind_complex   = kind(1.0)

! Maximum allowed rank
integer, parameter :: iotk_maxrank      = __IOTK_MAXRANK ! Controlled by cpp
integer, parameter :: iotk_maxrank_hard = 7       ! Controlled by sprep preprocessing

! Internal type dealing with io units
type iotk_unit
  integer                     :: unit  ! fortran unit
  character(iotk_namlenx)     :: root  ! name of the root tag
  logical                     :: skip_root ! if true, root tag is not written automatically
  logical                     :: raw   ! if true, the file is raw data
  integer                     :: level ! the hierarchical level inside the file
  logical                     :: close_at_end ! if true, the file has to be fortran-closed when iotk_close_* is called
  type (iotk_unit),   pointer :: son    ! a pointer to the son in the multi-file model
  type (iotk_unit),   pointer :: parent ! a pointer to the parent in the multi-file model
  type (iotk_unit),   pointer :: next   ! a pointer to the next unit in the linked list
end type iotk_unit

! Special type used to force optional argument labelling.
type iotk_dummytype
  integer :: dummy
end type iotk_dummytype

! Linked list of iotk_unit objects
logical                   :: iotk_units_init = .false.
type (iotk_unit), pointer :: iotk_units

type iotk_error
  character, pointer :: str(:)
end type iotk_error

integer, parameter :: iotk_error_linelength  = 120
integer, parameter :: iotk_error_pool_size   = 100

type(iotk_error) :: iotk_error_pool       (iotk_error_pool_size)
logical          :: iotk_error_pool_used  (iotk_error_pool_size) = .false.
integer          :: iotk_error_pool_order (iotk_error_pool_size) = 0

end module iotk_base
