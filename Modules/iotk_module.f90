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

module iotk_module
! This module is a wrapper for the entities in iotk_base which need to be public
  use iotk_base
  use iotk_interface
! All names are private ...
  private
! ... except the names listed below
  public :: iotk_open_write
  public :: iotk_open_read
  public :: iotk_close_write
  public :: iotk_close_read
  public :: iotk_write_begin
  public :: iotk_write_end
  public :: iotk_write_pi
  public :: iotk_write_comment
  public :: iotk_write_empty
  public :: iotk_write_dat
  public :: iotk_write_attr
  public :: iotk_scan_begin
  public :: iotk_scan_end
  public :: iotk_scan_pi
  public :: iotk_scan_empty
  public :: iotk_scan_dat
  public :: iotk_scan_attr
  public :: iotk_taglenx
  public :: iotk_attlenx
  public :: iotk_vallenx
  public :: iotk_namlenx
  public :: iotk_fillenx
  public :: iotk_index
  public :: iotk_version
  public :: iotk_binary_format
  public :: iotk_header_kind
  public :: iotk_copy_tag
  public :: iotk_unit_print
  public :: iotk_unit_get
  public :: iotk_free_unit
  public :: iotk_basefmt
  public :: iotk_defkind_character
  public :: iotk_defkind_logical
  public :: iotk_defkind_integer
  public :: iotk_defkind_real
  public :: iotk_defkind_complex
  public :: iotk_print_kinds
  public :: iotk_set_options
  public :: iotk_get_options
  public :: iotk_getline
  public :: iotk_phys_unit
  public :: iotk_link
  public :: iotk_read
  public :: iotk_copyfile
  public :: iotk_newline
  public :: iotk_eos
  public :: iotk_error_clear
  public :: iotk_error_print
end module iotk_module

