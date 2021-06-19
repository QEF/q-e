!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine print_mem_usage (out_unit, title)
  !! in human-readable format. Uses a standard f90 call.
  !
  USE wrappers,  ONLY : memstat, get_mem_usage, get_mem_avail
  implicit none
  integer, intent(in) :: out_unit
  character(len=*), intent(in) :: title
  !
  INTEGER mem_used, mem_avail
  !
  write(out_unit, '(/, A, A)') "---- Real-time Memory Report at ", trim(title)
  mem_used = get_mem_usage()
  write(out_unit, '(I14, A)') mem_used / 1024, " MiB given to the printing process from OS"
  call memstat(mem_used)
  write(out_unit, '(I14, A)') mem_used / 1024, " MiB allocation reported by mallinfo(arena+hblkhd)"
  mem_avail = get_mem_avail()
  write(out_unit, '(I14, A)') mem_avail / 1024, " MiB available memory on the node where the printing process lives"
  write(out_unit, '(A)') "------------------"
end subroutine print_mem_usage
