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
  use clib_wrappers, only : memstat, get_mem_usage, get_mem_avail
  implicit none
  integer, intent(in) :: out_unit
  character(len=*), intent(in) :: title
  !
  integer mem_used, mem_avail
  !
  write(out_unit, '(/, A, A)') "---- Real-time Memory Report at ", trim(title)
  mem_used = get_mem_usage()
  write(out_unit, '(I14, A)') mem_used / 1024, " MiB given to the printing process from OS"
  call memstat(mem_used)
  write(out_unit, '(I14, A)') mem_used / 1024, " MiB allocation reported by mallinfo(arena+hblkhd)"
  mem_avail = get_mem_avail()
  write(out_unit, '(I14, A)') mem_avail / 1024, " MiB available memory on the node where the printing process lives"
  call print_gpu_mem(out_unit)
  write(out_unit, '(A)') "------------------"
end subroutine print_mem_usage

subroutine print_gpu_mem(out_unit)
#if defined(__CUDA)
  use cudafor
  implicit none
  integer, intent(in) :: out_unit
  integer(kind=cuda_count_kind) :: freeMem, totalMem, usedMB, freeMB, totalMB
  integer :: istat
  istat = CudaMemGetInfo(freeMem, totalMem)
  totalMB = totalMem/(1024*1024)
  freeMB = freeMem/(1024*1024);
  usedMB = (totalMem-freeMem)/(1024*1024);
  write(out_unit,'(5X, "GPU memory used/free/total (MiB): ", I0, A3, I0,A3, I0)') &
              & usedMB, " / ", freeMB, " / ", totalMB
#else
  implicit none
  integer, intent(in) :: out_unit
#endif
end subroutine
