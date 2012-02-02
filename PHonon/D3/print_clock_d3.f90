!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

subroutine print_clock_d3
  USE io_global,  ONLY : stdout
  use d3com
  implicit none
  WRITE( stdout, * )
  call print_clock ('D3TOTEN')
  call print_clock ('d3_setup')

  call print_clock ('phq_init')
  WRITE( stdout, * )
  call print_clock ('solve_linter')
  call print_clock ('ortho')
  call print_clock ('cgsolve')
  call print_clock ('incdrhoscf')
  call print_clock ('dv_of_drho')
#ifdef __MPI
  call print_clock ('psymdvscf')
  call print_clock ('psymd0rho')
#else

  call print_clock ('symdvscf')
#endif
  WRITE( stdout, * )
  call print_clock ('cgsolve')

  call print_clock ('ch_psi')
  WRITE( stdout, * )
  call print_clock ('ch_psi')
  call print_clock ('h_psiq')

  call print_clock ('last')
  WRITE( stdout, * )
  call print_clock ('h_psiq')
  call print_clock ('firstfft')
  call print_clock ('product')

  call print_clock ('secondfft')
  WRITE( stdout, * )
  WRITE( stdout,  * ) '     General routines'
  call print_clock ('calbec')
  call print_clock ('fft')
  call print_clock ('ffts')
  call print_clock ('fftw')
  call print_clock ('cinterpolate')
  call print_clock ('davcio')
  WRITE( stdout, * )
#ifdef __MPI
  WRITE( stdout,  * ) '     Parallel routines'
  call print_clock ('reduce')
#endif
  return
end subroutine print_clock_d3
