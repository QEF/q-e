!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine print_clock_kcw
  !-----------------------------------------------------------------------
  !
  !! Print timing
  !
  USE io_global,  ONLY : stdout
  USE uspp,       ONLY : nlcc_any
  USE uspp_init,  ONLY : init_us_2
  implicit none
  !
  WRITE( stdout, * )
  call print_clock ('KCW')
  WRITE( stdout,  * ) '    INITIALIZATION: '
  call print_clock ('phq_setup')
  call print_clock ('phq_init')
  call print_clock ('map')
  call print_clock ('rho_of_q')
  WRITE( stdout, * )
  if (nlcc_any) call print_clock ('set_drhoc')
  call print_clock ('init_vloc')
  call print_clock ('init_us_1')
  !call print_clock ('init_us_2')
  call print_clock ('newd')
  call print_clock ('dvanqq')
  call print_clock ('drho')

  WRITE( stdout, * )
  call print_clock ('phqscf')
  call print_clock ('solve_linter')
  call print_clock ('kcw_run_nscf')
  WRITE( stdout, * )
  call print_clock ('solve_linter')
  call print_clock ('dvqpsi_us')
  call print_clock ('ortho')
  call print_clock ('cgsolve')
  call print_clock ('incdrhoscf')
  call print_clock ('addusddens')
  call print_clock ('vpsifft')
  call print_clock ('dv_of_drho')
  call print_clock ('mix_pot')
  call print_clock ('ef_shift')
  call print_clock ('localdos')
#if defined(__MPI)
  call print_clock ('psymdvscf')
#else
  call print_clock ('symdvscf')
#endif
  call print_clock ('newdq')
  call print_clock ('adddvscf')


  call print_clock ('drhodvus')
  WRITE( stdout, * )
  call print_clock ('dvqpsi_us')

  call print_clock ('dvqpsi_us_on')
  WRITE( stdout, * )
  call print_clock ('cgsolve')

  call print_clock ('ch_psi')
  WRITE( stdout, * )
  call print_clock ('ch_psi')
  call print_clock ('first')
  call print_clock ('h_psi')

  call print_clock ('last')
  WRITE( stdout, * )
  call print_clock ('h_psi')
  call print_clock ('firstfft')
  call print_clock ('product')
  call print_clock ('secondfft')

  call print_clock ('add_vuspsi')
  WRITE( stdout, * )
  call print_clock ('incdrhoscf')

  call print_clock ('addusdbec')
  WRITE( stdout, * )
  call print_clock ('drhodvus')

  call print_clock ('addusddort')
  WRITE( stdout, * )
  WRITE( stdout,  * ) '     General routines'
  call print_clock ('calbec')
  call print_clock ('fft')
  call print_clock ('ffts')
  call print_clock ('fftw')
  call print_clock ('cinterpolate')
  call print_clock ('davcio')
  call print_clock ('write_rec')
  WRITE( stdout, * )
  return
end subroutine print_clock_kcw
