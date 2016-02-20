!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE print_clock_vdw
  !-----------------------------------------------------------------------

  USE io_global,  ONLY : stdout
  USE uspp,       ONLY: okvan
  USE ramanm,     ONLY: lraman, elop
  IMPLICIT NONE
  !
  WRITE( stdout, * )
  CALL print_clock ('VdW')
  !
  WRITE( stdout, * )
  WRITE( stdout,  * ) '    INITIALIZATION: '
  CALL print_clock ('vdw_init')
  WRITE( stdout, * )
  !
  WRITE( stdout, * )
  WRITE( stdout,  * ) '    EFFECTIVE POTENTIAL:'
  CALL print_clock ('eff_pot')
  CALL print_clock ('ite_veff')
  !
  WRITE( stdout, * )
  WRITE( stdout,  * ) '    POLARIZABILITY:'
  CALL print_clock ('solve_e')
  CALL print_clock ('polariz')
  WRITE( stdout, * )
  CALL print_clock ('dvpsi_e')
  CALL print_clock ('adddvscf')
  CALL print_clock ('gmressolve')
  CALL print_clock ('incdrhoscf')
  CALL print_clock ('dv_of_drho')
  CALL print_clock ('mix_pot')
  !
  WRITE( stdout, * )
  WRITE( stdout,  * ) '    General routines'
  CALL print_clock ('calbec')
  CALL print_clock ('fft')
  CALL print_clock ('ffts')
  CALL print_clock ('fftw')
  CALL print_clock ('davcio')
  WRITE( stdout, * )
#ifdef __MPI
  WRITE( stdout,  * ) '    Parallel routines'
  CALL print_clock ('reduce')
  WRITE( stdout, * )
#endif
  RETURN
END SUBROUTINE print_clock_vdw
