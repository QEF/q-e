!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine print_clock_vdw
  !-----------------------------------------------------------------------

  USE io_global,  ONLY : stdout
  USE uspp, only: okvan
  USE control_ph
  USE ramanm, ONLY: lraman, elop
  USE nlcc_ph, ONLY: nlcc_any
  implicit none
  !
  WRITE( stdout, * )
  call print_clock ('VdW')
  !
  WRITE( stdout, * )
  WRITE( stdout,  * ) '    INITIALIZATION: '
  call print_clock ('vdw_init')
  WRITE( stdout, * )
  !
  WRITE( stdout, * )
  WRITE( stdout,  * ) '    EFFECTIVE POTENTIAL:'
  call print_clock ('eff_pot')
  call print_clock ('ite_veff')
  !
  WRITE( stdout, * )
  WRITE( stdout,  * ) '    POLARIZABILITY:'
  call print_clock ('solve_e')
  call print_clock ('polariz')
  WRITE( stdout, * )
  call print_clock ('dvpsi_e')
  call print_clock ('adddvscf')
  call print_clock ('gmressolve')
  call print_clock ('incdrhoscf')
  call print_clock ('dv_of_drho')
  call print_clock ('mix_pot')
  !
  WRITE( stdout, * )
  WRITE( stdout,  * ) '    General routines'
  call print_clock ('ccalbec')
  call print_clock ('cft3')
  call print_clock ('cft3s')
  call print_clock ('davcio')
  WRITE( stdout, * )
#ifdef __PARA
  WRITE( stdout,  * ) '    Parallel routines'
  call print_clock ('reduce')
  call print_clock ('poolreduce')
  WRITE( stdout, * )
#endif
  return
end subroutine print_clock_vdw
