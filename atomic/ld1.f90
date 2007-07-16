!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
program ld1
  !---------------------------------------------------------------
  !
  !     atomic self-consistent local-density program
  !     atomic rydberg units are used : e^2=2, m=1/2, hbar=1
  !     psi(r) = rR(r), where R(r) is the radial part of the wfct
  !     rho(r) = psi(r)^2 => rho(r) = (true charge density)*(4\pi r^2)
  !                       The same applies to the core charge
  !---------------------------------------------------------------
  !
  USE global_version,    ONLY : version_number
  USE io_files,          ONLY : nd_nmbr
  USE mp,                ONLY : mp_barrier, mp_end
  USE ld1inc
  !
  implicit none
  character :: day*9, hour*9
  CHARACTER (LEN=9) :: code = 'LD1'
  !
  !   write initialization information
  !
  call startup( nd_nmbr, code, version_number )
  !
  !    read input, possible pseudopotential and set the main variables
  !
  call ld1_readin
  call ld1_setup
  !
  !   three possible working mode:
  !
  if (iswitch.eq.1) then
     !
     !   all-electron calculation
     !
     call all_electron(.true.)
     !
  elseif (iswitch.eq.2) then
     !
     !   pseudopotential test
     !
     call run_test ( )
     call ld1_writeout ( )
     !
  elseif (iswitch.eq.3) then
     !
     !  pseudopotential generation and test
     !
     call all_electron(.false.)
     call gener_pseudo ( )
     call run_test ( )
     call ld1_writeout ( )
     !
  else
     call errore('ld1','iswitch not implemented',1)
  endif
  call mp_barrier()
  call mp_end()

end program ld1
