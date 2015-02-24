!
! Copyright (C) 2004-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
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
  USE mp_global,         ONLY : mp_startup, mp_global_end
  USE environment,       ONLY : environment_start
  USE ld1inc,            ONLY : iswitch, write_coulomb, grid, lgipaw_reconstruction
  USE radial_grids,      ONLY : deallocate_radial_grid
  USE command_line_options, ONLY: input_file_
  !
  implicit none
  CHARACTER (LEN=9) :: code = 'LD1'
  !
  !   write initialization information
  !
  call mp_startup( )
  call environment_start ( code )
  !
  !    read input, possible pseudopotential and set the main variables
  !
  call ld1_readin (input_file_)
  call ld1_setup ( )
  !
  !   four possible working mode:
  !
  if (iswitch.eq.1) then
     !
     !   all-electron calculation
     !
     call all_electron(.true.,1)
     if ( write_coulomb ) call write_ae_pseudo ( )
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
     call all_electron(.false.,1)
     call gener_pseudo ( )
     !if(.not. lgipaw_reconstruction) 
     call run_test ( )
     call ld1_writeout ( )
     !
  elseif (iswitch.eq.4) then
     !
     ! LDA-1/2 correction to the input pseudopotential 
     !
     call run_lda_half ( )
     call ld1_writeout ( )
     !
  else
     call errore('ld1','iswitch not implemented',1)
  endif
  call deallocate_radial_grid( grid )

  call mp_global_end()

end program ld1
!
