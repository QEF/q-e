!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!------------------------------------------------------------
SUBROUTINE hp_dealloc_1()
!------------------------------------------------------------
  !
  !  Deallocates some arrays
  !
  USE lr_symm_base,   ONLY : rtau
  USE start_k,        ONLY : xk_start, wk_start
  USE ldaU_hp,        ONLY : Rvect, dnsscf, dns0, dnsscf_tot, dns0_tot, &
                             x_q, comp_iq
  !
  IMPLICIT NONE
  !
  IF (ALLOCATED( x_q ))          DEALLOCATE( x_q )
  IF (ALLOCATED( comp_iq ))      DEALLOCATE( comp_iq )
  IF (ALLOCATED( Rvect ))        DEALLOCATE( Rvect )
  IF (ALLOCATED( dnsscf ))       DEALLOCATE( dnsscf )
  IF (ALLOCATED( dns0 ))         DEALLOCATE( dns0 )
  IF (ALLOCATED( dnsscf_tot ))   DEALLOCATE( dnsscf_tot )
  IF (ALLOCATED( dns0_tot ))     DEALLOCATE( dns0_tot )
  IF (ALLOCATED( rtau ))         DEALLOCATE( rtau )
  !
  ! Note that these two variables are allocated by read_file.
  ! They cannot be deallocated by clean_pw because the starting xk and wk
  ! points must be known at each q point.
  ! The logic of these two variables must be improved.
  !
  IF (ALLOCATED( xk_start ))  DEALLOCATE( xk_start )
  IF (ALLOCATED( wk_start ))  DEALLOCATE( wk_start )
  !
  RETURN
  !
END SUBROUTINE hp_dealloc_1
