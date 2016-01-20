!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine raman
  !-----------------------------------------------------------------------
  !
  USE kinds,    ONLY : DP
  USE klist,    ONLY : lgauss
  USE lsda_mod, ONLY : lsda
  USE control_flags, ONLY : gamma_only
  USE uspp, ONLY: okvan
  USE control_ph, ONLY : epsil, convt, rec_code_read
  USE ph_restart, ONLY : ph_writefile
  USE ramanm, ONLY: lraman, elop, done_lraman, done_elop

  USE control_lr, ONLY : lgamma

  implicit none

  INTEGER :: ierr

  if (okvan) &
      call errore ('raman','Ultrasoft pseudopotentials not implemented',1)
  if (lsda) call errore ('raman',' spin-polarized case not implemented',1)
  if (lgauss .or..not.lgamma) &
      call errore ('raman','called in the wrong case',1)
  if (epsil.and..not.convt) &
      call errore ('raman','epsil calcul. not converged',1)
  !
  ! Computes Pc [DH,Drho] |psi>
  !
  IF (rec_code_read == -10) THEN
     ! restart from a previous calculation
     write (6,'(/,5x,''Skipping computation of Pc [DH,Drho] |psi> '')')
  ELSE
     write (6,'(/,5x,''Computing Pc [DH,Drho] |psi> '')')
     call dhdrhopsi ( ) 
  END IF
  !
  ! Computes the electro-optic tensor
  !
  IF (elop.AND..NOT.done_elop) THEN
     call el_opt()
  ELSEIF (done_elop) THEN
     CALL summarize_elopt()
  ENDIF
     
  if (.not.lraman) return
  write (6,'(/,5x,''Computing Second order response '')')
  !
  ! Computes the potential that remains unchanged in the scf-cycle
  !
  call dvpsi_e2 ( ) 
  !
  ! Self-consistent cycle to compute the second derivative of the charge
  !
  call solve_e2 ( )
  !
  ! Computes and writes the Raman tensor
  !
  call raman_mat ( )

  done_lraman=.TRUE.
  CALL ph_writefile('tensors',0, 0, ierr)
  return
end subroutine raman
