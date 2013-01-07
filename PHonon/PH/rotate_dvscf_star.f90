!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine rotate_dvscf_star(iq_)
  !-----------------------------------------------------------------------
  !
  !
  ! Given dvscf or drho, this routine obtains dvscf or drho
  ! over the star{q} and writes it on a file. 
  !

  USE kinds,         ONLY : DP
  USE cell_base,     ONLY : at, bg
  USE ions_base,     ONLY : ntyp => nsp, ityp
  USE symm_base,     ONLY : s, sr, irt, nsym, time_reversal, invs
  USE qpoint,        ONLY : xq
  USE output,        ONLY : fildrho, fildvscf
  USE dfile_star,    ONLY : write_dfile_star, drho_star, dvscf_star !write_dfile_mq
  USE units_ph,      ONLY : iudrho, iudvscf
  USE modes,         ONLY : u

  INTEGER :: nq,  isq (48), imq, iq_
  LOGICAL :: opnd
  REAL(DP) :: sxq (3, 48)

  if(.not.drho_star%open.and..not.dvscf_star%open) return

  call start_clock('rotate_dvscf_star')
  !
  !   Generates the star of q
  !
  call star_q (xq, at, bg, nsym, s, invs, nq, sxq, isq, imq, .TRUE. )
  !
  ! Rotates and write drho_q* (to be improved)
  IF(drho_star%open) THEN
     INQUIRE (UNIT = iudrho, OPENED = opnd)
     IF (opnd) CLOSE(UNIT = iudrho, STATUS='keep')
     CALL write_dfile_star(drho_star, fildrho, nsym, xq, u, nq, sxq, isq, &
          s, sr, invs, irt, ntyp, ityp,(imq==0), -1 )
  ENDIF
  IF(dvscf_star%open) THEN
     INQUIRE (UNIT = iudvscf, OPENED = opnd)
     IF (opnd) CLOSE(UNIT = iudvscf, STATUS='keep')
     CALL write_dfile_star(dvscf_star, fildvscf, nsym, xq, u, nq, sxq, isq, &
          s, sr, invs, irt, ntyp, ityp,(imq==0), iq_ )
  ENDIF
  !
  call stop_clock('rotate_dvscf_star')
  !
end subroutine rotate_dvscf_star
