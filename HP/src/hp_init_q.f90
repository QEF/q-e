!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE hp_init_q()
  !----------------------------------------------------------------------------
  !
  ! This subroutine prepares several variables which are needed in the
  ! HP program at fixed q point:
  !
  ! 1) Compute the phase factor (for the USPP case)
  ! 2) Compute becp = <vkb|evc> for every k point
  ! 3) Calculate and write the S\phi for k and k+q
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : tpiba2
  USE ions_base,            ONLY : nat, tau
  USE becmod,               ONLY : calbec
  USE constants,            ONLY : eps8, tpi
  USE klist,                ONLY : xk, ngk, igk_k
  USE io_global,            ONLY : stdout
  USE buffers,              ONLY : get_buffer
  USE wavefunctions,        ONLY : evc
  USE noncollin_module,     ONLY : noncolin, npol
  USE uspp,                 ONLY : vkb, okvan
  USE eqv,                  ONLY : evq
  USE lrus,                 ONLY : becp1
  USE control_lr,           ONLY : lgamma
  USE units_lr,             ONLY : lrwfc, iuwfc
  USE qpoint,               ONLY : xq, nksq, eigqts, ikks, ikqs
  !
  IMPLICIT NONE
  !
  ! Local variables
  !
  INTEGER :: ik, ikk, ikq, ipol, na
    ! generic counter 
    ! counter on k points
    ! counter on k+q points
    ! counter on polarizations
    ! counter on atoms
  INTEGER :: npw
    ! number of plane waves at k
  REAL(DP) :: arg
    ! the argument of the phase
  !
  CALL start_clock( 'hp_init_q' )
  !
  ! 1) USPP: Compute the phase factor exp(-i q*\tau) 
  !
  IF (okvan) THEN
     DO na = 1, nat
        arg = ( xq(1) * tau(1,na) + &
                xq(2) * tau(2,na) + &
                xq(3) * tau(3,na) ) * tpi
        eigqts(na) = CMPLX( COS( arg ), - SIN( arg ) ,kind=DP)
     ENDDO
  ENDIF
  !
  DO ik = 1, nksq
     !
     ikk = ikks(ik)
     ikq = ikqs(ik)
     npw = ngk(ikk)
     !
     IF ( .NOT. lgamma ) THEN
        !
        IF ( ABS( xq(1) - ( xk(1,ikq) - xk(1,ikk) ) ) > eps8 .OR. &
             ABS( xq(2) - ( xk(2,ikq) - xk(2,ikk) ) ) > eps8 .OR. &
             ABS( xq(3) - ( xk(3,ikq) - xk(3,ikk) ) ) > eps8 ) THEN
           WRITE( stdout,'(/,5x,"k points #",i6," and ", &
                  & i6,5x," total number ",i6)') ikk, ikq, nksq
           WRITE( stdout, '(  5x,"Expected q ",3f10.7)')(xq(ipol), ipol=1,3)
           WRITE( stdout, '(  5x,"Found      ",3f10.7)')((xk(ipol,ikq) &
                                                -xk(ipol,ikk)), ipol = 1, 3)
           CALL errore( 'hp_init_q', 'wrong order of k points', 1 )
        ENDIF
        !
     ENDIF
     !
     ! 2) USPP: Compute the becp terms which are used in the rest of the code
     !
     IF (okvan) THEN
        !
        ! Compute the beta function vkb(k+G)
        ! 
        CALL init_us_2 (npw, igk_k(1,ikk), xk(1,ikk), vkb)
        !
        ! Read the wavefunctions evc at k
        !
        CALL get_buffer (evc, lrwfc, iuwfc, ikk)
        !
        ! becp1 = <vkb|evc>
        !
        CALL calbec (npw, vkb, evc, becp1(ik))
        !
     ENDIF
     !
  ENDDO
  !
  ! 3) Calculate and write to file S\phi for k and k+q
  !
  CALL lr_orthoUwfc (.FALSE.)
  !
  CALL stop_clock ( 'hp_init_q' )
  !
  RETURN
  !
END SUBROUTINE hp_init_q
