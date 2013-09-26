!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE drho_drc (iudrho_x, u_x, xq_x, drc_x, scalef)
  !-----------------------------------------------------------------------
  !  Reads the variation of the charge saved on a file and changes
  !  it according to the variation of the core_charge
  !  It is used by drho_cc. Have a look there for more explanation
  !
  USE ions_base,  ONLY : nat, ityp, ntyp => nsp, tau
  USE kinds,      ONLY : DP
  USE fft_base,   ONLY : dfftp
  USE fft_interfaces, ONLY : invfft
  USE pwcom
  USE phcom
  USE d3com
  USE uspp_param, ONLY : upf
  USE mp,         ONLY : mp_barrier
  USE mp_world,   ONLY : world_comm

  IMPLICIT NONE

  INTEGER :: iudrho_x
  !input: the unit containing the charge variation
  REAL (DP) :: xq_x (3), scalef
  !input: q point
  !input: drhocore will be added to the valence charge scaled by this factor
  COMPLEX (DP) :: u_x (3 * nat, 3 * nat), drc_x (ngm, ntyp)
  !input: the transformation modes patterns
  !input: contain the rhoc (without structu

  INTEGER :: ipert, na, mu, nt, ig, errcode, iudrho_tmp
  REAL (DP) :: gtau
  COMPLEX (DP) :: guexp
  COMPLEX (DP), ALLOCATABLE :: drhoc (:), drhov (:), uact (:)

  iudrho_tmp = iudrho_x
  iudrho_x = iudrho_tmp+1000 ! this must be already opened, see openfild3

  ALLOCATE  (drhoc( dfftp%nnr))
  ALLOCATE  (drhov( dfftp%nnr))
  ALLOCATE  (uact( 3 * nat))

  DO ipert = 1, 3 * nat
     drhoc(:) = (0.d0, 0.d0)
     uact(:) = u_x (:, ipert)
     DO na = 1, nat
        mu = 3 * (na - 1)
        IF (ABS (uact (mu + 1) ) + ABS (uact (mu + 2) ) + &
             ABS (uact (mu + 3) ) > 1.0d-12) THEN
           nt = ityp (na)
           IF (upf(nt)%nlcc) THEN
              DO ig = 1, ngm
                 gtau = tpi * ( (g (1, ig) + xq_x (1) ) * tau (1, na) &
                              + (g (2, ig) + xq_x (2) ) * tau (2, na) &
                              + (g (3, ig) + xq_x (3) ) * tau (3, na) )
                 guexp = tpiba * ( (g (1, ig) + xq_x (1) ) * uact (mu + 1) &
                                 + (g (2, ig) + xq_x (2) ) * uact (mu + 2) &
                                 + (g (3, ig) + xq_x (3) ) * uact (mu + 3) )&
                               * CMPLX(0.d0, - 1.d0,kind=DP) &
                               * CMPLX(COS (gtau), - SIN (gtau) ,kind=DP)
                 drhoc (nl (ig) ) = drhoc (nl (ig) ) + drc_x (ig, nt) * guexp
              ENDDO
           ENDIF
        ENDIF
     ENDDO

     CALL invfft ('Dense', drhoc, dfftp)
     CALL davcio_drho2 (drhov, lrdrho, iudrho_tmp, ipert, - 1)
     drhov(:) = drhov(:) + scalef * drhoc(:)
     CALL davcio_drho2 (drhov, lrdrho, iudrho_x, ipert, + 1)
  ENDDO

  CALL mp_barrier( world_comm )

  DEALLOCATE (drhoc)
  DEALLOCATE (drhov)
  DEALLOCATE (uact)
  RETURN
END SUBROUTINE drho_drc
