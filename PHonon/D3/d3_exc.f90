!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE d3_exc
  !-----------------------------------------------------------------------
  !
  !    Calculates the contribution to the derivative of the dynamical
  !    matrix due to the third derivative of the exchange and correlation
  !    energy
  !
  USE ions_base,  ONLY : nat
  USE kinds,      ONLY : DP
  USE pwcom
  USE scf, only : rho, rho_core
  USE fft_base, only : dfftp
  USE phcom
  USE d3com
  USE io_global, ONLY : ionode_id
  USE mp_global, ONLY : inter_pool_comm, my_pool_id, &
                        npool, intra_pool_comm
  USE mp,        ONLY : mp_bcast, mp_sum

  IMPLICIT NONE

  INTEGER :: errcode, ir, ipert, jpert, kpert, npert1, npert2
  REAL (DP) :: d2mxc, rhotot, xq0 (3)
  REAL (DP), ALLOCATABLE :: d2muxc (:)
  COMPLEX (DP) :: aux
  COMPLEX (DP), ALLOCATABLE :: work1 (:), work2 (:), &
                                      work3 (:), d3dyn1 (:,:,:)

  ALLOCATE (d2muxc( dfftp%nnr))
  ALLOCATE (work1 ( dfftp%nnr))
  ALLOCATE (work2 ( dfftp%nnr))
  ALLOCATE (work3 ( dfftp%nnr))
  ALLOCATE (d3dyn1( 3*nat, 3*nat, 3*nat))

!  IF ( my_pool_id == 0 ) THEN
     !
     ! Calculates third derivative of Exc
     !
     d2muxc(:) = 0.d0
     DO ir = 1, dfftp%nnr
        rhotot = rho%of_r (ir, 1) + rho_core (ir)
        IF (rhotot > 1.d-30) d2muxc (ir) = d2mxc (rhotot)
        IF (rhotot < - 1.d-30) d2muxc (ir) = - d2mxc ( - rhotot)
     ENDDO
     !
     ! Calculates the contribution to d3dyn
     !
     d3dyn1 (:,:,:) = (0.d0, 0.d0)
     DO ipert = 1, 3 * nat
        IF (q0mode (ipert) ) THEN
           CALL davcio_drho (work1, lrdrho, iud0rho, ipert, - 1)
           DO jpert = 1, 3 * nat
              CALL davcio_drho (work2, lrdrho, iudrho, jpert, - 1)
              DO kpert = 1, 3 * nat
                 CALL davcio_drho (work3, lrdrho, iudrho, kpert, - 1)
                 aux = CMPLX(0.d0, 0.d0,kind=DP)
                 DO ir = 1, dfftp%nnr
                    aux = aux + &
                          d2muxc (ir) * work1 (ir) * &
                          CONJG (work2 (ir) ) * work3 (ir)
                 ENDDO
                 !
                 CALL mp_sum ( aux, intra_pool_comm )
                 !
                 d3dyn1 (ipert, jpert, kpert) = omega * aux / (dfftp%nr1 * dfftp%nr2 * dfftp%nr3)
                 !
              ENDDO
           ENDDO
        ENDIF

     ENDDO
     !
!  END IF
  !
  IF ( npool /= 1 ) CALL mp_bcast( d3dyn1, ionode_id, inter_pool_comm )
  !
  d3dyn = d3dyn  + d3dyn1
  d3dyn_aux9 = d3dyn1
  !
  DEALLOCATE (d2muxc)
  DEALLOCATE (work1)
  DEALLOCATE (work2)
  DEALLOCATE (work3)
  DEALLOCATE (d3dyn1)
  !
  RETURN
  !
END SUBROUTINE d3_exc
