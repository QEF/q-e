!-----------------------------------------------------------------------
SUBROUTINE interp_atwfc_gpu ( npw, qg, nwfcm, chiq )
  !-----------------------------------------------------------------------
  !
  ! computes chiq: radial fourier transform of atomic orbitals chi
  !
  USE upf_kinds,  ONLY : dp
  USE uspp_param, ONLY : upf, nsp
  USE uspp_data,  ONLY : dq, tab_at
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)  :: npw
  INTEGER, INTENT(IN)  :: nwfcm
  REAL(dp), INTENT(IN) :: qg(npw)
  REAL(dp), INTENT(OUT):: chiq(npw,nwfcm,nsp)
!$acc declare deviceptr(qg,chiq)
  !
  INTEGER :: nt, nb, ig
  INTEGER :: i0, i1, i2, i3
  REAL(dp):: qgr, px, ux, vx, wx
  !
  !$acc data present_or_copyin(qg) present_or_copyout(chiq) present(tab_at)
  DO nt = 1, nsp
     DO nb = 1, upf(nt)%nwfc
        IF ( upf(nt)%oc(nb) >= 0.d0 ) THEN
           !$acc parallel loop
           DO ig = 1, npw
              qgr = qg(ig)
              px = qgr / dq - DBLE(INT(qgr/dq))
              ux = 1.d0 - px
              vx = 2.d0 - px
              wx = 3.d0 - px
              i0 = INT(qgr/dq) + 1
              i1 = i0 + 1
              i2 = i0 + 2
              i3 = i0 + 3
              chiq(ig,nb,nt) = &
                     tab_at(i0,nb,nt) * ux * vx * wx / 6.d0 + &
                     tab_at(i1,nb,nt) * px * vx * wx / 2.d0 - &
                     tab_at(i2,nb,nt) * px * ux * wx / 2.d0 + &
                     tab_at(i3,nb,nt) * px * ux * vx / 6.d0
           END DO
           !
        END IF
     END DO
  END DO
  !$acc end data
END SUBROUTINE interp_atwfc_gpu
