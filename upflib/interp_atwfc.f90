!-----------------------------------------------------------------------
SUBROUTINE interp_atwfc ( npw, qg, nwfcm, chiq )
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
  !
  INTEGER :: nt, nb, ig
  INTEGER :: i0, i1, i2, i3
  REAL(dp):: qgr, px, ux, vx, wx
  !
  DO nt = 1, nsp
     DO nb = 1, upf(nt)%nwfc
        IF ( upf(nt)%oc (nb) >= 0.d0) THEN
           DO ig = 1, npw
              qgr = qg(ig)
              px = qgr / dq - int (qgr/dq)
              ux = 1.d0 - px
              vx = 2.d0 - px
              wx = 3.d0 - px
              i0 = INT(qgr/dq) + 1
              i1 = i0 + 1
              i2 = i0 + 2
              i3 = i0 + 3
              chiq (ig, nb, nt) = &
                     tab_at (i0, nb, nt) * ux * vx * wx / 6.d0 + &
                     tab_at (i1, nb, nt) * px * vx * wx / 2.d0 - &
                     tab_at (i2, nb, nt) * px * ux * wx / 2.d0 + &
                     tab_at (i3, nb, nt) * px * ux * vx / 6.d0
           ENDDO
        ENDIF
     ENDDO
  ENDDO

END SUBROUTINE interp_atwfc

!-----------------------------------------------------------------------
SUBROUTINE interp_atdwfc ( npw, qg, nwfcm, dchiq )
  !-----------------------------------------------------------------------
  !
  ! computes dchi/dq
  !
  USE upf_kinds,  ONLY : dp
  USE uspp_param, ONLY : upf, nsp
  USE uspp_data,  ONLY : tab_at, dq
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)  :: npw
  INTEGER, INTENT(IN)  :: nwfcm
  REAL(dp), INTENT(IN) :: qg(npw)
  REAL(dp), INTENT(OUT):: dchiq(npw,nwfcm,nsp)
  !
  INTEGER :: nt, nb, ig
  INTEGER :: i0, i1, i2, i3
  REAL(dp):: px, ux, vx, wx
  !
  DO nt=1,nsp
     DO nb=1,upf(nt)%nwfc
        IF (upf(nt)%oc(nb) >= 0.d0) THEN
           DO ig = 1, npw
              px = qg(ig) / dq - INT(qg(ig)/dq)
              ux = 1.d0 - px
              vx = 2.d0 - px
              wx = 3.d0 - px
              i0 = qg(ig) / dq + 1
              i1 = i0 + 1
              i2 = i0 + 2
              i3 = i0 + 3
              dchiq(ig,nb,nt) = &
                   ( tab_at (i0, nb, nt) * (-vx*wx-ux*wx-ux*vx)/6.d0 + &
                   tab_at (i1, nb, nt) * (+vx*wx-px*wx-px*vx)/2.d0 - &
                   tab_at (i2, nb, nt) * (+ux*wx-px*wx-px*ux)/2.d0 + &
                   tab_at (i3, nb, nt) * (+ux*vx-px*vx-px*ux)/6.d0 )/dq
           ENDDO
        ENDIF
     ENDDO
  ENDDO

END SUBROUTINE interp_atdwfc
