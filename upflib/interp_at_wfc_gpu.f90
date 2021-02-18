!-----------------------------------------------------------------------
SUBROUTINE interp_at_wfc_gpu ( npw, qg_d, nwfcm, ntyp, chiq_d )
  !-----------------------------------------------------------------------
  !
  ! computes chiq: radial fourier transform of atomic orbitals chi
  !
  USE upf_kinds,  ONLY : dp
  USE uspp_param, ONLY : upf
  USE uspp_data,  ONLY : dq, tab_at, tab_at_d
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)  :: npw
  INTEGER, INTENT(IN)  :: nwfcm
  INTEGER, INTENT(IN)  :: ntyp
  REAL(dp), INTENT(IN) :: qg_d(npw)
  REAL(dp), INTENT(OUT):: chiq_d(npw,nwfcm,ntyp)
#if defined(__CUDA)
  attributes(DEVICE) :: qg_d, chiq_d
#endif 
  !
  INTEGER :: nt, nb, ig
  INTEGER :: i0, i1, i2, i3
  REAL(dp):: qgr, px, ux, vx, wx
  !
  !CALL using_tab_at(0)
  !CALL using_tab_at_d(0)
  !
  DO nt = 1, ntyp
     DO nb = 1, upf(nt)%nwfc
        IF ( upf(nt)%oc(nb) >= 0.d0 ) THEN
           !
           !$cuf kernel do (1) <<<*,*>>>
           DO ig = 1, npw
              qgr = qg_d(ig)
              px = qgr / dq - DBLE(INT(qgr/dq))
              ux = 1.d0 - px
              vx = 2.d0 - px
              wx = 3.d0 - px
              i0 = INT(qgr/dq) + 1
              i1 = i0 + 1
              i2 = i0 + 2
              i3 = i0 + 3
              chiq_d(ig,nb,nt) = &
                     tab_at_d(i0,nb,nt) * ux * vx * wx / 6.d0 + &
                     tab_at_d(i1,nb,nt) * px * vx * wx / 2.d0 - &
                     tab_at_d(i2,nb,nt) * px * ux * wx / 2.d0 + &
                     tab_at_d(i3,nb,nt) * px * ux * vx / 6.d0
           END DO
           !
        END IF
     END DO
  END DO

END SUBROUTINE interp_at_wfc_gpu
