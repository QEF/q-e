SUBROUTINE compute_w0n_unitcell(w0)
  ! this function calculates the wannier function
  ! in the unit cell (just because we have everything)
  ! Be careful, the wannier function is not periodic
  ! if we don't consider the supercell! so this is just 
  ! for checks of dmn
  USE kinds,                 ONLY : DP
  USE klist,                 ONLY : nks, xk, igk_k, ngk
  USE lsda_mod,              ONLY : lsda, isk, nspin, current_spin
  USE buffers,               ONLY : get_buffer
  USE constants,             ONLY : tpi
  USE control_kcw,           ONLY : rir, num_wann, evc0
  USE control_kcw,           ONLY : r, ir_end
  USE control_kcw,           ONLY : iuwfc_wann, spin_component
  USE fft_base,              ONLY : dffts
  USE wvfct,                 ONLY : npwx, current_k
  USE klist,                 ONLY:  xk, nkstot, nks
  !
  IMPLICIT NONE 
  !
  INTEGER                            :: iwann, ik !for loops
  INTEGER                            :: ik_eff, ir ! for loops
  INTEGER                            :: lrwfc, npw_k
  COMPLEX(DP)                        :: imag = (0.D0,1.D0)
  !TOBEREMOVED
  REAL(DP) :: gvect(3)
  ! ... the G vector in input 
  !
  !TOBEREMOVED
  COMPLEX(DP), ALLOCATABLE :: phase(:)
  !
  COMPLEX(DP)                        :: w0(dffts%nnr, num_wann)
  !wannier function on the unit cell
  !
  COMPLEX(DP), ALLOCATABLE :: evc_k_G_iwann(:)
  !wf at k, for a fixed wannier - G space
  COMPLEX(DP), ALLOCATABLE :: evc_k_r_iwann(:)
  !wf at k, for a fixed wannier - r space
  !
  ALLOCATE( evc_k_G_iwann(npwx) )
  ALLOCATE( evc_k_r_iwann(dffts%nnr) )
  !TOBEREMOVED
  ALLOCATE(phase(dffts%nnr))
  !
  w0 = 0
  !
  !TODO: Change this with a factorized function
  ! that only sets up ir_end and r
  !CALL calculate_r()
  CALL calculate_phase (gvect, phase)

  DO ik = 1, nks
    current_k = ik
    IF ( lsda ) current_spin = isk(ik)
    lrwfc = num_wann * npwx 
    IF ( lsda .AND. isk(ik) /= spin_component) CYCLE
    ik_eff = ik-(spin_component-1)*nkstot/nspin
    !
    ! get wf at current k point from buffer
    !
    CALL get_buffer ( evc0, lrwfc, iuwfc_wann, ik_eff )
    WRITE(*,*) "MAXVAL(evc0: )", MAXVAL(ABS(evc0))
    DO iwann = 1, num_wann
      !
      ! go from G to r for wf at k
      !
      npw_k = ngk(ik)
      evc_k_G_iwann(:) =  evc0(:,iwann)
      WRITE(*,*) "MAXVAL(evc_k_G_iwann: )", MAXVAL(ABS(evc_k_G_iwann))
      evc_k_r_iwann(:) = CMPLX(0.D0,0.D0,kind=DP)
      CALL invfft_wave (npw_k, igk_k (1,ik), evc_k_G_iwann , evc_k_r_iwann )
      WRITE(*,*) "MAXVAL(evc_k_r_iwann: )", MAXVAL(ABS(evc_k_r_iwann))
      DO ir = 1, ir_end
        !
        ! w0 = \sum_k Exp(i*k.r)*u_k(r)
        !
        w0( ir, iwann ) = w0( ir, iwann ) &
            + EXP(+IMAG*dot_product( xk(:, ik_eff), r(ir, :) ))*evc_k_r_iwann(ir) 
      END DO!ir
      WRITE(*,*) "MAXVAL w0", iwann, MAXVAL(ABS(w0(:, iwann)))
    END DO!iwann
    !WRITE(*,*) "iwann1 = ", iwann1, "SUM: ", SUM( ABS(gw0-gw0_) ) 

  END DO!ik 
END SUBROUTINE


SUBROUTINE Verify_dmn(dmatrix)
  USE kinds,                 ONLY : DP
  USE control_kcw,           ONLY : num_wann
  USE fft_base,              ONLY : dffts
  USE symm_base,             ONLY : nsym
  !
  IMPLICIT NONE 
  !
  INTEGER                            :: isym, iwann1, iwann2
  COMPLEX(DP), ALLOCATABLE           :: w0(:, :)
  COMPLEX(DP), ALLOCATABLE           :: gw0(:)
  COMPLEX(DP), ALLOCATABLE           :: gw0_(:)
  COMPLEX(DP)                        :: dmatrix( num_wann, num_wann, nsym ) 
  !wannier function on the unit cell
  !
  ALLOCATE(gw0(dffts%nnr))
  ALLOCATE(gw0_(dffts%nnr))
  ALLOCATE( w0(dffts%nnr, num_wann) )
  !
  CALL compute_w0n_unitcell(w0)
  WRITE(*,*) "MAXVAL(w0)", MAXVAL(ABS(w0))
  !
  DO isym = 1, nsym
    DO iwann1 = 1, num_wann
      !
      ! gw0_ is obtained rotating the space of the wannier function
      ! iwann1
      !
      gw0_(:) = w0(iwann1, :)
      CALL rotate_evc(isym, gw0_, gw0_)
      !
      ! gw0 on iwann1 using representation matrix
      !
      gw0 = 0
      DO iwann2 = 1, num_wann
        gw0(:) = gw0(:) + dmatrix(iwann2, iwann1, isym)*w0(:, iwann2)
      END DO !iwann2
      WRITE(*,*) "max(iwann1)", MAXVAL(ABS(gw0(:))), MAXVAL(ABS(gw0_(:)))
      WRITE(*,*) "iwann1 = ", iwann1, "SUM: ", SUM( ABS(gw0-gw0_) ) 
    END DO !iwann1
  END DO!isym
  !          
  !
  DEALLOCATE(gw0)
  DEALLOCATE(gw0_)
  DEALLOCATE(w0)
  !
END SUBROUTINE
