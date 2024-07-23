!
!Giovanni Cistaro
!
SUBROUTINE compute_dmn()
!For all symmetry operations, we compute how wannier function w_{n0} 
!transforms under them. With this we can construct 
!irreducible representations as well as point group of the wannier
!functions. 
!for details about dmn, check for example eq. 30 in 
!   Pizzi, JPCM 32 (2020) 165902
!The matrix we compute here is the Fourier transform at R=0 of that eq., i.e. 
!the sum over k.
!   d_{mn}(g) = \int dr w^*_{m, Rk}(r) w_{n, k}(R^{-1).r-f) = sum_k d_{mnk}(g)
!
  USE kinds,                 ONLY : DP
  USE constants,             ONLY : tpi
  USE klist,                 ONLY : nkstot, nks, xk, igk_k, ngk
  USE wvfct,                 ONLY : npwx, current_k
  USE buffers,               ONLY : get_buffer
  USE io_global,             ONLY : stdout
  USE lsda_mod,              ONLY : lsda, isk, nspin, current_spin
  USE symm_base,             ONLY : nsym, ft
  USE control_kcw,           ONLY : spin_component, iuwfc_wann, num_wann
  USE cell_base,             ONLY : omega, at
  USE fft_base,              ONLY : dffts
  USE fft_interfaces,        ONLY : fwfft, invfft
  !
  !
  IMPLICIT NONE 
  !
  COMPLEX(DP) :: imag = (0.D0,1.D0)
  INTEGER     :: isym              
  ! for loop over symmetry operations
  INTEGER     :: ik                
  ! for loop over k points
  INTEGER     :: iwann1, iwann2    
  ! for loop over wannier wf
  INTEGER     :: lrwfc             
  ! length of wf buffer
  INTEGER     :: ik_eff            
  ! right index in buffer for wannier wf
  !   
  INTEGER     :: iRk               
  ! index of rotated k point
  REAL(DP)    :: Gvector(3)        
  ! G vector to bring to FBZ the rotated k
  INTEGER    :: npw_k
  !number of plane waves at k
  INTEGER    :: npw_Rk
  !number of plane waves at Rk
  !
  COMPLEX(DP), ALLOCATABLE :: evc0_k_G(:, :) 
  !wf at k point ik
  COMPLEX(DP), ALLOCATABLE :: evc0_Rk_G(:, :) 
  !wf at rotated k point
  COMPLEX(DP), ALLOCATABLE :: evc_k_G_iwann1(:)
  !wf at k, for a fixed wannier - G space
  COMPLEX(DP), ALLOCATABLE :: evc_Rk_G_iwann2(:)
  !wf at rotated k, for a fixed wannier - G space
  COMPLEX(DP), ALLOCATABLE :: evc_k_r_iwann1(:)
  !wf at k, for a fixed wannier - r space
  COMPLEX(DP), ALLOCATABLE :: g_evc_k_r_iwann1(:)
  !wf at k, for a fixed wannier - r space rotated
  COMPLEX(DP), ALLOCATABLE :: evc_Rk_r_iwann2(:)
  !wf at rotated k, for a fixed wannier - r space
  COMPLEX(DP), ALLOCATABLE :: dmatrix_r(:, :, :, :)
  !dmatrix not integrated over space
  COMPLEX(DP), ALLOCATABLE :: dmatrix( :, :, : )
  ! It is the integral over space of dmatrix_r. It is the dmatrix
  ! we are looking for
  COMPLEX(DP), ALLOCATABLE :: aux(:)
  ! it will contain the dmatrix_r for fixed indices, 
  ! then Fourier tranformed
  COMPLEX(DP)              :: phase
  ! to compute exp(-i k.f)
  COMPLEX(DP), ALLOCATABLE :: phaseG(:)
  ! to compute exp(-i G r)
  REAL(DP)                 :: xk_cryst(3,1)
  !k vector in crystal coordinates
  !
  ALLOCATE( evc0_k_G(npwx, num_wann) ) 
  ALLOCATE( evc0_Rk_G(npwx, num_wann) ) 
  ALLOCATE( evc_k_G_iwann1(npwx) ) 
  ALLOCATE( evc_Rk_G_iwann2(npwx) ) 
  ALLOCATE( evc_k_r_iwann1(dffts%nnr) ) 
  ALLOCATE( g_evc_k_r_iwann1(dffts%nnr) )   
  ALLOCATE( evc_Rk_r_iwann2(dffts%nnr) ) 
  ALLOCATE( dmatrix_r(num_wann, num_wann, nsym, dffts%nnr) ) 
  ALLOCATE( phaseG(dffts%nnr))
  !
  !this must be a global variable in kcw 
  !Maybe better to allocate somewhere else?
  ALLOCATE( dmatrix( num_wann, num_wann, nsym ) ) 
  ALLOCATE( aux(dffts%nnr) )
  !
  !initialize dmatrix with 0, so we can accumulate sums over k.
  dmatrix = 0
  !
  ! construct rir
  CALL kcw_set_symm( dffts%nr1,  dffts%nr2,  dffts%nr3, &
  dffts%nr1x, dffts%nr2x, dffts%nr3x )  
  !
  !
  DO isym = 1, nsym
    DO ik = 1, nks
      current_k = ik
      IF ( lsda ) current_spin = isk(ik)
      lrwfc = num_wann * npwx 
      IF ( lsda .AND. isk(ik) /= spin_component) CYCLE
      ik_eff = ik-(spin_component-1)*nkstot/nspin
      WRITE(stdout, * ) "isym = ", isym, "ik = ", ik, "ik_eff = ", ik_eff
      !
      ! get wf at current k point from buffer
      !
      CALL get_buffer ( evc0_k_G, lrwfc, iuwfc_wann, ik_eff )
      !
      ! rotate k with symmetry isym
      !
      CALL rotate_xk(ik_eff, isym, iRk, Gvector)
      !iRk = ik
      WRITE(*,*) ik, ik_eff, isym, iRk, Gvector
      !
      ! get wf at rotated k point from buffer
      !
      CALL get_buffer ( evc0_Rk_G, lrwfc, iuwfc_wann, iRk )
      !
      DO iwann1 = 1, num_wann
        !
        ! go from G to r for wf at k
        !
        npw_k = ngk(ik)
        evc_k_G_iwann1(:) =  evc0_k_G(:,iwann1)
        evc_k_r_iwann1(:) = CMPLX(0.D0,0.D0,kind=DP)
        CALL invfft_wave (npw_k, igk_k (1,ik_eff), evc_k_G_iwann1 , evc_k_r_iwann1 )
        !
        ! rotate space of wf at k with symmetry isym. what we need is indeed 
        !                       u_{nk}(R^{-1).r-f)
        !
        CALL rotate_evc(isym, evc_k_r_iwann1, g_evc_k_r_iwann1)
        DO iwann2 = 1, num_wann
          !
          ! go from G to r for wf at rotated k
          !
          npw_Rk = ngk(iRk)
          evc_Rk_G_iwann2(:) =  evc0_Rk_G(:,iwann2)
          evc_Rk_r_iwann2(:) = CMPLX(0.D0,0.D0,kind=DP)
          CALL invfft_wave (npw_Rk, igk_k (1,iRk), evc_Rk_G_iwann2 , evc_Rk_r_iwann2 )
          !apply phase exp(-iGr)
          CALL calculate_phase (Gvector, phaseG)      
          evc_Rk_r_iwann2(:) = phaseG(:)*evc_Rk_r_iwann2(:)
          !
          !     Compute
          ! u^*_{i1, Rk}(r) * u{i2,k}(R-1r-f) * exp(-ikf)
          !
          xk_cryst(:,1)=xk(:,ik_eff)
          CALL cryst_to_cart(1,xk_cryst,at,-1)
          phase = EXP(-imag*tpi*dot_product(xk_cryst(:,1), ft(:, isym))) 
          dmatrix_r(iwann1, iwann2, isym, :) = conjg(evc_Rk_r_iwann2(:))*evc_k_r_iwann1(:)*phase  
          !
          ! We want to integrate dmatrix_r in r space. This is equivalent to go to G space
          ! and get 0-th component (ig=1)
          !
          aux(:) = dmatrix_r(iwann1, iwann2, isym, :)!/omega
          CALL fwfft ('Rho', aux, dffts)
          dmatrix( iwann1, iwann2, isym ) = dmatrix(iwann1, iwann2, isym) + aux(dffts%nl(1))
        END DO!iwann2
      END DO!iwann1 
    END DO!ik
  END DO!isym
  !
  ! normalize with respect to number of k points
  dmatrix = dmatrix/(nkstot/nspin)
  !
  !Verify the dmatrix is the correct one
  !
  open (unit = 637, file = "dmn_matrix.txt")
  DO isym=1, nsym
    WRITE(637, *) isym
    DO iwann1 = 1, num_wann
      DO iwann2 = 1, num_wann
        WRITE(637, '(1F10.6)', advance='no') REAL(dmatrix(iwann1, iwann2, isym))
        WRITE(637, '(1F10.6)', advance='no') AIMAG(dmatrix(iwann1, iwann2, isym))
      END DO
      WRITE(637,*) 
    END DO
  END DO
  CLOSE(unit=637)
  STOP
  !
  CALL verify_dmn(dmatrix)
  !
  DEALLOCATE( evc0_k_G ) 
  DEALLOCATE( evc0_Rk_G ) 
  DEALLOCATE( evc_k_G_iwann1 ) 
  DEALLOCATE( evc_Rk_G_iwann2 ) 
  DEALLOCATE( evc_k_r_iwann1 ) 
  DEALLOCATE( g_evc_k_r_iwann1 )   
  DEALLOCATE( evc_Rk_r_iwann2 ) 
  DEALLOCATE( dmatrix_r ) 
  DEALLOCATE( dmatrix ) 
  DEALLOCATE( aux )
  !
  STOP
!
!
END SUBROUTINE


