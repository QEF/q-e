
! Copyright (C) 2002-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
#if defined(__CUDA)



! Copyright (C) 2002-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE h_psi_gpu( lda, n, m, psi_d, hpsi_d )
  !----------------------------------------------------------------------------
  !
  ! ... This routine computes the product of the Hamiltonian
  ! ... matrix with m wavefunctions contained in psi
  !
  ! ... input:
  ! ...    lda   leading dimension of arrays psi, spsi, hpsi
  ! ...    n     true dimension of psi, spsi, hpsi
  ! ...    m     number of states psi
  ! ...    psi
  !
  ! ... output:
  ! ...    hpsi  H*psi
  !
  ! --- Wrapper routine: performs bgrp parallelization on non-distributed bands
  ! --- if suitable and required, calls old H\psi routine h_psi_
  !
  USE kinds,            ONLY : DP
  USE noncollin_module, ONLY : npol
  USE funct,            ONLY : exx_is_active
  USE mp_bands,         ONLY : use_bgrp_in_hpsi, inter_bgrp_comm
  USE mp,               ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)      :: lda, n, m
  COMPLEX(DP), DEVICE, INTENT(IN)  :: psi_d(lda*npol,m) 
  COMPLEX(DP), DEVICE, INTENT(OUT) :: hpsi_d(lda*npol,m)   
  !
  INTEGER     :: m_start, m_end
  !
  CALL start_clock( 'h_psi_bgrp' ); !write (*,*) 'start h_psi_bgrp'; FLUSH(6)

  ! band parallelization with non-distributed bands is performed if
  ! 1. enabled (variable use_bgrp_in_hpsi must be set to .T.)
  ! 2. exact exchange is not active (if it is, band parallelization is already
  !    used in exx routines called by Hpsi)
  ! 3. there is more than one band, otherwise there is nothing to parallelize
  !
  IF (use_bgrp_in_hpsi .AND. .NOT. exx_is_active() .AND. m > 1) THEN
     ! use band parallelization here
     hpsi_d(:,:) = (0.d0,0.d0)
     CALL divide(inter_bgrp_comm,m,m_start,m_end)
     ! Check if there at least one band in this band group
     IF (m_end >= m_start) &
        CALL h_psi__gpu( lda, n, m_end-m_start+1, psi_d(1,m_start), hpsi_d(1,m_start) )
     CALL mp_sum(hpsi_d,inter_bgrp_comm)
  ELSE
     ! don't use band parallelization here
     CALL h_psi__gpu( lda, n, m, psi_d, hpsi_d )
  END IF

  CALL stop_clock( 'h_psi_bgrp' )
  RETURN
  !
END SUBROUTINE h_psi_gpu
!
!----------------------------------------------------------------------------
SUBROUTINE h_psi__gpu( lda, n, m, psi_d, hpsi_d )
  !----------------------------------------------------------------------------
  !
  ! ... This routine computes the product of the Hamiltonian
  ! ... matrix with m wavefunctions contained in psi
  !
  ! ... input:
  ! ...    lda   leading dimension of arrays psi, spsi, hpsi
  ! ...    n     true dimension of psi, spsi, hpsi
  ! ...    m     number of states psi
  ! ...    psi
  !
  ! ... output:
  ! ...    hpsi  H*psi
  !
  USE kinds,    ONLY : DP
  USE bp,       ONLY : lelfield,l3dstring,gdir, efield, efield_cry
  USE becmod,   ONLY : bec_type, becp, calbec
  USE lsda_mod, ONLY : current_spin
  USE scf_gpum, ONLY : vrs_d, using_vrs_d
  USE uspp,     ONLY : vkb, nkb
  USE ldaU,     ONLY : lda_plus_u, U_projection
  USE gvect,    ONLY : gstart
  USE funct,    ONLY : dft_is_meta
  USE control_flags,    ONLY : gamma_only
  USE noncollin_module, ONLY: npol, noncolin
  USE realus,   ONLY : real_space, &
                       invfft_orbital_gamma, fwfft_orbital_gamma, calbec_rs_gamma, add_vuspsir_gamma, & 
                       invfft_orbital_k, fwfft_orbital_k, calbec_rs_k, add_vuspsir_k, & 
                       v_loc_psir_inplace
  USE fft_base, ONLY : dffts
  USE exx,      ONLY : use_ace, vexx, vexxace_gamma, vexxace_k
  USE funct,    ONLY : exx_is_active
  USE fft_helper_subroutines
  !
  USE wvfct_gpum,    ONLY : g2kin_d, using_g2kin_d
  USE uspp_gpum,     ONLY : using_vkb
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)     :: lda, n, m
  COMPLEX(DP), DEVICE, INTENT(IN)  :: psi_d(lda*npol,m) 
  COMPLEX(DP), DEVICE, INTENT(OUT) :: hpsi_d(lda*npol,m)   
  !
  COMPLEX(DP), ALLOCATABLE, PINNED :: psi_host(:,:)
  COMPLEX(DP), ALLOCATABLE, PINNED :: hpsi_host(:,:)
  !
  INTEGER     :: ipol, ibnd, incr, i
  REAL(dp)    :: ee
  !
  LOGICAL     :: need_host_copy
  !
  CALL start_clock( 'h_psi' ); !write (*,*) 'start h_psi';FLUSH(6)
  CALL using_g2kin_d(.false.)
  CALL using_vrs_d(.false.)
  CALL using_vkb(.false.)
  !
  hpsi_d (:, 1:m) = (0.0_dp, 0.0_dp)

  need_host_copy = ( real_space .and. nkb > 0  ) .OR. &
                    dft_is_meta() .OR. &
                    (lda_plus_u .AND. U_projection.NE."pseudo" ) .OR. &
                    ( nkb > 0 .AND. .NOT. real_space) .OR. &
                    exx_is_active() .OR. lelfield

  if (need_host_copy) then
      ALLOCATE(psi_host(lda*npol,m) , hpsi_host(lda*npol,m) )
      psi_host = psi_d
      hpsi_host = hpsi_d ! this is not needed
  end if

  CALL start_clock( 'h_psi:pot' ); !write (*,*) 'start h_pot';FLUSH(6)
  !
  ! ... Here the product with the local potential V_loc psi
  !
  IF ( gamma_only ) THEN
     ! 
     IF ( real_space .and. nkb > 0  ) then 
        !
        ! ... real-space algorithm
        ! ... fixme: real_space without beta functions does not make sense
        !
        IF ( dffts%has_task_groups ) then 
           incr = 2 * fftx_ntgrp(dffts)
        ELSE
           incr = 2
        ENDIF
        DO ibnd = 1, m, incr
           ! ... transform psi to real space -> psic 
           CALL invfft_orbital_gamma(psi_host,ibnd,m) 
           ! ... compute becp%r = < beta|psi> from psic in real space
     CALL start_clock( 'h_psi:calbec' ) 
           CALL calbec_rs_gamma(ibnd,m,becp%r) 
     CALL stop_clock( 'h_psi:calbec' )
           ! ... psic -> vrs * psic (psic overwritten will become hpsi)
           CALL v_loc_psir_inplace(ibnd,m) 
           ! ... psic (hpsi) -> psic + vusp
           CALL  add_vuspsir_gamma(ibnd,m)
           ! ... transform psic back in reciprocal space and assign it to hpsi
           CALL fwfft_orbital_gamma(hpsi_host,ibnd,m) 
        END DO
        hpsi_d = hpsi_host
        !
     ELSE
        ! ... usual reciprocal-space algorithm
        CALL vloc_psi_gamma_gpu ( lda, n, m, psi_d, vrs_d(1,current_spin), hpsi_d )
        !
     ENDIF 
     !
  ELSE IF ( noncolin ) THEN 
     !
     CALL vloc_psi_nc_gpu ( lda, n, m, psi_d, vrs_d, hpsi_d )
     !
  ELSE  
     ! 
     IF ( real_space .and. nkb > 0  ) then 
        !
        ! ... real-space algorithm
        ! ... fixme: real_space without beta functions does not make sense
        !
        IF ( dffts%has_task_groups ) then 
           incr = fftx_ntgrp(dffts)
        ELSE
           incr = 1
        ENDIF
        DO ibnd = 1, m
           ! ... transform psi to real space -> psic 
           CALL invfft_orbital_k(psi_host,ibnd,m) 
           ! ... compute becp%r = < beta|psi> from psic in real space
     CALL start_clock( 'h_psi:calbec' )
           CALL calbec_rs_k(ibnd,m) 
     CALL stop_clock( 'h_psi:calbec' )
           ! ... psic -> vrs * psic (psic overwritten will become hpsi)
           CALL v_loc_psir_inplace(ibnd,m) 
           ! ... psic (hpsi) -> psic + vusp
           CALL  add_vuspsir_k(ibnd,m)
           ! ... transform psic back in reciprocal space and assign it to hpsi
           CALL fwfft_orbital_k(hpsi_host,ibnd,m) 
        END DO
        IF (need_host_copy) hpsi_d = hpsi_host
        !
     ELSE
        !
        CALL vloc_psi_k_gpu ( lda, n, m, psi_d, vrs_d(1,current_spin), hpsi_d )
        !
     END IF
     !
  END IF  
  !
  ! ... Here the product with the non local potential V_NL psi
  ! ... (not in the real-space case: it is done together with V_loc)
  !
  IF ( nkb > 0 .AND. .NOT. real_space) THEN
     !
     CALL start_clock( 'h_psi:calbec' )
     CALL calbec ( n, vkb, psi_host, becp, m )
     CALL stop_clock( 'h_psi:calbec' )
     hpsi_host = hpsi_d
     CALL add_vuspsi( lda, n, m, hpsi_host )
     hpsi_d = hpsi_host
     !
  END IF
  !  
  CALL stop_clock( 'h_psi:pot' )
  !
  ! ... Here we add the kinetic energy (k+G)^2 psi
  !
  !$cuf kernel do(2)
  DO ibnd = 1, m
     DO i=1, n
        hpsi_d (i, ibnd) = hpsi_d(i, ibnd) + g2kin_d (i) * psi_d (i, ibnd)
        IF ( noncolin ) THEN
           hpsi_d (lda+i, ibnd) = hpsi_d(lda+i,ibnd) + g2kin_d (i) * psi_d (lda+i, ibnd)
        END IF
     END DO
  END DO
  !
  if (dft_is_meta()) then
     hpsi_host = hpsi_d
     call h_psi_meta (lda, n, m, psi_host, hpsi_host)
     hpsi_d = hpsi_host
  end if
  !
  ! ... Here we add the Hubbard potential times psi
  !
  IF ( lda_plus_u .AND. U_projection.NE."pseudo" ) THEN
     !
     hpsi_host = hpsi_d
     IF (noncolin) THEN
        CALL vhpsi_nc( lda, n, m, psi_host, hpsi_host )
     ELSE
        call vhpsi( lda, n, m, psi_host, hpsi_host )
     ENDIF
     hpsi_d = hpsi_host
     !
  ENDIF
  !
  ! ... Here the exact-exchange term Vxx psi
  !
  IF ( exx_is_active() ) THEN
     hpsi_host = hpsi_d
     IF ( use_ace) THEN
        IF (gamma_only) THEN
           CALL vexxace_gamma(lda,m,psi_host,ee,hpsi_host)
        ELSE
           CALL vexxace_k(lda,m,psi_host,ee,hpsi_host) 
        END IF
     ELSE
        CALL vexx( lda, n, m, psi_host, hpsi_host, becp )
     END IF
     hpsi_d = hpsi_host
  END IF
  !
  ! ... electric enthalpy if required
  !
  IF ( lelfield ) THEN
     !
     hpsi_host = hpsi_d
     IF ( .NOT.l3dstring ) THEN
        CALL h_epsi_her_apply( lda, n, m, psi_host, hpsi_host,gdir, efield )
     ELSE
        DO ipol=1,3
           CALL h_epsi_her_apply( lda, n, m, psi_host, hpsi_host,ipol,efield_cry(ipol) )
        END DO
     END IF
     hpsi_d = hpsi_host
     !
  END IF
  !
  ! ... With Gamma-only trick, Im(H*psi)(G=0) = 0 by definition,
  ! ... but it is convenient to explicitly set it to 0 to prevent trouble
  !
  IF ( gamma_only .AND. gstart == 2 ) then
      !$cuf kernel do(1)
      do i=1,m
         hpsi_d(1,i) = CMPLX( DBLE( hpsi_d(1,i) ), 0.D0 ,kind=DP)
      end do
  end if
  !
  if (need_host_copy) then
      DEALLOCATE(psi_host , hpsi_host )
  end if
  CALL stop_clock( 'h_psi' )
  !
  RETURN
  !
END SUBROUTINE h_psi__gpu



SUBROUTINE h_psi_gpu_compatibility( lda, n, m, psi_d, hpsi_d )
  USE kinds,            ONLY : DP
  USE noncollin_module, ONLY : npol
  USE cudafor
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)      :: lda, n, m
  COMPLEX(DP), DEVICE, INTENT(IN)  :: psi_d(lda*npol,m) 
  COMPLEX(DP), DEVICE, INTENT(OUT) :: hpsi_d(lda*npol,m)


  COMPLEX(DP), ALLOCATABLE :: psi(:,:)
  COMPLEX(DP), ALLOCATABLE :: hpsi(:,:)

  
  ALLOCATE(psi(lda*npol,m), hpsi(lda*npol,m) )
  
  psi = psi_d
  CALL h_psi( lda, n, m, psi, hpsi )
  hpsi_d = hpsi
  
  DEALLOCATE(psi, hpsi )
  
END SUBROUTINE h_psi_gpu_compatibility
#endif
