!
! Copyright (C) 2002-2022 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE h_psi_gpu( lda, n, m, psi, hpsi )
  !----------------------------------------------------------------------------
  !! This routine computes the product of the Hamiltonian matrix with m 
  !! wavefunctions contained in psi.
  !
  !! \(\textit{Wrapper routine}\): performs bgrp parallelization on 
  !! non-distributed bands. If suitable and required, calls old H\psi 
  !! routine h_psi_ .
  !
  USE kinds,              ONLY: DP
  USE noncollin_module,   ONLY: npol
  USE xc_lib,             ONLY: exx_is_active
  USE mp_bands,           ONLY: use_bgrp_in_hpsi, inter_bgrp_comm
  USE mp,                 ONLY: mp_allgather, mp_size, &
                                mp_type_create_column_section, mp_type_free
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of arrays psi, spsi, hpsi
  INTEGER, INTENT(IN) :: n
  !! true dimension of psi, spsi, hpsi
  INTEGER, INTENT(IN) :: m
  !! number of states psi
  COMPLEX(DP), INTENT(IN) :: psi(lda*npol,m)
  !! the wavefunction
  COMPLEX(DP), INTENT(OUT) :: hpsi(lda*npol,m)
  !! Hamiltonian dot psi
  !
  ! ... local variables
  !
  INTEGER :: m_start, m_end
  INTEGER :: column_type
  INTEGER, ALLOCATABLE :: recv_counts(:), displs(:)
  !
  !
  CALL start_clock( 'h_psi_bgrp' ); !write (*,*) 'start h_psi_bgrp'; FLUSH(6)
  !
  ! band parallelization with non-distributed bands is performed if
  ! 1. enabled (variable use_bgrp_in_hpsi must be set to .T.)
  ! 2. exact exchange is not active (if it is, band parallelization is already
  !    used in exx routines called by Hpsi)
  ! 3. there is more than one band, otherwise there is nothing to parallelize
  !
  IF (use_bgrp_in_hpsi .AND. .NOT. exx_is_active() .AND. m > 1) THEN
     !
     ! use band parallelization here
     ALLOCATE( recv_counts(mp_size(inter_bgrp_comm)), displs(mp_size(inter_bgrp_comm)) )
     CALL divide_all( inter_bgrp_comm, m, m_start, m_end, recv_counts,displs )
     !$acc host_data use_device(hpsi)
     CALL mp_type_create_column_section( hpsi(1,1), 0, lda*npol, lda*npol, column_type )
     !$acc end host_data
     !
     ! Check if there at least one band in this band group
     IF (m_end >= m_start) &
        CALL h_psi__gpu( lda, n, m_end-m_start+1, psi(1,m_start), hpsi(1,m_start) )
     !$acc host_data use_device(hpsi)
     CALL mp_allgather( hpsi, column_type, recv_counts, displs, inter_bgrp_comm)
     !$acc end host_data
     !
     CALL mp_type_free( column_type )
     DEALLOCATE( recv_counts )
     DEALLOCATE( displs )
     !
  ELSE
     ! don't use band parallelization here
     CALL h_psi__gpu( lda, n, m, psi, hpsi )
     !
  ENDIF
  !
  CALL stop_clock( 'h_psi_bgrp' )
  !
  !
  RETURN
  !
END SUBROUTINE h_psi_gpu
!
!----------------------------------------------------------------------------
SUBROUTINE h_psi__gpu( lda, n, m, psi, hpsi )
  !----------------------------------------------------------------------------
  !! This routine computes the product of the Hamiltonian matrix with m 
  !! wavefunctions contained in psi.
  !
#if defined(__CUDA)
  USE cudafor
  USE becmod,                  ONLY: calbec
#endif
  USE kinds,                   ONLY: DP
  USE bp,                      ONLY: lelfield, l3dstring, gdir, efield, efield_cry
  USE becmod,                  ONLY: bec_type, becp
  USE lsda_mod,                ONLY: current_spin
  USE uspp,                    ONLY: nkb, vkb
  USE ldaU,                    ONLY: lda_plus_u, lda_plus_u_kind, Hubbard_projectors
  USE gvect,                   ONLY: gstart
  USE control_flags,           ONLY: gamma_only, offload_type
  USE noncollin_module,        ONLY: npol, noncolin
  USE realus,                  ONLY: real_space, invfft_orbital_gamma, fwfft_orbital_gamma, &
                                     calbec_rs_gamma, add_vuspsir_gamma, invfft_orbital_k,  &
                                     fwfft_orbital_k, calbec_rs_k, add_vuspsir_k,           & 
                                     v_loc_psir_inplace
  USE fft_base,                ONLY: dffts
  USE exx,                     ONLY: use_ace, vexx, vexxace_gamma_gpu, vexxace_k_gpu
  USE xc_lib,                  ONLY: exx_is_active, xclib_dft_is
  USE fft_helper_subroutines
  USE device_memcpy_m,         ONLY: dev_memcpy, dev_memset
  !
  USE scf,                     ONLY: vrs  
  USE wvfct,                   ONLY: g2kin  
#if defined(__OSCDFT)
  USE plugin_flags,            ONLY : use_oscdft
  USE oscdft_base,             ONLY : oscdft_ctx
  USE oscdft_functions_gpu,    ONLY : oscdft_h_psi_gpu
#endif
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)     :: lda, n, m
  COMPLEX(DP), INTENT(IN)  :: psi(lda*npol,m)
  COMPLEX(DP), INTENT(OUT) :: hpsi(lda*npol,m)
  !
  INTEGER     :: ipol, ibnd, i
  REAL(dp)    :: ee
  !
  CALL start_clock( 'h_psi' ); !write (*,*) 'start h_psi';FLUSH(6)
  !
  ! ... Here we add the kinetic energy (k+G)^2 psi and clean up garbage
  !
  !$acc parallel loop collapse(2) present(g2kin, hpsi, psi)
  DO ibnd = 1, m
     DO i=1, lda
        IF (i <= n) THEN
           hpsi (i, ibnd) = g2kin (i) * psi (i, ibnd)
           IF ( noncolin ) THEN
              hpsi (lda+i, ibnd) = g2kin (i) * psi (lda+i, ibnd)
           END IF
        ELSE
           hpsi (i, ibnd) = (0.0_dp, 0.0_dp)
           IF ( noncolin ) THEN
              hpsi (lda+i, ibnd) = (0.0_dp, 0.0_dp)
           END IF
        END IF
     END DO
  END DO

  CALL start_clock( 'h_psi:pot' ); !write (*,*) 'start h_pot';FLUSH(6)
  !
  ! ... Here the product with the local potential V_loc psi
  !
  IF ( gamma_only ) THEN
     ! 
     IF ( real_space .AND. nkb > 0  ) THEN
        !
        ! ... real-space algorithm
        ! ... fixme: real_space without beta functions does not make sense
        !
        IF ( dffts%has_task_groups ) &
             CALL errore( 'h_psi', 'task_groups not implemented with real_space', 1 )
        DO ibnd = 1, m, 2
           ! ... transform psi to real space -> psic 
           CALL invfft_orbital_gamma(psi, ibnd, m )
           ! ... compute becp%r = < beta|psi> from psic in real space
           CALL start_clock( 'h_psi:calbec' )
           CALL calbec_rs_gamma( ibnd, m, becp%r )
           !$acc update device(becp%r)
           CALL stop_clock( 'h_psi:calbec' )
           ! ... psic -> vrs * psic (psic overwritten will become hpsi)
           CALL v_loc_psir_inplace( ibnd, m ) 
           ! ... psic (hpsi) -> psic + vusp
           CALL  add_vuspsir_gamma( ibnd, m )
           ! ... transform psic back in reciprocal space and add it to hpsi
           CALL fwfft_orbital_gamma( hpsi, ibnd, m, add_to_orbital=.TRUE. )
        ENDDO
        !
     ELSE
        ! ... usual reciprocal-space algorithm
        CALL vloc_psi_gamma_acc ( lda, n, m, psi, vrs(1,current_spin), hpsi )
        !
     ENDIF 
     !
  ELSE IF ( noncolin ) THEN 
     !
     CALL vloc_psi_nc_acc ( lda, n, m, psi, vrs, hpsi )
     !
  ELSE  
     ! 
     IF ( real_space .and. nkb > 0  ) then 
        !
        ! ... real-space algorithm
        ! ... fixme: real_space without beta functions does not make sense
        !
        IF ( dffts%has_task_groups ) &
             CALL errore( 'h_psi', 'task_groups not implemented with real_space', 1 )
        !
        DO ibnd = 1, m
           ! ... transform psi to real space -> psic 
           CALL invfft_orbital_k(psi, ibnd, m )
           ! ... compute becp%r = < beta|psi> from psic in real space
           CALL start_clock( 'h_psi:calbec' )
           CALL calbec_rs_k( ibnd, m )
           !$acc update device(becp%k)
           CALL stop_clock( 'h_psi:calbec' )
           ! ... psic -> vrs * psic (psic overwritten will become hpsi)
           CALL v_loc_psir_inplace( ibnd, m )
           ! ... psic (hpsi) -> psic + vusp
           CALL add_vuspsir_k( ibnd, m )
           ! ... transform psic back in reciprocal space and add it to hpsi
           CALL fwfft_orbital_k( hpsi, ibnd, m, add_to_orbital=.TRUE. )
           !
        ENDDO
        !
     ELSE
        !
        CALL vloc_psi_k_acc ( lda, n, m, psi, vrs(1,current_spin), hpsi )
        !
     ENDIF
     !
  ENDIF  
  !
  ! ... Here the product with the non local potential V_NL psi
  ! ... (not in the real-space case: it is done together with V_loc)
  !
  IF ( nkb > 0 .AND. .NOT. real_space) THEN
     !
     CALL start_clock( 'h_psi:calbec' )
#if defined(__CUDA)
     Call calbec(offload_type, n, vkb, psi, becp, m )
#endif
     CALL stop_clock( 'h_psi:calbec' )
     !$acc host_data use_device(hpsi)
     CALL add_vuspsi_gpu( lda, n, m, hpsi )
     !$acc end host_data
     !
  END IF
  !  
  CALL stop_clock( 'h_psi:pot' )
  !
  IF (xclib_dft_is('meta')) THEN
     call h_psi_meta (lda, n, m, psi, hpsi)
  end if
  !
  ! ... Here we add the Hubbard potential times psi
  !
  IF ( lda_plus_u .AND. Hubbard_projectors.NE."pseudo" ) THEN
     !
     IF ( noncolin ) THEN
        !FIXME: vhpsi_nc must be ported to GPU
        !$acc update host(psi, hpsi)
        CALL vhpsi_nc( lda, n, m, psi, hpsi )
        !$acc update device(psi, hpsi)
     ELSE
        IF ( lda_plus_u_kind.EQ.0 .OR. lda_plus_u_kind.EQ.1 ) THEN
          ! DFT + U
          CALL vhpsi_gpu( lda, n, m, psi, hpsi )  ! DFT+U
        ELSEIF ( lda_plus_u_kind.EQ.2 ) THEN
          ! DFT+U+V  FIXME: must be ported to GPU
          !$acc update host(psi, hpsi)
          CALL vhpsi( lda, n, m, psi, hpsi )
          !$acc update device(psi, hpsi)
        ENDIF
     ENDIF
     !
  ENDIF
  !
  ! ... Here the exact-exchange term Vxx psi
  !
  IF ( exx_is_active() ) THEN
     IF ( use_ace) THEN
        IF (gamma_only) THEN
           !$acc host_data use_device(psi, hpsi)
           CALL vexxace_gamma_gpu(lda,m,psi,ee,hpsi)
           !$acc end host_data
        ELSE
           !$acc host_data use_device(psi, hpsi)
           CALL vexxace_k_gpu(lda,m,psi,ee,hpsi)
           !$acc end host_data
        END IF
     ELSE
        !$acc update host (psi,hpsi)
        CALL vexx( lda, n, m, psi, hpsi, becp )
        !$acc update device(hpsi)
     END IF
  END IF
  !
  ! ... electric enthalpy if required
  !
  IF ( lelfield ) THEN
     !
     !$acc update host (psi,hpsi)
     IF ( .NOT.l3dstring ) THEN
        CALL h_epsi_her_apply( lda, n, m, psi, hpsi, gdir, efield )
     ELSE
        DO ipol=1,3
           CALL h_epsi_her_apply( lda, n, m, psi, hpsi, ipol, efield_cry(ipol) )
        END DO
     END IF
     !$acc update device(hpsi)
     !
  END IF
#if defined(__OSCDFT)
  IF ( use_oscdft ) THEN
     CALL oscdft_h_psi_gpu(oscdft_ctx, lda, n, m, psi, hpsi)
  END IF
#endif
  !
  ! ... With Gamma-only trick, Im(H*psi)(G=0) = 0 by definition,
  ! ... but it is convenient to explicitly set it to 0 to prevent trouble
  !
  IF ( gamma_only .AND. gstart == 2 ) then
      !$acc kernels 
      do i=1,m
         hpsi(1,i) = CMPLX( DBLE( hpsi(1,i) ), 0.D0 ,kind=DP)
      end do
      !$acc end kernels
  end if
  !
  CALL stop_clock( 'h_psi' )
  !
  RETURN
  !
END SUBROUTINE h_psi__gpu
