!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE s_psi_acc( lda, n, m, psi, spsi )
  !--------------------------------------------------------------------
  !! This routine applies the S matrix to m wavefunctions psi and puts 
  !! the results in spsi.
  !! Requires the products of psi with all beta functions in array 
  !! becp(nkb,m) (calculated in h_psi or by calbec).
  !
  !! \(\textit{Wrapper routine}\): performs bgrp parallelization on 
  !! non-distributed bands if suitable and required, calls old S\psi
  !! routine s_psi_ . See comments in h_psi.f90 about band 
  !! parallelization.
  !
  USE kinds,            ONLY : DP
  USE noncollin_module, ONLY : npol
  USE xc_lib,           ONLY : exx_is_active
  USE mp_bands,         ONLY : use_bgrp_in_hpsi, inter_bgrp_comm
  USE mp,               ONLY : mp_allgather, mp_size, &
                               mp_type_create_column_section, mp_type_free
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of arrays psi, spsi
  INTEGER, INTENT(IN) :: n
  !! true dimension of psi, spsi
  INTEGER, INTENT(IN) :: m
  !! number of states psi
  COMPLEX(DP), INTENT(IN) :: psi(lda*npol,m)
  !! the m wavefunctions
  COMPLEX(DP), INTENT(OUT)::spsi(lda*npol,m)
  !! S matrix dot wavefunctions psi
  !
  ! ... local variables
  !
  INTEGER :: m_start, m_end
  INTEGER :: column_type
  INTEGER, ALLOCATABLE :: recv_counts(:), displs(:)
  !
  CALL start_clock( 's_psi_bgrp' )
  !
  IF (use_bgrp_in_hpsi .AND. .NOT. exx_is_active() .AND. m > 1) THEN
     ! use band parallelization here
     ALLOCATE( recv_counts(mp_size(inter_bgrp_comm)), displs(mp_size(inter_bgrp_comm)) )
     CALL divide_all( inter_bgrp_comm,m,m_start,m_end, recv_counts,displs )
     !$acc host_data use_device(spsi)
     CALL mp_type_create_column_section( spsi(1,1), 0, lda*npol, lda*npol, column_type )
     !$acc end host_data
     !
     ! Check if there at least one band in this band group
     IF (m_end >= m_start) &
        CALL s_psi__acc( lda, n, m_end-m_start+1, psi(1,m_start), spsi(1,m_start) )
     !$acc host_data use_device(spsi)
     CALL mp_allgather( spsi, column_type, recv_counts, displs, inter_bgrp_comm )
     !$acc end host_data
     !
     CALL mp_type_free( column_type )
     DEALLOCATE( recv_counts )
     DEALLOCATE( displs )
  ELSE
     ! don't use band parallelization here
     CALL s_psi__acc( lda, n, m, psi, spsi )
  ENDIF
  !
  CALL stop_clock( 's_psi_bgrp' )
  !
  RETURN
  !
END SUBROUTINE s_psi_acc
!
!
!----------------------------------------------------------------------------
SUBROUTINE s_psi__acc( lda, n, m, psi, spsi )
  !----------------------------------------------------------------------------
  !! This routine applies the S matrix to m wavefunctions psi and puts 
  !! the results in spsi.
  !! Requires the products of psi with all beta functions in array 
  !! becp(nkb,m) (calculated in h_psi or by calbec).
  !
  USE kinds,            ONLY: DP
  USE becmod,           ONLY: becp
  USE uspp,             ONLY: vkb, nkb, okvan, qq_at, qq_so, ofsbeta
  USE uspp_param,       ONLY: upf, nh, nhm
  USE ions_base,        ONLY: nat, nsp, ityp
  USE control_flags,    ONLY: gamma_only 
  USE noncollin_module, ONLY: npol, noncolin, lspinorb
  USE realus,           ONLY: real_space, fwfft_orbital_gamma, s_psir_gamma, &
                              fwfft_orbital_k, s_psir_k
  USE wavefunctions,    ONLY: psic
  USE fft_base,         ONLY: dffts
#if defined (__CUDA)
  USE device_memcpy_m,  ONLY : dev_memcpy
#endif
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of arrays psi, spsi
  INTEGER, INTENT(IN) :: n
  !! true dimension of psi, spsi
  INTEGER, INTENT(IN) :: m
  !! number of states psi
  COMPLEX(DP), INTENT(IN) :: psi(lda*npol,m)
  !! the m wavefunctions
  COMPLEX(DP), INTENT(OUT)::spsi(lda*npol,m)
  !! S matrix dot wavefunctions psi
  !
  ! ... local variables
  !
  INTEGER :: ibnd
  !
  ! ... initialize  spsi
  !
#if defined(__CUDA)  
  !$acc host_data use_device(psi, spsi)
  CALL dev_memcpy( spsi , psi )
  !$acc end host_data
#else
  CALL threaded_memcpy( spsi, psi, lda*npol*m*2 )
#endif
  !
  IF ( nkb == 0 .OR. .NOT. okvan ) RETURN
  !
  CALL start_clock( 's_psi' )  
  !
  ! ... The product with the beta functions
  !
  IF ( gamma_only ) THEN
     !
     IF ( real_space ) THEN
        !
        DO ibnd = 1, m, 2
!SdG: the becp are already computed ! no need to invfft psi to real space.
!           CALL invfft_orbital_gamma( psi, ibnd, m ) 
!SdG: we just need to clean psic in real space ...
           CALL threaded_barrier_memset(psic, 0.D0, dffts%nnr*2)
!SdG: ... before computing the us-only contribution ...
           CALL s_psir_gamma( ibnd, m )
!SdG: ... and add it to spsi (already containing psi).
           CALL fwfft_orbital_gamma( spsi, ibnd, m, add_to_orbital=.TRUE. )
        ENDDO
        !
     ELSE
        !
        CALL s_psi_gamma_acc()
        !
     ENDIF
     !
  ELSEIF ( noncolin ) THEN
     !
     CALL s_psi_nc_acc()
     !
  ELSE 
     !
     IF ( real_space ) THEN
        !
        DO ibnd = 1, m
!SdG: the becp are already computed ! no need to invfft psi to real space.
!           CALL invfft_orbital_k( psi, ibnd, m )
!SdG: we just need to clean psic in real space ...
           CALL threaded_barrier_memset(psic, 0.D0, dffts%nnr*2)
!SdG: ... before computing the us-only contribution ...
           CALL s_psir_k( ibnd, m )
!SdG: ... and add it to spsi (already containing psi).
           CALL fwfft_orbital_k( spsi, ibnd, m, add_to_orbital=.TRUE. )
        ENDDO
        !
     ELSE
        !
        CALL s_psi_k_acc()
        !
     ENDIF    
     !
  ENDIF    
  !
  CALL stop_clock( 's_psi' )
  !
  RETURN
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     SUBROUTINE s_psi_gamma_acc()
       !---------------------------------------------------------------------
       !! Gamma version of \(\textrm{s_psi}\) routine.
       !
       USE mp,            ONLY : mp_get_comm_null, mp_circular_shift_left
       !
       IMPLICIT NONE  
       !
       ! ... local variables
       !
       INTEGER :: ikb, jkb, ih, jh, na, nt, ibnd, ierr
       ! counters
       REAL(DP), ALLOCATABLE :: ps(:,:), becpr(:,:)
       !$acc declare device_resident(ps, becpr)
       ! the product vkb and psi
       !
       !$acc data present( becp%r, vkb )
       !
       ALLOCATE( becpr( size(becp%r,1), size(becp%r,2)) )
       !$acc kernels 
       becpr(:,:) = becp%r(:,:)
       !$acc end  kernels
       !
       ! becp(l,i) = <beta_l|psi_i>, with vkb(n,l)=|beta_l>
       !
       ALLOCATE( ps( nkb, m ), STAT=ierr )
       IF( ierr /= 0 ) &
          CALL errore( ' s_psi_gamma ', ' cannot allocate memory (ps) ', ABS(ierr) )
       !    
       !$acc kernels
       ps(:,:) = 0.0_DP
       !$acc end kernels
       !
       !   In becp=<vkb_i|psi_j> terms corresponding to atom na of type nt
       !   run from index i=ofsbeta(na)+1 to i=ofsbeta(na)+nh(nt)
       !
       DO nt = 1, nsp
          IF ( upf(nt)%tvanp ) THEN
             DO na = 1, nat
                IF ( ityp(na) == nt ) THEN
                   !
                   ! Next operation computes ps(l',i)=\sum_m qq(l,m) becp(m',i)
                   ! (l'=l+ijkb0, m'=m+ijkb0, indices run from 1 to nh(nt))
                   !
                   !$acc host_data use_device(qq_at, becpr, ps)
                   CALL MYDGEMM('N', 'N', nh(nt), m, nh(nt), 1.0_dp, &
                                  qq_at(1,1,na), nhm, becpr(ofsbeta(na)+1,1),&
                                  nkb, 0.0_dp, ps(ofsbeta(na)+1,1), nkb )
                   !$acc end host_data
                ENDIF
             ENDDO
          ENDIF
       ENDDO
       !
       IF ( m == 1 ) THEN
          !$acc host_data use_device(vkb, ps, spsi)
          CALL MYDGEMV( 'N', 2 * n, nkb, 1.D0, vkb, &
               2 * lda, ps, 1, 1.D0, spsi, 1 )
          !$acc end host_data
       ELSE
          !$acc host_data use_device(vkb, ps, spsi)
          CALL MYDGEMM( 'N', 'N', 2 * n, m, nkb, 1.D0, vkb, &
               2 * lda, ps, nkb, 1.D0, spsi, 2*lda )
          !$acc end host_data
       ENDIF
       !
       DEALLOCATE( ps, becpr ) 
       !
       !$acc end data 
       !
       RETURN
       !
     END SUBROUTINE s_psi_gamma_acc
     !
     !-----------------------------------------------------------------------
     SUBROUTINE s_psi_k_acc()
       !-----------------------------------------------------------------------
       !! k-points version of \(\textrm{s_psi}\) routine.
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       INTEGER :: ikb, jkb, ih, jh, na, nt, ibnd, ierr
       ! counters
       COMPLEX(DP), ALLOCATABLE :: ps(:,:), qqc(:,:), becpk(:,:)
       !$acc declare device_resident(ps, qqc, becpk)
       ! ps = product vkb and psi ; qqc = complex version of qq
       !
       !$acc data present( becp%k, vkb )
       !
       ALLOCATE( becpk( size(becp%k,1), size(becp%k,2)) )
       !$acc kernels 
       becpk(:,:) = becp%k(:,:)
       !$acc end  kernels
       !
       ALLOCATE( ps( nkb, m ), STAT=ierr )
       !
       IF( ierr /= 0 ) &
          CALL errore( ' s_psi_k ', ' cannot allocate memory (ps) ', ABS(ierr) )
       !
       DO nt = 1, nsp
          !
          IF ( upf(nt)%tvanp ) THEN
             ! qq is real:  copy it into a complex variable to perform
             ! a zgemm - simple but sub-optimal solution
             ALLOCATE( qqc(nh(nt),nh(nt)) )
             DO na = 1, nat
                IF ( ityp(na) == nt ) THEN
                   !$acc kernels
                   qqc(:,:) = CMPLX ( qq_at(1:nh(nt),1:nh(nt),na), 0.0_DP, KIND=DP )
                   !$acc end kernels
                   !$acc host_data use_device( qqc, becpk, ps)
                   CALL MYZGEMM('N','N', nh(nt), m, nh(nt), (1.0_DP,0.0_DP), &
                        qqc, nh(nt), becpk(ofsbeta(na)+1,1), nkb, &
                        (0.0_DP,0.0_DP), ps(ofsbeta(na)+1,1), nkb )
                   !$acc end host_data
                ENDIF
             ENDDO
             DEALLOCATE( qqc )
             !
          ELSE
             !
             IF (nh(nt)>0) THEN
                !$acc kernels present_or_copyin(ityp, ofsbeta, nh)
                DO na = 1, nat
                   IF ( ityp(na) == nt ) THEN
                      ps(ofsbeta(na)+1:ofsbeta(na)+nh(nt),1:m) = (0.0_DP,0.0_DP)
                   ENDIF
                ENDDO
                !$acc end kernels
             ENDIF
             !
          ENDIF
          !
       ENDDO
       !
       IF ( m == 1 ) THEN
          !
          !$acc host_data use_device(vkb, ps, spsi)
          CALL MYZGEMV( 'N', n, nkb, ( 1.D0, 0.D0 ), vkb, &
                      lda, ps, 1, ( 1.D0, 0.D0 ), spsi, 1 )
          !$acc end host_data
          !
       ELSE
          !
          !$acc host_data use_device(vkb, ps, spsi)
          CALL MYZGEMM( 'N', 'N', n, m, nkb, ( 1.D0, 0.D0 ), vkb, &
                      lda, ps, nkb, ( 1.D0, 0.D0 ), spsi, lda )
          !$acc end host_data
          !
       ENDIF
       !
       DEALLOCATE( ps, becpk )
       !
       !$acc end data
       !
       RETURN
       !
     END SUBROUTINE s_psi_k_acc
     !
     !
     !-----------------------------------------------------------------------
       SUBROUTINE s_psi_nc_acc ( )
       !-----------------------------------------------------------------------
       !! k-points noncolinear/spinorbit version of \(\textrm{s_psi}\) routine.
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       INTEGER :: ikb, jkb, ih, jh, na, nt, ibnd, ipol, ierr, nh_nt, ofs_na
       ! counters
       COMPLEX (DP), ALLOCATABLE :: ps(:,:,:), becpnc(:,:,:) 
       !$acc declare device_resident(ps, becpnc)
       ! the product vkb and psi
       !
       !$acc data present(vkb, becp%nc)
       !
       ALLOCATE( becpnc(size(becp%nc,1),size(becp%nc,2),size(becp%nc,3)) )
       !$acc kernels
       becpnc(:,:,:) = becp%nc(:,:,:)
       !$acc end kernels
       !
       ALLOCATE( ps(nkb,npol,m), STAT=ierr )
       IF( ierr /= 0 ) &
          CALL errore( ' s_psi_nc ', ' cannot allocate memory (ps) ', ABS(ierr) )
       !$acc kernels
       ps(:,:,:) = (0.D0,0.D0)
       !$acc end kernels
       !
       DO nt = 1, nsp
          !
          IF ( upf(nt)%tvanp ) THEN
             !
             nh_nt = nh(nt)
             !
             !$acc parallel present_or_copyin(ityp, nh_nt, nt, ofsbeta) 
             !$acc loop gang private(ofs_na, ikb, jkb) 
             DO na = 1, nat
                IF ( ityp(na) == nt ) THEN
                   ofs_na = ofsbeta(na)
                   DO ih = 1, nh_nt
                      ikb = ofs_na + ih
                      DO jh = 1, nh_nt
                         jkb = ofs_na + jh
                         IF ( .NOT. lspinorb ) THEN
                            !$acc loop vector collapse(2)
                            DO ipol = 1, npol
                               DO ibnd = 1, m
                                  ps(ikb,ipol,ibnd) = ps(ikb,ipol,ibnd) + &
                                       qq_at(ih,jh,na)*becpnc(jkb,ipol,ibnd)
                               ENDDO
                            ENDDO
                         ELSE
                            !$acc loop vector 
                            DO ibnd = 1, m
                               ps(ikb,1,ibnd) = ps(ikb,1,ibnd) + &
                                    qq_so(ih,jh,1,nt)*becpnc(jkb,1,ibnd)+ &
                                    qq_so(ih,jh,2,nt)*becpnc(jkb,2,ibnd)
                               ps(ikb,2,ibnd) = ps(ikb,2,ibnd) + &
                                    qq_so(ih,jh,3,nt)*becpnc(jkb,1,ibnd)+ &
                                    qq_so(ih,jh,4,nt)*becpnc(jkb,2,ibnd)
                            ENDDO
                         ENDIF
                      ENDDO
                   ENDDO
                ENDIF
             ENDDO
             !$acc end parallel 
             !
          ENDIF
          !
       ENDDO
       !
       !$acc host_data use_device(vkb, ps, spsi)
       CALL MYZGEMM ( 'N', 'N', n, m*npol, nkb, (1.d0,0.d0) , vkb, &
                    lda, ps, nkb, (1.d0,0.d0) , spsi(1,1), lda )
       !$acc end host_data
       !
       DEALLOCATE( ps, becpnc )
       !
       !$acc end data
       !
       RETURN
       !
    END SUBROUTINE s_psi_nc_acc
    !
END SUBROUTINE s_psi__acc
