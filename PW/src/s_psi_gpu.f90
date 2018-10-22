!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
#if ! defined (__CUDA)
#define cublasDGEMM DGEMM
#endif
!@njs: s_psi, psi, spsi, s_psi_, s_psi_k, s_psi_gamma, s_psi_nc
SUBROUTINE s_psi_gpu( lda, n, m, psi_d, spsi_d )
  !----------------------------------------------------------------------------
  !
  ! ... This routine applies the S matrix to m wavefunctions psi
  ! ... and puts the results in spsi.
  ! ... Requires the products of psi with all beta functions
  ! ... in array becp(nkb,m) (calculated in h_psi or by calbec)
  !
  ! ... input:
  !
  ! ...    lda   leading dimension of arrays psi, spsi
  ! ...    n     true dimension of psi, spsi
  ! ...    m     number of states psi
  ! ...    psi
  !
  ! ... output:
  !
  ! ...    spsi  S*psi
  !
  ! --- Wrapper routine: performs bgrp parallelization on non-distributed bands
  ! --- if suitable and required, calls old S\psi routine s_psi_
  ! --- See comments in h_psi.f90 about band parallelization
  !
#if defined(__CUDA)
  USE cudafor
#endif
  USE kinds,            ONLY : DP
  USE noncollin_module, ONLY : npol
  USE funct,            ONLY : exx_is_active
  USE mp_bands,         ONLY : use_bgrp_in_hpsi, inter_bgrp_comm
  USE mp,               ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda, n, m
  COMPLEX(DP), INTENT(IN) :: psi_d(lda*npol,m)
  COMPLEX(DP), INTENT(OUT)::spsi_d(lda*npol,m)
#if defined(__CUDA)
  attributes(DEVICE) :: psi_d, spsi_d
#endif
  !
  INTEGER     :: m_start, m_end, i
  !
  CALL start_clock( 's_psi_bgrp' )

  IF (use_bgrp_in_hpsi .AND. .NOT. exx_is_active() .AND. m > 1) THEN
     ! use band parallelization here
     spsi_d(:,:) = (0.d0,0.d0)
     CALL divide(inter_bgrp_comm,m,m_start,m_end)
     !write(6,*) m, m_start,m_end
     ! Check if there at least one band in this band group
     IF (m_end >= m_start) &
        CALL s_psi__gpu( lda, n, m_end-m_start+1, psi_d(1,m_start), spsi_d(1,m_start) )
     CALL mp_sum(spsi_d,inter_bgrp_comm)
  ELSE
     ! don't use band parallelization here
     CALL s_psi__gpu( lda, n, m, psi_d, spsi_d )
  END IF

  CALL stop_clock( 's_psi_bgrp' )
  RETURN

END SUBROUTINE s_psi_gpu
!
!----------------------------------------------------------------------------
SUBROUTINE s_psi__gpu( lda, n, m, psi_d, spsi_d )
  !----------------------------------------------------------------------------
  !
  ! ... This routine applies the S matrix to m wavefunctions psi
  ! ... and puts the results in spsi.
  ! ... Requires the products of psi with all beta functions
  ! ... in array becp(nkb,m) (calculated in h_psi or by calbec)
  !
  ! ... input:
  !
  ! ...    lda   leading dimension of arrays psi, spsi
  ! ...    n     true dimension of psi, spsi
  ! ...    m     number of states psi
  ! ...    psi
  !
  ! ... output:
  !
  ! ...    spsi  S*psi
  !
#if defined (__CUDA)
  USE cudafor
  USE cublas
#endif
  USE kinds,      ONLY : DP
  USE becmod_gpum,ONLY : becp_d
  USE uspp,       ONLY : nkb, okvan, indv_ijkb0
  USE spin_orb,   ONLY : lspinorb
  USE uspp_param, ONLY : upf, nh, nhm
  USE ions_base,  ONLY : nat, nsp, ityp
  USE control_flags,    ONLY: gamma_only 
  USE noncollin_module, ONLY: npol, noncolin
  USE realus,     ONLY :  real_space, &
                  invfft_orbital_gamma, fwfft_orbital_gamma, calbec_rs_gamma, s_psir_gamma, &
                  invfft_orbital_k, fwfft_orbital_k, calbec_rs_k, s_psir_k
  !
  USE uspp_gpum,  ONLY : vkb_d, using_vkb_d, using_indv_ijkb0
  USE becmod_gpum, ONLY : using_becp_r_d, using_becp_k_d, using_becp_nc_d
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda, n, m
  COMPLEX(DP), INTENT(IN) :: psi_d(lda*npol,m)
  COMPLEX(DP), INTENT(OUT)::spsi_d(lda*npol,m)
#if defined(__CUDA)
  attributes(DEVICE) :: psi_d, spsi_d
  !
  COMPLEX(DP), PINNED, ALLOCATABLE :: psi_host(:,:)
  COMPLEX(DP), PINNED, ALLOCATABLE ::spsi_host(:,:)
#else
  COMPLEX(DP), ALLOCATABLE :: psi_host(:,:)
  COMPLEX(DP), ALLOCATABLE ::spsi_host(:,:)
#endif
  !
  INTEGER :: ibnd
  !
  LOGICAL :: need_host_copy
  !
  ! ... initialize  spsi
  !
  spsi_d = psi_d
  !
  IF ( nkb == 0 .OR. .NOT. okvan ) RETURN
  !
  need_host_copy = real_space
  IF (need_host_copy) THEN
      ALLOCATE(psi_host(lda*npol,m), spsi_host(lda*npol,m))
      psi_host  = psi_d
      spsi_host = spsi_d
  END IF
  !
  CALL start_clock( 's_psi' )  
  !
  ! ... The product with the beta functions
  !
  IF ( gamma_only ) THEN
     !
     IF (real_space ) THEN
        !
        DO ibnd = 1, m, 2
           !   transform the orbital to real space
           CALL invfft_orbital_gamma(psi_host,ibnd,m) 
           CALL s_psir_gamma(ibnd,m)
           CALL fwfft_orbital_gamma(spsi_host,ibnd,m)
        END DO
        spsi_d = spsi_host
        !
     ELSE
        !
        CALL s_psi_gamma_gpu()
        !
     END IF
     !
  ELSE IF ( noncolin ) THEN
     !
     CALL s_psi_nc_gpu()
     !
  ELSE 
     !
     IF (real_space ) THEN
        !
        DO ibnd = 1, m
           !   transform the orbital to real space
           CALL invfft_orbital_k(psi_host,ibnd,m) 
           CALL s_psir_k(ibnd,m)
           CALL fwfft_orbital_k(spsi_host,ibnd,m)
           spsi_d = spsi_host
        END DO
        !
     ELSE
        !
        CALL s_psi_k_gpu()
        !
     END IF    
     !
  END IF    
  !
  CALL stop_clock( 's_psi' )
  !
  RETURN
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     SUBROUTINE s_psi_gamma_gpu()
       !-----------------------------------------------------------------------
       ! 
       ! ... gamma version
       !
       USE mp, ONLY: mp_get_comm_null, mp_circular_shift_left
       USE gbuffers,   ONLY : dev_buf
       USE uspp_gpum,  ONLY : qq_at_d, using_qq_at_d
       !
       IMPLICIT NONE  
       !
       ! ... here the local variables
       !
       INTEGER :: ikb, jkb, ih, jh, na, nt, ibnd, ierr
         ! counters
       INTEGER :: nproc, mype, m_loc, m_begin, ibnd_loc, icyc, icur_blk, m_max
         ! data distribution indexes
       INTEGER, EXTERNAL :: ldim_block, gind_block
         ! data distribution functions
       REAL(DP), POINTER :: ps_d(:,:)
#if defined(__CUDA)
       attributes(DEVICE) :: ps_d
#endif
         ! the product vkb and psi
       !
       CALL using_vkb_d(0)
       CALL using_becp_r_d(0)
       CALL using_qq_at_d(0)
       CALL using_indv_ijkb0(0)
       !
       IF( becp_d%comm == mp_get_comm_null() ) THEN
          nproc   = 1
          mype    = 0
          m_loc   = m
          m_begin = 1
          m_max   = m
       ELSE
          !
          ! becp(l,i) = <beta_l|psi_i>, with vkb(n,l)=|beta_l>
          ! in this case becp(l,i) are distributed (index i is)
          !
          nproc   = becp_d%nproc
          mype    = becp_d%mype
          m_loc   = becp_d%nbnd_loc
          m_begin = becp_d%ibnd_begin
          m_max   = SIZE( becp_d%r_d, 2 )
          IF( ( m_begin + m_loc - 1 ) > m ) m_loc = m - m_begin + 1
       END IF
       !
       CALL dev_buf%lock_buffer(ps_d, (/ nkb, m_max /), ierr)
       
       IF( ierr /= 0 .and. ierr /= -1 ) &
          CALL errore( ' s_psi_gamma_gpu ', ' cannot allocate buffer (ps_d) ', ABS(ierr) )
       !
       ps_d(1:nkb,1:m_max) = 0.D0
       !
       !   In becp=<vkb_i|psi_j> terms corresponding to atom na of type nt
       !   run from index i=indv_ijkb0(na)+1 to i=indv_ijkb0(na)+nh(nt)
       !
       DO nt = 1, nsp
          IF ( upf(nt)%tvanp ) THEN
             DO na = 1, nat
                IF ( ityp(na) == nt ) THEN
                   !
                   ! Next operation computes ps(l',i)=\sum_m qq(l,m) becp(m',i)
                   ! (l'=l+ijkb0, m'=m+ijkb0, indices run from 1 to nh(nt))
                   !
                   IF ( m_loc > 0 ) THEN
                      CALL DGEMM('N', 'N', nh(nt), m_loc, nh(nt), 1.0_dp, &
                                  qq_at_d(1,1,na), nhm, becp_d%r_d(indv_ijkb0(na)+1,1),&
                                  nkb, 0.0_dp, ps_d(indv_ijkb0(na)+1,1), nkb )
                   END IF
                END IF
             END DO
          END IF
       END DO
       !
       IF( becp_d%comm == mp_get_comm_null() ) THEN
          IF ( m == 1 ) THEN
             CALL cudaDGEMV( 'N', 2 * n, nkb, 1.D0, vkb_d, &
                  2 * lda, ps_d, 1, 1.D0, spsi_d, 1 )
          ELSE
             CALL cublasDGEMM( 'N', 'N', 2 * n, m, nkb, 1.D0, vkb_d, &
                  2 * lda, ps_d, nkb, 1.D0, spsi_d, 2 * lda )
          END IF
       ELSE
          !
          ! parallel block multiplication of vkb and ps
          !
          icur_blk = mype
          !
          DO icyc = 0, nproc - 1

             m_loc   = ldim_block( becp_d%nbnd , nproc, icur_blk )
             m_begin = gind_block( 1,  becp_d%nbnd, nproc, icur_blk )

             IF( ( m_begin + m_loc - 1 ) > m ) m_loc = m - m_begin + 1

             IF( m_loc > 0 ) THEN
                CALL cublasDGEMM( 'N', 'N', 2 * n, m_loc, nkb, 1.D0, vkb_d, &
                            2 * lda, ps_d, nkb, 1.D0, spsi_d( 1, m_begin ), 2 * lda )
             END IF

             ! block rotation
             !
             CALL mp_circular_shift_left( ps_d, icyc, becp_d%comm )

             icur_blk = icur_blk + 1
             IF( icur_blk == nproc ) icur_blk = 0

          END DO
          !
       END IF
       !
       CALL dev_buf%release_buffer(ps_d, ierr)
       !
       RETURN
       !
     END SUBROUTINE s_psi_gamma_gpu
     !
     !-----------------------------------------------------------------------
     SUBROUTINE s_psi_k_gpu()
       !-----------------------------------------------------------------------
       !
       ! ... k-points version
       !
       USE gbuffers,   ONLY : dev_buf
       USE uspp_gpum,  ONLY : qq_at_d, using_qq_at_d

       IMPLICIT NONE
       !
       ! ... local variables
       !
       INTEGER :: ikb, jkb, ih, jh, na, nt, ibnd, ierr
         ! counters
       COMPLEX(DP), POINTER :: qqc_d(:,:,:)
       COMPLEX(DP), POINTER :: ps_d(:,:)
#if defined(__CUDA)
       attributes(DEVICE) :: ps_d, qqc_d
#endif
         ! ps = product vkb and psi ; qqc = complex version of qq
       !
       CALL dev_buf%lock_buffer(ps_d, (/ nkb, m /), ierr)
       !
       IF( ierr /= 0 .and. ierr /= -1 ) &
          CALL errore( ' s_psi_k_gpu ', ' cannot allocate buffer (ps_d) ', ABS(ierr) )

       ! sync vkb if needed
       CALL using_vkb_d(0)
       CALL using_becp_k_d(0)
       CALL using_qq_at_d(0)
       !
       ps_d(1:nkb,1:m) = ( 0.D0, 0.D0 )
       !
       ! qq is real:  copy it into a complex variable to perform
       ! a zgemm - simple but sub-optimal solution
       !
       ! here we need to use qq_at_d instead of qq_nt_d otherwise real space augmentation brakes!
       !  qq_nt_d would be much faster and works for calculations without real space augmentation
       CALL dev_buf%lock_buffer( qqc_d, (/ nhm, nhm, nat/), ierr )
       IF( ierr /= 0 .and. ierr /= -1 ) &
          CALL errore( ' s_psi_k_gpu ', ' cannot allocate buffer (qqc_d) ', ABS(ierr) )

!$cuf kernel do(3) <<<*,*>>>
       DO na = 1, nat
          DO jh = 1, nhm
             DO ih = 1, nhm
                qqc_d(ih,jh, na) = CMPLX ( qq_at_d(ih,jh, na), 0.0_dp, KIND=dp )
             END DO
          END DO
       END DO

       DO nt = 1, nsp
          IF ( upf(nt)%tvanp ) THEN
             DO na = 1, nat
                IF ( ityp(na) == nt ) THEN
                   CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                        qqc_d(1,1,na), nhm, becp_d%k_d(indv_ijkb0(na)+1,1), nkb, &
                        (0.0_dp,0.0_dp), ps_d(indv_ijkb0(na)+1,1), nkb )
                   !
                END IF
             END DO
          END IF
       END DO
       CALL dev_buf%release_buffer(qqc_d, ierr)
       !
       IF ( m == 1 ) THEN
          !
          CALL ZGEMV( 'N', n, nkb, ( 1.D0, 0.D0 ), vkb_d, &
                      lda, ps_d, 1, ( 1.D0, 0.D0 ), spsi_d, 1 )
          !
       ELSE
          !
          CALL ZGEMM( 'N', 'N', n, m, nkb, ( 1.D0, 0.D0 ), vkb_d, &
                      lda, ps_d, nkb, ( 1.D0, 0.D0 ), spsi_d, lda )
          !
       END IF
       !
       CALL dev_buf%release_buffer(ps_d, ierr)
       !
       RETURN
       !
     END SUBROUTINE s_psi_k_gpu     
     !
     !
     !-----------------------------------------------------------------------
      SUBROUTINE s_psi_nc_gpu ( )
     !-----------------------------------------------------------------------
       !
       ! ... k-points noncolinear/spinorbit version
       !
       USE gbuffers,   ONLY : dev_buf
       USE uspp_gpum,  ONLY : qq_at_d, using_qq_at_d, qq_so_d, using_qq_so_d

       IMPLICIT NONE
       !
       !    here the local variables
       !
       INTEGER :: ikb, jkb, ih, jh, na, nt, ibnd, ipol, ierr
       ! counters
       COMPLEX(DP), POINTER :: ps_d(:,:,:)
       COMPLEX(DP), POINTER :: qqc_d(:,:,:)
#if defined(__CUDA)
       attributes(DEVICE) :: ps_d, qqc_d
#endif
       ! the product vkb and psi
       !
       ! sync if needed
       CALL using_vkb_d(0)
       CALL using_becp_nc_d(0)
       CALL using_indv_ijkb0(0)
       IF ( .not. lspinorb ) CALL using_qq_at_d(0)
       IF ( lspinorb ) CALL using_qq_so_d(0)

       CALL dev_buf%lock_buffer(ps_d, (/ nkb, npol, m /), ierr)
       IF( ierr /= 0 .and. ierr /= -1 ) &
          CALL errore( ' s_psi_nc_gpu ', ' cannot allocate buffer (ps_d) ', ABS(ierr) )

       ps_d(1:nkb,1:npol,1:m) = (0.D0,0.D0)
       !
       IF ( .NOT. lspinorb ) THEN
          CALL dev_buf%lock_buffer( qqc_d, (/ nhm, nhm, nat /), ierr )
          IF( ierr /= 0 .and. ierr /= -1 ) &
             CALL errore( ' s_psi_nc_gpu ', ' cannot allocate buffer (qqc_d) ', ABS(ierr) )
          ! Possibly convert only what's needed??
!$cuf kernel do(3) <<<*,*>>>
          DO na = 1, nat
             DO jh = 1, nhm
                DO ih = 1, nhm
                   qqc_d(ih, jh, na) = CMPLX ( qq_at_d(ih,jh, na), 0.0_dp, KIND=dp )
                END DO
             END DO
          END DO
       END IF
       !
       DO nt = 1, nsp
          !
          IF ( upf(nt)%tvanp ) THEN
             !
             IF ( .NOT. lspinorb ) THEN
                DO na = 1, nat
                   IF ( ityp(na) == nt ) THEN
                      DO ipol=1,npol
                         CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                              qqc_d(1,1, na), nhm, becp_d%nc_d(indv_ijkb0(na)+1,ipol,1), nkb*npol, &
                              (0.0_dp,0.0_dp), ps_d(indv_ijkb0(na)+1,ipol,1), nkb*npol )
                       END DO
                    END IF
                END DO
             ELSE
                DO na = 1, nat
                   IF ( ityp(na) == nt ) THEN
                      CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                           qq_so_d(1,1,1,nt), nhm, becp_d%nc_d(indv_ijkb0(na)+1,1,1), nkb*npol, &
                           (0.0_dp,0.0_dp), ps_d(indv_ijkb0(na)+1,1,1), nkb*npol )
                      CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                           qq_so_d(1,1,2,nt), nhm, becp_d%nc_d(indv_ijkb0(na)+1,2,1), nkb*npol, &
                           (1.0_dp,0.0_dp), ps_d(indv_ijkb0(na)+1,1,1), nkb*npol )
                      !
                      CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                           qq_so_d(1,1,3,nt), nhm, becp_d%nc_d(indv_ijkb0(na)+1,1,1), nkb*npol, &
                           (0.0_dp,0.0_dp), ps_d(indv_ijkb0(na)+1,2,1), nkb*npol )
                      CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                           qq_so_d(1,1,4,nt), nhm, becp_d%nc_d(indv_ijkb0(na)+1,2,1), nkb*npol, &
                           (1.0_dp,0.0_dp), ps_d(indv_ijkb0(na)+1,2,1), nkb*npol )
                    END IF
                END DO
             END IF
          END IF
       END DO
       IF ( .NOT. lspinorb ) CALL dev_buf%release_buffer(qqc_d, ierr)

       call ZGEMM ('N', 'N', n, m*npol, nkb, (1.d0, 0.d0) , vkb_d, &
          lda, ps_d, nkb, (1.d0, 0.d0) , spsi_d(1,1), lda)

       CALL dev_buf%release_buffer(ps_d, ierr)

       RETURN

    END SUBROUTINE s_psi_nc_gpu

END SUBROUTINE s_psi__gpu

