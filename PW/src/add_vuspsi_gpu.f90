!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE add_vuspsi_gpu( lda, n, m, hpsi_d )
  !----------------------------------------------------------------------------
  !! This routine applies the Ultra-Soft Hamiltonian to a
  !! vector psi and puts the result in hpsi. 
  !! It requires the products of psi with all beta functions
  !! in array becp(nkb,m) (calculated by calbec).
  !
  USE kinds,           ONLY: DP
  USE ions_base,       ONLY: nat, ntyp => nsp, ityp
  USE lsda_mod,        ONLY: current_spin
  USE control_flags,   ONLY: gamma_only
  USE noncollin_module
  USE uspp,            ONLY: ofsbeta, nkb, vkb, deeq, deeq_nc
  USE uspp_param,      ONLY: nh, nhm
  USE becmod,          ONLY: becp
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of arrays psi, spsi
  INTEGER, INTENT(IN) :: n
  !! true dimension of psi, spsi
  INTEGER, INTENT(IN) :: m
  !! number of states psi
  COMPLEX(DP), INTENT(INOUT) :: hpsi_d(lda*npol,m)
  !! V_US|psi> is added to hpsi
  !
  ! ... here the local variables
  !
#if defined(__CUDA)
  attributes(DEVICE) :: hpsi_d
  !
  !
  CALL start_clock_gpu( 'add_vuspsi' )  
  !
  IF ( gamma_only ) THEN
     !
     CALL add_vuspsi_gamma_gpu()
     !
  ELSE IF ( noncolin) THEN
     !
     CALL add_vuspsi_nc_gpu ()
     !
  ELSE
     !
     CALL add_vuspsi_k_gpu()
     !
  END IF
  !
  CALL stop_clock_gpu( 'add_vuspsi' )  
  !
  RETURN
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     SUBROUTINE add_vuspsi_gamma_gpu()
       !-----------------------------------------------------------------------
       !! See comments inside
       !
       USE mp, ONLY: mp_get_comm_null, mp_circular_shift_left
       USE device_fbuff_m, ONLY : dev_buf
       !
#if defined(__CUDA)
       USE cudafor
       USE cublas
#endif
       !
       IMPLICIT NONE
       INTEGER, EXTERNAL :: ldim_block, gind_block
       REAL(DP), POINTER :: ps_d (:,:)
       INTEGER :: ierr
       INTEGER :: nproc, mype, m_loc, m_begin, ibnd_loc, icyc, icur_blk, m_max
       ! counters
       INTEGER :: jkb, ikb, ih, jh, na, nt, ibnd
       !
#if defined(__CUDA)
       attributes(device) :: ps_d
#endif
       REAL(DP), ALLOCATABLE :: becp_r(:,:)
       !$acc declare device_resident(becp_r)
       !
       IF ( nkb == 0 ) RETURN
       !
       IF( becp%comm == mp_get_comm_null() ) THEN
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
          nproc   = becp%nproc
          mype    = becp%mype
          m_loc   = becp%nbnd_loc
          m_begin = becp%ibnd_begin
          m_max   = SIZE( becp%r, 2 )
          IF( ( m_begin + m_loc - 1 ) > m ) m_loc = m - m_begin + 1
       END IF
       !
       CALL dev_buf%lock_buffer(ps_d, (/ nkb, m_max /), ierr ) !ALLOCATE (ps_d (nkb,m_max), STAT=ierr )
       !
       IF( ierr /= 0 ) &
          CALL errore( ' add_vuspsi_gamma ', ' cannot allocate ps_d ', ABS(ierr) )
       !
       ps_d(:,:) = 0.D0
       !
       ALLOCATE( becp_r(size(becp%r,1),size(becp%r,2)), stat=ierr )
       IF( ierr /= 0 ) &
          CALL errore( ' add_vuspsi_gamma ', ' cannot allocate becp_r', ABS(ierr) )
       !$acc kernels
       becp_r = becp%r
       !$acc end kernels
       !
       !   In becp=<vkb_i|psi_j> terms corresponding to atom na of type nt
       !   run from index i=ofsbeta(na)+1 to i=ofsbeta(na)+nh(nt)
       !
       DO nt = 1, ntyp
          !
          IF ( nh(nt) == 0 ) CYCLE
          DO na = 1, nat
             !
             IF ( ityp(na) == nt ) THEN
                !
                ! Next operation computes ps(l',i) = \sum_m deeq(l,m) becp(m',i)
                ! (l'=l+ijkb0, m'=m+ijkb0, indices run from 1 to nh(nt))
                !
                IF ( m_loc > 0 ) THEN
                  !$acc host_data use_device(deeq,becp_r)
                  CALL DGEMM('N', 'N', nh(nt), m_loc, nh(nt), 1.0_dp, &
                           deeq(1,1,na,current_spin), nhm, &
                           becp_r(ofsbeta(na)+1,1), nkb, 0.0_dp, &
                               ps_d(ofsbeta(na)+1,1), nkb )
                  !$acc end host_data
                END IF
                !
             END IF
             !
          END DO
          !
       END DO
       !
       IF( becp%comm == mp_get_comm_null() ) THEN
          !
          ! Normal case: hpsi(n,i) = \sum_l beta(n,l) ps(l,i) 
          ! (l runs from 1 to nkb)
          !
          !$acc host_data use_device(vkb)
          CALL DGEMM( 'N', 'N', ( 2 * n ), m, nkb, 1.D0, vkb, &
                   ( 2 * lda ), ps_d, nkb, 1.D0, hpsi_d, ( 2 * lda ) )
          !$acc end host_data
       ELSE
          !
          ! parallel block multiplication of vkb and ps
          !
          icur_blk = mype
          !
          DO icyc = 0, nproc - 1

             m_loc   = ldim_block( becp%nbnd , nproc, icur_blk )
             m_begin = gind_block( 1,  becp%nbnd, nproc, icur_blk )

             IF( ( m_begin + m_loc - 1 ) > m ) m_loc = m - m_begin + 1

             IF( m_loc > 0 ) THEN
                !$acc host_data use_device(vkb)
                CALL DGEMM( 'N', 'N', ( 2 * n ), m_loc, nkb, 1.D0, vkb, &
                   ( 2 * lda ), ps_d, nkb, 1.D0, hpsi_d( 1, m_begin ), ( 2 * lda ) )
                !$acc end host_data
             ENDIF

             ! block rotation
             !
             CALL mp_circular_shift_left( ps_d, icyc, becp%comm )

             icur_blk = icur_blk + 1
             IF( icur_blk == nproc ) icur_blk = 0

          ENDDO
       ENDIF
       !
       CALL dev_buf%release_buffer(ps_d, ierr) ! DEALLOCATE (ps_d)
       !
       DEALLOCATE( becp_r )
       !
       RETURN
       !
     END SUBROUTINE add_vuspsi_gamma_gpu
     !
     !-----------------------------------------------------------------------
     SUBROUTINE add_vuspsi_k_gpu()
       !-----------------------------------------------------------------------
       !! See add_vuspsi_gamma for comments
       !
#if defined(__CUDA)
       USE cudafor
       USE cublas
#endif
       USE device_fbuff_m, ONLY : dev_buf
       !
       IMPLICIT NONE
       COMPLEX(DP), POINTER :: ps_d (:,:), deeaux_d (:,:)
       INTEGER :: ierr
       ! counters
       INTEGER :: i, j, k, jkb, ikb, ih, jh, na, nt, ibnd, nhnt

#if defined(__CUDA)
       ATTRIBUTES( DEVICE ) :: ps_d, deeaux_d
#endif
       COMPLEX(DP), ALLOCATABLE :: becp_k(:,:)
       !$acc declare device_resident(becp_k)
       !
       IF ( nkb == 0 ) RETURN
       !
       CALL dev_buf%lock_buffer(ps_d, (/ nkb,m /), ierr ) ! ALLOCATE (ps_d (nkb,m), STAT=ierr )
       IF( ierr /= 0 ) &
          CALL errore( ' add_vuspsi_k ', ' cannot allocate ps_d ', ABS( ierr ) )
       !
       CALL dev_buf%lock_buffer(deeaux_d, (/ nhm, nhm /), ierr ) !ALLOCATE ( deeaux_d(nhm, nhm) )
       IF( ierr /= 0 ) &
          CALL errore( ' add_vuspsi_k ', ' cannot allocate deeaux_d ', ABS( ierr ) )
       !
       ALLOCATE( becp_k(size(becp%k,1), size(becp%k,2) ), stat=ierr )
       IF( ierr /= 0 ) &
          CALL errore( ' add_vuspsi_k ', ' cannot allocate becp_k ', ABS( ierr ) )
       !$acc kernels
       becp_k = becp%k
       !$acc end kernels
       !
       DO nt = 1, ntyp
          !
          IF ( nh(nt) == 0 ) CYCLE
          !
          nhnt = nh(nt)
          !
          DO na = 1, nat
             !
             IF ( ityp(na) == nt ) THEN
                !
                ! deeq is real: copy it into a complex variable to perform
                ! a zgemm - simple but sub-optimal solution
                !
                !deeaux_d(:,:) = CMPLX(deeq(1:nh(nt),1:nh(nt),na,current_spin), 0.0_dp, KIND=dp )
                !
                !$acc parallel loop collapse(2) present(deeq)
                DO j = 1, nhnt
                   DO k = 1, nhnt
                      deeaux_d(k,j) = CMPLX(deeq(k,j,na,current_spin), 0.0_dp, KIND=DP )
                   END DO
                END DO
                !
                !$acc host_data use_device(becp_k)
                CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                           deeaux_d, nhm, becp_k(ofsbeta(na)+1,1), nkb, &
                          (0.0_dp, 0.0_dp), ps_d(ofsbeta(na)+1,1), nkb )
                !$acc end host_data
                !
             END IF
             !
          END DO
          !
       END DO
       CALL dev_buf%release_buffer(deeaux_d, ierr) ! DEALLOCATE (deeaux_d)
       !
       !$acc host_data use_device(vkb)
       CALL ZGEMM( 'N', 'N', n, m, nkb, ( 1.D0, 0.D0 ) , vkb, &
                   lda, ps_d, nkb, ( 1.D0, 0.D0 ) , hpsi_d, lda )
       !$acc end host_data
       !
       CALL dev_buf%release_buffer(ps_d, ierr) !DEALLOCATE (ps_d)
       !
       DEALLOCATE( becp_k )
       !
       RETURN
       !
     END SUBROUTINE add_vuspsi_k_gpu
     !  
     !-----------------------------------------------------------------------
     SUBROUTINE add_vuspsi_nc_gpu()
       !-----------------------------------------------------------------------
       !! See add_vuspsi_k for comments
       !
#if defined(__CUDA)
       USE cudafor
       USE cublas
#endif
       USE device_fbuff_m,      ONLY : dev_buf
       IMPLICIT NONE
       COMPLEX(DP), POINTER :: ps_d (:,:,:)
       INTEGER :: ierr
       ! counters
       INTEGER :: na, nt, ibnd

#if defined(__CUDA)
       ATTRIBUTES( DEVICE ) :: ps_d
#endif
       COMPLEX(DP), ALLOCATABLE :: becp_nc(:,:,:)
       !$acc declare device_resident(becp_nc)
       !
       IF ( nkb == 0 ) RETURN
       !
       ! ALLOCATE (ps_d( nkb, npol, m), STAT=ierr )
       CALL dev_buf%lock_buffer(ps_d, (/ nkb, npol, m /), ierr )
       !
       IF( ierr /= 0 ) &
          CALL errore( ' add_vuspsi_nc ', ' error allocating ps_d ', ABS( ierr ) )
       !
       ALLOCATE( becp_nc(size(becp%nc,1),size(becp%nc,2),size(becp%nc,3)), stat=ierr )
       IF( ierr /= 0 ) &
          CALL errore( ' add_vuspsi_nc ', ' error allocating becp_nc ', ABS( ierr ) )
       !$acc kernels
       becp_nc = becp%nc
       !$acc end kernels
       !
       !  OPTIMIZE HERE: possibly streamline
       !
       DO nt = 1, ntyp
          !
          IF ( nh(nt) == 0 ) CYCLE
          DO na = 1, nat
             !
             IF ( ityp(na) == nt ) THEN
                !
                !$acc host_data use_device(deeq_nc,becp_nc)
                CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                           deeq_nc(1,1,na,1), nhm, becp_nc(ofsbeta(na)+1,1,1), 2*nkb, &
                          (0.0_dp, 0.0_dp), ps_d(ofsbeta(na)+1,1,1), 2*nkb )

                CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                           deeq_nc(1,1,na,2), nhm, becp_nc(ofsbeta(na)+1,2,1), 2*nkb, &
                          (1.0_dp, 0.0_dp), ps_d(ofsbeta(na)+1,1,1), 2*nkb )


                CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                           deeq_nc(1,1,na,3), nhm, becp_nc(ofsbeta(na)+1,1,1), 2*nkb, &
                          (0.0_dp, 0.0_dp), ps_d(ofsbeta(na)+1,2,1), 2*nkb )


                CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                           deeq_nc(1,1,na,4), nhm, becp_nc(ofsbeta(na)+1,2,1), 2*nkb, &
                          (1.0_dp, 0.0_dp), ps_d(ofsbeta(na)+1,2,1), 2*nkb )
                !$acc end host_data
                !
!                DO ibnd = 1, m
!                   !
!                   DO jh = 1, nh(nt)
!                      !
!!$acc parallel loop present(deeq_nc)
!                      DO ih = 1, nh(nt)
!                         !
!                         ikb = ijkb0 + ih
!                         jkb = ijkb0 + jh   
!                         becpup_jkb = becp_nc_d(jkb,1,ibnd)
!                         becpdn_jkb = becp_nc_d(jkb,2,ibnd)
!                         !
!                         ps_d(ikb,1,ibnd) = ps_d(ikb,1,ibnd) +   & 
!                              deeq_nc(ih,jh,na,1)*becpup_jkb + & 
!                              deeq_nc(ih,jh,na,2)*becpdn_jkb
!                         ps_d(ikb,2,ibnd) = ps_d(ikb,2,ibnd) +   & 
!                              deeq_nc(ih,jh,na,3)*becpup_jkb + &
!                              deeq_nc(ih,jh,na,4)*becpdn_jkb
!                         !
!                      END DO
!                      !
!                   END DO
!                   !
!                END DO
                !
             END IF
             !
          END DO
          !
       END DO
       !
       !$acc host_data use_device(vkb)
       call ZGEMM ('N', 'N', n, m*npol, nkb, ( 1.D0, 0.D0 ) , vkb, &
                   lda, ps_d, nkb, ( 1.D0, 0.D0 ) , hpsi_d, lda )
       !$acc end host_data
       !
       CALL dev_buf%release_buffer(ps_d, ierr ) ! DEALLOCATE (ps_d)
       !
       DEALLOCATE( becp_nc )
       !
       RETURN
       !
     END SUBROUTINE add_vuspsi_nc_gpu
#endif
     !
END SUBROUTINE add_vuspsi_gpu
