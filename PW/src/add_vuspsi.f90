!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE add_vuspsi( lda, n, m, hpsi )
  !----------------------------------------------------------------------------
  !! This routine applies the Ultra-Soft Hamiltonian to a
  !! vector psi and puts the result in hpsi. 
  !! It requires the products of psi with all beta functions
  !! in array becp(nkb,m) (calculated by calbec).
  !
  USE kinds,            ONLY: DP
  USE ions_base,        ONLY: nat, ntyp => nsp, ityp
  USE lsda_mod,         ONLY: current_spin
  USE control_flags,    ONLY: gamma_only
  USE noncollin_module
  USE uspp,             ONLY: vkb, nkb, deeq, deeq_nc, ofsbeta
  USE uspp_param,       ONLY: nh, nhm
  USE becmod,           ONLY: becp
  !
  IMPLICIT NONE
  !
  ! This is a unique interface to offload add_vuspsi with openmp
  ! INOUT  : hpsi is assumed on device
  ! Global : becp, deeq, deeq_nc are mapped to
  !          vkb is mapped in init_us_2
  !
  ! ... I/O variables
  !
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of arrays psi, spsi
  INTEGER, INTENT(IN) :: n
  !! true dimension of psi, spsi
  INTEGER, INTENT(IN) :: m
  !! number of states psi
  COMPLEX(DP), INTENT(INOUT) :: hpsi(lda*npol,m)
  !! V_US|psi> is added to hpsi
  !
  ! ... here the local variables
  !
  INTEGER :: jkb, ikb, ih, jh, na, nt, ibnd
  ! counters
  !
  !
  CALL start_clock( 'add_vuspsi' )  
  !
  IF ( gamma_only ) THEN
     !
     CALL add_vuspsi_gamma()
     !
  ELSE IF ( noncolin) THEN
     !
     CALL add_vuspsi_nc()
     !
  ELSE
     !
     CALL add_vuspsi_k()
     !
  ENDIF
  !
  CALL stop_clock( 'add_vuspsi' )  
  !
  RETURN
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     SUBROUTINE add_vuspsi_gamma()
       !-----------------------------------------------------------------------
       !! See comments inside
       !
       USE mp, ONLY: mp_get_comm_null, mp_circular_shift_left
       !
       IMPLICIT NONE
       !
       INTEGER, EXTERNAL :: ldim_block, gind_block
       REAL(DP), ALLOCATABLE :: ps (:,:)
       INTEGER :: ierr, i, j
       INTEGER :: nproc, mype, m_loc, m_begin, ibnd_loc, icyc, icur_blk, m_max
       !
       IF ( nkb == 0 ) RETURN
       !
       IF ( becp%comm == mp_get_comm_null() ) THEN
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
       ENDIF
       !
       ALLOCATE( ps (nkb,m_max), STAT=ierr )
       IF ( ierr /= 0 ) &
          CALL errore( ' add_vuspsi_gamma ', ' cannot allocate ps ', ABS(ierr) )
       !
#if defined(__OPENMP_GPU)
       !$omp target data map(to:becp%r,deeq) map(alloc:ps)
       !
       !$omp target teams distribute parallel do collapse(2)
       DO j = 1, m_max 
          DO i = 1, nkb
             ps(i,j) = 0.D0
          END DO
       END DO
#else
       ps(:,:) = 0.D0
#endif
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
                  CALL MYDGEMM2('N', 'N', nh(nt), m_loc, nh(nt), 1.0_dp, &
                           deeq(1,1,na,current_spin), nhm, &
                           becp%r(ofsbeta(na)+1,1), nkb, 0.0_dp, &
                               ps(ofsbeta(na)+1,1), nkb, .true. )
                ENDIF
                !
             ENDIF
             !
          ENDDO
          !
       ENDDO
       !
       IF( becp%comm == mp_get_comm_null() ) THEN
          !
          ! Normal case: hpsi(n,i) = \sum_l beta(n,l) ps(l,i) 
          ! (l runs from 1 to nkb)
          !
          CALL MYDGEMM2( 'N', 'N', ( 2 * n ), m, nkb, 1.D0, vkb, &
                   ( 2 * lda ), ps, nkb, 1.D0, hpsi, ( 2 * lda ), .true. )
       ELSE
          !
          ! parallel block multiplication of vkb and ps
          !
          icur_blk = mype
          !
          DO icyc = 0, nproc - 1
             !
             m_loc   = ldim_block( becp%nbnd , nproc, icur_blk )
             m_begin = gind_block( 1,  becp%nbnd, nproc, icur_blk )
             !
             IF( ( m_begin + m_loc - 1 ) > m ) m_loc = m - m_begin + 1
             !
             IF( m_loc > 0 ) THEN
                CALL MYDGEMM2( 'N', 'N', ( 2 * n ), m_loc, nkb, 1.D0, vkb, &
                   ( 2 * lda ), ps, nkb, 1.D0, hpsi( 1, m_begin ), ( 2 * lda ), .true. )
             ENDIF
             !
             ! block rotation
             !
             !$omp target update from(ps)
             CALL mp_circular_shift_left( ps, icyc, becp%comm )
             !$omp target update to(ps)
             !
             icur_blk = icur_blk + 1
             IF( icur_blk == nproc ) icur_blk = 0
             !
          ENDDO
       ENDIF
#if defined(__OPENMP_GPU)
       !$omp end target data
#endif
       !
       DEALLOCATE( ps )
       !
       RETURN
       !
     END SUBROUTINE add_vuspsi_gamma
     !
     !-----------------------------------------------------------------------
     SUBROUTINE add_vuspsi_k()
       !-----------------------------------------------------------------------
       !! See add_vuspsi_gamma for comments
       !
       IMPLICIT NONE
       !
       COMPLEX(DP), ALLOCATABLE :: ps(:,:), deeaux(:,:)
       INTEGER :: ierr, j, k, nhnt
       !
       IF ( nkb == 0 ) RETURN
       !
       ALLOCATE( ps(nkb,m), STAT=ierr )
       IF( ierr /= 0 ) &
          CALL errore( ' add_vuspsi_k ', ' cannot allocate ps ', ABS( ierr ) )
       !
#if defined(__OPENMP_GPU)
       ALLOCATE ( deeaux(nhm, nhm), STAT=ierr )
       IF( ierr /= 0 ) &
          CALL errore( ' add_vuspsi_k ', ' cannot allocate deeaux ', ABS( ierr ) )
       !
       !$omp target data map(alloc:ps,deeaux) map(to:deeq,becp%k)
#endif
       DO nt = 1, ntyp
          !
          IF ( nh(nt) == 0 ) CYCLE
#if !defined(__OPENMP_GPU)
          ALLOCATE( deeaux(nh(nt),nh(nt)) )
#else
          !
          nhnt = nh(nt)
#endif
          !
          DO na = 1, nat
             !
             IF ( ityp(na) == nt ) THEN
                !
                ! deeq is real: copy it into a complex variable to perform
                ! a zgemm - simple but sub-optimal solution
                !
#if !defined(__OPENMP_GPU)
                deeaux(:,:) = CMPLX(deeq(1:nh(nt),1:nh(nt),na,current_spin),&
                                    0.0_dp, KIND=dp )
                CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                           deeaux, nh(nt), becp%k(ofsbeta(na)+1,1), nkb, &
                          (0.0_dp, 0.0_dp), ps(ofsbeta(na)+1,1), nkb )
#else
                !$omp target teams distribute parallel do collapse(2)
                DO j = 1, nhnt
                   DO k = 1, nhnt
                      deeaux(k,j) = CMPLX(deeq(k,j,na,current_spin), 0.0_dp, KIND=DP )
                   END DO
                END DO
                !
                CALL MYZGEMM2('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                           deeaux, nhm, becp%k(ofsbeta(na)+1,1), nkb, &
                          (0.0_dp, 0.0_dp), ps(ofsbeta(na)+1,1), nkb, .true. )
#endif
                !
             ENDIF
             !
          ENDDO
#if !defined(__OPENMP_GPU)
          DEALLOCATE( deeaux )
#endif
          !
       ENDDO
       !

       CALL MYZGEMM2( 'N', 'N', n, m, nkb, ( 1.D0, 0.D0 ) , vkb, &
                   lda, ps, nkb, ( 1.D0, 0.D0 ) , hpsi, lda, .true. )
#if defined(__OPENMP_GPU)
       !$omp end target data
       DEALLOCATE (deeaux)
#endif
       DEALLOCATE( ps )
       !
       RETURN
       !
     END SUBROUTINE add_vuspsi_k     
     !  
     !-----------------------------------------------------------------------
     SUBROUTINE add_vuspsi_nc()
       !-----------------------------------------------------------------------
       !! See add_vuspsi_k for comments
       !
       IMPLICIT NONE
       !
       COMPLEX(DP), ALLOCATABLE :: ps(:,:,:)
       INTEGER :: ierr, ijkb0, i, j, k 
       !
       IF ( nkb == 0 ) RETURN
       !
       ALLOCATE( ps(  nkb,npol, m), STAT=ierr )
       IF( ierr /= 0 ) &
          CALL errore( ' add_vuspsi_nc ', ' error allocating ps ', ABS( ierr ) )
       !
#if defined(__OPENMP_GPU)
       !$omp target data map(alloc:ps) map(to:becp%nc,deeq_nc)
       !$omp target teams distribute parallel do collapse(3)
       DO k = 1, m
         DO j = 1, npol
           DO i = 1, nkb
              ps(i,j,k) = (0.d0, 0.d0)
           END DO
         END DO
       END DO
#else
       ps(:,:,:) = (0.d0, 0.d0)
#endif
       !
       DO nt = 1, ntyp
          !
          IF ( nh(nt) == 0 ) CYCLE
          DO na = 1, nat
             !
             IF ( ityp(na) == nt ) THEN
                !
#if defined(__OPENMP_GPU)
                CALL MYZGEMM2('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                           deeq_nc(1,1,na,1), nhm, becp%nc(ofsbeta(na)+1,1,1), 2*nkb, &
                          (0.0_dp, 0.0_dp), ps(ofsbeta(na)+1,1,1), 2*nkb, .true. )

                CALL MYZGEMM2('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                           deeq_nc(1,1,na,2), nhm, becp%nc(ofsbeta(na)+1,2,1), 2*nkb, &
                          (1.0_dp, 0.0_dp), ps(ofsbeta(na)+1,1,1), 2*nkb, .true. )


                CALL MYZGEMM2('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                           deeq_nc(1,1,na,3), nhm, becp%nc(ofsbeta(na)+1,1,1), 2*nkb, &
                          (0.0_dp, 0.0_dp), ps(ofsbeta(na)+1,2,1), 2*nkb, .true. )


                CALL MYZGEMM2('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                           deeq_nc(1,1,na,4), nhm, becp%nc(ofsbeta(na)+1,2,1), 2*nkb, &
                          (1.0_dp, 0.0_dp), ps(ofsbeta(na)+1,2,1), 2*nkb, .true. )
#else
                ! The following is left to easily merge the two routines in the future
                DO ibnd = 1, m
                   !
                   DO jh = 1, nh(nt)
                      !
                      jkb = ofsbeta(na) + jh
                      !
                      DO ih = 1, nh(nt)
                         !
                         ikb = ofsbeta(na) + ih
                         !
                         ps(ikb,1,ibnd) = ps(ikb,1,ibnd) +    & 
                              deeq_nc(ih,jh,na,1)*becp%nc(jkb,1,ibnd)+ & 
                              deeq_nc(ih,jh,na,2)*becp%nc(jkb,2,ibnd) 
                         ps(ikb,2,ibnd) = ps(ikb,2,ibnd)  +   & 
                              deeq_nc(ih,jh,na,3)*becp%nc(jkb,1,ibnd)+&
                              deeq_nc(ih,jh,na,4)*becp%nc(jkb,2,ibnd) 
                         !
                      ENDDO
                      !
                   ENDDO
                   !
                ENDDO
                ! The aboce is left to easily merge the two routines in the future
#endif
                !
             ENDIF
             !
          ENDDO
          !
       ENDDO
       !
       CALL MYZGEMM2('N', 'N', n, m*npol, nkb, ( 1.D0, 0.D0 ) , vkb, &
                   lda, ps, nkb, ( 1.D0, 0.D0 ) , hpsi, lda, .true. )
       !
#if defined(__OPENMP_GPU)
       !$omp end target data
#endif
       DEALLOCATE( ps )
       !
       RETURN
       !
     END SUBROUTINE add_vuspsi_nc
     !  
     !  
END SUBROUTINE add_vuspsi
