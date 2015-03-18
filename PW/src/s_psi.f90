!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE s_psi( lda, n, m, psi, spsi )
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
  USE kinds,      ONLY : DP
  USE becmod,     ONLY : becp
  USE uspp,       ONLY : vkb, nkb, okvan, qq, qq_so, indv_ijkb0
  USE spin_orb,   ONLY : lspinorb
  USE uspp_param, ONLY : upf, nh, nhm
  USE ions_base,  ONLY : nat, nsp, ityp
  USE control_flags,    ONLY: gamma_only 
  USE noncollin_module, ONLY: npol, noncolin
  USE realus,     ONLY :  real_space, fft_orbital_gamma, initialisation_level,&
                          bfft_orbital_gamma, calbec_rs_gamma, s_psir_gamma
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda, n, m
  COMPLEX(DP), INTENT(IN) :: psi(lda*npol,m)
  COMPLEX(DP), INTENT(OUT)::spsi(lda*npol,m)
  !
  INTEGER :: ibnd
  !
  ! ... initialize  spsi
  !
  spsi = psi
  !
  IF ( nkb == 0 .OR. .NOT. okvan ) RETURN
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
           CALL  fft_orbital_gamma(psi,ibnd,m) 
           CALL s_psir_gamma(ibnd,m)
           CALL bfft_orbital_gamma(spsi,ibnd,m)
        END DO
        !
     ELSE
        !
        CALL s_psi_gamma()
        !
     END IF
     !
  ELSE IF ( noncolin ) THEN
     !
     CALL s_psi_nc()
     !
  ELSE 
     !
     CALL s_psi_k()
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
     SUBROUTINE s_psi_gamma()
       !-----------------------------------------------------------------------
       ! 
       ! ... gamma version
       !
       USE mp, ONLY: mp_get_comm_null, mp_circular_shift_left
       !
       IMPLICIT NONE  
       !
       ! ... here the local variables
       !
       INTEGER :: ikb, jkb, ih, jh, na, nt, ibnd, ierr
         ! counters
       INTEGER :: nproc, mype, m_loc, m_begin, ibnd_loc, icyc, icur_blk, m_max
         ! data distribution indexes
       INTEGER, EXTERNAL :: ldim_block, lind_block, gind_block
         ! data distribution functions
       REAL(DP), ALLOCATABLE :: ps(:,:)
         ! the product vkb and psi
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
       ALLOCATE( ps( nkb, m_max ), STAT=ierr )
       IF( ierr /= 0 ) &
          CALL errore( ' s_psi_gamma ', ' cannot allocate memory (ps) ', ABS(ierr) )
       !    
       ps(:,:) = 0.D0
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
                   CALL DGEMM('N', 'N', nh(nt), m_loc, nh(nt), 1.0_dp, &
                        qq(1,1,nt), nhm, &
                        becp%r(indv_ijkb0(na)+1,1), nkb, 0.0_dp, &
                            ps(indv_ijkb0(na)+1,1), nkb )
                END IF
             END DO
          END IF
       END DO
       !
       IF( becp%comm == mp_get_comm_null() ) THEN
          IF ( m == 1 ) THEN
             CALL DGEMV( 'N', 2 * n, nkb, 1.D0, vkb, &
                  2 * lda, ps, 1, 1.D0, spsi, 1 )
          ELSE
             CALL DGEMM( 'N', 'N', 2 * n, m, nkb, 1.D0, vkb, &
                  2 * lda, ps, nkb, 1.D0, spsi, 2 * lda )
          END IF
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
                CALL DGEMM( 'N', 'N', 2 * n, m_loc, nkb, 1.D0, vkb, &
                            2 * lda, ps, nkb, 1.D0, spsi( 1, m_begin ), 2 * lda )
             END IF

             ! block rotation
             !
             CALL mp_circular_shift_left( ps, icyc, becp%comm )

             icur_blk = icur_blk + 1
             IF( icur_blk == nproc ) icur_blk = 0

          END DO
          !
       END IF
       !
       DEALLOCATE( ps ) 
       !
       RETURN
       !
     END SUBROUTINE s_psi_gamma
     !
     !-----------------------------------------------------------------------
     SUBROUTINE s_psi_k()
       !-----------------------------------------------------------------------
       !
       ! ... k-points version
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       INTEGER :: ikb, jkb, ih, jh, na, nt, ibnd, ierr
         ! counters
       COMPLEX(DP), ALLOCATABLE :: ps(:,:), qqc(:,:)
         ! ps = product vkb and psi ; qqc = complex version of qq
       !
       ALLOCATE( ps( nkb, m ), STAT=ierr )    
       IF( ierr /= 0 ) &
          CALL errore( ' s_psi_k ', ' cannot allocate memory (ps) ', ABS(ierr) )
       !
       ps(:,:) = ( 0.D0, 0.D0 )
       !
       DO nt = 1, nsp
          IF ( upf(nt)%tvanp ) THEN
             ! qq is real:  copy it into a complex variable to perform
             ! a zgemm - simple but sub-optimal solution
             ALLOCATE( qqc(nh(nt),nh(nt)) )
             qqc(:,:) = CMPLX ( qq(1:nh(nt),1:nh(nt),nt), 0.0_dp, KIND=dp )
             DO na = 1, nat
                IF ( ityp(na) == nt ) THEN
                   CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                        qqc, nh(nt), becp%k(indv_ijkb0(na)+1,1), nkb, &
                        (0.0_dp,0.0_dp), ps(indv_ijkb0(na)+1,1), nkb )
                   !
                END IF
             END DO
             DEALLOCATE (qqc)
          END IF
       END DO
       !
       IF ( m == 1 ) THEN
          !
          CALL ZGEMV( 'N', n, nkb, ( 1.D0, 0.D0 ), vkb, &
                      lda, ps, 1, ( 1.D0, 0.D0 ), spsi, 1 )
          !
       ELSE
          !
          CALL ZGEMM( 'N', 'N', n, m, nkb, ( 1.D0, 0.D0 ), vkb, &
                      lda, ps, nkb, ( 1.D0, 0.D0 ), spsi, lda )
          !
       END IF
       !
       DEALLOCATE( ps )
       !
       RETURN
       !
     END SUBROUTINE s_psi_k     
     !
     !
     !-----------------------------------------------------------------------
      SUBROUTINE s_psi_nc ( )
     !-----------------------------------------------------------------------
       !
       ! ... k-points noncolinear/spinorbit version
       !
       IMPLICIT NONE
       !
       !    here the local variables
       !
       INTEGER :: ikb, jkb, ih, jh, na, nt, ijkb0, ibnd, ipol, ierr
       ! counters
       COMPLEX (DP), ALLOCATABLE :: ps (:,:,:), qqc(:,:)
       ! the product vkb and psi
       !
       ALLOCATE (ps(nkb,npol,m),STAT=ierr)    
       IF( ierr /= 0 ) &
          CALL errore( ' s_psi_nc ', ' cannot allocate memory (ps) ', ABS(ierr) )
       ps(:,:,:) = (0.D0,0.D0)
       !
       ijkb0 = 1
       do nt = 1, nsp
          if ( upf(nt)%tvanp ) then
             !
             IF ( .NOT. lspinorb ) THEN
                ! qq is real:  copy it into a complex variable to perform
                ! a zgemm - simple but sub-optimal solution
                ALLOCATE( qqc(nh(nt),nh(nt)) )
                qqc(:,:) = CMPLX ( qq(1:nh(nt),1:nh(nt),nt), 0.0_dp, KIND=dp )
                DO na = 1, nat
                   IF ( ityp(na) == nt ) THEN
                      !  ps(ikb,ipol,ibnd) = ps(ikb,ipol,ibnd) + &
                      !     qq(ih,jh,nt)*becp%nc(jkb,ipol,ibnd)
                      DO ipol=1,npol
                         CALL ZGEMM('N','N', nh(nt),m,nh(nt), (1.0_dp,0.0_dp),&
                              qqc, nh(nt), becp%nc(indv_ijkb0(na)+1,ipol,1),  &
                              2*nkb,(0.0_dp,0.0_dp), &
                              ps(indv_ijkb0(na)+1,ipol,1), 2*nkb )
                      END DO
                      !
                   END IF
                END DO
                DEALLOCATE (qqc)
             ELSE
                DO na = 1, nat
                   IF (ityp (na) == nt) THEN
                      !
                      ! The quantities to be computed are:
                      ! ps(ikb,1,ibnd)=ps(ikb,1,ibnd) + &
                      !     qq_so(ih,jh,1,nt)*becp%nc(jkb,1,ibnd)+ &
                      !     qq_so(ih,jh,2,nt)*becp%nc(jkb,2,ibnd)
                      ! ps(ikb,2,ibnd)=ps(ikb,2,ibnd) + &
                      !      qq_so(ih,jh,3,nt)*becp%nc(jkb,1,ibnd)+ &
                      !      qq_so(ih,jh,4,nt)*becp%nc(jkb,2,ibnd)
                      ! Not sure this is the best solution to use zgemm, though:
                      !
                      ijkb0 = indv_ijkb0(na)+1
                      CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                           qq_so(1,1,1,nt), nh(nt), becp%nc(ijkb0,1,1),2*nkb,&
                           (0.0_dp,0.0_dp), ps(ijkb0,1,1), 2*nkb )
                      CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                           qq_so(1,1,2,nt), nh(nt), becp%nc(ijkb0,2,1),2*nkb,&
                           (1.0_dp,0.0_dp), ps(ijkb0,1,1), 2*nkb )
                      CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                           qq_so(1,1,3,nt), nh(nt), becp%nc(ijkb0,1,1),2*nkb,&
                           (0.0_dp,0.0_dp), ps(ijkb0,2,1), 2*nkb )
                      CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                           qq_so(1,1,4,nt), nh(nt), becp%nc(ijkb0,2,1),2*nkb,&
                           (1.0_dp,0.0_dp), ps(ijkb0,2,1), 2*nkb )
                   ENDIF
                END DO
             END IF
          ENDIF
       END DO

       call ZGEMM ('N', 'N', n, m*npol, nkb, (1.d0, 0.d0) , vkb, &
          lda, ps, nkb, (1.d0, 0.d0) , spsi(1,1), lda)

       DEALLOCATE(ps)

       RETURN

    END SUBROUTINE s_psi_nc

END SUBROUTINE s_psi
