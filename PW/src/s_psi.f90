!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
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
  USE uspp,       ONLY : vkb, nkb, qq, okvan
  USE uspp_param, ONLY : upf, nh 
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
       USE becmod, ONLY : bec_type, becp
       USE mp, ONLY: mp_get_comm_null, mp_circular_shift_left
       !
       IMPLICIT NONE  
       !
       ! ... here the local variables
       !
       INTEGER :: ikb, jkb, ih, jh, na, nt, ijkb0, ibnd, ierr
         ! counters
       INTEGER :: nproc, mype, m_loc, m_begin, ibnd_loc, icyc, icur_blk, m_max
         ! data distribution indexes
       INTEGER, EXTERNAL :: ldim_block, lind_block, gind_block
         ! data distribution functions
       REAL(DP), ALLOCATABLE :: ps(:,:)
         ! the product vkb and psi
       !
       m_loc   = m
       m_begin = 1
       m_max   = m
       nproc   = 1
       mype    = 0
       !
       IF( becp%comm /= mp_get_comm_null() ) THEN
          nproc   = becp%nproc
          mype    = becp%mype
          m_loc   = becp%nbnd_loc
          m_begin = becp%ibnd_begin
          m_max   = SIZE(becp%r,2)
          IF( ( m_begin + m_loc - 1 ) > m ) m_loc = m - m_begin + 1
       END IF
       !
       ALLOCATE( ps( nkb, m_max ), STAT=ierr )
       IF( ierr /= 0 ) &
          CALL errore( ' s_psi_gamma ', ' cannot allocate memory (ps) ', ABS(ierr) )
       !    
       ps(:,:) = 0.D0
       !
       ijkb0 = 0
       DO nt = 1, nsp
          IF ( upf(nt)%tvanp ) THEN
             DO na = 1, nat
                IF ( ityp(na) == nt ) THEN
                   DO ibnd_loc = 1, m_loc
                      DO jh = 1, nh(nt)
                         jkb = ijkb0 + jh
                         DO ih = 1, nh(nt)
                            ikb = ijkb0 + ih
                            ps(ikb,ibnd_loc) = ps(ikb,ibnd_loc) + &
                                           qq(ih,jh,nt) * becp%r(jkb,ibnd_loc)
                         END DO
                      END DO
                   END DO
                   ijkb0 = ijkb0 + nh(nt)
                END IF
             END DO
          ELSE
             DO na = 1, nat
                IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
             END DO
          END IF
       END DO
       !
       IF( becp%comm /= mp_get_comm_null() ) THEN
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
       ELSE IF ( m == 1 ) THEN
          !
          CALL DGEMV( 'N', 2 * n, nkb, 1.D0, vkb, &
                      2 * lda, ps, 1, 1.D0, spsi, 1 )
          !
       ELSE
          !
          CALL DGEMM( 'N', 'N', 2 * n, m, nkb, 1.D0, vkb, &
                      2 * lda, ps, nkb, 1.D0, spsi, 2 * lda )
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
       USE becmod,  ONLY : becp
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       INTEGER :: ikb, jkb, ih, jh, na, nt, ijkb0, ibnd, ierr
         ! counters
       COMPLEX(DP), ALLOCATABLE :: ps(:,:)
         ! the product vkb and psi
       !
       ALLOCATE( ps( nkb, m ), STAT=ierr )    
       IF( ierr /= 0 ) &
          CALL errore( ' s_psi_k ', ' cannot allocate memory (ps) ', ABS(ierr) )
       !
       ps(:,:) = ( 0.D0, 0.D0 )
       !
       ijkb0 = 0
       DO nt = 1, nsp
          IF ( upf(nt)%tvanp ) THEN
             DO na = 1, nat
                IF ( ityp(na) == nt ) THEN
                   DO ibnd = 1, m
                      DO jh = 1, nh(nt)
                         jkb = ijkb0 + jh
                         DO ih = 1, nh(nt)
                            ikb = ijkb0 + ih
                            ps(ikb,ibnd) = ps(ikb,ibnd) + &
                                           qq(ih,jh,nt) * becp%k(jkb,ibnd)
                         END DO
                      END DO
                   END DO
                   ijkb0 = ijkb0 + nh(nt)
                END IF
             END DO
          ELSE
             DO na = 1, nat
                IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
             END DO
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
       USE uspp,   ONLY: qq_so
       USE becmod, ONLY: bec_type, becp
       USE spin_orb, ONLY: lspinorb
       IMPLICIT NONE
       !
       !    here the local variables
       !
       INTEGER :: ikb, jkb, ih, jh, na, nt, ijkb0, ibnd, ipol, ierr
       ! counters
       COMPLEX (DP), ALLOCATABLE :: ps (:,:,:)
       ! the product vkb and psi
       !
       ALLOCATE (ps(nkb,npol,m),STAT=ierr)    
       IF( ierr /= 0 ) &
          CALL errore( ' s_psi_nc ', ' cannot allocate memory (ps) ', ABS(ierr) )
       ps(:,:,:) = (0.D0,0.D0)
       !
       ijkb0 = 0
       do nt = 1, nsp
          if ( upf(nt)%tvanp ) then
             do na = 1, nat
                if (ityp (na) == nt) then
                   do ih = 1,nh(nt)
                      ikb = ijkb0 + ih
                      do ibnd = 1, m
                         do jh = 1, nh (nt)
                            jkb = ijkb0 + jh
                            if (lspinorb) then
                               ps(ikb,1,ibnd)=ps(ikb,1,ibnd) + &
                                 qq_so(ih,jh,1,nt)*becp%nc(jkb,1,ibnd)+ &
                                 qq_so(ih,jh,2,nt)*becp%nc(jkb,2,ibnd)
                               ps(ikb,2,ibnd)=ps(ikb,2,ibnd) + &
                                 qq_so(ih,jh,3,nt)*becp%nc(jkb,1,ibnd)+ &
                                 qq_so(ih,jh,4,nt)*becp%nc(jkb,2,ibnd)
                            else
                               do ipol=1,npol
                                  ps(ikb,ipol,ibnd)=ps(ikb,ipol,ibnd) + &
                                        qq(ih,jh,nt)*becp%nc(jkb,ipol,ibnd)
                               enddo
                            endif
                         enddo
                      enddo
                   enddo
                   ijkb0 = ijkb0 + nh (nt)
                endif
             enddo
          else
             do na = 1, nat
                if (ityp (na) == nt) ijkb0 = ijkb0 + nh (nt)
             enddo
          endif
       enddo

       call ZGEMM ('N', 'N', n, m*npol, nkb, (1.d0, 0.d0) , vkb, &
          lda, ps, nkb, (1.d0, 0.d0) , spsi(1,1), lda)

       DEALLOCATE(ps)

       RETURN

    END SUBROUTINE s_psi_nc

END SUBROUTINE s_psi
