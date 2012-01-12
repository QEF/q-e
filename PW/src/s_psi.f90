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
  USE ldaU,       ONLY : lda_plus_u
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
       !
       IMPLICIT NONE  
       !
       ! ... here the local variables
       !
       INTEGER :: ikb, jkb, ih, jh, na, nt, ijkb0, ibnd
         ! counters
       REAL(DP), ALLOCATABLE :: ps(:,:)
         ! the product vkb and psi
       !
       !
       ALLOCATE( ps( nkb, m ) )
       !    
       ps(:,:) = 0.D0
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
                                           qq(ih,jh,nt) * becp%r(jkb,ibnd)
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
       INTEGER :: ikb, jkb, ih, jh, na, nt, ijkb0, ibnd
         ! counters
       COMPLEX(DP), ALLOCATABLE :: ps(:,:)
         ! the product vkb and psi
       !
       ALLOCATE( ps( nkb, m ) )    
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
       INTEGER :: ikb, jkb, ih, jh, na, nt, ijkb0, ibnd, ipol
       ! counters
       COMPLEX (DP), ALLOCATABLE :: ps (:,:,:)
       ! the product vkb and psi
       !
       ALLOCATE (ps(nkb,npol,m))    
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
