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
       INTEGER :: ierr
       !
       IF ( nkb == 0 ) RETURN
       !
       ALLOCATE( ps (nkb,m), STAT=ierr )
       IF ( ierr /= 0 ) &
          CALL errore( ' add_vuspsi_gamma ', ' cannot allocate ps ', ABS(ierr) )
       !
       ps(:,:) = 0.D0
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
                CALL DGEMM('N', 'N', nh(nt), m, nh(nt), 1.0_dp, &
                           deeq(1,1,na,current_spin), nhm, &
                           becp%r(ofsbeta(na)+1,1), nkb, 0.0_dp, &
                               ps(ofsbeta(na)+1,1), nkb )
                !
             ENDIF
             !
          ENDDO
          !
       ENDDO
       !
       ! Normal case: hpsi(n,i) = \sum_l beta(n,l) ps(l,i) 
       ! (l runs from 1 to nkb)
       !
       CALL DGEMM( 'N', 'N', ( 2 * n ), m, nkb, 1.D0, vkb, &
                   ( 2 * lda ), ps, nkb, 1.D0, hpsi, ( 2 * lda ) )
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
       INTEGER :: ierr
       !
       IF ( nkb == 0 ) RETURN
       !
       ALLOCATE( ps(nkb,m), STAT=ierr )
       IF( ierr /= 0 ) &
          CALL errore( ' add_vuspsi_k ', ' cannot allocate ps ', ABS( ierr ) )
       !
       DO nt = 1, ntyp
          !
          IF ( nh(nt) == 0 ) CYCLE
          ALLOCATE( deeaux(nh(nt),nh(nt)) )
          DO na = 1, nat
             !
             IF ( ityp(na) == nt ) THEN
                !
                ! deeq is real: copy it into a complex variable to perform
                ! a zgemm - simple but sub-optimal solution
                !
                deeaux(:,:) = CMPLX(deeq(1:nh(nt),1:nh(nt),na,current_spin),&
                                    0.0_dp, KIND=dp )
                CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                           deeaux, nh(nt), becp%k(ofsbeta(na)+1,1), nkb, &
                          (0.0_dp, 0.0_dp), ps(ofsbeta(na)+1,1), nkb )
                !
             ENDIF
             !
          ENDDO
          DEALLOCATE( deeaux )
          !
       ENDDO
       !
       CALL ZGEMM( 'N', 'N', n, m, nkb, ( 1.D0, 0.D0 ) , vkb, &
                   lda, ps, nkb, ( 1.D0, 0.D0 ) , hpsi, lda )
       !
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
       INTEGER :: ierr, ijkb0
       !
       IF ( nkb == 0 ) RETURN
       !
       ALLOCATE( ps(  nkb,npol, m), STAT=ierr )
       IF( ierr /= 0 ) &
          CALL errore( ' add_vuspsi_nc ', ' error allocating ps ', ABS( ierr ) )
       !
       ps(:,:,:) = (0.d0, 0.d0)
       !
       DO nt = 1, ntyp
          !
          IF ( nh(nt) == 0 ) CYCLE
          DO na = 1, nat
             !
             IF ( ityp(na) == nt ) THEN
                !
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
                !
             ENDIF
             !
          ENDDO
          !
       ENDDO
       !
       CALL ZGEMM('N', 'N', n, m*npol, nkb, ( 1.D0, 0.D0 ) , vkb, &
                   lda, ps, nkb, ( 1.D0, 0.D0 ) , hpsi, lda )
       !
       DEALLOCATE( ps )
       !
       RETURN
       !
     END SUBROUTINE add_vuspsi_nc
     !  
     !  
END SUBROUTINE add_vuspsi
