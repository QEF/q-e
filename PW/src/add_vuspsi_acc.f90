!
! Copyright (C) 2001-2025 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE add_vuspsi_acc( lda, n, m, hpsi )
  !----------------------------------------------------------------------------
  !! This routine applies the nonlocal potential V=sum_lm |beta_l>D_lm<beta_m|
  !! to a set of vectors |psi_i> and puts the result in |hpsi_i>. 
  !! It requires on input the products becp(m,i) = <beta_m|psi_i> between psi
  !! and the beta functions (in array becp(nkb,m), calculated by calbec).
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
  COMPLEX(DP), INTENT(INOUT) :: hpsi(lda*npol,m)
  !! V_US|psi> is added to hpsi
  !
  CALL start_clock( 'add_vuspsi' )  
  !
  IF ( gamma_only ) THEN
     !
     CALL add_vuspsi_gamma_acc()
     !
  ELSE IF ( noncolin) THEN
     !
     CALL add_vuspsi_nc_acc ()
     !
  ELSE
     !
     CALL add_vuspsi_k_acc()
     !
  END IF
  !
  CALL stop_clock( 'add_vuspsi' )  
  !
  RETURN
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     SUBROUTINE add_vuspsi_gamma_acc()
       !-----------------------------------------------------------------------
       !! Gamma-only version
       !
       IMPLICIT NONE
       REAL(DP), ALLOCATABLE :: ps (:,:)
       INTEGER :: ierr
       ! counters
       INTEGER :: jkb, ikb, ih, jh, na, nt, ibnd
       !
       REAL(DP), ALLOCATABLE :: becp_r(:,:)
       !$acc declare device_resident(ps,becp_r)
       !
       IF ( nkb == 0 ) RETURN
       !
       ! becp(l,i) = <beta_l|psi_i>, with vkb(n,l)=|beta_l>
       !
       ALLOCATE (ps (nkb,m), STAT=ierr )
       !
       IF( ierr /= 0 ) &
          CALL errore( ' add_vuspsi_gamma ', ' cannot allocate ps ', ABS(ierr) )
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
                !$acc host_data use_device(deeq,becp_r,ps)
                CALL MYDGEMM('N', 'N', nh(nt), m, nh(nt), 1.0_dp, &
                           deeq(1,1,na,current_spin), nhm, &
                           becp_r(ofsbeta(na)+1,1), nkb, 0.0_dp, &
                               ps(ofsbeta(na)+1,1), nkb )
                !$acc end host_data
                !
             END IF
             !
          END DO
          !
       END DO
       !
       ! hpsi(n,i) = hpis(n,i) + \sum_l beta(n,l) ps(l,i) 
       ! (l runs from 1 to nkb)
       !
       !$acc host_data use_device(vkb,ps,hpsi)
       CALL MYDGEMM( 'N', 'N', ( 2 * n ), m, nkb, 1.D0, vkb, &
                   ( 2 * lda ), ps, nkb, 1.D0, hpsi, ( 2 * lda ) )
       !$acc end host_data
       !
       DEALLOCATE (ps)
       DEALLOCATE( becp_r )
       !
       RETURN
       !
     END SUBROUTINE add_vuspsi_gamma_acc
     !
     !-----------------------------------------------------------------------
     SUBROUTINE add_vuspsi_k_acc()
       !-----------------------------------------------------------------------
       !! k-point version, see add_vuspsi_gamma for comments
       !
       IMPLICIT NONE
       COMPLEX(DP), ALLOCATABLE :: ps (:,:), deeaux (:,:)
       INTEGER :: ierr
       ! counters
       INTEGER :: i, j, k, jkb, ikb, ih, jh, na, nt, ibnd, nhnt

       COMPLEX(DP), ALLOCATABLE :: becp_k(:,:)
       !$acc declare device_resident(ps, deeaux,becp_k)
       !
       IF ( nkb == 0 ) RETURN
       !
       ALLOCATE (ps (nkb,m), STAT=ierr )
       IF( ierr /= 0 ) &
          CALL errore( ' add_vuspsi_k ', ' cannot allocate ps ', ABS( ierr ) )
       !
       ALLOCATE ( deeaux(nhm, nhm) )
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
       !$acc data present(deeq)
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
                !$acc parallel loop collapse(2)
                DO j = 1, nhnt
                   DO k = 1, nhnt
                      deeaux(k,j) = CMPLX(deeq(k,j,na,current_spin), 0.0_dp, KIND=DP )
                   END DO
                END DO
                !
                !$acc host_data use_device(deeaux,becp_k,ps)
                CALL MYZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                           deeaux, nhm, becp_k(ofsbeta(na)+1,1), nkb, &
                          (0.0_dp, 0.0_dp), ps(ofsbeta(na)+1,1), nkb )
                !$acc end host_data
                !
             END IF
             !
          END DO
          !
       END DO
       !$acc end data
       !
       !$acc host_data use_device(vkb,ps,hpsi)
       CALL MYZGEMM( 'N', 'N', n, m, nkb, ( 1.D0, 0.D0 ) , vkb, &
                   lda, ps, nkb, ( 1.D0, 0.D0 ) , hpsi, lda )
       !$acc end host_data
       !
       DEALLOCATE (ps)
       DEALLOCATE (deeaux)
       DEALLOCATE( becp_k )
       !
       RETURN
       !
     END SUBROUTINE add_vuspsi_k_acc
     !  
     !-----------------------------------------------------------------------
     SUBROUTINE add_vuspsi_nc_acc()
       !-----------------------------------------------------------------------
       !! NOn-colinear/spinorbit case, see add_vuspsi_gamma for comments
       !
       IMPLICIT NONE
       COMPLEX(DP), ALLOCATABLE :: ps (:,:,:)
       INTEGER :: ierr
       ! counters
       INTEGER :: na, nt, ibnd
       COMPLEX(DP), ALLOCATABLE :: becp_nc(:,:,:)
       !$acc declare device_resident(ps,becp_nc)
       !
       IF ( nkb == 0 ) RETURN
       !
       ALLOCATE (ps( nkb, npol, m), STAT=ierr )
       IF( ierr /= 0 ) &
          CALL errore( ' add_vuspsi_nc ', ' error allocating ps ', ABS( ierr ) )
       !
       ALLOCATE( becp_nc(size(becp%nc,1),size(becp%nc,2),size(becp%nc,3)), stat=ierr )
       IF( ierr /= 0 ) &
          CALL errore( ' add_vuspsi_nc ', ' error allocating becp_nc ', ABS( ierr ) )
       !$acc kernels
       becp_nc = becp%nc
       !$acc end kernels
       !
       DO nt = 1, ntyp
          !
          IF ( nh(nt) == 0 ) CYCLE
          DO na = 1, nat
             !
             IF ( ityp(na) == nt ) THEN
                !
                !$acc host_data use_device(deeq_nc,becp_nc,ps)
                CALL MYZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                           deeq_nc(1,1,na,1), nhm, becp_nc(ofsbeta(na)+1,1,1), 2*nkb, &
                          (0.0_dp, 0.0_dp), ps(ofsbeta(na)+1,1,1), 2*nkb )

                CALL MYZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                           deeq_nc(1,1,na,2), nhm, becp_nc(ofsbeta(na)+1,2,1), 2*nkb, &
                          (1.0_dp, 0.0_dp), ps(ofsbeta(na)+1,1,1), 2*nkb )


                CALL MYZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                           deeq_nc(1,1,na,3), nhm, becp_nc(ofsbeta(na)+1,1,1), 2*nkb, &
                          (0.0_dp, 0.0_dp), ps(ofsbeta(na)+1,2,1), 2*nkb )


                CALL MYZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                           deeq_nc(1,1,na,4), nhm, becp_nc(ofsbeta(na)+1,2,1), 2*nkb, &
                          (1.0_dp, 0.0_dp), ps(ofsbeta(na)+1,2,1), 2*nkb )
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
       !$acc host_data use_device(vkb,ps,hpsi)
       call MYZGEMM ('N', 'N', n, m*npol, nkb, ( 1.D0, 0.D0 ) , vkb, &
                   lda, ps, nkb, ( 1.D0, 0.D0 ) , hpsi, lda )
       !$acc end host_data
       !
       DEALLOCATE (ps)
       DEALLOCATE( becp_nc )
       !
       RETURN
       !
     END SUBROUTINE add_vuspsi_nc_acc
     !
END SUBROUTINE add_vuspsi_acc
