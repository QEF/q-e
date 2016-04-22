!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE force_us( forcenl )
  !----------------------------------------------------------------------------
  !
  ! ... nonlocal potential contribution to forces
  !
  USE kinds,                ONLY : DP
  USE control_flags,        ONLY : gamma_only
  USE cell_base,            ONLY : at, bg, tpiba
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE klist,                ONLY : nks, xk, ngk, igk_k
  USE gvect,                ONLY : g
  USE uspp,                 ONLY : nkb, vkb, qq, deeq, qq_so, deeq_nc, indv_ijkb0
  USE uspp_param,           ONLY : upf, nh, newpseudo, nhm
  USE wvfct,                ONLY : nbnd, npwx, wg, et
  USE lsda_mod,             ONLY : lsda, current_spin, isk, nspin
  USE symme,                ONLY : symvector
  USE wavefunctions_module, ONLY : evc
  USE noncollin_module,     ONLY : npol, noncolin
  USE spin_orb,             ONLY : lspinorb
  USE io_files,             ONLY : iunwfc, nwordwfc
  USE buffers,              ONLY : get_buffer
  USE becmod,               ONLY : calbec, becp, bec_type, allocate_bec_type, &
                                   deallocate_bec_type
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum, mp_get_comm_null
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: forcenl(3,nat) ! the nonlocal contribution
  !
  COMPLEX(DP), ALLOCATABLE :: vkb1(:,:)   ! contains g*|beta>
  COMPLEX(DP), ALLOCATABLE :: deff_nc(:,:,:,:)
  REAL(DP), ALLOCATABLE :: deff(:,:,:)
  TYPE(bec_type) :: dbecp                 ! contains <dbeta|psi>
  INTEGER    :: npw, ik, ipol, ig, jkb
  !
  forcenl(:,:) = 0.D0
  !
  CALL allocate_bec_type ( nkb, nbnd, becp, intra_bgrp_comm )   
  CALL allocate_bec_type ( nkb, nbnd, dbecp, intra_bgrp_comm )   
  ALLOCATE( vkb1( npwx, nkb ) )   
  IF (noncolin) then
     ALLOCATE( deff_nc(nhm,nhm,nat,nspin) )
  ELSE IF (.NOT. gamma_only ) THEN
     ALLOCATE( deff(nhm,nhm,nat) )
  ENDIF
  !
  ! ... the forces are a sum over the K points and over the bands
  !   
  DO ik = 1, nks
     !
     IF ( lsda ) current_spin = isk(ik)
     npw = ngk (ik)

     IF ( nks > 1 ) THEN
        CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
        IF ( nkb > 0 ) &
             CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb )
     END IF
     !
     CALL calbec ( npw, vkb, evc, becp )
     !
     DO ipol = 1, 3
        DO jkb = 1, nkb
!$omp parallel do default(shared) private(ig)
           do ig = 1, npw
              vkb1(ig,jkb) = vkb(ig,jkb) * (0.D0,-1.D0) * g(ipol,igk_k(ig,ik))
           END DO
!$omp end parallel do
        END DO
        !
        CALL calbec ( npw, vkb1, evc, dbecp )
        !
        IF ( gamma_only ) THEN
           !
           CALL force_us_gamma( forcenl )
           !
        ELSE
           !
           CALL force_us_k( forcenl )
           !
        END IF
     END DO
  END DO
  !
  ! ... if sums over bands are parallelized over the band group
  !
  IF( becp%comm /= mp_get_comm_null() ) CALL mp_sum( forcenl, becp%comm )
  !
  IF (noncolin) THEN
     DEALLOCATE( deff_nc )
  ELSE IF ( .NOT. GAMMA_ONLY) THEN
     DEALLOCATE( deff )
  ENDIF
  DEALLOCATE( vkb1 )
  CALL deallocate_bec_type ( dbecp )   
  CALL deallocate_bec_type ( becp )   
  !
  ! ... The total D matrix depends on the ionic position via the
  ! ... augmentation part \int V_eff Q dr, the term deriving from the 
  ! ... derivative of Q is added in the routine addusforce
  !
  CALL addusforce( forcenl )
  !
  ! ... collect contributions across pools from all k-points
  !
  CALL mp_sum( forcenl, inter_pool_comm )
  !
  ! ... Since our summation over k points was only on the irreducible 
  ! ... BZ we have to symmetrize the forces.
  !
  CALL symvector ( nat, forcenl )
  !
  RETURN
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     SUBROUTINE force_us_gamma( forcenl )
       !-----------------------------------------------------------------------
       !
       ! ... calculation at gamma
       !
       IMPLICIT NONE
       !
       REAL(DP) :: forcenl(3,nat)
       REAL(DP), ALLOCATABLE :: aux(:,:)
       INTEGER ::  nt, na, ibnd, ibnd_loc, ih, jh, ijkb0 ! counters
       !
       ! ... Important notice about parallelization over the band group of processors:
       ! ... 1) internally, "calbec" parallelises on plane waves over the band group
       ! ... 2) the results of "calbec" are distributed across processors of the band
       ! ...    group: the band index of becp, dbecp is distributed
       ! ... 3) the band group is subsequently used to parallelize over bands
       !
       !
       DO nt = 1, ntyp
          IF ( nh(nt) == 0 ) CYCLE
          ALLOCATE ( aux(nh(nt),becp%nbnd_loc) )
          DO na = 1, nat
             IF ( ityp(na) == nt ) THEN
                ijkb0 = indv_ijkb0(na)
                ! this is \sum_j q_{ij} <beta_j|psi>
                CALL DGEMM ('N','N', nh(nt), becp%nbnd_loc, nh(nt), &
                     1.0_dp, qq(1,1,nt), nhm, becp%r(ijkb0+1,1),&
                     nkb, 0.0_dp, aux, nh(nt) )
                ! multiply by -\epsilon_n
!$omp parallel do default(shared) private(ibnd_loc,ibnd,ih)
                DO ih = 1, nh(nt)
                   DO ibnd_loc = 1, becp%nbnd_loc
                      ibnd = ibnd_loc + becp%ibnd_begin - 1
                      aux(ih,ibnd_loc) = - et(ibnd,ik) * aux(ih,ibnd_loc)
                   END DO
                END DO
!$omp end parallel do
                ! add  \sum_j d_{ij} <beta_j|psi>
                CALL DGEMM ('N','N', nh(nt), becp%nbnd_loc, nh(nt), &
                     1.0_dp, deeq(1,1,na,current_spin), nhm, &
                     becp%r(ijkb0+1,1), nkb, 1.0_dp, aux, nh(nt) )
!$omp parallel do default(shared) private(ibnd_loc,ibnd,ih) reduction(-:forcenl)
                DO ih = 1, nh(nt)
                   DO ibnd_loc = 1, becp%nbnd_loc
                      ibnd = ibnd_loc + becp%ibnd_begin - 1
                      forcenl(ipol,na) = forcenl(ipol,na) - &
                           2.0_dp * tpiba * aux(ih,ibnd_loc) * &
                           dbecp%r(ijkb0+ih,ibnd_loc) * wg(ibnd,ik)
                   END DO
                END DO
!$omp end parallel do
                !
             END IF
          END DO
          DEALLOCATE (aux)
       END DO
       !
     END SUBROUTINE force_us_gamma
     !     
     !-----------------------------------------------------------------------
     SUBROUTINE force_us_k( forcenl )
       !-----------------------------------------------------------------------
       !  
       IMPLICIT NONE
       !
       REAL(DP) :: forcenl(3,nat)
       !
       REAL(DP) :: fac
       INTEGER  :: ibnd, ih, jh, na, nt, ikb, jkb, ijkb0, is, js, ijs !counters
       !
       DO ibnd = 1, nbnd
          IF (noncolin) THEN
             CALL compute_deff_nc(deff_nc,et(ibnd,ik))
          ELSE
             CALL compute_deff(deff,et(ibnd,ik))
          ENDIF
          fac=wg(ibnd,ik)*tpiba
          DO nt = 1, ntyp
             DO na = 1, nat
                ijkb0 = indv_ijkb0(na)
                IF ( ityp(na) == nt ) THEN
                   DO ih = 1, nh(nt)
                      ikb = ijkb0 + ih
                      IF (noncolin) THEN
                         ijs=0
                         DO is=1,npol
                            DO js=1,npol
                               ijs=ijs+1
                               forcenl(ipol,na) = forcenl(ipol,na)- &
                                    deff_nc(ih,ih,na,ijs)*fac*( &
                                    CONJG(dbecp%nc(ikb,is,ibnd))* &
                                    becp%nc(ikb,js,ibnd)+ &
                                    CONJG(becp%nc(ikb,is,ibnd))* &
                                    dbecp%nc(ikb,js,ibnd) )
                            END DO
                         END DO
                      ELSE
                         forcenl(ipol,na) = forcenl(ipol,na) - &
                              2.D0 * fac * deff(ih,ih,na)*&
                              DBLE( CONJG( dbecp%k(ikb,ibnd) ) * &
                              becp%k(ikb,ibnd) )
                      END IF
                   END DO
                   !
                   IF ( upf(nt)%tvanp .OR. newpseudo(nt) ) THEN
                      DO ih = 1, nh(nt)
                         ikb = ijkb0 + ih
                         !
                         ! ... in US case there is a contribution for jh<>ih. 
                         ! ... We use here the symmetry in the interchange 
                         ! ... of ih and jh
                         !
                         DO jh = ( ih + 1 ), nh(nt)
                            jkb = ijkb0 + jh
                            IF (noncolin) THEN
                               ijs=0
                               DO is=1,npol
                                  DO js=1,npol
                                     ijs=ijs+1
                                     forcenl(ipol,na)=forcenl(ipol,na)- &
                                          deff_nc(ih,jh,na,ijs)*fac*( &
                                          CONJG(dbecp%nc(ikb,is,ibnd))* &
                                          becp%nc(jkb,js,ibnd)+ &
                                          CONJG(becp%nc(ikb,is,ibnd))* &
                                          dbecp%nc(jkb,js,ibnd))- &
                                          deff_nc(jh,ih,na,ijs)*fac*( &
                                          CONJG(dbecp%nc(jkb,is,ibnd))* &
                                          becp%nc(ikb,js,ibnd)+ &
                                          CONJG(becp%nc(jkb,is,ibnd))* &
                                          dbecp%nc(ikb,js,ibnd) )
                                  END DO
                               END DO
                            ELSE
                               forcenl(ipol,na) = forcenl (ipol,na) - &
                                    2.D0 * fac * deff(ih,jh,na)* &
                                    DBLE( CONJG( dbecp%k(ikb,ibnd) ) * &
                                    becp%k(jkb,ibnd) +       &
                                    dbecp%k(jkb,ibnd) * &
                                    CONJG( becp%k(ikb,ibnd) ) )
                            END IF
                         END DO !jh
                      END DO !ih
                   END IF ! tvanp
                END IF ! ityp(na) == nt
             END DO ! nat
          END DO ! ntyp
       END DO ! nbnd

       !
     END SUBROUTINE force_us_k
     !     
END SUBROUTINE force_us
