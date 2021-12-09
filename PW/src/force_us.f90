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
  !! The nonlocal potential contribution to forces.
  !
  USE kinds,                ONLY : DP
  USE control_flags,        ONLY : gamma_only
  USE cell_base,            ONLY : at, bg, tpiba
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE klist,                ONLY : nks, xk, ngk, igk_k
  USE gvect,                ONLY : g
  USE uspp,                 ONLY : nkb, vkb, qq_at, deeq, qq_so, deeq_nc, ofsbeta
  USE uspp_param,           ONLY : upf, nh, nhm
  USE wvfct,                ONLY : nbnd, npwx, wg, et
  USE lsda_mod,             ONLY : lsda, current_spin, isk, nspin
  USE symme,                ONLY : symvector
  USE wavefunctions,        ONLY : evc
  USE noncollin_module,     ONLY : npol, noncolin
  USE io_files,             ONLY : iunwfc, nwordwfc
  USE buffers,              ONLY : get_buffer
  USE becmod,               ONLY : calbec, becp, bec_type, allocate_bec_type, &
                                   deallocate_bec_type
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum, mp_get_comm_null
  USE wavefunctions_gpum,   ONLY : using_evc
  USE wvfct_gpum,           ONLY : using_et
  USE becmod_subs_gpum,     ONLY : using_becp_auto
  USE uspp_init,            ONLY : init_us_2
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: forcenl(3,nat)
  !! the nonlocal contribution
  !
  ! ... local variables
  !
  COMPLEX(DP), ALLOCATABLE :: vkb1(:,:)   ! contains g*|beta>
  COMPLEX(DP), ALLOCATABLE :: deff_nc(:,:,:,:)
  REAL(DP), ALLOCATABLE :: deff(:,:,:)
  TYPE(bec_type) :: dbecp                 ! contains <dbeta|psi>
  INTEGER    :: npw, ik, ipol, ig, jkb
  !
  !
  forcenl(:,:) = 0.D0
  !
  CALL allocate_bec_type( nkb, nbnd, becp, intra_bgrp_comm )   
  CALL using_becp_auto(2)
  CALL allocate_bec_type( nkb, nbnd, dbecp, intra_bgrp_comm )   
  !
  ALLOCATE( vkb1( npwx, nkb ) )   
  !
  IF (noncolin) THEN
     ALLOCATE( deff_nc(nhm,nhm,nat,nspin) )
  ELSEIF (.NOT. gamma_only ) THEN
     ALLOCATE( deff(nhm,nhm,nat) )
  ENDIF
  !
  ! ... the forces are a sum over the K points and over the bands
  !   
  CALL using_evc(0)
  !
  DO ik = 1, nks
     !
     IF ( lsda ) current_spin = isk(ik)
     npw = ngk (ik)

     IF ( nks > 1 ) THEN
        CALL get_buffer( evc, nwordwfc, iunwfc, ik )
        CALL using_evc(1)
        IF ( nkb > 0 ) CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb )
     ENDIF
     !
     CALL using_becp_auto(2)
     CALL calbec( npw, vkb, evc, becp )
     !
     DO ipol = 1, 3
!$omp parallel do collapse(2) private(ig)
        DO jkb = 1, nkb
           DO ig = 1, npw
              vkb1(ig,jkb) = vkb(ig,jkb) * (0.D0,-1.D0) * g(ipol,igk_k(ig,ik))
           ENDDO
        ENDDO
!$omp end parallel do
        !
        CALL calbec( npw, vkb1, evc, dbecp )
        !
        IF ( gamma_only ) THEN
           !
           CALL force_us_gamma( forcenl )
           !
        ELSE
           !
           CALL force_us_k( forcenl )
           !
        ENDIF
     ENDDO
  ENDDO
  !
  ! ... if sums over bands are parallelized over the band group
  !
  CALL using_becp_auto(1)
  IF ( becp%comm /= mp_get_comm_null() ) CALL mp_sum( forcenl, becp%comm )
  !
  IF (noncolin) THEN
     DEALLOCATE( deff_nc )
  ELSEIF ( .NOT. GAMMA_ONLY) THEN
     DEALLOCATE( deff )
  ENDIF
  !
  DEALLOCATE( vkb1 )
  !
  CALL deallocate_bec_type( dbecp )
  CALL deallocate_bec_type( becp )
  CALL using_becp_auto(2)
  !
  ! ... collect contributions across pools from all k-points
  !
  CALL mp_sum( forcenl, inter_pool_comm )
  !
  ! ... The total D matrix depends on the ionic position via the
  ! ... augmentation part \int V_eff Q dr, the term deriving from the 
  ! ... derivative of Q is added in the routine addusforce
  !
  CALL addusforce( forcenl )
  !
  ! ... Since our summation over k points was only on the irreducible 
  ! ... BZ we have to symmetrize the forces.
  !
  CALL symvector( nat, forcenl )
  !
  RETURN
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     SUBROUTINE force_us_gamma( forcenl )
       !-----------------------------------------------------------------------
       !! Nonlocal contributiuon. Calculation at gamma.
       !
       IMPLICIT NONE
       !
       REAL(DP) :: forcenl(3,nat)
       !! the nonlocal contribution
       !
       ! ... local variables
       !
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
       CALL using_et(0)
       !
       DO nt = 1, ntyp
          IF ( nh(nt) == 0 ) CYCLE
          ALLOCATE( aux(nh(nt),becp%nbnd_loc) )
          DO na = 1, nat
             IF ( ityp(na) == nt ) THEN
                ijkb0 = ofsbeta(na)
                ! this is \sum_j q_{ij} <beta_j|psi>
                CALL DGEMM( 'N','N', nh(nt), becp%nbnd_loc, nh(nt),        &
                            1.0_dp, qq_at(1,1,na), nhm, becp%r(ijkb0+1,1), &
                            nkb, 0.0_dp, aux, nh(nt) )
                ! multiply by -\epsilon_n
                !
                !$omp parallel do default(shared) private(ibnd_loc,ibnd,ih)
                DO ih = 1, nh(nt)
                   DO ibnd_loc = 1, becp%nbnd_loc
                      ibnd = ibnd_loc + becp%ibnd_begin - 1
                      aux(ih,ibnd_loc) = - et(ibnd,ik) * aux(ih,ibnd_loc)
                   ENDDO
                ENDDO
                !$omp end parallel do
                !
                ! add  \sum_j d_{ij} <beta_j|psi>
                CALL DGEMM( 'N','N', nh(nt), becp%nbnd_loc, nh(nt), &
                            1.0_dp, deeq(1,1,na,current_spin), nhm, &
                            becp%r(ijkb0+1,1), nkb, 1.0_dp, aux, nh(nt) )
                !$omp parallel do default(shared) private(ibnd_loc,ibnd,ih) reduction(-:forcenl)
                DO ih = 1, nh(nt)
                   DO ibnd_loc = 1, becp%nbnd_loc
                      ibnd = ibnd_loc + becp%ibnd_begin - 1
                      forcenl(ipol,na) = forcenl(ipol,na) -    &
                           2.0_dp * tpiba * aux(ih,ibnd_loc) * &
                           dbecp%r(ijkb0+ih,ibnd_loc) * wg(ibnd,ik)
                   ENDDO
                ENDDO
                !$omp end parallel do
                !
             ENDIF
          ENDDO
          DEALLOCATE( aux )
       ENDDO
       !
     END SUBROUTINE force_us_gamma
     !     
     !-----------------------------------------------------------------------
     SUBROUTINE force_us_k( forcenl )
       !-----------------------------------------------------------------------
       !! Nonlocal contributiuon. Calculation for k-points.
       !
       IMPLICIT NONE
       !
       REAL(DP) :: forcenl(3,nat)
       !! the nonlocal contribution
       !
       ! ... local variables
       !
       REAL(DP) :: fac
       INTEGER  :: ibnd, ih, jh, na, nt, ikb, jkb, ijkb0, is, js, ijs !counters
       !
       CALL using_et(0)
       !
       DO ibnd = 1, nbnd
          !
          IF (noncolin) THEN
             CALL compute_deff_nc( deff_nc, et(ibnd,ik) )
          ELSE
             CALL compute_deff( deff, et(ibnd,ik) )
          ENDIF
          !
          fac = wg(ibnd,ik)*tpiba
          !
          DO nt = 1, ntyp
             DO na = 1, nat
                ijkb0 = ofsbeta(na)
                IF ( ityp(na) == nt ) THEN
                   DO ih = 1, nh(nt)
                      ikb = ijkb0 + ih
                      IF (noncolin) THEN
                         ijs=0
                         DO is = 1, npol
                            DO js = 1, npol
                               ijs=ijs+1
                               forcenl(ipol,na) = forcenl(ipol,na)- &
                                    deff_nc(ih,ih,na,ijs)*fac*(     &
                                    CONJG(dbecp%nc(ikb,is,ibnd))*   &
                                    becp%nc(ikb,js,ibnd)+           &
                                    CONJG(becp%nc(ikb,is,ibnd))*    &
                                    dbecp%nc(ikb,js,ibnd) )
                            ENDDO
                         ENDDO
                      ELSE
                         forcenl(ipol,na) = forcenl(ipol,na) -   &
                              2.D0 * fac * deff(ih,ih,na)*       &
                              DBLE( CONJG( dbecp%k(ikb,ibnd) ) * &
                              becp%k(ikb,ibnd) )
                      ENDIF
                   ENDDO
                   !
                   IF ( upf(nt)%tvanp .OR. upf(nt)%is_multiproj ) THEN
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
                               DO is = 1, npol
                                  DO js = 1, npol
                                     ijs = ijs + 1
                                     forcenl(ipol,na) = forcenl(ipol,na)- &
                                          deff_nc(ih,jh,na,ijs)*fac*(     &
                                          CONJG(dbecp%nc(ikb,is,ibnd))*   &
                                          becp%nc(jkb,js,ibnd)+           &
                                          CONJG(becp%nc(ikb,is,ibnd))*    &
                                          dbecp%nc(jkb,js,ibnd))-         &
                                          deff_nc(jh,ih,na,ijs)*fac*(     &
                                          CONJG(dbecp%nc(jkb,is,ibnd))*   &
                                          becp%nc(ikb,js,ibnd)+           &
                                          CONJG(becp%nc(jkb,is,ibnd))*    &
                                          dbecp%nc(ikb,js,ibnd) )
                                  ENDDO
                               ENDDO
                            ELSE
                               forcenl(ipol,na) = forcenl(ipol,na) -     &
                                    2.D0 * fac * deff(ih,jh,na) *        &
                                    DBLE( CONJG( dbecp%k(ikb,ibnd) ) *   &
                                    becp%k(jkb,ibnd) + dbecp%k(jkb,ibnd) &
                                    * CONJG( becp%k(ikb,ibnd) ) )
                            ENDIF
                         ENDDO !jh
                      ENDDO !ih
                   ENDIF ! tvanp
                   !
                ENDIF ! ityp(na) == nt
             ENDDO ! nat
          ENDDO ! ntyp
       ENDDO ! nbnd
       !
       !
     END SUBROUTINE force_us_k
     !
     !
END SUBROUTINE force_us
