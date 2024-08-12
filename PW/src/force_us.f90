!
! Copyright (C) 2001-2024 Quantum ESPRESSO group
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
  USE control_flags,        ONLY : gamma_only, offload_type
  USE cell_base,            ONLY : tpiba
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE klist,                ONLY : nks, xk, ngk, igk_k
  USE gvect,                ONLY : g
  USE uspp,                 ONLY : nkb, vkb, qq_at, deeq, qq_so, deeq_nc, ofsbeta
  USE uspp_param,           ONLY : upf, nh, nhm
  USE wvfct,                ONLY : nbnd, npwx, wg, et
  USE lsda_mod,             ONLY : lsda, current_spin, isk, nspin
  USE symme,                ONLY : symvector
  USE wavefunctions,        ONLY : evc
  USE noncollin_module,     ONLY : npol, noncolin, lspinorb
  USE io_files,             ONLY : iunwfc, nwordwfc
  USE buffers,              ONLY : get_buffer
  USE becmod,               ONLY : calbec, becp, bec_type, &
                                   allocate_bec_type, deallocate_bec_type, &
                                   allocate_bec_type_acc, deallocate_bec_type_acc
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp_bands,             ONLY : intra_bgrp_comm, me_bgrp, nproc_bgrp
  USE mp,                   ONLY : mp_sum
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
  TYPE(bec_type) :: becd                  ! contains <dbeta|psi>
  COMPLEX(DP) :: deff_nc
  REAL(DP) :: deff, fnl
  INTEGER :: npw, ik, ipol, ig, na, na_s, na_e, mykey
  INTEGER :: nt, ibnd, nhnt, ih, jh, ijkb0, ikb, jkb, is, js, ijs
  !
  forcenl(:,:) = 0.D0
  !
  CALL allocate_bec_type_acc( nkb, nbnd, becp, intra_bgrp_comm )
  CALL allocate_bec_type_acc( nkb, nbnd, becd, intra_bgrp_comm )
  ALLOCATE( vkb1(npwx,nkb) )
  !$acc data create(vkb1)
  ! 
  ! ... the forces are summed over K-points
  !
  DO ik = 1, nks
     !
     IF ( lsda ) current_spin = isk(ik)
     npw = ngk(ik)
     !
     IF ( nks > 1 ) THEN
        CALL get_buffer( evc, nwordwfc, iunwfc, ik )
        !$acc update device( evc )
     ENDIF
     !
     IF ( nkb > 0 ) CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb, .TRUE. )
     !$acc data present (evc, vkb, becp)
     CALL calbec( offload_type, npw, vkb, evc, becp )
     !$acc end data
     !
     DO ipol = 1, 3
        !
#if defined(_OPENACC)
        !$acc parallel loop collapse(2) present(vkb, g, igk_k) 
#else
        !$omp parallel do collapse(2) private(ig)
#endif
        DO jkb = 1, nkb
           DO ig = 1, npw
              vkb1(ig,jkb) = vkb(ig,jkb) * (0.D0,-1.D0) * g(ipol,igk_k(ig,ik))
           ENDDO
        ENDDO
        !$acc data present (evc, becd)
        CALL calbec( offload_type, npw, vkb1, evc, becd )
        !$acc end data
        !
        ! becp = <beta|psi>, becd = <dbeta/dG_ipol|psi>
        ! Now sum over bands and over projectors belonging to each atom
        !
        ! ... NOTE: calls to calbec are parallelized over the bgrp group
        ! ... The rest of the calculation is parallelized by subdividing 
        ! ... the atoms over the bgrp group
        !
        CALL block_distribute( nat, me_bgrp, nproc_bgrp, na_s, na_e, mykey )
        !
        IF ( mykey /= 0 ) CYCLE
        !
        !$acc data present(becp, becd, deeq, qq_at, deeq_nc, qq_so, et) copyin(wg) 
        DO na = na_s, na_e
           fnl = 0.0_dp
           nt = ityp(na)
           nhnt = nh(nt)
           ijkb0 = ofsbeta(na)
           IF ( gamma_only ) THEN
              !$acc parallel loop collapse(3) present(becp%r,becd%r) reduction(+:fnl)
              DO ibnd = 1, nbnd
                 DO ih = 1, nhnt
                    DO jh = 1, nhnt
                       ikb = ijkb0 + ih
                       jkb = ijkb0 + jh
                       deff = deeq(ih,jh,na,current_spin) - &
                            et(ibnd,ik) * qq_at(ih,jh,na)
                       fnl = fnl + wg(ibnd,ik) * deff *  &
                            becd%r(ikb,ibnd) * becp%r(jkb,ibnd)
                    END DO
                 END DO
              END DO
           ELSE IF ( .NOT. noncolin ) THEN
              !$acc parallel loop collapse(3) present(becp%k,becd%k) reduction(+:fnl)
              DO ibnd = 1, nbnd
                 DO ih = 1, nhnt
                    DO jh = 1, nhnt
                       ikb = ijkb0 + ih
                       jkb = ijkb0 + jh
                       deff = deeq(ih,jh,na,current_spin) - et(ibnd,ik) * qq_at(ih,jh,na)
                       fnl = fnl + wg(ibnd,ik) * deff *  &
                            DBLE(CONJG(becp%k(ikb,ibnd)) * becd%k(jkb,ibnd) )
                    END DO
                 END DO
              END DO
           ELSE
              !$acc parallel loop collapse(3) present(becp%nc,becd%nc) reduction(+:fnl)
              DO ibnd = 1, nbnd
                 DO ih = 1, nhnt
                    DO jh = 1, nhnt
                       ikb = ijkb0 + ih
                       jkb = ijkb0 + jh
                       !$acc loop seq collapse(2)
                       DO is = 1, npol
                          DO js = 1, npol
                             ijs = (is-1)*npol + js
                             deff_nc = deeq_nc(ih,jh,na,ijs)
                             IF ( lspinorb ) THEN
                                deff_nc = deff_nc - et(ibnd,ik) * qq_so(ih,jh,ijs,nt)
                             ELSE IF ( is == js ) THEN
                                deff_nc = deff_nc - et(ibnd,ik) * qq_at(ih,jh,na)
                             END IF
                             fnl = fnl + wg(ibnd,ik) * DBLE ( &
                                  deff_nc * CONJG(becp%nc(ikb,is,ibnd)) * &
                                  becd%nc(jkb,js,ibnd) )
                          END DO
                       END DO
                    END DO
                 END DO
              END DO
           END IF
           ! factor 2 from Ry a.u. (e^2=2)? tpiba from k+G, minus sign
           forcenl(ipol,na) = forcenl(ipol,na) - 2.0_dp * tpiba* fnl
        END DO
        !$acc end data
        !
     ENDDO
  ENDDO
  !
  !$acc end data
  DEALLOCATE( vkb1 )
  CALL deallocate_bec_type_acc( becd )
  CALL deallocate_bec_type_acc( becp )
  !
  ! ... collect contributions across processors and pools from all k-points
  !
  CALL mp_sum( forcenl, intra_bgrp_comm )
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
END SUBROUTINE force_us
