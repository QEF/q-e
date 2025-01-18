!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE addusforce( forcenl )
  !----------------------------------------------------------------------
  !! Wrapper to \(\texttt{addusforce_g}\) or \(\texttt{addusforce_r}\).
  !
  USE kinds,         ONLY : dp
  USE ions_base,     ONLY : nat
  USE control_flags, ONLY : tqr
  USE realus,        ONLY : addusforce_r
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: forcenl(3,nat)
  !! the non-local contribution to the force
  !
  IF ( tqr ) THEN
     CALL addusforce_r( forcenl )
  ELSE
     CALL addusforce_g( forcenl )
  ENDIF
  !
END SUBROUTINE addusforce
!
!----------------------------------------------------------------------
SUBROUTINE addusforce_g( forcenl )
  !----------------------------------------------------------------------
  !! This routine computes the contribution to atomic forces due
  !! to the dependence of the Q function on the atomic position.
  !! \[ F_{j,\text{at}} = \sum_G \sum_{lm} iG_j\ \text{exp}(-iG*R_\text{at})
  !!    V^*(G)\ Q_{lm}(G)\ \text{becsum}(lm,\text{at}) \]
  !! where:
  !! \[ \text{becsum}(lm,\text{at}) = \sum_i \langle \psi_i|\beta_l\rangle
  !!    w_i\langle \beta_m|\psi_i\rangle \]
  !! On output: the contribution is added to \(\text{forcenl}\).
  !
  USE kinds,              ONLY : DP
  USE ions_base,          ONLY : nat, ntyp => nsp, ityp
  USE cell_base,          ONLY : omega, tpiba
  USE fft_base,           ONLY : dfftp
  USE fft_rho,            ONLY : rho_r2g
  USE gvect,              ONLY : ngm, gg, g, eigts1, eigts2, eigts3,&
                                 mill
  USE noncollin_module,   ONLY : nspin_mag
  USE scf,                ONLY : v, vltot
  USE uspp,               ONLY : becsum, okvan
  USE uspp_param,         ONLY : upf, lmaxq, nh, nhm
  USE mp_bands,           ONLY : intra_bgrp_comm
  USE mp_pools,           ONLY : inter_pool_comm
  USE mp,                 ONLY : mp_sum
  USE control_flags,      ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: forcenl(3,nat)
  !! the non-local contribution to the force
  !
  ! ... local variables
  !
  INTEGER :: ngm_s, ngm_e, ngm_l
  INTEGER :: ig, nt, ih, jh, ijh, nij, ipol, is, na, nb, nab, ir
  REAL(DP), ALLOCATABLE :: forceq(:,:)
  REAL(DP), ALLOCATABLE :: ddeeq(:,:,:,:), qmod(:), ylmk0(:,:)
  REAL(DP) :: fact, forceqx, forceqy, forceqz
  COMPLEX(DP) :: cfac
  COMPLEX(DP), ALLOCATABLE :: aux1(:,:,:), vg(:,:), qgm(:,:)
  !
  IF (.NOT.okvan) RETURN
  !
  ALLOCATE( forceq(3,nat) )
  forceq(:,:) = 0._DP
  IF ( gamma_only ) THEN
     fact = 2._DP*omega
  ELSE
     fact = 1._DP*omega
  ENDIF
  !
  ! ... Fourier transform of the total effective potential
  !
  ALLOCATE( vg(ngm,nspin_mag) )
  !$acc data copyin( v, vltot ) create( vg ) 
  !$acc data copyin( v%of_r )
  !
  DO is = 1, nspin_mag
    IF ( nspin_mag == 4 .and. is /= 1 ) THEN
       CALL rho_r2g( dfftp, v%of_r(:,is), vg(:,is:is) )
    ELSE
       CALL rho_r2g( dfftp, v%of_r(:,is), vg(:,is:is), v=vltot )
    ENDIF
    ! ... Note the factors -i and 2pi/a *units of G) here in V(G) !
    !$acc kernels
    vg(:,is) = vg(:,is) * tpiba * (0.d0,-1.d0)
    !$acc end kernels
  ENDDO
  !
  ! ... With k-point parallelization, distribute G-vectors across processors
  ! ... ngm_s = index of first G-vector for this processor
  ! ... ngm_e = index of last  G-vector for this processor
  ! ... ngm_l = local number of G-vectors
  !
  CALL divide( inter_pool_comm, ngm, ngm_s, ngm_e )
  ngm_l = ngm_e-ngm_s+1
  !
  ! ... for the extraordinary unlikely case of more processors than G-vectors
  IF ( ngm_l <= 0 ) GO TO 10
  !
  ALLOCATE( ylmk0(ngm_l,lmaxq*lmaxq), qmod(ngm_l) )
  !$acc data create(ylmk0,qmod)
  !
  CALL ylmr2( lmaxq*lmaxq, ngm_l, g(1,ngm_s), gg(ngm_s), ylmk0 )
  !
  !$acc parallel loop
  DO ig = 1, ngm_l
     qmod(ig) = SQRT( gg(ngm_s+ig-1) )*tpiba
  ENDDO
  !
  DO nt = 1, ntyp
     IF ( upf(nt)%tvanp ) THEN
        !
        ! ... nij = max number of (ih,jh) pairs per atom type nt
        ! ... qgm contains the Q functions in G space
        !
        nij = nh(nt)*(nh(nt)+1)/2
        !
        ALLOCATE( qgm(ngm_l,nij) )
        !$acc data create(qgm)
        !
        ijh = 0
        DO ih = 1, nh(nt)
           DO jh = ih, nh(nt)
              ijh = ijh + 1
              CALL qvan2( ngm_l, ih, jh, nt, qmod, qgm(1,ijh), ylmk0 )
           ENDDO
        ENDDO
        !
        ! ... nab = number of atoms of type nt
        !
        nab = 0
        DO na = 1, nat
           IF ( ityp(na)==nt ) nab = nab + 1
        ENDDO
        !
        ALLOCATE( aux1(ngm_l,nab,3), ddeeq(nij,nab,3,nspin_mag) )
        !$acc data create(ddeeq)
        !$acc data create(aux1)
        !
        DO is = 1, nspin_mag
           nb = 0
           DO na = 1, nat
              IF (ityp(na) == nt) THEN
                 nb = nb + 1
                 !
                 ! ... aux1 = product of potential, structure factor and iG
                 !
#if defined(_OPENACC)
!$acc parallel loop present(eigts1,eigts2,eigts3,mill,g)
#else
!$omp parallel do default(shared) private(ig,cfac)
#endif
                 DO ig = 1, ngm_l
                    cfac = vg(ngm_s+ig-1,is) * &
                           CONJG(eigts1(mill(1,ngm_s+ig-1),na) * &
                                 eigts2(mill(2,ngm_s+ig-1),na) * &
                                 eigts3(mill(3,ngm_s+ig-1),na) )
                    aux1(ig,nb,1) = g(1,ngm_s+ig-1) * cfac
                    aux1(ig,nb,2) = g(2,ngm_s+ig-1) * cfac
                    aux1(ig,nb,3) = g(3,ngm_s+ig-1) * cfac
                 ENDDO
#if !defined(_OPENACC)
!$omp end parallel do
#endif
                 !
              ENDIF
           ENDDO
           !
           ! ... ddeeq = dot product of aux1 with the Q functions
           ! ... No need for special treatment of the G=0 term (is zero)
           !
           !$acc host_data use_device(qgm,aux1,ddeeq)
           DO ipol = 1, 3
              CALL MYDGEMM( 'C', 'N', nij, nab, 2*ngm_l, fact, qgm, 2*ngm_l, &
                            aux1(1,1,ipol), 2*ngm_l, 0.0_DP, ddeeq(1,1,ipol,is), nij )
           ENDDO
           !$acc end host_data
           !
        ENDDO
        !
        !$acc end data
        DEALLOCATE( aux1 )
        !
        DO is = 1, nspin_mag
           nb = 0
           DO na = 1, nat
              IF (ityp(na) == nt) THEN
                 nb = nb + 1
                 forceqx = 0
                 forceqy = 0
                 forceqz = 0
                 !$acc parallel loop reduction(+:forceqx,forceqy,forceqz) present(becsum)
                 DO ijh = 1, nij
                   forceqx = forceqx + ddeeq(ijh,nb,1,is) * becsum(ijh,na,is)
                   forceqy = forceqy + ddeeq(ijh,nb,2,is) * becsum(ijh,na,is)
                   forceqz = forceqz + ddeeq(ijh,nb,3,is) * becsum(ijh,na,is)
                 ENDDO
                 forceq(1,na) = forceq(1,na) + forceqx
                 forceq(2,na) = forceq(2,na) + forceqy
                 forceq(3,na) = forceq(3,na) + forceqz
              ENDIF
           ENDDO
        ENDDO
        !$acc end data
        DEALLOCATE( ddeeq )
        !$acc end data
        DEALLOCATE( qgm )
     ENDIF
  ENDDO
  !
  !$acc end data
  DEALLOCATE( ylmk0, qmod )
  !
  10 CONTINUE
  CALL mp_sum( forceq, inter_pool_comm )
  CALL mp_sum( forceq, intra_bgrp_comm )
  !
  forcenl(:,:) = forcenl(:,:) + forceq(:,:)
  !
  !$acc end data
  !$acc end data
  DEALLOCATE( vg )
  DEALLOCATE( forceq )
  !
  RETURN
  !
END SUBROUTINE addusforce_g
