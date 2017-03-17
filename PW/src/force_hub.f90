!
! Copyright (C) 2002-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE force_hub(forceh)
   !----------------------------------------------------------------------
   !
   ! This routine computes the Hubbard contribution to the force. It gives
   ! in output the product (dE_{hub}/dn_{ij}^{alpha})(dn_{ij}^{alpha}
   ! /du(alpha,ipol)) which is the force acting on the atom at tau_{alpha}
   ! (in the unit cell) along the direction ipol.
   !
   USE kinds,                ONLY : DP
   USE ions_base,            ONLY : nat, ntyp => nsp, ityp
   USE cell_base,            ONLY : at, bg
   USE ldaU,                 ONLY : hubbard_lmax, hubbard_l, U_projection, &
                                    nwfcU, wfcU, is_hubbard, lda_plus_u_kind, &
                                    copy_U_wfc, offsetU
   USE basis,                ONLY : natomwfc
   USE symme,                ONLY : symvector
   USE io_files,             ONLY : prefix
   USE wvfct,                ONLY : nbnd, npwx
   USE control_flags,        ONLY : gamma_only
   USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
   USE scf,                  ONLY : v
   USE mp_pools,             ONLY : inter_pool_comm, intra_pool_comm, me_pool, &
                                    nproc_pool
   USE mp,                   ONLY : mp_sum
   USE becmod,               ONLY : bec_type, becp, calbec, allocate_bec_type, &
                                    deallocate_bec_type
   USE uspp,                 ONLY : nkb, vkb, indv_ijkb0
   USE uspp_param,           ONLY : nh
   USE wavefunctions_module, ONLY : evc
   USE klist,                ONLY : nks, xk, ngk, igk_k
   USE io_files,             ONLY : nwordwfc, iunwfc
   USE buffers,              ONLY : get_buffer

   IMPLICIT NONE
   REAL (DP) :: forceh(3,nat)  ! output: the Hubbard forces

   type (bec_type) :: proj     ! proj(nwfcU,nbnd)
   COMPLEX (DP), ALLOCATABLE :: spsi(:,:), wfcatom(:,:) 
   REAL (DP), ALLOCATABLE :: dns(:,:,:,:)
   !       dns(ldim,ldim,nspin,nat) ! the derivative of the atomic occupations
   INTEGER :: npw, alpha, na, nt, is, m1, m2, ipol, ldim, ik, ijkb0
   INTEGER :: nb_s, nb_e, mykey

   IF (U_projection .NE. "atomic") CALL errore("force_hub", &
                   " forces for this U_projection_type not implemented",1)
   IF (lda_plus_u_kind == 1) CALL errore("force_hub", &
                   " forces in full LDA+U scheme are not yet implemented",1)

   call start_clock('force_hub')
   ldim= 2 * Hubbard_lmax + 1
   ALLOCATE ( dns(ldim,ldim,nspin,nat) )
   ALLOCATE ( spsi(npwx,nbnd) ) 
   ALLOCATE ( wfcatom (npwx,natomwfc) ) 
   call allocate_bec_type ( nkb, nbnd, becp) 
   call allocate_bec_type ( nwfcU, nbnd, proj )
   !
   ! poor-man parallelization over bands
   ! - if nproc_pool=1   : nb_s=1, nb_e=nbnd, mykey=0
   ! - if nproc_pool<=nbnd:each processor calculates band nb_s to nb_e; mykey=0
   ! - if nproc_pool>nbnd :each processor takes care of band na_s=nb_e;
   !   mykey labels how many times each band appears (mykey=0 first time etc.)
   !
   CALL block_distribute( nbnd, me_pool, nproc_pool, nb_s, nb_e, mykey )
   !
   forceh(:,:) = 0.d0
   !
   !    we start a loop on k points
   !
   DO ik = 1, nks
      !
      IF (lsda) current_spin = isk(ik)
      npw = ngk (ik)

      IF (nks > 1) &
         CALL get_buffer (evc, nwordwfc, iunwfc, ik)

      CALL init_us_2 (npw,igk_k(1,ik),xk(1,ik),vkb)
      CALL calbec( npw, vkb, evc, becp )
      CALL s_psi  (npwx, npw, nbnd, evc, spsi )

      ! re-calculate atomic wfc - wfcatom is used here as work space

      CALL atomic_wfc (ik, wfcatom)
      call copy_U_wfc (wfcatom)

      ! wfcU contains Hubbard-U atomic wavefunctions
      ! proj=<wfcU|S|evc> - no need to read S*wfcU from buffer

      CALL calbec( npw, wfcU, spsi, proj )

      ! now we need the first derivative of proj with respect to tau(alpha,ipol)

      DO alpha = 1,nat  ! forces are calculated for atom alpha ...
         !
         ijkb0 = indv_ijkb0(alpha) ! positions of beta functions for atom alpha
         DO ipol = 1,3  ! forces are calculated for coordinate ipol ...
            !
            IF ( gamma_only ) THEN
               CALL dndtau_gamma ( ldim, proj%r, spsi, alpha, ijkb0, ipol, ik, &
                                   nb_s, nb_e, mykey, dns )
            ELSE
               CALL dndtau_k     ( ldim, proj%k, spsi, alpha, ijkb0, ipol, ik, &
                                   nb_s, nb_e, mykey, dns )
            ENDIF
!!omp parallel do default(shared) private(na,nt,m1,m2,is)
            DO na = 1,nat                 ! the Hubbard atom
               nt = ityp(na)
               IF ( is_hubbard(nt) ) THEN
                  DO is = 1,nspin
                     DO m2 = 1,ldim
                        DO m1 = 1,ldim
                           forceh(ipol,alpha) = forceh(ipol,alpha) -    &
                              v%ns(m2,m1,is,na) * dns(m1,m2,is,na)
                        END DO
                     END DO
                  END DO
               END IF
            END DO
!!omp end parallel do
         END DO
      END DO
   END DO
   !
   CALL mp_sum( forceh, inter_pool_comm )
   !
   call deallocate_bec_type (becp)
   call deallocate_bec_type (proj)
   DEALLOCATE( wfcatom  ) 
   DEALLOCATE( spsi  ) 
   DEALLOCATE( dns ) 
   
   IF (nspin == 1) forceh(:,:) = 2.d0 * forceh(:,:)
   !
   ! ...symmetrize...
   !
   CALL symvector ( nat, forceh )
#if defined(__DEBUG)
   write(66,'("Hubbard contribution Begin")')
   write(66,'(3f12.6)') forceh(:,:)
   write(66,'("Hubbard contribution End")')
#endif
   !
   call stop_clock('force_hub')
   !
   RETURN
END SUBROUTINE force_hub
!
!-----------------------------------------------------------------------
SUBROUTINE dndtau_k &
     (ldim, proj, spsi, alpha, jkb0, ipol, ik, nb_s, nb_e, mykey, dns)
   !-----------------------------------------------------------------------
   !
   ! This routine computes the derivative of the ns with respect to the ionic
   ! displacement u(alpha,ipol) used to obtain the Hubbard contribution to the
   ! atomic forces.
   !
   USE kinds,                ONLY : DP
   USE ions_base,            ONLY : nat, ityp
   USE lsda_mod,             ONLY : nspin, current_spin
   USE ldaU,                 ONLY : is_hubbard, Hubbard_l, nwfcU, offsetU
   USE wvfct,                ONLY : nbnd, npwx, wg
   USE mp_pools,             ONLY : intra_pool_comm, me_pool, nproc_pool
   USE mp,                   ONLY : mp_sum
   
   IMPLICIT NONE

   INTEGER, INTENT(IN) :: alpha, jkb0, ipol, ik, ldim
   INTEGER, INTENT(IN) :: nb_s, nb_e, mykey
   COMPLEX (DP), INTENT(IN) :: proj(nwfcU,nbnd), spsi(npwx,nbnd)
   REAL (DP), INTENT (OUT) :: dns(ldim,ldim,nspin,nat)
   !
   INTEGER ::  ibnd, is, na, nt, m1, m2
   COMPLEX (DP), ALLOCATABLE :: dproj(:,:)
   !
   !
   CALL start_clock('dndtau')
   !
   ALLOCATE ( dproj(nwfcU,nb_s:nb_e) )
   CALL dprojdtau_k ( spsi, alpha, jkb0, ipol, ik, nb_s, nb_e, mykey, dproj )
   !
   ! compute the derivative of occupation numbers (the quantities dn(m1,m2))
   ! of the atomic orbitals. They are real quantities as well as n(m1,m2)
   !
   dns(:,:,:,:) = 0.d0
   ! band parallelization. If each band appears more than once
   ! compute its contribution only once (i.e. when mykey=0)
   IF ( mykey /= 0 ) GO TO 10
!!omp parallel do default(shared) private(na,nt,m1,m2,ibnd)
   DO na = 1, nat
      nt = ityp(na)
      IF ( is_hubbard(nt) ) THEN
         DO m1 = 1, 2*Hubbard_l(nt)+1
            DO m2 = m1, 2*Hubbard_l(nt)+1
               DO ibnd = nb_s, nb_e
                  dns(m1,m2,current_spin,na) = dns(m1,m2,current_spin,na) + &
                                          wg(ibnd,ik) *            &
                              DBLE( proj(offsetU(na)+m1,ibnd)  *   &
                             CONJG(dproj(offsetU(na)+m2,ibnd))  +   &
                                   dproj(offsetU(na)+m1,ibnd)  *   &
                             CONJG( proj(offsetU(na)+m2,ibnd)) )
               END DO
            END DO
         END DO
      END IF
   END DO
!!omp end parallel do
10   DEALLOCATE ( dproj ) 
   !
   CALL mp_sum(dns, intra_pool_comm)
   !
   ! In nspin.eq.1 k-point weight wg is normalized to 2 el/band 
   ! in the whole BZ but we are interested in dns of one spin component
   !
   IF (nspin == 1) dns = 0.5d0 * dns
   !
   ! impose hermiticity of dn_{m1,m2}
   !
!!omp parallel do default(shared) private(na,is,m1,m2)
   DO na = 1,nat
      DO is = 1,nspin
         DO m1 = 1,ldim
            DO m2 = m1+1,ldim
               dns(m2,m1,is,na) = dns(m1,m2,is,na)
            END DO
         END DO
      END DO
   END DO
!!omp end parallel do

   CALL stop_clock('dndtau')
   RETURN
END SUBROUTINE dndtau_k
!
!-----------------------------------------------------------------------
SUBROUTINE dndtau_gamma &
     (ldim, rproj, spsi, alpha, jkb0, ipol, ik, nb_s, nb_e, mykey, dns)
  !-----------------------------------------------------------------------
   !
   ! This routine computes the derivative of the ns with respect to the ionic
   ! displacement u(alpha,ipol) used to obtain the Hubbard contribution to the
   ! atomic forces.
   !
   USE kinds,                ONLY : DP
   USE ions_base,            ONLY : nat, ityp
   USE lsda_mod,             ONLY : nspin, current_spin
   USE ldaU,                 ONLY : is_hubbard, Hubbard_l, nwfcU, offsetU
   USE wvfct,                ONLY : nbnd, npwx, wg
   USE mp_pools,             ONLY : intra_pool_comm, me_pool, nproc_pool
   USE mp,                   ONLY : mp_sum

   IMPLICIT NONE

   INTEGER, INTENT(IN) ::  alpha, jkb0, ipol, ik, ldim
   COMPLEX (DP), INTENT(IN) :: spsi(npwx,nbnd)
   REAL(DP), INTENT (IN) ::  rproj(nwfcU,nbnd)
   REAL (DP), INTENT (OUT) :: dns(ldim,ldim,nspin,nat)
   INTEGER, INTENT(IN) :: nb_s, nb_e, mykey
   !
   INTEGER ::  ibnd, is, na, nt, m1, m2
   REAL (DP), ALLOCATABLE :: dproj(:,:)
   !
   !
   CALL start_clock('dndtau')
   !
   ALLOCATE ( dproj(nwfcU,nb_s:nb_e) )
   CALL dprojdtau_gamma ( spsi, alpha, jkb0, ipol, ik, nb_s, nb_e, mykey, dproj )
   !
   ! compute the derivative of occupation numbers (the quantities dn(m1,m2))
   ! of the atomic orbitals. They are real quantities as well as n(m1,m2)
   !
   dns(:,:,:,:) = 0.d0
   ! band parallelization. If each band appears more than once
   ! compute its contribution only once (i.e. when mykey=0)
   IF ( mykey /= 0 ) GO TO 10
!!omp parallel do default(shared) private(na,nt,m1,m2,is)
   DO na = 1, nat
      nt = ityp(na)
      IF ( is_hubbard(nt) ) THEN
         DO m1 = 1, 2*Hubbard_l(nt)+1
            DO m2 = m1, 2*Hubbard_l(nt)+1
               DO ibnd = nb_s, nb_e
                  dns(m1,m2,current_spin,na) = dns(m1,m2,current_spin,na) + &
                                          wg(ibnd,ik) * (   &
                              rproj(offsetU(na)+m1,ibnd)  *   &
                              dproj(offsetU(na)+m2,ibnd)  +   &
                              dproj(offsetU(na)+m1,ibnd)  *   &
                              rproj(offsetU(na)+m2,ibnd) )
               END DO
            END DO
         END DO
      END IF
   END DO
!!omp end parallel do
10   DEALLOCATE ( dproj ) 
   !
   CALL mp_sum(dns, intra_pool_comm)
   !
   ! In nspin.eq.1 k-point weight wg is normalized to 2 el/band 
   ! in the whole BZ but we are interested in dns of one spin component
   !
   IF (nspin == 1) dns = 0.5d0 * dns
   !
   ! impose hermiticity of dn_{m1,m2}
   !
!!omp parallel do default(shared) private(na,is,m1,m2)
   DO na = 1,nat
      DO is = 1,nspin
         DO m1 = 1,ldim
            DO m2 = m1+1,ldim
               dns(m2,m1,is,na) = dns(m1,m2,is,na)
            END DO
         END DO
      END DO
   END DO
!!omp end parallel do

   CALL stop_clock('dndtau')
   RETURN
END SUBROUTINE dndtau_gamma
!
!-----------------------------------------------------------------------
SUBROUTINE dprojdtau_k (spsi, alpha, ijkb0, ipol, ik, nb_s, nb_e, mykey, dproj)
   !-----------------------------------------------------------------------
   !
   ! This routine computes the first derivative of the projection
   ! <\fi^{at}_{I,m1}|S|\psi_{k,v,s}> with respect to the atomic displacement
   ! u(alpha,ipol) (we remember that ns_{I,s,m1,m2} = \sum_{k,v}
   ! f_{kv} <\fi^{at}_{I,m1}|S|\psi_{k,v,s}><\psi_{k,v,s}|S|\fi^{at}_{I,m2}>)
   !
   USE kinds,                ONLY : DP
   USE ions_base,            ONLY : nat, ntyp => nsp, ityp
   USE cell_base,            ONLY : tpiba
   USE gvect,                ONLY : g
   USE klist,                ONLY : nks, xk, ngk, igk_k
   USE ldaU,                 ONLY : is_hubbard, Hubbard_l, nwfcU, wfcU, offsetU
   USE wvfct,                ONLY : nbnd, npwx, wg
   USE uspp,                 ONLY : nkb, vkb, qq
   USE uspp_param,           ONLY : nh
   USE wavefunctions_module, ONLY : evc
   USE becmod,               ONLY : bec_type, becp, calbec
   USE mp_bands,             ONLY : intra_bgrp_comm
   USE mp,                   ONLY : mp_sum
   
   IMPLICIT NONE
   INTEGER, INTENT (IN) :: ik,      &! k-point index
                           alpha,   &! the displaced atom
                           ipol,    &! the component of displacement
                           ijkb0     ! position of beta functions for atom alpha
   INTEGER, INTENT (IN) :: nb_s, nb_e, mykey       ! band parallelization
   COMPLEX (DP), INTENT (IN) :: spsi(npwx,nbnd)    ! S|evc>
   COMPLEX (DP), INTENT (OUT) :: dproj(nwfcU,nb_s:nb_e) ! derivative of projection
   !
   INTEGER :: npw, nt, ig, na_, m1, ibnd, iwf, nt_, ih, jh, ldim
   REAL (DP) :: gvec
   COMPLEX (DP), ALLOCATABLE :: dproj0(:,:), dwfc(:,:), dbeta(:,:), &
                                betapsi(:,:), dbetapsi(:,:), &
                                wfatbeta(:,:), wfatdbeta(:,:)
   !      dwfc(npwx,ldim),       ! the derivative of the atomic d wfc
   !      dbeta(npwx,nhm),       ! the derivative of the beta function
   !      betapsi(nhm,nbnd),     ! <beta|evc>
   !      dbetapsi(nhm,nbnd),    ! <dbeta|evc>
   !      wfatbeta(nwfcU,nhm),   ! <wfc|beta>
   !      wfatdbeta(nwfcU,nhm)   ! <wfc|dbeta>

   call start_clock('dprojdtau')
   nt = ityp(alpha)
   npw= ngk(ik)
   ldim = 2 * Hubbard_l(nt) + 1

   dproj(:,:) = (0.d0, 0.d0)
   !
   ! First the derivatives of the atomic wfc and the beta are computed
   ! Note: parallelization here is over plane waves, not over bands!
   !
   IF ( is_hubbard(nt) ) THEN
      ALLOCATE ( dproj0(ldim,nbnd) )
      ALLOCATE ( dwfc(npwx,ldim) )
!!omp parallel do default(shared) private(ig,gvec,m1)
      DO ig = 1,npw
         gvec = g(ipol,igk_k(ig,ik)) * tpiba

         ! in the expression of dwfc we don't need (k+G) but just G; k always
         ! multiplies the underived quantity and gives an opposite contribution
         ! in c.c. term because the sign of the imaginary unit.
   
         DO m1 = 1, ldim
            dwfc(ig,m1) = (0.d0,-1.d0) * gvec * wfcU(ig,offsetU(alpha)+m1)
         END DO
      END DO
!!omp end parallel do

      CALL ZGEMM('C','N',ldim, nbnd, npw, (1.d0,0.d0), &
                  dwfc, npwx, spsi, npwx, (0.d0,0.d0), &
                  dproj0, ldim)

      DEALLOCATE ( dwfc ) 
      CALL mp_sum( dproj0, intra_bgrp_comm )
      ! copy to dproj results for the bands treated by this processor
      dproj( offsetU(alpha)+1:offsetU(alpha)+ldim, :) = dproj0(:, nb_s:nb_e)
      DEALLOCATE ( dproj0 ) 
      !
   END IF
   !
   ALLOCATE (dbetapsi(nh(nt),nbnd) ) 
   ALLOCATE (wfatdbeta(nwfcU,nh(nt)) )
   ALLOCATE ( wfatbeta(nwfcU,nh(nt)) )
   ALLOCATE ( dbeta(npwx,nh(nt)) )
!!omp parallel do default(shared) private(ig,ih)
   DO ih=1,nh(nt)
      DO ig = 1, npw
         dbeta(ig,ih) = vkb(ig,ijkb0+ih)
      END DO
   END DO
!!omp end parallel do
   CALL calbec ( npw, wfcU, dbeta, wfatbeta ) 
!!omp parallel do default(shared) private(ig,ih)
   DO ih=1,nh(nt)
      DO ig = 1, npw
         gvec = g(ipol,igk_k(ig,ik)) * tpiba
         dbeta(ig,ih) = (0.d0,-1.d0) * dbeta(ig,ih) * gvec
      END DO
   END DO
!!omp end parallel do
   CALL calbec ( npw, dbeta, evc, dbetapsi ) 
   CALL calbec ( npw, wfcU, dbeta, wfatdbeta ) 
   DEALLOCATE ( dbeta )
   ! calculate \sum_j qq(i,j)*dbetapsi(j)
   ! betapsi is used here as work space 
   ALLOCATE ( betapsi(nh(nt), nbnd) ) 
   betapsi(:,:) = (0.0_dp, 0.0_dp)
   ! here starts band parallelization
!!omp parallel do default(shared) private(ih,ibnd,jh)
   DO ih=1,nh(nt)
      DO ibnd=nb_s, nb_e
         DO jh=1,nh(nt)
            betapsi(ih,ibnd) = betapsi(ih,ibnd) + &
                               qq(ih,jh,nt) * dbetapsi(jh,ibnd)
         END DO
      END DO
   END DO
!!omp end parallel do
   dbetapsi(:,:) = betapsi(:,:)
   ! calculate \sum_j qq(i,j)*betapsi(j)
   betapsi(:,:) = (0.0_dp, 0.0_dp)
!!omp parallel do default(shared) private(ih,ibnd,jh)
   DO ih=1,nh(nt)
      DO ibnd=nb_s, nb_e
         DO jh=1,nh(nt)
            betapsi(ih,ibnd) = betapsi(ih,ibnd) + &
                               qq(ih,jh,nt) * becp%k(ijkb0+jh,ibnd)
         END DO
      END DO
   END DO
!!omp end parallel do
   !
   ! dproj(iwf,ibnd) = \sum_ih wfatdbeta(iwf,ih)*betapsi(ih,ibnd) +
   !                           wfatbeta(iwf,ih)*dbetapsi(ih,ibnd) 
   !
   IF ( mykey == 0 .AND. nh(nt) > 0 ) THEN
      CALL ZGEMM('N','N',nwfcU, nb_e-nb_s+1, nh(nt), (1.0_dp,0.0_dp), &
           wfatdbeta, nwfcU, betapsi(1,nb_s), nh(nt),(1.0_dp,0.0_dp), &
           dproj(1,nb_s), nwfcU)
      CALL ZGEMM('N','N',nwfcU,nb_e-nb_s+1, nh(nt), (1.0_dp,0.0_dp),  &
           wfatbeta, nwfcU, dbetapsi(1,nb_s), nh(nt),(1.0_dp,0.0_dp), &
           dproj(1,nb_s), nwfcU)
   END IF
   ! end band parallelization - only dproj(1,nb_s:nb_e) are calculated
   DEALLOCATE ( betapsi )
   DEALLOCATE ( wfatbeta ) 
   DEALLOCATE (wfatdbeta )
   DEALLOCATE (dbetapsi )
   !
   call stop_clock('dprojdtau')

   RETURN
END SUBROUTINE dprojdtau_k
!
!-----------------------------------------------------------------------
SUBROUTINE dprojdtau_gamma (spsi, alpha, ijkb0, ipol, ik, nb_s, nb_e, mykey, dproj)
   !-----------------------------------------------------------------------
   !
   ! This routine computes the first derivative of the projection
   ! <\fi^{at}_{I,m1}|S|\psi_{k,v,s}> with respect to the atomic displacement
   ! u(alpha,ipol) (we remember that ns_{I,s,m1,m2} = \sum_{k,v}
   ! f_{kv} <\fi^{at}_{I,m1}|S|\psi_{k,v,s}><\psi_{k,v,s}|S|\fi^{at}_{I,m2}>)
   !
   USE kinds,                ONLY : DP
   USE ions_base,            ONLY : nat, ntyp => nsp, ityp
   USE cell_base,            ONLY : tpiba
   USE gvect,                ONLY : g
   USE klist,                ONLY : nks, xk, ngk, igk_k
   USE ldaU,                 ONLY : is_hubbard, Hubbard_l, nwfcU, wfcU, offsetU
   USE wvfct,                ONLY : nbnd, npwx,  wg
   USE uspp,                 ONLY : nkb, vkb, qq
   USE uspp_param,           ONLY : nh
   USE wavefunctions_module, ONLY : evc
   USE becmod,               ONLY : bec_type, becp, calbec
   USE mp_bands,             ONLY : intra_bgrp_comm
   USE mp_pools,             ONLY : intra_pool_comm, me_pool, nproc_pool
   USE mp,                   ONLY : mp_sum
   
   IMPLICIT NONE

   INTEGER, INTENT (IN) :: ik,      &! k-point index
                           alpha,   &! the displaced atom
                           ipol,    &! the component of displacement
                           ijkb0     ! position of beta functions for atom alpha
   INTEGER, INTENT (IN) :: nb_s, nb_e, mykey       ! band parallelization
   COMPLEX (DP), INTENT (IN) :: spsi(npwx,nbnd)   ! S|evc>
   REAL (DP), INTENT (OUT) ::  dproj(nwfcU,nb_s:nb_e) ! derivative of projection
   !
   INTEGER :: npw, nt, ig, na_, m1, ibnd, iwf, nt_, ih, jh, ldim
   REAL (DP) :: gvec
   COMPLEX (DP), ALLOCATABLE :: dwfc(:,:), dbeta(:,:)
   REAL (DP), ALLOCATABLE ::    dproj0(:,:), betapsi(:,:), dbetapsi(:,:), &
                                wfatbeta(:,:), wfatdbeta(:,:), bproj(:,:)
   !      dwfc(npwx,ldim),       ! the derivative of the atomic d wfc
   !      dbeta(npwx,nhm),       ! the derivative of the beta function
   !      betapsi(nhm,nbnd),     ! <beta|evc>
   !      dbetapsi(nhm,nbnd),    ! <dbeta|evc>
   !      wfatbeta(nwfcU,nhm),   ! <wfcU|beta>
   !      wfatdbeta(nwfcU,nhm)   ! <wfcU|dbeta>

   call start_clock('dprojdtau')
   nt = ityp(alpha)
   npw=ngk(ik)
   ldim = 2 * Hubbard_l(nt) + 1
   !
   ! At first the derivatives of the atomic wfc and the beta are computed
   ! Note: parallelization here is over plane waves, not over bands!
   !
   dproj(:,:) = 0.0_dp
   IF (is_hubbard(nt) ) THEN
      ALLOCATE ( dproj0(ldim,nbnd) )
      ALLOCATE ( dwfc(npwx,ldim) )
!!omp parallel do default(shared) private(ig,m1,gvec)
      DO ig = 1,npw
         gvec = g(ipol,igk_k(ig,ik)) * tpiba
         ! in the expression of dwfc we don't need (k+G) but just G; k always
         ! multiplies the underived quantity and gives an opposite contribution
         ! in c.c. term because the sign of the imaginary unit.
         DO m1 = 1, ldim   
            dwfc(ig,m1) = (0.d0,-1.d0) * gvec * wfcU(ig,offsetU(alpha)+m1)
         END DO
      END DO
!!omp end parallel do
      ! there is no G=0 term
      CALL DGEMM('T','N',ldim, nbnd, 2*npw, 2.0_dp,  &
                  dwfc, 2*npwx, spsi, 2*npwx, 0.0_dp,&
                  dproj0, ldim)
      DEALLOCATE ( dwfc ) 
      CALL mp_sum( dproj0, intra_bgrp_comm )
      ! copy to dproj results for the bands treated by this processor
      dproj( offsetU(alpha)+1:offsetU(alpha)+ldim, :) = dproj0(:, nb_s:nb_e)
      DEALLOCATE ( dproj0 ) 
      !
   END IF
   !
   ALLOCATE (dbetapsi(nh(nt),nbnd) ) 
   ALLOCATE (wfatdbeta(nwfcU,nh(nt)) )
   ALLOCATE ( wfatbeta(nwfcU,nh(nt)) )
   ALLOCATE ( dbeta(npwx,nh(nt)) )
!!omp parallel do default(shared) private(ih,ig)
   DO ih=1,nh(nt)
      DO ig = 1, npw
         dbeta(ig,ih) = vkb(ig,ijkb0+ih)
      END DO
   END DO
!!omp end parallel do
   CALL calbec ( npw, wfcU, dbeta, wfatbeta ) 
!!omp parallel do default(shared) private(ih,ig,gvec)
   DO ih=1,nh(nt)
      DO ig = 1, npw
         gvec = g(ipol,igk_k(ig,ik)) * tpiba
         dbeta(ig,ih) = (0.d0,-1.d0) * dbeta(ig,ih) * gvec
      END DO
   END DO
!!omp end parallel do
   CALL calbec ( npw, dbeta, evc, dbetapsi ) 
   CALL calbec ( npw, wfcU, dbeta, wfatdbeta ) 
   DEALLOCATE ( dbeta )
   !
   ! calculate \sum_j qq(i,j)*dbetapsi(j)
   ! betapsi is used here as work space 
   ALLOCATE ( betapsi(nh(nt), nbnd) ) 
   betapsi(:,:) = (0.0_dp, 0.0_dp)
   ! here starts band parallelization
!!omp parallel do default(shared) private(ih,ibnd,jh)
   DO ih=1,nh(nt)
      DO ibnd=nb_s,nb_e
         DO jh=1,nh(nt)
            betapsi(ih,ibnd) = betapsi(ih,ibnd) + &
                               qq(ih,jh,nt) * dbetapsi(jh,ibnd)
         END DO
      END DO
   END DO
!!omp end parallel do
   dbetapsi(:,:) = betapsi(:,:)
   ! calculate \sum_j qq(i,j)*betapsi(j)
   betapsi(:,:) = (0.0_dp, 0.0_dp)
!!omp parallel do default(shared) private(ih,ibnd,jh)
   DO ih=1,nh(nt)
      DO ibnd=nb_s,nb_e
         DO jh=1,nh(nt)
            betapsi(ih,ibnd) = betapsi(ih,ibnd) + &
                               qq(ih,jh,nt) * becp%r(ijkb0+jh,ibnd)
         END DO
      END DO
   END DO
!!omp end parallel do
   !
   ! dproj(iwf,ibnd) = \sum_ih wfatdbeta(iwf,ih)*betapsi(ih,ibnd) +
   !                           wfatbeta(iwf,ih)*dbetapsi(ih,ibnd) 
   !
   IF ( mykey == 0 .AND. nh(nt) > 0 ) THEN
      CALL DGEMM('N','N',nwfcU, nb_e-nb_s+1, nh(nt), 1.0_dp,  &
           wfatdbeta, nwfcU, betapsi(1,nb_s), nh(nt), 1.0_dp,&
           dproj(1,nb_s), nwfcU)
      CALL DGEMM('N','N',nwfcU, nb_e-nb_s+1, nh(nt), 1.0_dp,  &
           wfatbeta, nwfcU, dbetapsi(1,nb_s), nh(nt), 1.0_dp,&
           dproj(1,nb_s), nwfcU)
   END IF
   ! end band parallelization - only dproj(1,nb_s:nb_e) are calculated
   DEALLOCATE ( betapsi )
   DEALLOCATE ( wfatbeta ) 
   DEALLOCATE (wfatdbeta )
   DEALLOCATE (dbetapsi )
   !
   call stop_clock('dprojdtau')

   RETURN
END SUBROUTINE dprojdtau_gamma
