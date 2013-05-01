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
   USE ions_base,            ONLY : nat, ityp
   USE cell_base,            ONLY : at, bg
   USE ldaU,                 ONLY : hubbard_lmax, hubbard_l, U_projection, &
                                    nwfcU, wfcU, is_hubbard, lda_plus_u_kind, &
                                    oatwfc, copy_U_wfc, offsetU
   USE basis,                ONLY : swfcatom
   USE symme,                ONLY : symvector
   USE io_files,             ONLY : prefix
   USE wvfct,                ONLY : nbnd, npwx, npw, igk
   USE control_flags,        ONLY : gamma_only
   USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
   USE scf,                  ONLY : v
   USE mp_global,            ONLY : inter_pool_comm
   USE mp,                   ONLY : mp_sum
   USE becmod,               ONLY : bec_type, becp, calbec, allocate_bec_type, deallocate_bec_type
   USE uspp,                 ONLY : nkb, vkb
   USE wavefunctions_module, ONLY : evc
   USE klist,                ONLY : nks, xk, ngk
   USE io_files,             ONLY : iunigk, nwordwfc, iunwfc, &
                                    iunhub, nwordwfcU, nwordatwfc
   USE buffers,              ONLY : get_buffer

   IMPLICIT NONE
   REAL (DP) :: forceh(3,nat)  ! output: the Hubbard forces

   type (bec_type) :: proj     ! proj(nwfcU,nbnd)
   COMPLEX (DP), ALLOCATABLE :: spsi(:,:)
   !                            spsi(npwx,nbnd)
   REAL (DP), ALLOCATABLE :: dns(:,:,:,:)
   !       dns(ldim,ldim,nspin,nat) ! the derivative of the atomic occupations
   INTEGER :: alpha, na, nt, is, m1, m2, ipol, ldim, ik

   IF (U_projection .NE. "atomic") CALL errore("force_hub", &
                   " forces for this U_projection_type not implemented",1)
   IF (lda_plus_u_kind == 1) CALL errore("force_hub", &
                   " forces in full LDA+U scheme are not yet implemented",1)

   call start_clock('force_hub')
   ldim= 2 * Hubbard_lmax + 1
   ALLOCATE ( dns(ldim,ldim,nspin,nat), spsi(npwx,nbnd) )
   call allocate_bec_type ( nkb, nbnd, becp) 
   call allocate_bec_type ( nwfcU, nbnd, proj )

   forceh(:,:) = 0.d0

   ! Offset of atomic wavefunctions initialized in setup and stored in oatwfc
   !
   !    we start a loop on k points
   !
   IF (nks > 1) REWIND (iunigk)
   DO ik = 1, nks
      IF (lsda) current_spin = isk(ik)
      !
      ! now we need the first derivative of proj with respect to tau(alpha,ipol)
      !
      npw = ngk (ik)
      IF (nks > 1) THEN
         READ (iunigk) igk
         CALL get_buffer (evc, nwordwfc, iunwfc, ik)
      END IF
      CALL get_buffer (wfcU, nwordwfcU, iunhub, ik)
      CALL init_us_2 (npw,igk,xk(1,ik),vkb)
      CALL calbec( npw, wfcU, evc, proj )
      CALL calbec( npw, vkb, evc, becp )
      CALL s_psi  (npwx, npw, nbnd, evc, spsi )

! re-calculate atomic wfc - swfcatom is used here as work space
! (has to be modified for noncolinear case)

      CALL atomic_wfc (ik, swfcatom)
      call copy_U_wfc (swfcatom)

      DO ipol = 1,3
         DO alpha = 1,nat                 ! the displaced atom
            IF ( gamma_only ) THEN
               CALL dndtau_gamma(ldim,proj%r,nwfcU,wfcU,offsetU,spsi,alpha,ipol,ik,dns)
            ELSE
               CALL dndtau_k (ldim,proj%k,nwfcU,wfcU,offsetU,spsi,alpha,ipol,ik,dns)
            ENDIF
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
         END DO
      END DO
   END DO
   !
   CALL mp_sum( forceh, inter_pool_comm )
   !
   DEALLOCATE(dns, spsi)
   call deallocate_bec_type (proj)
   call deallocate_bec_type (becp)
   
   IF (nspin == 1) forceh(:,:) = 2.d0 * forceh(:,:)
   !
   ! ...symmetrize...
   !
   CALL symvector ( nat, forceh )
!write(66,'("Hubbard contribution Begin")')
!write(66,'(3f12.6)') forceh(:,:)
!write(66,'("Hubbard contribution End")')
   !
   call stop_clock('force_hub')
   !!!
   call print_clock('force_hub')
   call print_clock('atomic_wfc')
   call print_clock('dndtau')
   call print_clock('dprojdtau')
   call print_clock('dprojdtau:1')
   call print_clock('dprojdtau:2')
   call print_clock('dprojdtau:3')
   call print_clock('dprojdtau:4')
   !
   RETURN
END SUBROUTINE force_hub
!
!-----------------------------------------------------------------------
SUBROUTINE dndtau_k (ldim, proj, nwfcU, wfcU, offsetU, &
     spsi, alpha, ipol, ik, dns)
   !-----------------------------------------------------------------------
   !
   ! This routine computes the derivative of the ns with respect to the ionic
   ! displacement u(alpha,ipol) used to obtain the Hubbard contribution to the
   ! atomic forces.
   !
   USE kinds,                ONLY : DP
   USE ions_base,            ONLY : nat, ityp
   USE lsda_mod,             ONLY : nspin, current_spin
   USE ldaU,                 ONLY : is_hubbard, Hubbard_l
   USE wvfct,                ONLY : nbnd, npwx, npw, wg
   USE mp_pools,             ONLY : intra_pool_comm, me_pool, nproc_pool
   USE mp,                   ONLY : mp_sum
   
   IMPLICIT NONE

   INTEGER, INTENT(IN) :: nwfcU, alpha, ipol, ik, ldim, offsetU(nat)
   ! offsetU(nat): offset of d electrons 
   COMPLEX (DP), INTENT(IN) :: &
             proj(nwfcU,nbnd), wfcU(npwx,nwfcU), spsi(npwx,nbnd)
   REAL (DP), INTENT (OUT) :: dns(ldim,ldim,nspin,nat)
   !
   INTEGER ::  ibnd, is, na, nt, m1, m2, na_s, na_e, mykey
   COMPLEX (DP), ALLOCATABLE :: dproj(:,:)
   !
   !
   CALL start_clock('dndtau')
   !
   ALLOCATE ( dproj(nwfcU,nbnd) )
   CALL dprojdtau_k ( nwfcU, wfcU, offsetU, spsi, alpha, ipol, dproj )
   !
   ! compute the derivative of occupation numbers (the quantities dn(m1,m2))
   ! of the atomic orbitals. They are real quantities as well as n(m1,m2)
   !
   ! poor-man parallelization over atoms
   ! - if nproc_pool=1   : na_s=1, na_e=nat, mykey=0
   ! - if nproc_pool<=nat: each processor calculates atom na_s to na_e; mykey=0
   ! - if nproc_pool> nat: each processor takes care of atom na_s=na_e;
   !   mykey labels how many times each atom appears (mykey=0 first time etc.)
   !
   CALL block_distribute( nat, me_pool, nproc_pool, na_s, na_e, mykey )
   dns(:,:,:,:) = 0.d0
   DO na = na_s, na_e
      nt = ityp(na)
      IF ( is_hubbard(nt) .AND. mykey == 0 ) THEN
         DO m1 = 1, 2*Hubbard_l(nt)+1
            DO m2 = m1, 2*Hubbard_l(nt)+1
!$omp parallel do default(shared) private(ibnd)
               DO ibnd = 1,nbnd
                  dns(m1,m2,current_spin,na) = dns(m1,m2,current_spin,na) + &
                                          wg(ibnd,ik) *            &
                              DBLE( proj(offsetU(na)+m1,ibnd)  *   &
                             CONJG(dproj(offsetU(na)+m2,ibnd))  +   &
                                   dproj(offsetU(na)+m1,ibnd)  *   &
                             CONJG( proj(offsetU(na)+m2,ibnd)) )
               END DO
!$omp end parallel do
            END DO
         END DO
      END IF
   END DO
   DEALLOCATE ( dproj ) 
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
   DO na = 1,nat
      DO is = 1,nspin
         DO m1 = 1,ldim
            DO m2 = m1+1,ldim
               dns(m2,m1,is,na) = dns(m1,m2,is,na)
            END DO
         END DO
      END DO
   END DO

   CALL stop_clock('dndtau')
   RETURN
END SUBROUTINE dndtau_k
!
!-----------------------------------------------------------------------
SUBROUTINE dndtau_gamma (ldim, rproj, nwfcU, wfcU, offsetU, &
     spsi, alpha, ipol, ik, dns)
  !-----------------------------------------------------------------------
   !
   ! This routine computes the derivative of the ns with respect to the ionic
   ! displacement u(alpha,ipol) used to obtain the Hubbard contribution to the
   ! atomic forces.
   !
   USE kinds,                ONLY : DP
   USE ions_base,            ONLY : nat, ityp
   USE lsda_mod,             ONLY : nspin, current_spin
   USE ldaU,                 ONLY : is_hubbard, Hubbard_l
   USE wvfct,                ONLY : nbnd, npwx, npw, wg
   USE mp_pools,             ONLY : intra_pool_comm, me_pool, nproc_pool
   USE mp,                   ONLY : mp_sum

   IMPLICIT NONE

   INTEGER, INTENT(IN) ::  alpha, ipol, ik, ldim, nwfcU, offsetU(nat)
   ! offsetU(nat): offset of d electrons of atom 
   COMPLEX (DP), INTENT(IN) :: wfcU(npwx,nwfcU), spsi(npwx,nbnd)
   REAL(DP), INTENT (IN) ::  rproj(nwfcU,nbnd)
   REAL (DP), INTENT (OUT) :: dns(ldim,ldim,nspin,nat)
   !
   INTEGER ::  ibnd, is, na, nt, m1, m2, na_s, na_e, mykey
   REAL (DP), ALLOCATABLE :: dproj(:,:)
   !
   !
   CALL start_clock('dndtau')
   !
   ALLOCATE ( dproj(nwfcU,nbnd) )
   CALL dprojdtau_gamma ( nwfcU, wfcU, offsetU, spsi, alpha, ipol, dproj )
   !
   ! compute the derivative of occupation numbers (the quantities dn(m1,m2))
   ! of the atomic orbitals. They are real quantities as well as n(m1,m2)
   !
   ! poor-man parallelization over atoms
   ! - if nproc_pool=1   : na_s=1, na_e=nat, mykey=0
   ! - if nproc_pool<=nat: each processor calculates atom na_s to na_e; mykey=0
   ! - if nproc_pool> nat: each processor takes care of atom na_s=na_e;
   !   mykey labels how many times each atom appears (mykey=0 first time etc.)
   !
   CALL block_distribute( nat, me_pool, nproc_pool, na_s, na_e, mykey )
   dns(:,:,:,:) = 0.d0
   DO na = na_s, na_e
      nt = ityp(na)
      IF (is_hubbard(nt) .AND. mykey == 0 ) THEN
         DO m1 = 1, 2*Hubbard_l(nt)+1
            DO m2 = m1, 2*Hubbard_l(nt)+1
!$omp parallel do default(shared) private(ibnd)
               DO ibnd = 1,nbnd
                  dns(m1,m2,current_spin,na) = dns(m1,m2,current_spin,na) + &
                                          wg(ibnd,ik) * (   &
                              rproj(offsetU(na)+m1,ibnd)  *   &
                              dproj(offsetU(na)+m2,ibnd)  +   &
                              dproj(offsetU(na)+m1,ibnd)  *   &
                              rproj(offsetU(na)+m2,ibnd) )
               END DO
!$omp end parallel do
            END DO
         END DO
      END IF
   END DO
   DEALLOCATE ( dproj ) 
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
   DO na = 1,nat
      DO is = 1,nspin
         DO m1 = 1,ldim
            DO m2 = m1+1,ldim
               dns(m2,m1,is,na) = dns(m1,m2,is,na)
            END DO
         END DO
      END DO
   END DO

   CALL stop_clock('dndtau')
   RETURN
END SUBROUTINE dndtau_gamma
!
!-----------------------------------------------------------------------
SUBROUTINE dprojdtau_k (nwfcU, wfcU, offsetU, spsi, alpha, ipol, dproj)
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
   USE klist,                ONLY : nks, xk
   USE ldaU,                 ONLY : is_hubbard, Hubbard_l
   USE wvfct,                ONLY : nbnd, npwx, npw, igk, wg
   USE uspp,                 ONLY : nkb, vkb, qq
   USE uspp_param,           ONLY : nhm, nh
   USE wavefunctions_module, ONLY : evc
   USE becmod,               ONLY : bec_type, becp, calbec
   USE mp_global,            ONLY : intra_bgrp_comm
   USE mp,                   ONLY : mp_sum
   
   IMPLICIT NONE
   INTEGER, INTENT (IN) :: &
              alpha,   &! the displaced atom
              ipol,    &! the component of displacement
              nwfcU, offsetU(nat)
   COMPLEX (DP), INTENT (IN) :: &
           wfcU(npwx,nwfcU), &! the atomic wfc
           spsi(npwx,nbnd)          ! S|evc>
   COMPLEX (DP), INTENT (OUT) :: &
           dproj(nwfcU,nbnd)     ! output: the derivative of the projection
   !
   INTEGER :: nt, ig, ijkb0, na_, m1, ibnd, iwf, nt_, ih, jh, ldim
   REAL (DP) :: gvec
   COMPLEX (DP), ALLOCATABLE :: dwfc(:,:), dbeta(:,:), &
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

   ldim = 2 * Hubbard_l(nt) + 1

   dproj(:,:) = (0.d0, 0.d0)
   !
   ! At first the derivatives of the atomic wfc and the beta are computed
   !
   call start_clock('dprojdtau:1')
   IF ( is_hubbard(nt) ) THEN
      ALLOCATE ( dwfc(npwx,ldim) )
      DO ig = 1,npw
         gvec = g(ipol,igk(ig)) * tpiba

         ! in the expression of dwfc we don't need (k+G) but just G; k always
         ! multiplies the underived quantity and gives an opposite contribution
         ! in c.c. term because the sign of the imaginary unit.
   
         DO m1 = 1, ldim
            dwfc(ig,m1) = (0.d0,-1.d0) * gvec * wfcU(ig,offsetU(alpha)+m1)
         END DO
      END DO

      CALL ZGEMM('C','N',ldim, nbnd, npw, (1.d0,0.d0), &
                  dwfc, npwx, spsi, npwx, (0.d0,0.d0), &
                  dproj(offsetU(alpha)+1,1), nwfcU)

      DEALLOCATE ( dwfc ) 
      CALL mp_sum( dproj, intra_bgrp_comm )
   END IF
   call stop_clock('dprojdtau:1')
   !
   ! FIXME: ijkb0 should be calculated and stored once for all
   ijkb0 = 0
   DO nt_=1,ntyp
      DO na_=1,nat
         IF ( ityp(na_) .EQ. nt_ ) THEN
            IF ( na_ == alpha ) GO TO 10
            ijkb0 = ijkb0 + nh(nt_)
         END IF
      END DO
   END DO
10 CONTINUE
   ! 
   ! ijkb0 points now to the beta functions of atom alpha
   !
   call start_clock('dprojdtau:2')
   ALLOCATE (dbetapsi(nh(nt),nbnd) ) 
   ALLOCATE (wfatdbeta(nwfcU,nh(nt)) )
   ALLOCATE ( wfatbeta(nwfcU,nh(nt)) )
   ALLOCATE ( dbeta(npwx,nh(nt)) )
   DO ih=1,nh(nt)
      DO ig = 1, npw
         dbeta(ig,ih) = vkb(ig,ijkb0+ih)
      END DO
   END DO
   CALL calbec ( npw, wfcU, dbeta, wfatbeta ) 
   DO ih=1,nh(nt)
      DO ig = 1, npw
         gvec = g(ipol,igk(ig)) * tpiba
         dbeta(ig,ih) = (0.d0,-1.d0) * dbeta(ig,ih) * gvec
      END DO
   END DO
   CALL calbec ( npw, dbeta, evc, dbetapsi ) 
   CALL calbec ( npw, wfcU, dbeta, wfatdbeta ) 
   DEALLOCATE ( dbeta )
   call stop_clock('dprojdtau:2')
   call start_clock('dprojdtau:3')
   ! calculate \sum_j qq(i,j)*dbetapsi(j)
   ! betapsi is used here as work space 
   ALLOCATE ( betapsi(nh(nt), nbnd) ) 
   betapsi(:,:) = (0.0_dp, 0.0_dp)
   DO ih=1,nh(nt)
      DO ibnd=1,nbnd
         DO jh=1,nh(nt)
            betapsi(ih,ibnd) = betapsi(ih,ibnd) + &
                               qq(ih,jh,nt) * dbetapsi(jh,ibnd)
         END DO
      END DO
   END DO
   dbetapsi(:,:) = betapsi(:,:)
   ! calculate \sum_j qq(i,j)*betapsi(j)
   betapsi(:,:) = (0.0_dp, 0.0_dp)
   DO ih=1,nh(nt)
      DO ibnd=1,nbnd
         DO jh=1,nh(nt)
            betapsi(ih,ibnd) = betapsi(ih,ibnd) + &
                               qq(ih,jh,nt) * becp%k(ijkb0+jh,ibnd)
         END DO
      END DO
   END DO
   call stop_clock('dprojdtau:3')
   !
   ! dproj(iwf,ibnd) = \sum_ih wfatdbeta(iwf,ih)*betapsi(ih,ibnd) +
   !                           wfatbeta(iwf,ih)*dbetapsi(ih,ibnd) 
   !
   call start_clock('dprojdtau:4')
   CALL ZGEMM('N','N',nwfcU, nbnd, nh(nt), 1.0_dp,  &
        wfatdbeta, nwfcU, betapsi, nh(nt), 1.0_dp,&
        dproj, nwfcU)
   CALL ZGEMM('N','N',nwfcU, nbnd, nh(nt), 1.0_dp,  &
        wfatbeta, nwfcU, dbetapsi, nh(nt), 1.0_dp,&
        dproj, nwfcU)
   !
   DEALLOCATE ( betapsi )
   DEALLOCATE ( wfatbeta ) 
   DEALLOCATE (wfatdbeta )
   DEALLOCATE (dbetapsi )
   call stop_clock('dprojdtau:4')

   call stop_clock('dprojdtau')

   RETURN
END SUBROUTINE dprojdtau_k
!
!-----------------------------------------------------------------------
SUBROUTINE dprojdtau_gamma (nwfcU, wfcU, offsetU, spsi, alpha, ipol, dproj)
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
   USE klist,                ONLY : nks, xk
   USE ldaU,                 ONLY : is_hubbard, Hubbard_l
   USE wvfct,                ONLY : nbnd, npwx, npw, igk, wg
   USE uspp,                 ONLY : nkb, vkb, qq
   USE uspp_param,           ONLY : nhm, nh
   USE wavefunctions_module, ONLY : evc
   USE becmod,               ONLY : bec_type, becp, calbec
   USE mp_global,            ONLY : intra_bgrp_comm
   USE mp,                   ONLY : mp_sum
   
   IMPLICIT NONE

   INTEGER, INTENT (IN) :: &
              alpha,   &! the displaced atom
              ipol,    &! the component of displacement
              nwfcU, offsetU(nat)
   COMPLEX (DP), INTENT (IN) :: &
           wfcU(npwx,nwfcU),     &! the atomic wfc
           spsi(npwx,nbnd)         ! S|evc>
   REAL (DP), INTENT (OUT) :: &
           dproj(nwfcU,nbnd)     ! output: the derivative of the projection
   !
   INTEGER :: nt, ig, ijkb0, na_, m1, ibnd, iwf, nt_, ih, jh, ldim
   REAL (DP) :: gvec
   COMPLEX (DP), ALLOCATABLE :: dwfc(:,:), dbeta(:,:)
   REAL (DP), ALLOCATABLE ::    betapsi(:,:), dbetapsi(:,:), &
                                wfatbeta(:,:), wfatdbeta(:,:)
   !      dwfc(npwx,ldim),       ! the derivative of the atomic d wfc
   !      dbeta(npwx,nhm),       ! the derivative of the beta function
   !      betapsi(nhm,nbnd),     ! <beta|evc>
   !      dbetapsi(nhm,nbnd),    ! <dbeta|evc>
   !      wfatbeta(nwfcU,nhm),   ! <wfcU|beta>
   !      wfatdbeta(nwfcU,nhm)   ! <wfcU|dbeta>

   call start_clock('dprojdtau')
   nt = ityp(alpha)
   ldim = 2 * Hubbard_l(nt) + 1

   dproj(:,:) = 0.0_dp
   !
   ! At first the derivatives of the atomic wfc and the beta are computed
   !
   call start_clock('dprojdtau:1')
   IF (is_hubbard(nt) ) THEN
      ALLOCATE ( dwfc(npwx,ldim) )
!$omp parallel do default(shared) private(m1,ig,gvec)
      DO m1 = 1, ldim
         DO ig = 1,npw
            gvec = g(ipol,igk(ig)) * tpiba

         ! in the expression of dwfc we don't need (k+G) but just G; k always
         ! multiplies the underived quantity and gives an opposite contribution
         ! in c.c. term because the sign of the imaginary unit.
   
            dwfc(ig,m1) = (0.d0,-1.d0) * gvec * wfcU(ig,offsetU(alpha)+m1)
         END DO
      END DO
!$omp end parallel do
      ! there is no G=0 term
      CALL DGEMM('T','N',ldim, nbnd, 2*npw, 2.0_dp,  &
                  dwfc, 2*npwx, spsi, 2*npwx, 0.0_dp,&
                  dproj(offsetU(alpha)+1,1), nwfcU)
      DEALLOCATE ( dwfc ) 
      CALL mp_sum( dproj, intra_bgrp_comm )
   END IF
   call stop_clock('dprojdtau:1')
   !
   ! FIXME: ijkb0 should be calculated and stored once for all
   ijkb0 = 0
   DO nt_=1,ntyp
      DO na_=1,nat
         IF ( ityp(na_) .EQ. nt_ ) THEN
            IF ( na_ == alpha ) GO TO 10
            ijkb0 = ijkb0 + nh(nt_)
         END IF
      END DO
   END DO
10 CONTINUE
   ! 
   ! ijkb0 points now to the beta functions of atom alpha
   !
   call start_clock('dprojdtau:2')
   ALLOCATE (dbetapsi(nh(nt),nbnd) ) 
   ALLOCATE (wfatdbeta(nwfcU,nh(nt)) )
   ALLOCATE ( wfatbeta(nwfcU,nh(nt)) )
   ALLOCATE ( dbeta(npwx,nh(nt)) )
!$omp parallel do default(shared) private(ih,ig)
   DO ih=1,nh(nt)
      DO ig = 1, npw
         dbeta(ig,ih) = vkb(ig,ijkb0+ih)
      END DO
   END DO
!$omp end parallel do
   CALL calbec ( npw, wfcU, dbeta, wfatbeta ) 
!$omp parallel do default(shared) private(ih,ig,gvec)
   DO ih=1,nh(nt)
      DO ig = 1, npw
         gvec = g(ipol,igk(ig)) * tpiba
         dbeta(ig,ih) = (0.d0,-1.d0) * dbeta(ig,ih) * gvec
      END DO
   END DO
!$omp end parallel do
   CALL calbec ( npw, dbeta, evc, dbetapsi ) 
   CALL calbec ( npw, wfcU, dbeta, wfatdbeta ) 
   DEALLOCATE ( dbeta )
   call stop_clock('dprojdtau:2')
   ! calculate \sum_j qq(i,j)*dbetapsi(j)
   ! betapsi is used here as work space 
   call start_clock('dprojdtau:3')
   ALLOCATE ( betapsi(nh(nt), nbnd) ) 
   betapsi(:,:) = (0.0_dp, 0.0_dp)
!$omp parallel do default(shared) private(ih,ibnd,jh)
   DO ih=1,nh(nt)
      DO ibnd=1,nbnd
         DO jh=1,nh(nt)
            betapsi(ih,ibnd) = betapsi(ih,ibnd) + &
                               qq(ih,jh,nt) * dbetapsi(jh,ibnd)
         END DO
      END DO
   END DO
!$omp end parallel do
   dbetapsi(:,:) = betapsi(:,:)
   ! calculate \sum_j qq(i,j)*betapsi(j)
   betapsi(:,:) = (0.0_dp, 0.0_dp)
!$omp parallel do default(shared) private(ih,ibnd,jh)
   DO ih=1,nh(nt)
      DO ibnd=1,nbnd
         DO jh=1,nh(nt)
            betapsi(ih,ibnd) = betapsi(ih,ibnd) + &
                               qq(ih,jh,nt) * becp%r(ijkb0+jh,ibnd)
         END DO
      END DO
   END DO
!$omp end parallel do
   call stop_clock('dprojdtau:3')
   !
   ! dproj(iwf,ibnd) = \sum_ih wfatdbeta(iwf,ih)*betapsi(ih,ibnd) +
   !                           wfatbeta(iwf,ih)*dbetapsi(ih,ibnd) 
   !
   call start_clock('dprojdtau:4')
   CALL DGEMM('N','N',nwfcU, nbnd, nh(nt), 1.0_dp,  &
        wfatdbeta, nwfcU, betapsi, nh(nt), 1.0_dp,&
        dproj, nwfcU)
   CALL DGEMM('N','N',nwfcU, nbnd, nh(nt), 1.0_dp,  &
        wfatbeta, nwfcU, dbetapsi, nh(nt), 1.0_dp,&
        dproj, nwfcU)
   !
   DEALLOCATE ( betapsi )
   DEALLOCATE ( wfatbeta ) 
   DEALLOCATE (wfatdbeta )
   DEALLOCATE (dbetapsi )
   call stop_clock('dprojdtau:4')

   call stop_clock('dprojdtau')

   RETURN
END SUBROUTINE dprojdtau_gamma
