!
! Copyright (C) 2002-2009 Quantum ESPRESSO group
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
   USE ldaU,                 ONLY : hubbard_lmax, hubbard_l, hubbard_u, &
                                    hubbard_alpha, U_projection, &
                                    swfcatom, oatwfc
   USE symme,                ONLY : symvector
   USE io_files,             ONLY : prefix, iunocc
   USE wvfct,                ONLY : nbnd, npwx, npw, igk
   USE control_flags,        ONLY : gamma_only
   USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
   USE scf,                  ONLY : v
   USE mp_global,            ONLY : me_pool, inter_pool_comm, intra_pool_comm
   USE mp,                   ONLY : mp_sum
   USE basis,                ONLY : natomwfc
   USE becmod,               ONLY : bec_type, becp, calbec, allocate_bec_type, deallocate_bec_type
   USE uspp,                 ONLY : nkb, vkb
   USE wavefunctions_module, ONLY : evc
   USE klist,                ONLY : nks, xk, ngk
   USE io_files,             ONLY : iunigk, nwordwfc, iunwfc, &
                                    iunat, iunsat, nwordatwfc
   USE buffers,              ONLY : get_buffer

   IMPLICIT NONE
   REAL (DP) :: forceh(3,nat)  ! output: the Hubbard forces

   type (bec_type) :: proj     ! proj(natomwfc,nbnd)
   COMPLEX (DP), ALLOCATABLE :: spsi(:,:)
   !                            spsi(npwx,nbnd)
   REAL (DP), ALLOCATABLE :: dns(:,:,:,:)
   !       dns(ldim,ldim,nspin,nat) ! the derivative of the atomic occupations

   COMPLEX (DP) :: c_one, c_zero

   INTEGER :: alpha, na, nt, is, m1, m2, ipol, ldim, ik

   IF (U_projection .NE. "atomic") CALL errore("force_hub", &
                   " forces for this U_projection_type not implemented",1)

   call start_clock('force_hub')
   ldim= 2 * Hubbard_lmax + 1
   ALLOCATE ( dns(ldim,ldim,nspin,nat), spsi(npwx,nbnd) )
   call allocate_bec_type ( nkb, nbnd, becp) 
   call allocate_bec_type ( natomwfc, nbnd, proj )

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
      IF (nks > 1) READ (iunigk) igk

      CALL get_buffer (evc, nwordwfc, iunwfc, ik)
      CALL davcio(swfcatom,nwordatwfc,iunsat,ik,-1)
      CALL init_us_2 (npw,igk,xk(1,ik),vkb)
      CALL calbec( npw, swfcatom, evc, proj )
      CALL calbec( npw, vkb, evc, becp )
      CALL s_psi  (npwx, npw, nbnd, evc, spsi )

! read atomic wfc - swfcatom is used here as work space
      CALL davcio(swfcatom,nwordatwfc,iunat,ik,-1)
   
      DO ipol = 1,3
         DO alpha = 1,nat                 ! the displaced atom
            IF ( gamma_only ) THEN
               CALL dndtau_gamma(ldim,oatwfc,proj%r,swfcatom,spsi,alpha,ipol,ik,dns)
            ELSE
               CALL dndtau_k (ldim,oatwfc,proj%k,swfcatom,spsi,alpha,ipol,ik,dns)
            ENDIF
            DO na = 1,nat                 ! the Hubbard atom
               nt = ityp(na)
               IF (Hubbard_U(nt).NE.0.d0.OR. Hubbard_alpha(nt).NE.0.d0) THEN
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

#ifdef __MPI
   CALL mp_sum( forceh, inter_pool_comm )
#endif
   DEALLOCATE(dns, spsi)
   call deallocate_bec_type (proj)
   call deallocate_bec_type (becp)
   
   IF (nspin.EQ.1) forceh(:,:) = 2.d0 * forceh(:,:)
   !
   ! ...symmetrize...
   !
   CALL symvector ( nat, forceh )
   !
   RETURN
END SUBROUTINE force_hub
!
!-----------------------------------------------------------------------
SUBROUTINE dndtau_k (ldim, offset, proj, wfcatom, spsi, alpha, ipol, ik, dns)
   !-----------------------------------------------------------------------
   !
   ! This routine computes the derivative of the ns with respect to the ionic
   ! displacement u(alpha,ipol) used to obtain the Hubbard contribution to the
   ! atomic forces.
   !
   USE kinds,                ONLY : DP
   USE ions_base,            ONLY : nat, ityp
   USE basis,                ONLY : natomwfc
   USE lsda_mod,             ONLY : nspin, current_spin
   USE ldaU,                 ONLY : Hubbard_U, Hubbard_alpha, Hubbard_l
   USE wvfct,                ONLY : nbnd, npwx, npw, wg
   
   IMPLICIT NONE

   INTEGER, INTENT(IN) ::  alpha, ipol, ik, ldim, offset(nat)
   ! offset(nat): offset of d electrons of atom d in the natomwfc ordering
   COMPLEX (DP), INTENT(IN) :: &
             proj(natomwfc,nbnd), wfcatom(npwx,natomwfc), spsi(npwx,nbnd)
   REAL (DP), INTENT (OUT) :: dns(ldim,ldim,nspin,nat)
   !
   INTEGER ::  ibnd, is, na, nt, m1, m2
   COMPLEX (DP), ALLOCATABLE :: dproj(:,:)
   !
   !
   CALL start_clock('dndtau')
   !
   ALLOCATE ( dproj(natomwfc,nbnd) )
   CALL dprojdtau_k ( wfcatom, spsi, alpha, ipol, offset(alpha), dproj )
   !
   ! compute the derivative of occupation numbers (the quantities dn(m1,m2))
   ! of the atomic orbitals. They are real quantities as well as n(m1,m2)
   !
   dns(:,:,:,:) = 0.d0
   DO na = 1,nat
      nt = ityp(na)
      IF ( Hubbard_U(nt) /= 0.d0 .OR. Hubbard_alpha(nt) /= 0.d0) THEN
         DO m1 = 1, 2*Hubbard_l(nt)+1
            DO m2 = m1, 2*Hubbard_l(nt)+1
               DO ibnd = 1,nbnd
                  dns(m1,m2,current_spin,na) = dns(m1,m2,current_spin,na) + &
                                          wg(ibnd,ik) *            &
                              DBLE(  proj(offset(na)+m1,ibnd)  *   &
                             CONJG(dproj(offset(na)+m2,ibnd))  +   &
                                    dproj(offset(na)+m1,ibnd)  *   &
                             CONJG(proj(offset(na)+m2,ibnd)) )
               END DO
            END DO
         END DO
      END IF
   END DO
   DEALLOCATE ( dproj ) 
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
SUBROUTINE dndtau_gamma (ldim, offset, rproj, wfcatom, spsi, alpha, ipol, ik, dns)
   !-----------------------------------------------------------------------
   !
   ! This routine computes the derivative of the ns with respect to the ionic
   ! displacement u(alpha,ipol) used to obtain the Hubbard contribution to the
   ! atomic forces.
   !
   USE kinds,                ONLY : DP
   USE ions_base,            ONLY : nat, ityp
   USE basis,                ONLY : natomwfc
   USE lsda_mod,             ONLY : nspin, current_spin
   USE ldaU,                 ONLY : Hubbard_U, Hubbard_alpha, Hubbard_l
   USE wvfct,                ONLY : nbnd, npwx, npw, wg
   
   IMPLICIT NONE

   INTEGER, INTENT(IN) ::  alpha, ipol, ik, ldim, offset(nat)
   ! offset(nat): offset of d electrons of atom d in the natomwfc ordering
   COMPLEX (DP), INTENT(IN) ::  wfcatom(npwx,natomwfc), spsi(npwx,nbnd)
   REAL(DP), INTENT (IN) ::  rproj(natomwfc,nbnd)
   REAL (DP), INTENT (OUT) :: dns(ldim,ldim,nspin,nat)
   !
   INTEGER ::  ibnd, is, na, nt, m1, m2
   REAL (DP), ALLOCATABLE :: dproj(:,:)
   !
   !
   CALL start_clock('dndtau')
   !
   ALLOCATE ( dproj(natomwfc,nbnd) )
   CALL dprojdtau_gamma ( wfcatom, spsi, alpha, ipol, offset(alpha), dproj )
   !
   ! compute the derivative of occupation numbers (the quantities dn(m1,m2))
   ! of the atomic orbitals. They are real quantities as well as n(m1,m2)
   !
   dns(:,:,:,:) = 0.d0
   DO na = 1,nat
      nt = ityp(na)
      IF (Hubbard_U(nt) /= 0.d0 .OR. Hubbard_alpha(nt) /= 0.d0) THEN
         DO m1 = 1, 2*Hubbard_l(nt)+1
            DO m2 = m1, 2*Hubbard_l(nt)+1
               DO ibnd = 1,nbnd
                  dns(m1,m2,current_spin,na) = dns(m1,m2,current_spin,na) + &
                                          wg(ibnd,ik) * (   &
                              rproj(offset(na)+m1,ibnd)  *   &
                              dproj(offset(na)+m2,ibnd)  +   &
                              dproj(offset(na)+m1,ibnd)  *   &
                              rproj(offset(na)+m2,ibnd) )
               END DO
            END DO
         END DO
      END IF
   END DO
   DEALLOCATE ( dproj ) 
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
SUBROUTINE dprojdtau_k (wfcatom, spsi, alpha, ipol, offset, dproj)
   !-----------------------------------------------------------------------
   !
   ! This routine computes the first derivative of the projection
   ! <\fi^{at}_{I,m1}|S|\psi_{k,v,s}> with respect to the atomic displacement
   ! u(alpha,ipol) (we remember that ns_{I,s,m1,m2} = \sum_{k,v}
   ! f_{kv} <\fi^{at}_{I,m1}|S|\psi_{k,v,s}><\psi_{k,v,s}|S|\fi^{at}_{I,m2}>)
   !
   USE kinds,                ONLY : DP
   USE ions_base,            ONLY : nat, ntyp => nsp, ityp
   USE basis,                ONLY : natomwfc
   USE cell_base,            ONLY : tpiba
   USE gvect,                ONLY : g
   USE klist,                ONLY : nks, xk
   USE ldaU,                 ONLY : Hubbard_l, Hubbard_U, Hubbard_alpha
   USE wvfct,                ONLY : nbnd, npwx, npw, igk, wg
   USE uspp,                 ONLY : nkb, vkb, qq
   USE uspp_param,           ONLY : nhm, nh
   USE wavefunctions_module, ONLY : evc
   USE becmod,               ONLY : bec_type, becp, calbec
   USE mp_global,            ONLY : intra_pool_comm
   USE mp,                   ONLY : mp_sum
   
   IMPLICIT NONE
   INTEGER, INTENT (IN) :: &
              alpha,   &! the displaced atom
              ipol,    &! the component of displacement
              offset    ! the offset of the wfcs of the atom "alpha"
   COMPLEX (DP), INTENT (IN) :: &
           wfcatom(npwx,natomwfc), &! the atomic wfc
           spsi(npwx,nbnd)          ! S|evc>
   COMPLEX (DP), INTENT (OUT) :: &
           dproj(natomwfc,nbnd)     ! output: the derivative of the projection
   !
   INTEGER :: ig, jkb2, na, m1, ibnd, iwf, nt, ih, jh, ldim
   REAL (DP) :: gvec
   COMPLEX (DP), ALLOCATABLE :: dwfc(:,:), dbeta(:,:), &
                                betapsi(:,:), dbetapsi(:,:), &
                                wfatbeta(:,:), wfatdbeta(:,:)
   !      dwfc(npwx,ldim),       ! the derivative of the atomic d wfc
   !      dbeta(npwx,nhm),       ! the derivative of the beta function
   !      betapsi(nhm,nbnd),     ! <beta|evc>
   !      dbetapsi(nhm,nbnd),    ! <dbeta|evc>
   !      wfatbeta(natomwfc,nhm),! <wfc|beta>
   !      wfatdbeta(natomwfc,nhm)! <wfc|dbeta>

   call start_clock('dprojdtau')
   nt = ityp(alpha)

   ldim = 2 * Hubbard_l(nt) + 1

   dproj(:,:) = (0.d0, 0.d0)
   !
   ! At first the derivatives of the atomic wfc and the beta are computed
   !
   IF (Hubbard_U(nt) /= 0.d0 .OR. Hubbard_alpha(nt) /= 0.d0) THEN
      ALLOCATE ( dwfc(npwx,ldim) )
      DO ig = 1,npw
         gvec = g(ipol,igk(ig)) * tpiba

         ! in the expression of dwfc we don't need (k+G) but just G; k always
         ! multiplies the underived quantity and gives an opposite contribution
         ! in c.c. term because the sign of the imaginary unit.
   
         DO m1 = 1, ldim
            dwfc(ig,m1) = (0.d0,-1.d0) * gvec * wfcatom(ig,offset+m1)
         END DO
      END DO

      CALL ZGEMM('C','N',ldim, nbnd, npw, (1.d0,0.d0), &
                  dwfc, npwx, spsi, npwx, (0.d0,0.d0), &
                  dproj(offset+1,1), natomwfc)

      DEALLOCATE ( dwfc ) 
   END IF
#ifdef __MPI
   CALL mp_sum( dproj, intra_pool_comm )
#endif

   jkb2 = 0
   DO nt=1,ntyp
      DO na=1,nat
         IF ( ityp(na) .EQ. nt ) THEN
         IF ( na.EQ.alpha ) THEN
            ALLOCATE (dbetapsi(nh(nt),nbnd) ) 
            ALLOCATE (wfatdbeta(natomwfc,nh(nt)) )
            ALLOCATE ( wfatbeta(natomwfc,nh(nt)) )
            ALLOCATE ( dbeta(npwx,nh(nt)) )
            DO ih=1,nh(nt)
               DO ig = 1, npw
                  gvec = g(ipol,igk(ig)) * tpiba
                  dbeta(ig,ih) = (0.d0,-1.d0) * vkb(ig,jkb2+ih) * gvec
               END DO
            END DO
            CALL calbec ( npw, dbeta, evc, dbetapsi ) 
            CALL calbec ( npw, wfcatom, dbeta, wfatdbeta ) 
            DO ih=1,nh(nt)
               DO ig = 1, npw
                  dbeta(ig,ih) = vkb(ig,jkb2+ih)
               END DO
            END DO
            CALL calbec ( npw, wfcatom, dbeta, wfatbeta ) 
            DEALLOCATE ( dbeta )
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
                                         qq(ih,jh,nt) * becp%k(jkb2+jh,ibnd)
                  END DO
               END DO
            END DO
            !
            DO ibnd=1,nbnd
               DO ih=1,nh(nt)
                  DO iwf=1,natomwfc
                      dproj(iwf,ibnd) = dproj(iwf,ibnd) +            &
                            ( wfatdbeta(iwf,ih)*betapsi(ih,ibnd) +   &
                               wfatbeta(iwf,ih)*dbetapsi(ih,ibnd) )
                  END DO
               END DO
            END DO
            DEALLOCATE ( betapsi )
            DEALLOCATE ( wfatbeta ) 
            DEALLOCATE (wfatdbeta )
            DEALLOCATE (dbetapsi )
         END IF
         jkb2 = jkb2 + nh(nt)
         END IF
      END DO
   END DO

   call stop_clock('dprojdtau')

   RETURN
END SUBROUTINE dprojdtau_k
!
!-----------------------------------------------------------------------
SUBROUTINE dprojdtau_gamma (wfcatom, spsi, alpha, ipol, offset, dproj)
   !-----------------------------------------------------------------------
   !
   ! This routine computes the first derivative of the projection
   ! <\fi^{at}_{I,m1}|S|\psi_{k,v,s}> with respect to the atomic displacement
   ! u(alpha,ipol) (we remember that ns_{I,s,m1,m2} = \sum_{k,v}
   ! f_{kv} <\fi^{at}_{I,m1}|S|\psi_{k,v,s}><\psi_{k,v,s}|S|\fi^{at}_{I,m2}>)
   !
   USE kinds,                ONLY : DP
   USE ions_base,            ONLY : nat, ntyp => nsp, ityp
   USE basis,                ONLY : natomwfc
   USE cell_base,            ONLY : tpiba
   USE gvect,                ONLY : g, gstart
   USE klist,                ONLY : nks, xk
   USE ldaU,                 ONLY : Hubbard_l, Hubbard_U, Hubbard_alpha
   USE wvfct,                ONLY : nbnd, npwx, npw, igk, wg
   USE uspp,                 ONLY : nkb, vkb, qq
   USE uspp_param,           ONLY : nhm, nh
   USE wavefunctions_module, ONLY : evc
   USE becmod,               ONLY : bec_type, becp, calbec
   USE mp_global,            ONLY : intra_pool_comm
   USE mp,                   ONLY : mp_sum
   
   IMPLICIT NONE
   INTEGER, INTENT (IN) :: &
              alpha,   &! the displaced atom
              ipol,    &! the component of displacement
              offset    ! the offset of the wfcs of the atom "alpha"
   COMPLEX (DP), INTENT (IN) :: &
           wfcatom(npwx,natomwfc), &! the atomic wfc
           spsi(npwx,nbnd)          ! S|evc>
   REAL (DP), INTENT (OUT) :: &
           dproj(natomwfc,nbnd)     ! output: the derivative of the projection
   !
   INTEGER :: ig, jkb2, na, m1, ibnd, iwf, nt, ih, jh, ldim
   REAL (DP) :: gvec
   COMPLEX (DP), ALLOCATABLE :: dwfc(:,:), dbeta(:,:)
   REAL (DP), ALLOCATABLE ::    betapsi(:,:), dbetapsi(:,:), &
                                wfatbeta(:,:), wfatdbeta(:,:)
   !      dwfc(npwx,ldim),       ! the derivative of the atomic d wfc
   !      dbeta(npwx,nhm),       ! the derivative of the beta function
   !      betapsi(nhm,nbnd),     ! <beta|evc>
   !      dbetapsi(nhm,nbnd),    ! <dbeta|evc>
   !      wfatbeta(natomwfc,nhm),! <wfc|beta>
   !      wfatdbeta(natomwfc,nhm)! <wfc|dbeta>

   call start_clock('dprojdtau')
   nt = ityp(alpha)

   ldim = 2 * Hubbard_l(nt) + 1

   ALLOCATE ( dwfc(npwx,ldim) )

   dproj(:,:) = (0.d0, 0.d0)
   !
   ! At first the derivatives of the atomic wfc and the beta are computed
   !
   IF (Hubbard_U(nt) /= 0.d0 .OR. Hubbard_alpha(nt) /= 0.d0) THEN
      DO ig = 1,npw
         gvec = g(ipol,igk(ig)) * tpiba

         ! in the expression of dwfc we don't need (k+G) but just G; k always
         ! multiplies the underived quantity and gives an opposite contribution
         ! in c.c. term because the sign of the imaginary unit.
   
         DO m1 = 1, ldim
            dwfc(ig,m1) = (0.d0,-1.d0) * gvec * wfcatom(ig,offset+m1)
         END DO
      END DO
      ! there is no G=0 term
      CALL DGEMM('T','N',ldim, nbnd, 2*npw, 2.0_dp,  &
                  dwfc, 2*npwx, spsi, 2*npwx, 0.0_dp,&
                  dproj(offset+1,1), natomwfc)
      DEALLOCATE ( dwfc ) 
   END IF

#ifdef __MPI
   CALL mp_sum( dproj, intra_pool_comm )
#endif
   jkb2 = 0
   DO nt=1,ntyp
      DO na=1,nat
         IF ( ityp(na) .EQ. nt ) THEN
         IF ( na.EQ.alpha ) THEN
            ALLOCATE (dbetapsi(nh(nt),nbnd) ) 
            ALLOCATE (wfatdbeta(natomwfc,nh(nt)) )
            ALLOCATE ( wfatbeta(natomwfc,nh(nt)) )
            ALLOCATE ( dbeta(npwx,nh(nt)) )
            DO ih=1,nh(nt)
               DO ig = 1, npw
                  gvec = g(ipol,igk(ig)) * tpiba
                  dbeta(ig,ih) = (0.d0,-1.d0) * vkb(ig,jkb2+ih) * gvec
               END DO
            END DO
            CALL calbec ( npw, dbeta, evc, dbetapsi ) 
            CALL calbec ( npw, wfcatom, dbeta, wfatdbeta ) 
            DO ih=1,nh(nt)
               DO ig = 1, npw
                  dbeta(ig,ih) = vkb(ig,jkb2+ih)
               END DO
            END DO
            CALL calbec ( npw, wfcatom, dbeta, wfatbeta ) 
            DEALLOCATE ( dbeta )
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
                                         qq(ih,jh,nt) * becp%r(jkb2+jh,ibnd)
                  END DO
               END DO
            END DO
            !
            DO ibnd=1,nbnd
               DO ih=1,nh(nt)
                  DO iwf=1,natomwfc
                        dproj(iwf,ibnd) = dproj(iwf,ibnd) +           &
                             ( wfatdbeta(iwf,ih)*betapsi(ih,ibnd) +   &
                                wfatbeta(iwf,ih)*dbetapsi(ih,ibnd) )
                  END DO
               END DO
            END DO
            DEALLOCATE ( betapsi )
            DEALLOCATE ( wfatbeta ) 
            DEALLOCATE (wfatdbeta )
            DEALLOCATE (dbetapsi )
         END IF
         jkb2 = jkb2 + nh(nt)
         END IF
      END DO
   END DO

   call stop_clock('dprojdtau')

   RETURN
END SUBROUTINE dprojdtau_gamma
