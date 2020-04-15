!
! Copyright (C) 2002-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE stres_hub ( sigmah )
   !----------------------------------------------------------------------
   !! This routines computes the Hubbard contribution to the internal stress
   !! tensor. It gives in output the array sigmah(i,j) which corresponds to
   !! the quantity \( -(1/\omega)dE_h/\epsilon_{i,j} \)
   !
   USE kinds,         ONLY : DP
   USE wavefunctions, ONLY : evc
   USE ions_base,     ONLY : nat, ityp, ntyp => nsp
   USE cell_base,     ONLY : omega, at, bg
   USE wvfct,         ONLY : nbnd, npwx
   USE ldaU,          ONLY : Hubbard_lmax, Hubbard_l, is_hubbard, &
                             lda_plus_u_kind, U_projection, is_hubbard_back, &
                             ldim_back, ldmx_b, nsg, v_nsg, max_num_neighbors, &
                             ldim_u, Hubbard_V, at_sc, neighood, ldmx_tot, &
                             wfcU, nwfcU, copy_U_wfc
   USE becmod,        ONLY : bec_type, becp, calbec, allocate_bec_type, deallocate_bec_type
   USE lsda_mod,      ONLY : lsda, nspin, current_spin, isk
   USE uspp,          ONLY : nkb, vkb
   USE klist,         ONLY : nks, xk, ngk, igk_k
   USE basis,         ONLY : natomwfc
   USE io_files,      ONLY : nwordwfc, iunwfc
   USE buffers,       ONLY : get_buffer
   USE scf,           ONLY : v, rho
   USE symme,         ONLY : symmatrix
   USE io_global,     ONLY : stdout
   USE mp_pools,      ONLY : inter_pool_comm, me_pool, nproc_pool
   USE mp,            ONLY : mp_sum
   USE control_flags, ONLY : gamma_only
   USE mp_bands,      ONLY : use_bgrp_in_hpsi
   !
   IMPLICIT NONE
   !
   REAL(DP), INTENT(OUT) :: sigmah(3,3) 
   !! the Hubbard contribution to stresses
   !
   ! ... local variables
   !
   INTEGER :: ipol, jpol, na, nt, is, m1, m2, na1, nt1, na2, nt2, viz, ik, npw
   INTEGER :: ldim, ldim1, ldim2, ldimb, equiv_na2, i_type, nb_s, nb_e, mykey
   REAL(DP), ALLOCATABLE :: dns(:,:,:,:), dnsb(:,:,:,:)
   COMPLEX(DP), ALLOCATABLE ::  dnsg(:,:,:,:,:)
   !! the derivative of the atomic occupations
   COMPLEX(DP), ALLOCATABLE :: spsi(:,:), wfcatom(:,:)
   TYPE (bec_type) :: proj
   INTEGER, EXTERNAL :: type_interaction
   LOGICAL :: lhubb
   LOGICAL :: save_flag
   !
   CALL start_clock( 'stres_hub' )
   !
   save_flag = use_bgrp_in_hpsi ; use_bgrp_in_hpsi = .false.
   !
   IF (U_projection .NE. "atomic") CALL errore("stres_hub", &
                   " stress for this U_projection_type not implemented",1)
   IF (lda_plus_u_kind.EQ.1) CALL errore("stres_hub", &
                   " stress in non collinear LDA+U scheme is not yet implemented",1)
   !
   sigmah(:,:) = 0.d0
   !
   ALLOCATE ( wfcatom (npwx,natomwfc) )
   ALLOCATE ( spsi(npwx,nbnd) )
   !
   CALL allocate_bec_type( nkb,   nbnd, becp )
   CALL allocate_bec_type( nwfcU, nbnd, proj )
   !
   ! poor-man parallelization over bands
   ! - if nproc_pool=1   : nb_s=1, nb_e=nbnd, mykey=0
   ! - if nproc_pool<=nbnd:each processor calculates band nb_s to nb_e; mykey=0
   ! - if nproc_pool>nbnd :each processor takes care of band na_s=nb_e;
   !   mykey labels how many times each band appears (mykey=0 first time etc.)
   !
   CALL block_distribute( nbnd, me_pool, nproc_pool, nb_s, nb_e, mykey )
   !
   IF (lda_plus_u_kind.EQ.0) THEN
      ldim = 2 * Hubbard_lmax + 1
      ALLOCATE (dns(ldim,ldim,nspin,nat))
      lhubb = .FALSE.
      DO nt = 1, ntyp
         IF (is_hubbard_back(nt)) lhubb = .TRUE.
      ENDDO
      IF (lhubb) THEN
         ldimb = ldmx_b
         ALLOCATE ( dnsb(ldimb,ldimb,nspin,nat) )
      ENDIF
   ELSEIF (lda_plus_u_kind.EQ.2) THEN
      ldim = ldmx_tot
      ALLOCATE(dnsg(ldim,ldim,max_num_neighbors,nat,nspin))
   ENDIF
   !
#ifdef DEBUG
   IF (lda_plus_u_kind.EQ.0) THEN
      DO na = 1, nat
         DO is = 1, nspin
            nt = ityp(na)
            IF ( is_hubbard(nt) ) THEN
               WRITE( stdout,'(a,2i3)') 'VNS(NA,IS) ', na, is
               DO m1 = 1, ldim
                  WRITE( stdout,'(7f10.4)') (v%ns(m1,m2,is,na),m2=1,ldim)
               ENDDO
            ENDIF
            IF ( is_hubbard_back(nt) ) THEN
               WRITE( stdout,'(a,2i3)') 'VNSB(NA,IS) ', na, is
               DO m1 = 1, ldimb
                  WRITE( stdout,'(7f10.4)') (v%nsb(m1,m2,is,na),m2=1,ldimb)
               ENDDO
            ENDIF
         ENDDO
      ENDDO
   ENDIF
#endif
   !
   ! We start a loop on k points
   !
   DO ik = 1, nks
      !
      IF (lsda) current_spin = isk(ik)
      npw = ngk(ik)
      !
      IF (nks > 1) CALL get_buffer (evc, nwordwfc, iunwfc, ik)
      !
      CALL init_us_2 (npw, igk_k(1,ik), xk(1,ik), vkb)
      CALL calbec (npw, vkb, evc, becp)
      CALL s_psi  (npwx, npw, nbnd, evc, spsi)
      !
      ! Re-calculate atomic wfc - wfcatom is used here as work space
      !
      CALL atomic_wfc (ik, wfcatom)
      CALL copy_U_wfc (wfcatom)
      !
      ! wfcU contains Hubbard-U atomic wavefunctions
      ! proj=<wfcU|S|evc> - no need to read S*wfcU from buffer
      !
      CALL calbec ( npw, wfcU, spsi, proj)
      !
      ! NB: both ipol and jpol must run from 1 to 3 because this stress 
      !     contribution is not in general symmetric when computed only 
      !     from k-points in the irreducible wedge of the BZ. 
      !     It is (must be) symmetric after symmetrization but this requires 
      !     the full stress tensor not only its upper triangular part.
      !
      DO ipol = 1, 3
         DO jpol = 1, 3
            !
            IF (lda_plus_u_kind.EQ.0) THEN
               !
               ! The DFT+U case
               !
               ! Compute the derivative of the occupation matrix w.r.t epsilon
               !
               IF (gamma_only) THEN
                  CALL dndepsilon_gamma(ipol,jpol,ldim,proj%r,spsi,ik,nb_s,nb_e,mykey,1,dns)
               ELSE
                  CALL dndepsilon_k(ipol,jpol,ldim,proj%k,spsi,ik,nb_s,nb_e,mykey,1,dns)
               ENDIF  
               !
               DO na = 1, nat                 
                  nt = ityp(na)
                  IF ( is_hubbard(nt) ) THEN
                     DO is = 1, nspin
                        DO m2 = 1, 2 * Hubbard_l(nt) + 1
                           DO m1 = 1, 2 * Hubbard_l(nt) + 1
                              sigmah(ipol,jpol) = sigmah(ipol,jpol) - &
                                 v%ns(m2,m1,is,na) * dns(m1,m2,is,na) 
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDIF
               ENDDO
               !
               ! The background part
               !
               IF (lhubb) THEN
                  !
                  ! Compute the derivative of the occupation matrix w.r.t epsilon
                  !
                  IF (gamma_only) THEN
                     CALL dndepsilon_gamma(ipol,jpol,ldimb,proj%r,spsi,ik,nb_s,nb_e,mykey,2,dnsb)
                  ELSE
                     CALL dndepsilon_k(ipol,jpol,ldimb,proj%k,spsi,ik,nb_s,nb_e,mykey,2,dnsb)
                  ENDIF
                  !
                  DO na = 1, nat                 
                     nt = ityp(na)
                     IF ( is_hubbard_back(nt) ) THEN
                        DO is = 1, nspin
                           DO m2 = 1, ldim_back(nt)
                              DO m1 = 1, ldim_back(nt)
                                 sigmah(ipol,jpol) = sigmah(ipol,jpol) - &
                                    v%nsb(m2,m1,is,na) * dnsb(m1,m2,is,na)
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDIF
                  ENDDO
                  ! 
               ENDIF
               !
            ELSEIF (lda_plus_u_kind.EQ.2) THEN
               !
               ! The DFT+U+V case
               !
               ! Compute the derivative of the occupation matrix w.r.t epsilon
               !
               IF (gamma_only) THEN
                  CALL dngdepsilon_gamma(ipol,jpol,ldim,proj%r,spsi,ik,nb_s,nb_e,mykey,dnsg)
               ELSE
                  CALL dngdepsilon_k(ipol,jpol,ldim,proj%k,spsi,ik,nb_s,nb_e,mykey,dnsg)
               ENDIF
               !
               DO is = 1, nspin
                  DO na1 = 1, nat
                     nt1 = ityp(na1)
                     IF ( is_hubbard(nt1) ) THEN
                        ldim1 = ldim_u(nt1)
                        DO viz = 1, neighood(na1)%num_neigh
                           na2 = neighood(na1)%neigh(viz)
                           equiv_na2 = at_sc(na2)%at
                           nt2 = ityp (equiv_na2)
                           ldim2 = ldim_u(nt2)
                           IF (Hubbard_V(na1,na2,1).NE.0.d0 .OR. &
                               Hubbard_V(na1,na2,2).NE.0.d0 .OR. &
                               Hubbard_V(na1,na2,3).NE.0.d0 .OR. &
                               Hubbard_V(na1,na2,4).NE.0.d0 ) THEN
                               DO m1 = 1, ldim1
                                  DO m2 = 1, ldim2
                                     i_type = type_interaction(na1, m1, equiv_na2, m2)
                                     sigmah(ipol,jpol) = sigmah(ipol,jpol) - &
                                           DBLE(v_nsg(m2,m1,viz,na1,is) * dnsg(m2,m1,viz,na1,is)) 
                                  ENDDO 
                               ENDDO 
                           ENDIF
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO 
               ! 
            ENDIF
            !
         ENDDO ! jpol
      ENDDO ! ipol
      !
   ENDDO ! ik
   !
   CALL mp_sum( sigmah, inter_pool_comm )
   !
   sigmah(:,:) = sigmah(:,:) / omega
   !
   IF (nspin.EQ.1) sigmah(:,:) = 2.d0 * sigmah(:,:)
   !
   ! Symmetrization
   ! 
   CALL symmatrix ( sigmah )
   !
   ! Impose symmetry s(i,j) = s(j,i) to the stress tensor.
   ! It should NOT be needed, but let us do it for safety.
   !
   DO ipol = 1, 3
      DO jpol = ipol, 3
         IF ( ABS( sigmah(ipol,jpol)-sigmah(jpol,ipol) ) > 1.d-6 ) THEN
             WRITE(stdout,'(2i3,2f12.7)') ipol,jpol,sigmah(ipol,jpol), &
                                                     sigmah(jpol,ipol)
             CALL errore('stres_hub',' non-symmetric stress contribution',1)
         ENDIF
         sigmah(ipol,jpol) = 0.5d0 * ( sigmah(ipol,jpol) + sigmah(jpol,ipol) )
         sigmah(jpol,ipol) = sigmah(ipol,jpol)
      ENDDO
   ENDDO
   !
   CALL deallocate_bec_type ( becp )
   CALL deallocate_bec_type ( proj )
   IF (ALLOCATED(dns))  DEALLOCATE (dns)
   IF (ALLOCATED(dnsb)) DEALLOCATE (dnsb)
   IF (ALLOCATED(dnsg)) DEALLOCATE (dnsg)
   DEALLOCATE (wfcatom)
   DEALLOCATE (spsi)
   !
   use_bgrp_in_hpsi = save_flag
   !
   CALL stop_clock( 'stres_hub' )
   !
   RETURN
   !
END  SUBROUTINE stres_hub
!
!--------------------------------------------------------------------------------
SUBROUTINE dndepsilon_k ( ipol,jpol,ldim,proj,spsi,ik,nb_s,nb_e,mykey,lpuk,dns )
   !-----------------------------------------------------------------------------
   !! This routine computes the derivative of the ns atomic occupations with
   !! respect to the strain epsilon(ipol,jpol) used to obtain the Hubbard
   !! contribution to the internal stres tensor.
   !
   USE kinds,             ONLY : DP
   USE ions_base,         ONLY : nat, ityp
   USE klist,             ONLY : ngk
   USE lsda_mod,          ONLY : nspin, current_spin
   USE wvfct,             ONLY : nbnd, npwx, wg
   USE becmod,            ONLY : bec_type, allocate_bec_type, deallocate_bec_type
   USE mp_pools,          ONLY : intra_pool_comm
   USE mp,                ONLY : mp_sum
   USE ldaU,              ONLY : nwfcU, offsetU, Hubbard_l, is_hubbard,  &
                                 ldim_back, offsetU_back, offsetU_back1, &
                                 is_hubbard_back, Hubbard_l_back, backall

   IMPLICIT NONE
   !
   ! I/O variables 
   !
   INTEGER, INTENT(IN) :: ipol
   !! component index 1 to 3
   INTEGER, INTENT(IN) :: jpol
   !! component index 1 to 3
   INTEGER, INTENT(IN) :: ldim
   !! number of m-values: ldim = 2*Hubbard_lmax+1
   COMPLEX (DP), INTENT(IN) :: proj(nwfcU,nbnd)
   !! projection
   COMPLEX (DP), INTENT(IN) :: spsi(npwx,nbnd)
   !! \(S|\ \text{evc}\rangle\)
   INTEGER, INTENT(IN) :: ik
   !! k-point index
   INTEGER, INTENT(IN) :: nb_s
   !! starting band number (for band parallelization)
   INTEGER, INTENT(IN) :: nb_e
   !! ending band number (for band parallelization)
   INTEGER, INTENT(IN) :: mykey
   !! If each band appears more than once
   !! compute its contribution only once (i.e. when mykey=0) 
   INTEGER, INTENT(IN) :: lpuk
   !! index to control the standard (lpuk=1) or 
   !! background (lpuk=2) contribution to the force
   REAL(DP), INTENT(OUT) :: dns(ldim,ldim,nspin,nat)
   !! the derivative of the ns atomic occupations
   !
   ! ... local variables
   !
   INTEGER :: ibnd,  & ! count on bands
              is,    & ! count on spins
              npw,   & ! number of plane waves
              na,    & ! atomic index
              nt,    & ! index of the atomic type
              m1, m2, m11, m22, & ! indices of magnetic quantum numbers
              off1, off2, off22   ! indices for the offsets
   TYPE (bec_type) :: dproj
   !
   CALL allocate_bec_type ( nwfcU,nbnd, dproj )
   !
   ! D_Sl for l=1 and l=2 are already initialized, for l=0 D_S0 is 1
   !
   ! Offset of atomic wavefunctions initialized in setup and stored in offsetU
   ! 
   dns(:,:,:,:) = 0.d0
   !
   npw = ngk(ik)
   !
   ! Calculate the first derivative of proj with respect to epsilon(ipol,jpol)
   !
   CALL dprojdepsilon_k (spsi, ik, ipol, jpol, nb_s, nb_e, mykey, dproj%k)
   !
   ! Band parallelization. If each band appears more than once
   ! compute its contribution only once (i.e. when mykey=0)
   !
   IF ( mykey /= 0 ) GO TO 10
   !
   DO na = 1, nat
      !
      nt = ityp(na)
      !
      IF ( is_hubbard(nt) .AND. lpuk.EQ.1 ) THEN        
         DO m1 = 1, 2 * Hubbard_l(nt) + 1
            DO m2 = m1, 2 * Hubbard_l(nt) + 1
               DO ibnd = nb_s, nb_e
                  dns(m1,m2,current_spin,na) = &
                     dns(m1,m2,current_spin,na) + wg(ibnd,ik) * &
                     DBLE(proj(offsetU(na)+m1,ibnd) *           &
                     CONJG(dproj%k(offsetU(na)+m2,ibnd) ) +     &
                     dproj%k(offsetU(na)+m1,ibnd)*              &
                     CONJG(proj(offsetU(na)+m2,ibnd) ) )
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      !
      ! The background part
      !
      IF ( is_hubbard_back(nt) .AND. lpuk.EQ.2 ) THEN        
         DO m1 = 1, ldim_back(nt)
            off1 = offsetU_back(na)
            m11 = m1
            IF (backall(nt) .AND. m1.GT.2*Hubbard_l_back(nt)+1) THEN
               off1 = offsetU_back1(na)
               m11 = m1 - 2*Hubbard_l_back(nt) - 1
            ENDIF
            DO m2 = m1, ldim_back(nt)
               off2 = offsetU_back(na)
               m22 = m2
               IF (backall(nt) .AND. m2.GT.2*Hubbard_l_back(nt)+1) THEN
                  off2 = offsetU_back1(na)
                  m22 = m2 - 2*Hubbard_l_back(nt) - 1
               ENDIF
               DO ibnd = nb_s, nb_e
                  dns(m1,m2,current_spin,na) = &
                     dns(m1,m2,current_spin,na) + wg(ibnd,ik) * &
                     DBLE( proj(off1+m11,ibnd) *                &
                     CONJG(dproj%k(off2+m22,ibnd) ) +           &
                     dproj%k(off1+m11,ibnd)*                    &
                     CONJG(proj(off2+m22,ibnd) ) )
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      !
   ENDDO
   !
10 CALL mp_sum(dns, intra_pool_comm)
   !
   ! In nspin=1 k-point weight wg is normalized to 2 el/band 
   ! in the whole BZ but we are interested in dns of one spin component
   !
   IF (nspin.EQ.1) dns = 0.5d0 * dns
   !
   ! Impose hermiticity of dns_{m1,m2}
   !
   DO na = 1, nat
      nt = ityp(na)
      DO is = 1, nspin
         DO m1 = 1, ldim 
            DO m2 = m1+1, ldim 
               dns(m2,m1,is,na) = dns(m1,m2,is,na)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   !
   CALL deallocate_bec_type (dproj)
   !
   RETURN
   !
END SUBROUTINE dndepsilon_k
!
!------------------------------------------------------------------------------------
SUBROUTINE dndepsilon_gamma ( ipol,jpol,ldim,proj,spsi,ik,nb_s,nb_e,mykey,lpuk,dns )
   !---------------------------------------------------------------------------------
   !! This routine computes the derivative of the ns atomic occupations with
   !! respect to the strain epsilon(ipol,jpol) used to obtain the Hubbard
   !! contribution to the internal stres tensor.
   !
   USE kinds,             ONLY : DP
   USE ions_base,         ONLY : nat, ityp
   USE klist,             ONLY : ngk
   USE lsda_mod,          ONLY : nspin, current_spin
   USE wvfct,             ONLY : nbnd, npwx, wg
   USE becmod,            ONLY : bec_type, allocate_bec_type, deallocate_bec_type
   USE mp_pools,          ONLY : intra_pool_comm
   USE mp,                ONLY : mp_sum
   USE ldaU,              ONLY : nwfcU, offsetU, Hubbard_l, is_hubbard,  &
                                 ldim_back, offsetU_back, offsetU_back1, &
                                 is_hubbard_back, Hubbard_l_back, backall
 
   IMPLICIT NONE
   !
   ! I/O variables 
   !
   INTEGER, INTENT(IN) :: ipol
   !! component index 1 to 3
   INTEGER, INTENT(IN) :: jpol
   !! component index 1 to 3
   INTEGER, INTENT(IN) :: ldim
   !! number of m-values: ldim = 2*Hubbard_lmax+1
   REAL (DP), INTENT(IN) :: proj(nwfcU,nbnd)
   !! projection
   COMPLEX (DP), INTENT(IN) :: spsi(npwx,nbnd)
   !! \(S|\ \text{evc}\rangle\)
   INTEGER, INTENT(IN) :: ik
   !! k-point index
   INTEGER, INTENT(IN) :: nb_s
   !! starting band number (for band parallelization)
   INTEGER, INTENT(IN) :: nb_e
   !! ending band number (for band parallelization)
   INTEGER, INTENT(IN) :: mykey
   !! If each band appears more than once
   !! compute its contribution only once (i.e. when mykey=0) 
   INTEGER, INTENT(IN) :: lpuk
   !! index to control the standard (lpuk=1) or 
   !! background (lpuk=2) contribution to the force
   REAL(DP), INTENT(OUT) :: dns(ldim,ldim,nspin,nat)
   !! the derivative of the ns atomic occupations
   !
   ! ... local variables
   !
   INTEGER :: ibnd,  & ! count on bands
              is,    & ! count on spins
              npw,   & ! number of plane waves
              na,    & ! atomic index
              nt,    & ! index of the atomic type
              m1, m2, m11, m22, & ! indices of magnetic quantum numbers
              off1, off2, off22   ! indices for the offsets
   TYPE (bec_type) :: dproj
   ! 
   CALL allocate_bec_type ( nwfcU,nbnd, dproj )
   !
   ! D_Sl for l=1 and l=2 are already initialized, for l=0 D_S0 is 1
   !
   ! Offset of atomic wavefunctions initialized in setup and stored in offsetU
  
   dns(:,:,:,:) = 0.d0
   !
   npw = ngk(ik)
   !
   ! Calculate the first derivative of proj with respect to epsilon(ipol,jpol)
   !
   CALL dprojdepsilon_gamma (spsi, ik, ipol, jpol, nb_s, nb_e, mykey, dproj%r)
   !
   ! Band parallelization. If each band appears more than once
   ! compute its contribution only once (i.e. when mykey=0)
   !
   IF ( mykey /= 0 ) GO TO 10
   !
   DO na = 1, nat
      !
      nt = ityp(na)
      !
      IF ( is_hubbard(nt) .AND. lpuk.EQ.1 ) THEN        
         DO m1 = 1, 2 * Hubbard_l(nt) + 1
            DO m2 = m1, 2 * Hubbard_l(nt) + 1
               DO ibnd = nb_s, nb_e
                  dns(m1,m2,current_spin,na) = &
                     dns(m1,m2,current_spin,na) + wg(ibnd,ik) * &
                     ( proj(offsetU(na)+m1,ibnd) *              &
                       dproj%r(offsetU(na)+m2,ibnd) +           &
                       dproj%r(offsetU(na)+m1,ibnd) *           &
                       proj(offsetU(na)+m2,ibnd))
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      !
      ! The background part
      !
      IF ( is_hubbard_back(nt) .AND. lpuk.EQ.2 ) THEN        
         DO m1 = 1, ldim_back(nt)
            off1 = offsetU_back(na)
            m11 = m1
            IF (backall(nt) .AND. m1.GT.2*Hubbard_l_back(nt)+1) THEN
               off1 = offsetU_back1(na)
               m11 = m1 - 2*Hubbard_l_back(nt) - 1
            ENDIF
            DO m2 = m1, ldim_back(nt)
               off2 = offsetU_back(na)
               m22 = m2
               IF (backall(nt) .AND. m2.GT.2*Hubbard_l_back(nt)+1) THEN
                  off2 = offsetU_back1(na)
                  m22 = m2 - 2*Hubbard_l_back(nt)-1
               ENDIF
               DO ibnd = nb_s, nb_e
                  dns(m1,m2,current_spin,na) = &
                     dns(m1,m2,current_spin,na) + wg(ibnd,ik) * &
                     ( proj(off1+m11,ibnd) *                    &
                       dproj%r(off2+m22,ibnd) +                 &
                       dproj%r(off1+m11,ibnd) *                 &
                       proj(off2+m22,ibnd))
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      !
   ENDDO
   !
10 CALL mp_sum(dns, intra_pool_comm)
   !
   ! In nspin=1 k-point weight wg is normalized to 2 el/band 
   ! in the whole BZ but we are interested in dns of one spin component
   !
   IF (nspin.EQ.1) dns = 0.5d0 * dns
   !
   ! Impose hermiticity of dns_{m1,m2}
   !
   DO na = 1, nat
      nt = ityp(na)
      DO is = 1, nspin
         DO m1 = 1, ldim
            DO m2 = m1+1, ldim 
               dns(m2,m1,is,na) = dns(m1,m2,is,na)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   !
   CALL deallocate_bec_type (dproj)
   !
   RETURN
   !
END SUBROUTINE dndepsilon_gamma
!
!-----------------------------------------------------------------------
SUBROUTINE dngdepsilon_k ( ipol,jpol,ldim,proj,spsi,ik,nb_s,nb_e,mykey,dnsg )
   !-----------------------------------------------------------------------
   !! This routine computes the derivative of the nsg atomic occupations with
   !! respect to the strain epsilon(ipol,jpol) used to obtain the generalized 
   !! Hubbard contribution to the internal stres tensor.
   !
   USE kinds,             ONLY : DP
   USE ions_base,         ONLY : nat, ityp
   USE klist,             ONLY : ngk
   USE lsda_mod,          ONLY : nspin, current_spin
   USE wvfct,             ONLY : nbnd, npwx, wg
   USE becmod,            ONLY : bec_type, allocate_bec_type, deallocate_bec_type
   USE mp_pools,          ONLY : intra_pool_comm
   USE mp,                ONLY : mp_sum
   USE ldaU,              ONLY : nwfcU, Hubbard_l, is_hubbard, Hubbard_l_back, &
                                 ldim_u, at_sc, neighood, max_num_neighbors,   &
                                 phase_fac, backall, offsetU, offsetU_back,    &
                                 offsetU_back1
   
   IMPLICIT NONE
   !
   ! I/O variables 
   !
   INTEGER, INTENT(IN) :: ipol
   !! component index 1 to 3
   INTEGER, INTENT(IN) :: jpol
   !! component index 1 to 3
   INTEGER, INTENT(IN) :: ldim
   !! number of m-values: ldim = 2*Hubbard_lmax+1
   COMPLEX (DP), INTENT(IN) :: proj(nwfcU,nbnd)
   !! projection
   COMPLEX (DP), INTENT(IN) :: spsi(npwx,nbnd)
   !! \(S|\ \text{evc}\rangle\)
   INTEGER, INTENT(IN) :: ik
   !! k-point index
   INTEGER, INTENT(IN) :: nb_s
   !! starting band number (for band parallelization)
   INTEGER, INTENT(IN) :: nb_e
   !! ending band number (for band parallelization)
   INTEGER, INTENT(IN) :: mykey
   !! If each band appears more than once
   !! compute its contribution only once (i.e. when mykey=0) 
   COMPLEX (DP), INTENT (OUT) :: dnsg(ldim,ldim,max_num_neighbors,nat,nspin)
   !! the derivative of the nsg atomic occupations
   !
   ! ... local variables
   !
   INTEGER :: ibnd,              & ! count on bands
              is,                & ! count on spins
              npw,               & ! number of plane waves
              na1, na2,          & ! atomic indices
              nt1, nt2,          & ! indices of the atomic type
              m1, m2, m11, m22,  & ! indices of magnetic quantum numbers
              off1, off2, off22, & ! indices for the offsets
              ldim1, ldim2,      & ! dimensions of the Hubbard manifold
              eq_na2, viz          ! variables to control neighbors
   TYPE (bec_type) :: dproj
   INTEGER, EXTERNAL :: find_viz
   !
   CALL allocate_bec_type ( nwfcU, nbnd, dproj )
   !
   ! D_Sl for l=1 and l=2 are already initialized, for l=0 D_S0 is 1
   !
   ! Offset of atomic wavefunctions initialized in setup and stored in oatwfc
   !
   dnsg(:,:,:,:,:) = (0.0_dp, 0.0_dp)
   !
   npw = ngk(ik)
   !
   ! Calculate the first derivative of proj with respect to epsilon(ipol,jpol)
   !
   CALL dprojdepsilon_k (spsi, ik, ipol, jpol, nb_s, nb_e, mykey, dproj%k)
   !
   ! Set the phases for each atom at this ik
   !
   CALL phase_factor(ik)
   !
   ! Band parallelization. If each band appears more than once
   ! compute its contribution only once (i.e. when mykey=0)
   !
   IF ( mykey /= 0 ) GO TO 10
   !
   DO na1 = 1, nat
      nt1 = ityp(na1)
      IF ( is_hubbard(nt1) ) THEN
         ldim1 = ldim_u(nt1)
         DO viz = 1, neighood(na1)%num_neigh
            na2 = neighood(na1)%neigh(viz)
            eq_na2 = at_sc(na2)%at
            nt2 = ityp(eq_na2)
            ldim2 = ldim_u(nt2)
            IF (na1.GT.na2) THEN 
               DO m1 = 1, ldim1
                  DO m2 = 1, ldim2
                     dnsg(m2,m1,viz,na1,current_spin) = &
                     CONJG(dnsg(m1,m2,find_viz(na2,na1), na2, current_spin))
                  ENDDO
               ENDDO
            ELSE
               DO m1 = 1, ldim1
                  off1 = offsetU(na1) + m1
                  IF (m1.GT.2*Hubbard_l(nt1)+1) &
                     off1 = offsetU_back(na1) + m1 - 2*Hubbard_l(nt1) - 1
                  IF (backall(nt1) .AND. &
                     m1.GT.2*(Hubbard_l(nt1)+Hubbard_l_back(nt1)+1) ) &
                     off1 = offsetU_back1(na1) + m1 - &
                            2*(Hubbard_l(nt1)+Hubbard_l_back(nt1)+1)
                  DO m2 = 1, ldim2
                     off2 = offsetU(eq_na2) + m2
                     IF (m2.GT.2*Hubbard_l(nt2)+1) & 
                        off2 = offsetU_back(eq_na2)+ m2 - 2*Hubbard_l(nt2) - 1
                     IF (backall(nt2) .AND. &
                        m2.GT.2*(Hubbard_l(nt2)+Hubbard_l_back(nt2)+1) ) &
                        off2 = offsetU_back1(eq_na2) + m2 - &
                               2*(Hubbard_l(nt2)+Hubbard_l_back(nt2)+1)
                     DO ibnd = nb_s, nb_e
                        dnsg(m2,m1,viz,na1,current_spin) =                  &
                            dnsg(m2,m1,viz,na1,current_spin) +              &
                            wg(ibnd,ik) * CONJG(phase_fac(na2)) *           &
                            (proj(off1,ibnd) * CONJG(dproj%k(off2,ibnd))  + &
                             dproj%k(off1,ibnd) * CONJG(proj(off2,ibnd)) )
                     ENDDO
                  ENDDO 
               ENDDO  
            ENDIF 
         ENDDO           
      ENDIF 
   ENDDO 
   !
10 CALL mp_sum(dnsg, intra_pool_comm)
   !
   ! In nspin=1 k-point weight wg is normalized to 2 el/band 
   ! in the whole BZ but we are interested in dnsg of one spin component
   !
   IF (nspin.EQ.1) dnsg = 0.5d0 * dnsg
   !
   CALL deallocate_bec_type (dproj)
   !
   RETURN
   !
END SUBROUTINE dngdepsilon_k
!
!-----------------------------------------------------------------------
SUBROUTINE dngdepsilon_gamma ( ipol,jpol,ldim,proj,spsi,ik,nb_s,nb_e,mykey,dnsg )
   !-----------------------------------------------------------------------
   !! This routine computes the derivative of the nsg atomic occupations with
   !! respect to the strain epsilon(ipol,jpol) used to obtain the generalized 
   !! Hubbard contribution to the internal stres tensor.
   !
   USE kinds,             ONLY : DP
   USE ions_base,         ONLY : nat, ityp
   USE klist,             ONLY : ngk
   USE lsda_mod,          ONLY : nspin, current_spin
   USE wvfct,             ONLY : nbnd, npwx, wg
   USE becmod,            ONLY : bec_type, allocate_bec_type, deallocate_bec_type
   USE mp_pools,          ONLY : intra_pool_comm
   USE mp,                ONLY : mp_sum
   USE ldaU,              ONLY : nwfcU, Hubbard_l, is_hubbard, Hubbard_l_back, &
                                 ldim_u, at_sc, neighood, max_num_neighbors,   &
                                 phase_fac, backall, offsetU, offsetU_back,    &
                                 offsetU_back1

   IMPLICIT NONE
   !
   ! I/O variables 
   !
   INTEGER, INTENT(IN) :: ipol
   !! component index 1 to 3
   INTEGER, INTENT(IN) :: jpol
   !! component index 1 to 3
   INTEGER, INTENT(IN) :: ldim
   !! number of m-values: ldim = 2*Hubbard_lmax+1
   REAL (DP), INTENT(IN) :: proj(nwfcU,nbnd)
   !! projection
   COMPLEX (DP), INTENT(IN) :: spsi(npwx,nbnd)
   !! \(S|\ \text{evc}\rangle\)
   INTEGER, INTENT(IN) :: ik
   !! k-point index
   INTEGER, INTENT(IN) :: nb_s
   !! starting band number (for band parallelization)
   INTEGER, INTENT(IN) :: nb_e
   !! ending band number (for band parallelization)
   INTEGER, INTENT(IN) :: mykey
   !! If each band appears more than once
   !! compute its contribution only once (i.e. when mykey=0) 
   COMPLEX (DP), INTENT (OUT) :: dnsg(ldim,ldim,max_num_neighbors,nat,nspin)
   !! the derivative of the nsg atomic occupations
   !
   ! ... local variables
   !
   INTEGER :: ibnd,              & ! count on bands
              is,                & ! count on spins
              npw,               & ! number of plane waves
              na1, na2,          & ! atomic indices
              nt1, nt2,          & ! indices of the atomic type
              m1, m2, m11, m22,  & ! indices of magnetic quantum numbers
              off1, off2, off22, & ! indices for the offsets
              ldim1, ldim2,      & ! dimensions of the Hubbard manifold
              eq_na2, viz          ! variables to control neighbors
   TYPE (bec_type) :: dproj
   INTEGER, EXTERNAL :: find_viz
   !
   CALL allocate_bec_type ( nwfcU, nbnd, dproj )
   !
   ! D_Sl for l=1 and l=2 are already initialized, for l=0 D_S0 is 1
   !
   ! Offset of atomic wavefunctions initialized in setup and stored in oatwfc
   !
   dnsg(:,:,:,:,:) = (0.0_dp, 0.0_dp)
   !
   npw = ngk(ik)
   !
   ! Calculate the first derivative of proj with respect to epsilon(ipol,jpol)
   !
   CALL dprojdepsilon_gamma (spsi, ik, ipol, jpol, nb_s, nb_e, mykey, dproj%r)
   !
   ! Set the phases for each atom at this ik
   !
   CALL phase_factor(ik)
   !
   ! Band parallelization. If each band appears more than once
   ! compute its contribution only once (i.e. when mykey=0)
   !
   IF ( mykey /= 0 ) GO TO 10
   !
   DO na1 = 1, nat
      nt1 = ityp(na1)
      IF ( is_hubbard(nt1) ) THEN
         ldim1 = ldim_u(nt1)
         DO viz = 1, neighood(na1)%num_neigh
            na2 = neighood(na1)%neigh(viz)
            eq_na2 = at_sc(na2)%at
            nt2 = ityp(eq_na2)
            ldim2 = ldim_u(nt2)
            IF (na1.GT.na2) THEN 
               DO m1 = 1, ldim1
                  DO m2 = 1, ldim2
                     dnsg(m2,m1,viz,na1,current_spin) = &
                     CONJG(dnsg(m1,m2,find_viz(na2,na1), na2, current_spin))
                  ENDDO
               ENDDO
            ELSE
               DO m1 = 1, ldim1
                  off1 = offsetU(na1) + m1
                  IF (m1.GT.2*Hubbard_l(nt1)+1) &
                       off1 = offsetU_back(na1) + m1 - 2*Hubbard_l(nt1) - 1
                  IF (backall(nt1) .AND. &
                     m1.GT.2*(Hubbard_l(nt1)+Hubbard_l_back(nt1)+1) ) &
                     off1 = offsetU_back1(na1) + m1 - &
                            2*(Hubbard_l(nt1)+Hubbard_l_back(nt1)+1)
                  DO m2 = 1, ldim2
                     off2 = offsetU(eq_na2) + m2
                     IF (m2.GT.2*Hubbard_l(nt2)+1) & 
                        off2 = offsetU_back(eq_na2) + m2 - 2*Hubbard_l(nt2) - 1
                     IF (backall(nt2) .AND. &
                        m2.GT.2*(Hubbard_l(nt2)+Hubbard_l_back(nt2)+1) ) &
                        off2 = offsetU_back1(eq_na2) + m2 - &
                               2*(Hubbard_l(nt2)+Hubbard_l_back(nt2)+1)
                     DO ibnd = nb_s, nb_e
                        dnsg(m2,m1,viz,na1,current_spin) =              &
                            dnsg(m2,m1,viz,na1,current_spin) +          &
                            wg(ibnd,ik) * DBLE( CONJG(phase_fac(na2)) * &
                            (proj(off1,ibnd) * dproj%r(off2,ibnd) +     &
                             dproj%r(off1,ibnd) * proj(off2,ibnd) ) )
                     ENDDO
                  ENDDO 
               ENDDO  
            ENDIF 
         ENDDO           
      ENDIF 
   ENDDO 
   !
10 CALL mp_sum(dnsg, intra_pool_comm)
   !
   ! In nspin=1 k-point weight wg is normalized to 2 el/band 
   ! in the whole BZ but we are interested in dnsg of one spin component
   !
   IF (nspin.EQ.1) dnsg = 0.5d0 * dnsg
   !
   CALL deallocate_bec_type (dproj)
   !
   RETURN
   !
END SUBROUTINE dngdepsilon_gamma
!
!-----------------------------------------------------------------------
SUBROUTINE dprojdepsilon_k ( spsi, ik, ipol, jpol, nb_s, nb_e, mykey, dproj )
   !-----------------------------------------------------------------------
   !! This routine computes the first derivative of the projection
   !! \( \langle\phi^{at}_{I,m1}|S|\psi_{k,v,s}\rangle \)
   !! with respect to the strain \( \epsilon(i,j) \).
   !
   !! We remember that: $$ \text{ns}_{I,s,m_1,m_2} = \sum_{k,v}
   !! f_{kv} \langle\phi^{at}_{I,m1}|S|\psi_{k,v,s}\rangle
   !! \langle\psi_{k,v,s}|S|\phi^{at}_{I,m2}\rangle $$
   !
   USE kinds,                ONLY : DP
   USE cell_base,            ONLY : tpiba
   USE ions_base,            ONLY : nat, ntyp => nsp, ityp
   USE gvect,                ONLY : g
   USE klist,                ONLY : nks, xk, igk_k, ngk
   USE ldaU,                 ONLY : Hubbard_l, is_hubbard, nwfcU, wfcU
   USE lsda_mod,             ONLY : lsda, nspin, isk
   USE wvfct,                ONLY : nbnd, npwx, wg
   USE uspp,                 ONLY : nkb, vkb, qq_at
   USE uspp_param,           ONLY : upf, nhm, nh
   USE wavefunctions,        ONLY : evc
   USE becmod,               ONLY : becp, calbec

   IMPLICIT NONE
   !
   ! I/O variables 
   !
   COMPLEX(DP), INTENT(IN)  :: spsi(npwx,nbnd)
   !! S|evc>
   INTEGER, INTENT(IN) :: ik
   !! k-point index
   INTEGER, INTENT(IN) :: ipol
   !! component index 1 to 3
   INTEGER, INTENT(IN) :: jpol
   !! component index 1 to 3
   INTEGER, INTENT(IN) :: nb_s
   !! band number start
   INTEGER, INTENT(IN) :: nb_e
   !! band number end
   INTEGER, INTENT(IN) :: mykey
   !! labels how many times each band appears (mykey=0 first time etc.)  
   COMPLEX(DP), INTENT(OUT) :: dproj(nwfcU,nbnd)     
   !! the derivative of the projection
   !
   ! ... local variables
   !
   INTEGER :: i, ig, ijkb0, na, ibnd, iwf, nt, ih, jh, npw
   REAL (DP) :: xyz(3,3), q, a1, a2
   REAL (DP), PARAMETER :: eps=1.0d-8
   COMPLEX (DP), ALLOCATABLE :: &
           dwfc(:,:), aux(:,:), dbeta(:,:), aux0(:,:), aux1(:,:), &
           betapsi(:,:), dbetapsi(:,:), wfatbeta(:,:), wfatdbeta(:,:)
   !
   !       dwfc(npwx,nwfcU),   ! the derivative of the atomic d wfc
   !       aux(npwx,nwfcU),    ! auxiliary array
   !       dbeta(npwx,nkb),    ! the derivative of the beta function
   !       aux0,aux1(npwx,nkb),! auxiliary arrays
   !       betapsi(nhm,nbnd),  ! <beta|evc>
   !       dbetapsi(nhm,nbnd), ! <dbeta|evc>
   !       wfatbeta(nwfcU,nhm),! <wfc|beta>
   !       wfatdbeta(nwfcU,nhm)! <wfc|dbeta>
   !
   REAL (DP), ALLOCATABLE :: gk(:,:), qm1(:)
   !       gk(3,npwx),
   !       qm1(npwx)
   !
   ! xyz are the three unit vectors in the x,y,z directions
   xyz(:,:) = 0.d0
   DO i=1,3
      xyz(i,i) = 1.d0
   END DO
   !
   dproj(:,:) = (0.d0, 0.d0)
   !
   ! At first the derivatives of the atomic wfcs: we compute the term
   ! <d\fi^{at}_{I,m1}/d\epsilon(ipol,jpol)|S|\psi_{k,v,s}>
   !
   ALLOCATE ( qm1(npwx), gk(3,npwx) )
   ALLOCATE ( dwfc(npwx,nwfcU), aux(npwx,nwfcU) )
   !
   ! The derivative of the Bessel function
   !
   CALL gen_at_dj ( ik, nwfcU, is_hubbard, Hubbard_l, dwfc )
   !
   ! The derivative of the spherical harmonic
   !
   CALL gen_at_dy ( ik, nwfcU, is_hubbard, Hubbard_l, xyz(1,ipol), aux )
   !
   ! Number of plane waves at the k point with the index ik
   !
   npw = ngk(ik)
   !
   DO ig = 1, npw
      !
      gk(1,ig) = (xk(1,ik) + g(1,igk_k(ig,ik))) * tpiba
      gk(2,ig) = (xk(2,ik) + g(2,igk_k(ig,ik))) * tpiba
      gk(3,ig) = (xk(3,ik) + g(3,igk_k(ig,ik))) * tpiba
      !
      q = SQRT(gk(1,ig)**2 + gk(2,ig)**2 + gk(3,ig)**2)
      !
      IF (q.GT.eps) THEN
         qm1(ig) = 1.d0/q
      ELSE
         qm1(ig) = 0.d0
      ENDIF
      !
      a1 = -gk(jpol,ig)
      a2 = -gk(ipol,ig)*gk(jpol,ig)*qm1(ig)
      ! 
      DO iwf = 1, nwfcU
         dwfc(ig,iwf) = aux(ig,iwf)*a1 + dwfc(ig,iwf)*a2
      ENDDO
      !
   ENDDO
   !
   IF (ipol.EQ.jpol) dwfc(1:npw,:) = dwfc(1:npw,:) - wfcU(1:npw,:)*0.5d0
   !
   CALL calbec ( npw, dwfc, spsi, dproj )
   !
   DEALLOCATE ( dwfc, aux )
   !
   ! Now the derivatives of the beta functions: we compute the term
   ! <\fi^{at}_{I,m1}|dS/d\epsilon(ipol,jpol)|\psi_{k,v,s}>
   !
   ALLOCATE (aux0(npwx,nkb), aux1(npwx,nkb) )
   !
   ! The derivative of the Bessel function
   !
   CALL gen_us_dj (ik, aux0)
   !
   ! The derivative of the spherical harmonic
   !
   CALL gen_us_dy (ik, xyz(1,ipol), aux1)
   !
   ijkb0 = 0
   !
   DO nt = 1, ntyp
      !
      ALLOCATE (dbeta(npwx,nh(nt)), dbetapsi(nh(nt),nbnd), betapsi(nh(nt),nbnd), &
                wfatbeta(nwfcU,nh(nt)), wfatdbeta(nwfcU,nh(nt)) )
      !
      DO na = 1, nat
         !
         IF ( ityp(na).EQ.nt ) THEN
            !
            DO ih = 1, nh(nt)
               ! now we compute the true dbeta function
               DO ig = 1, npw
                  dbeta(ig,ih) = - aux1(ig,ijkb0+ih)*gk(jpol,ig) - &
                       aux0(ig,ijkb0+ih) * gk(ipol,ig) * gk(jpol,ig) * qm1(ig)
                  IF (ipol.EQ.jpol) &
                       dbeta(ig,ih) = dbeta(ig,ih) - vkb(ig,ijkb0+ih)*0.5d0
               ENDDO
            ENDDO
            !
            CALL calbec(npw, dbeta, evc,  dbetapsi )
            CALL calbec(npw, wfcU, dbeta, wfatdbeta )
            !
            ! dbeta is now used as work space to store vkb
            DO ih = 1, nh(nt)
               DO ig = 1, npw
                  dbeta(ig,ih) = vkb(ig,ijkb0+ih)
               ENDDO
            ENDDO
            !
            CALL calbec(npw, wfcU, dbeta, wfatbeta )
            !
            ! here starts band parallelization
            ! beta is used here as a work space to calculate dbetapsi
            !
            betapsi(:,:) = (0.0_dp, 0.0_dp)
            !
            DO ih = 1, nh(nt)
               DO ibnd = nb_s, nb_e
                  DO jh = 1, nh(nt)
                     betapsi(ih,ibnd) = betapsi(ih,ibnd) + &
                          qq_at(ih,jh,na)  * dbetapsi(jh,ibnd)
                  ENDDO
               ENDDO
            ENDDO
            !
            dbetapsi(:,:) = betapsi(:,:)
            !
            DO ih = 1, nh(nt)
               DO ibnd = nb_s, nb_e
                  betapsi(ih,ibnd) = (0.0_dp, 0.0_dp)
                  DO jh = 1, nh(nt)
                     betapsi(ih,ibnd) = betapsi(ih,ibnd) + &
                          qq_at(ih,jh,na) * becp%k(ijkb0+jh,ibnd)
                  ENDDO
               ENDDO
            ENDDO
            !
            ijkb0 = ijkb0 + nh(nt)
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
            ENDIF
            ! end band parallelization - only dproj(1,nb_s:nb_e) are calculated
         ENDIF
      ENDDO
      DEALLOCATE (dbeta, dbetapsi, betapsi, wfatbeta, wfatdbeta )
   ENDDO
   !
   DEALLOCATE ( aux0, aux1 )
   DEALLOCATE ( qm1, gk )
   !
   RETURN
   !
END SUBROUTINE dprojdepsilon_k
!
!-----------------------------------------------------------------------
SUBROUTINE dprojdepsilon_gamma ( spsi, ik, ipol, jpol, nb_s, nb_e, mykey, dproj )
   !-----------------------------------------------------------------------
   !! This routine computes the first derivative of the projection
   !! \( \langle\phi^{at}_{I,m1}|S|\psi_{k,v,s}\rangle \)
   !! with respect to the strain \( \epsilon(i,j) \). Gamma-only case.
   !
   !! We remember that: $$ \text{ns}_{I,s,m_1,m_2} = \sum_{k,v}
   !! f_{kv} \langle\phi^{at}_{I,m1}|S|\psi_{k,v,s}\rangle
   !! \langle\psi_{k,v,s}|S|\phi^{at}_{I,m2}\rangle $$
   !
   USE kinds,                ONLY : DP
   USE cell_base,            ONLY : tpiba
   USE ions_base,            ONLY : nat, ntyp => nsp, ityp
   USE gvect,                ONLY : g, gstart
   USE klist,                ONLY : nks, xk, igk_k, ngk
   USE ldaU,                 ONLY : Hubbard_l, is_hubbard, nwfcU, wfcU
   USE lsda_mod,             ONLY : lsda, nspin, isk
   USE wvfct,                ONLY : nbnd, npwx, wg
   USE uspp,                 ONLY : nkb, vkb, qq_at
   USE uspp_param,           ONLY : upf, nhm, nh
   USE wavefunctions,        ONLY : evc
   USE becmod,               ONLY : becp, calbec

   IMPLICIT NONE
   !
   ! I/O variables
   !
   COMPLEX(DP), INTENT(IN)  :: spsi(npwx,nbnd)
   !! S|evc>
   INTEGER, INTENT(IN) :: ik
   !! k-point index
   INTEGER, INTENT(IN) :: ipol
   !! component index 1 to 3
   INTEGER, INTENT(IN) :: jpol
   !! component index 1 to 3
   INTEGER, INTENT(IN) :: nb_s
   !! band number start
   INTEGER, INTENT(IN) :: nb_e
   !! band number end
   INTEGER, INTENT(IN) :: mykey
   !! labels how many times each band appears (mykey=0 first time etc.)
   REAL(DP), INTENT(OUT) :: dproj(nwfcU,nbnd)
   !! the derivative of the projection
   !
   ! ... local variables
   !
   INTEGER :: i, ig, ijkb0, na, ibnd, iwf, nt, ih, jh, npw
   REAL (DP) :: xyz(3,3), q, a1, a2
   REAL (DP), PARAMETER :: eps=1.0d-8
   COMPLEX (DP), ALLOCATABLE :: &
           dwfc(:,:), aux(:,:), dbeta(:,:), aux0(:,:), aux1(:,:)
   !       dwfc(npwx,nwfcU),   ! the derivative of the atomic d wfc
   !       aux(npwx,nwfcU),    ! auxiliary array
   !       dbeta(npwx,nkb),    ! the derivative of the beta function
   !       aux0,aux1(npwx,nkb) ! auxiliary arrays
   REAL (DP), ALLOCATABLE :: &
           betapsi(:,:), dbetapsi(:,:), wfatbeta(:,:), wfatdbeta(:,:)
   !       betapsi(nhm,nbnd),     ! <beta|evc>
   !       dbetapsi(nhm,nbnd),    ! <dbeta|evc>
   !       wfatbeta(nwfcU,nhm),! <wfc|beta>
   !       wfatdbeta(nwfcU,nhm)! <wfc|dbeta>

   REAL (DP), ALLOCATABLE :: gk(:,:), qm1(:)
   !       gk(3,npwx),
   !       qm1(npwx)
   !
   ! xyz are the three unit vectors in the x,y,z directions
   xyz(:,:) = 0.d0
   DO i=1,3
      xyz(i,i) = 1.d0
   END DO

   dproj(:,:) = 0.d0
   !
   ! At first the derivatives of the atomic wfcs: we compute the term
   ! <d\fi^{at}_{I,m1}/d\epsilon(ipol,jpol)|S|\psi_{k,v,s}>
   !
   ALLOCATE ( qm1(npwx), gk(3,npwx) )
   ALLOCATE ( dwfc(npwx,nwfcU), aux(npwx,nwfcU) )
   !
   ! The derivative of the Bessel function
   !
   CALL gen_at_dj ( ik, nwfcU, is_hubbard, Hubbard_l, dwfc )
   !
   ! The derivative of the spherical harmonic
   !
   CALL gen_at_dy ( ik, nwfcU, is_hubbard, Hubbard_l, xyz(1,ipol), aux)
   !
   ! Number of plane waves at the k point with the index ik
   !
   npw = ngk(ik)
   !
   DO ig = 1, npw
      !
      gk(1,ig) = (xk(1,ik) + g(1,igk_k(ig,ik))) * tpiba
      gk(2,ig) = (xk(2,ik) + g(2,igk_k(ig,ik))) * tpiba
      gk(3,ig) = (xk(3,ik) + g(3,igk_k(ig,ik))) * tpiba
      !
      q = SQRT(gk(1,ig)**2+gk(2,ig)**2+gk(3,ig)**2)
      !
      IF (q.GT.eps) THEN
         qm1(ig) = 1.d0/q
      ELSE
         qm1(ig) = 0.d0
      ENDIF
      !
      a1 = -gk(jpol,ig)
      a2 = -gk(ipol,ig)*gk(jpol,ig)*qm1(ig)
      !
      DO iwf = 1, nwfcU
         dwfc(ig,iwf) = aux(ig,iwf)*a1 + dwfc(ig,iwf)*a2
      ENDDO
      !
   ENDDO
   !
   IF (ipol.EQ.jpol) dwfc(1:npw,:) = dwfc(1:npw,:) - wfcU(1:npw,:)*0.5d0
   !
   CALL calbec ( npw, dwfc, spsi, dproj )
   !
   DEALLOCATE ( dwfc, aux )
   !
   ! Now the derivatives of the beta functions: we compute the term
   ! <\fi^{at}_{I,m1}|dS/d\epsilon(ipol,jpol)|\psi_{k,v,s}>
   !
   ALLOCATE (aux0(npwx,nkb), aux1(npwx,nkb) )
   !
   ! The derivative of the Bessel function
   !
   CALL gen_us_dj (ik, aux0)
   !
   ! The derivative of the spherical harmonic
   !
   CALL gen_us_dy (ik, xyz(1,ipol), aux1)
   !
   ijkb0 = 0
   !
   DO nt = 1, ntyp
      !
      ALLOCATE (dbeta(npwx,nh(nt)), dbetapsi(nh(nt),nbnd), betapsi(nh(nt),nbnd), &
                wfatbeta(nwfcU,nh(nt)), wfatdbeta(nwfcU,nh(nt)) )
      !
      DO na = 1, nat
         !
         IF ( ityp(na).EQ.nt ) THEN
            !
            DO ih = 1, nh(nt)
               ! now we compute the true dbeta function
               DO ig = 1, npw
                  dbeta(ig,ih) = - aux1(ig,ijkb0+ih)*gk(jpol,ig) - &
                       aux0(ig,ijkb0+ih) * gk(ipol,ig) * gk(jpol,ig) * qm1(ig)
                  IF (ipol.EQ.jpol) &
                       dbeta(ig,ih) = dbeta(ig,ih) - vkb(ig,ijkb0+ih)*0.5d0
               ENDDO
            ENDDO
            !
            CALL calbec(npw, dbeta, evc, dbetapsi )
            CALL calbec(npw, wfcU, dbeta, wfatdbeta )
            !
            ! dbeta is now used as work space to store vkb
            DO ih = 1, nh(nt)
               DO ig = 1, npw
                  dbeta(ig,ih) = vkb(ig,ijkb0+ih)
               ENDDO
            ENDDO
            !
            CALL calbec(npw, wfcU, dbeta, wfatbeta )
            !
            ! here starts band parallelization
            ! beta is here used as work space to calculate dbetapsi

            betapsi(:,:) = 0.0_dp
            DO ih = 1, nh(nt)
               DO ibnd = nb_s, nb_e
                  DO jh = 1, nh(nt)
                     betapsi(ih,ibnd) = betapsi(ih,ibnd) + &
                          qq_at(ih,jh,na)  * dbetapsi(jh,ibnd)
                  ENDDO
               ENDDO
            ENDDO
            !
            dbetapsi(:,:) = betapsi(:,:)
            !
            DO ih = 1, nh(nt)
               DO ibnd = nb_s, nb_e
                  betapsi(ih,ibnd) = 0.0_dp
                  DO jh = 1, nh(nt)
                     betapsi(ih,ibnd) = betapsi(ih,ibnd) + &
                          qq_at(ih,jh,na) * becp%r(ijkb0+jh,ibnd)
                  ENDDO
               ENDDO
            ENDDO
            !
            ijkb0 = ijkb0 + nh(nt)
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
            ENDIF
            ! end band parallelization - only dproj(1,nb_s:nb_e) are calculated
         ENDIF
      ENDDO
      DEALLOCATE (dbeta, dbetapsi, betapsi, wfatbeta, wfatdbeta )
   ENDDO
   !
   DEALLOCATE ( aux0, aux1 )
   DEALLOCATE ( qm1, gk )
   !
   RETURN
   !
END SUBROUTINE dprojdepsilon_gamma
