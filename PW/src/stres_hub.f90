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
   USE uspp,          ONLY : nkb, vkb, okvan
   USE klist,         ONLY : nks, xk, ngk, igk_k
   USE basis,         ONLY : natomwfc, wfcatom, swfcatom
   USE io_files,      ONLY : nwordwfc, iunwfc, nwordwfcU
   USE buffers,       ONLY : get_buffer
   USE scf,           ONLY : v, rho
   USE symme,         ONLY : symmatrix
   USE io_global,     ONLY : stdout
   USE mp_pools,      ONLY : inter_pool_comm, me_pool, nproc_pool
   USE mp,            ONLY : mp_sum
   USE control_flags, ONLY : gamma_only
   USE mp_bands,      ONLY : use_bgrp_in_hpsi
   USE noncollin_module, ONLY : noncolin
   USE force_mod,        ONLY : eigenval, eigenvect, overlap_inv, at_dy, at_dj, &
                                us_dy, us_dj
   !
   IMPLICIT NONE
   !
   REAL(DP), INTENT(OUT) :: sigmah(3,3) 
   !! the Hubbard contribution to stresses
   !
   ! ... local variables
   !
   INTEGER :: ipol, jpol, na, nt, is, m1, m2, na1, nt1, na2, nt2, viz, ik, npw
   INTEGER :: ldim, ldim1, ldim2, ldimb, equiv_na2, i_type, nb_s, nb_e, mykey, i
   REAL(DP), ALLOCATABLE :: dns(:,:,:,:), dnsb(:,:,:,:)
   REAL(DP) :: xyz(3,3)
   COMPLEX(DP), ALLOCATABLE ::  dnsg(:,:,:,:,:)
   !! the derivative of the atomic occupations
   COMPLEX(DP), ALLOCATABLE :: spsi(:,:)
   TYPE (bec_type) :: proj
   INTEGER, EXTERNAL :: type_interaction
   LOGICAL :: lhubb
   LOGICAL :: save_flag
   !
   CALL start_clock( 'stres_hub' )
   !
   save_flag = use_bgrp_in_hpsi ; use_bgrp_in_hpsi = .false.
   !
   IF (.NOT.((U_projection.EQ."atomic") .OR. (U_projection.EQ."ortho-atomic"))) &
      CALL errore("force_hub", &
                   " forces for this U_projection_type not implemented",1)
   !
   IF (lda_plus_u_kind.EQ.1) CALL errore("stres_hub", &
                   " stress in non collinear LDA+U scheme is not yet implemented",1)
   !
   IF (noncolin) CALL errore ("forceh","Noncollinear case is not supported",1)
   !
   sigmah(:,:) = 0.d0
   !
   ALLOCATE (spsi(npwx,nbnd))
   ALLOCATE (wfcatom(npwx,natomwfc))
   ALLOCATE (at_dy(npwx,natomwfc), at_dj(npwx,natomwfc))
   IF (okvan) ALLOCATE (us_dy(npwx,nkb), us_dj(npwx,nkb))
   IF (U_projection.EQ."ortho-atomic") THEN
      ALLOCATE (swfcatom(npwx,natomwfc))
      ALLOCATE (eigenval(natomwfc))
      ALLOCATE (eigenvect(natomwfc,natomwfc))
      ALLOCATE (overlap_inv(natomwfc,natomwfc))
   ENDIF
   !
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
      ! Compute spsi = S * psi
      CALL allocate_bec_type ( nkb, nbnd, becp)
      CALL calbec (npw, vkb, evc, becp)
      CALL s_psi  (npwx, npw, nbnd, evc, spsi)
      CALL deallocate_bec_type (becp)
      !
      ! Set up various quantities, in particular wfcU which 
      ! contains Hubbard-U (ortho-)atomic wavefunctions (without ultrasoft S)
      CALL orthoUwfc2 (ik)
      !
      ! proj=<wfcU|S|evc>
      CALL calbec ( npw, wfcU, spsi, proj)
      !
      ! Compute derivatives of spherical harmonics and spherical Bessel functions
      !
      ! xyz are the three unit vectors in the x,y,z directions
      xyz(:,:) = 0.d0
      DO i=1,3
         xyz(i,i) = 1.d0
      ENDDO 
      ! The derivative of spherical Bessel functions (for atomic functions)
      CALL gen_at_dj (ik, at_dj)
      ! The derivative of spherical Bessel functions (for beta functions)
      IF (okvan) CALL gen_us_dj (ik, us_dj)
      !
      ! NB: both ipol and jpol must run from 1 to 3 because this stress 
      !     contribution is not in general symmetric when computed only 
      !     from k-points in the irreducible wedge of the BZ. 
      !     It is (must be) symmetric after symmetrization but this requires 
      !     the full stress tensor not only its upper triangular part.
      !
      DO ipol = 1, 3
         !
         ! The derivative of spherical harmonics (for atomic functions)
         CALL gen_at_dy (ik, xyz(1,ipol), at_dy) 
         ! The derivative of spherical harmonics (for beta functions)
         IF (okvan) CALL gen_us_dy (ik, xyz(1,ipol), us_dy)
         !
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
   CALL deallocate_bec_type ( proj )
   IF (ALLOCATED(dns))  DEALLOCATE (dns)
   IF (ALLOCATED(dnsb)) DEALLOCATE (dnsb)
   IF (ALLOCATED(dnsg)) DEALLOCATE (dnsg)
   DEALLOCATE (spsi)
   DEALLOCATE (wfcatom)
   DEALLOCATE (at_dy, at_dj)
   IF (okvan) DEALLOCATE (us_dy, us_dj)
   IF (U_projection.EQ."ortho-atomic") THEN
      DEALLOCATE (swfcatom)
      DEALLOCATE (eigenval)
      DEALLOCATE (eigenvect)
      DEALLOCATE (overlap_inv)
   ENDIF
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
   USE ldaU,                 ONLY : nwfcU, wfcU, is_hubbard, is_hubbard_back,  &
                                    offsetU, offsetU_back, offsetU_back1,      &
                                    oatwfc, oatwfc_back, oatwfc_back1, ldim_u, &
                                    U_projection, Hubbard_l, Hubbard_l_back,   &
                                    backall
   USE lsda_mod,             ONLY : lsda, nspin, isk
   USE wvfct,                ONLY : nbnd, npwx, wg
   USE uspp,                 ONLY : nkb, vkb, okvan
   USE uspp_param,           ONLY : upf, nhm, nh
   USE wavefunctions,        ONLY : evc
   USE becmod,               ONLY : becp, calbec
   USE basis,                ONLY : natomwfc, wfcatom, swfcatom
   USE force_mod,            ONLY : eigenval, eigenvect, overlap_inv, at_dy, at_dj
   USE mp_bands,             ONLY : intra_bgrp_comm
   USE mp,                   ONLY : mp_sum

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
   INTEGER :: i, ig, ijkb0, na, ibnd, iwf, nt, ih, jh, npw, offpm, offpmU, &
              m1, m2, ldim_std
   REAL (DP) :: q, a1, a2
   REAL (DP), PARAMETER :: eps = 1.0d-8
   COMPLEX (DP), ALLOCATABLE :: &
   dproj0(:,:),       & ! derivative of the projector
   dproj_us(:,:),     & ! USPP contribution to dproj0
   dwfc(:,:),         & ! the derivative of the (ortho-atomic) wavefunction
   doverlap(:,:),     & ! derivative of the overlap matrix  
   doverlap_us(:,:),  & ! USPP contribution to doverlap
   doverlap_inv(:,:)    ! derivative of (O^{-1/2})_JI (note the transposition)   
   REAL (DP), ALLOCATABLE :: &
   gk(:,:), & ! k+G
   qm1(:)     ! 1/|k+G|
   !
   CALL start_clock('dprojdepsilon')
   ! 
   ! Number of plane waves at the k point with the index ik
   npw = ngk(ik)
   !
   dproj(:,:) = (0.d0, 0.d0)
   !
   ! At first the derivatives of the atomic wfcs: we compute the term
   ! <d\fi^{at}_{I,m1}/d\epsilon(ipol,jpol)|S|\psi_{k,v,s}>
   !
   ALLOCATE ( qm1(npwx), gk(3,npwx) )
   ALLOCATE ( dwfc(npwx,nwfcU) )
   dwfc(:,:) = (0.d0, 0.d0)
   !
   ! 1. Derivative of the atomic wavefunctions
   !    (and then multiplied by (O^{-1/2})_JI in the ortho-atomic case)
   !
   DO ig = 1, npw
      !
      ! Compute k+G and 1/|k+G|
      DO i = 1, 3
         gk(i,ig) = (xk(i,ik) + g(i,igk_k(ig,ik))) * tpiba
      ENDDO
      q = SQRT(gk(1,ig)**2 + gk(2,ig)**2 + gk(3,ig)**2)
      IF (q.GT.eps) THEN
         qm1(ig)=1.d0/q
      ELSE
         qm1(ig)=0.d0
      ENDIF
      !
      ! - (k+G)_jpol
      a1 = -gk(jpol,ig)
      !
      ! - (k+G)_ipol * (k+G)_jpol / |k+G|
      a2 = -gk(ipol,ig) * gk(jpol,ig) * qm1(ig)
      !
      DO na = 1, nat
         nt = ityp(na)
         ldim_std = 2*Hubbard_l(nt)+1
         IF (is_hubbard(nt) .OR. is_hubbard_back(nt)) THEN
            DO m1 = 1, ldim_u(nt)
               IF (m1.LE.ldim_std) THEN
                  offpmU = offsetU(na)
                  offpm  = oatwfc(na)
               ELSE
                  offpmU = offsetU_back(na) - ldim_std
                  offpm  = oatwfc_back(na)  - ldim_std
                  IF (backall(nt) .AND. m1.GT.ldim_std+2*Hubbard_l_back(nt)+1) THEN
                     offpmU = offsetU_back1(na) - ldim_std - 2*Hubbard_l_back(nt) - 1
                     offpm  = oatwfc_back1(na)  - ldim_std - 2*Hubbard_l_back(nt) - 1
                  ENDIF
               ENDIF
               IF (U_projection.EQ."atomic") THEN
                  dwfc(ig,offpmU+m1) = at_dy(ig,offpm+m1) * a1 + at_dj(ig,offpm+m1) * a2
               ELSEIF (U_projection.EQ."ortho-atomic") THEN
                  IF (m1>ldim_std) CALL errore("dprojdtau_k", &
                        " Stress with background and ortho-atomic is not supported",1)
                  DO m2 = 1, natomwfc
                     dwfc(ig,offpmU+m1) = dwfc(ig,offpmU+m1) + &
                         overlap_inv(offpm+m1,m2) * ( at_dy(ig,m2) * a1 + at_dj(ig,m2) * a2 )
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
      ENDDO
      !
   ENDDO
   !
   ! The diagonal term
   IF (ipol.EQ.jpol) dwfc(1:npw,:) = dwfc(1:npw,:) - wfcU(1:npw,:)*0.5d0
   !
   ! 2. Contribution due to the derivative of (O^{-1/2})_JI which
   !    is multiplied by atomic wavefunctions (only for ortho-atomic case)
   !
   IF (U_projection.EQ."ortho-atomic") THEN
      !
      ! Compute the derivative dO_IJ/d\epsilon(ipol,jpol)
      !
      ALLOCATE (doverlap(natomwfc,natomwfc))
      ALLOCATE (doverlap_inv(natomwfc,natomwfc))
      doverlap(:,:) = (0.0d0, 0.0d0)
      doverlap_inv(:,:) = (0.0d0, 0.0d0)
      !
      ! Calculate:
      ! doverlap = < dphi_I/d\epsilon(ipol,jpol) | S | phi_J > 
      !            + < phi_I | S | dphi_J/d\epsilon(ipol,jpol) >
      !
      DO ig = 1, npw
         !
         ! - (k+G)_jpol
         a1 = -gk(jpol,ig)
         !
         ! - (k+G)_ipol * (k+G)_jpol / |k+G|
         a2 = -gk(ipol,ig) * gk(jpol,ig) * qm1(ig)
         !
         DO m1 = 1, natomwfc
            DO m2 = 1, natomwfc
               doverlap(m1,m2) = doverlap(m1,m2) &
                       + CONJG((at_dy(ig,m1)*a1 + at_dj(ig,m1)*a2)) * swfcatom(ig,m2) &
                       + CONJG(swfcatom(ig,m1)) * (at_dy(ig,m2)*a1 + at_dj(ig,m2)*a2)
               IF (ipol.EQ.jpol) THEN
                  doverlap(m1,m2) = doverlap(m1,m2) &
                       + CONJG((-wfcatom(ig,m1)*0.5d0)) * swfcatom(ig,m2) &
                       + CONJG(swfcatom(ig,m1)) * (-wfcatom(ig,m2)*0.5d0)
               ENDIF
            ENDDO
         ENDDO
         !
      ENDDO
      ! Sum over G vectors
      CALL mp_sum( doverlap, intra_bgrp_comm )
      !
      ! USPP term in dO_IJ/d\epsilon(ipol,jpol)
      !
      IF (okvan) THEN
         ! Calculate doverlap_us = < phi_I | dS/d\epsilon(ipol,jpol) | phi_J >
         ALLOCATE (doverlap_us(natomwfc,natomwfc))
         CALL matrix_element_of_dSdepsilon (ik, ipol, jpol, &
              natomwfc, wfcatom, natomwfc, wfcatom, doverlap_us, 1, natomwfc, 0)
         ! Sum up the "normal" and ultrasoft terms
         DO m1 = 1, natomwfc
            DO m2 = 1, natomwfc
               doverlap(m1,m2) = doverlap(m1,m2) + doverlap_us(m1,m2)
            ENDDO
         ENDDO
         DEALLOCATE (doverlap_us)
      ENDIF
      !
      ! Now compute dO^{-1/2}_JI/d\epsilon(ipol,jpol) using dO_IJ/d\epsilon(ipol,jpol)
      ! Note the transposition!
      ! 
      CALL calculate_doverlap_inv (natomwfc, eigenval, eigenvect, &
                                     doverlap, doverlap_inv)
      !
      ! Now compute \sum_J dO^{-1/2}_JI/d\epsilon(ipol,jpol) \phi_J
      ! and add it to another term (see above)
      !
      DO na = 1, nat
         nt = ityp(na)
         IF (is_hubbard(nt) .OR. is_hubbard_back(nt)) THEN
            offpmU = offsetU(na)
            offpm  = oatwfc(na)
            DO m1 = 1, ldim_u(nt)
               DO m2 = 1, natomwfc
                  DO ig = 1, npw
                     dwfc(ig,offpmU+m1) = dwfc(ig,offpmU+m1) + &
                                   doverlap_inv(offpm+m1,m2) * wfcatom(ig,m2)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ENDDO
      !
      DEALLOCATE (doverlap)
      DEALLOCATE (doverlap_inv)
      !
   ENDIF
   !
   ! Compute dproj = <dwfc|S|psi> = <dwfc|spsi>
   CALL calbec ( npw, dwfc, spsi, dproj )
   !
   DEALLOCATE ( dwfc, qm1, gk)
   !
   ! Now the derivatives of the beta functions: we compute the term
   ! <\phi^{at}_{I,m1}|dS/d\epsilon(ipol,jpol)|\psi_{k,v,s}>
   !
   IF (okvan) THEN
      ALLOCATE(dproj_us(nwfcU,nb_s:nb_e))
      CALL matrix_element_of_dSdepsilon (ik, ipol, jpol, &
                         nwfcU, wfcU, nbnd, evc, dproj_us, nb_s, nb_e, mykey)
      ! dproj + dproj_us
      IF (mykey == 0) THEN
         DO m1 = 1, nwfcU
            dproj(m1,nb_s:nb_e) = dproj(m1,nb_s:nb_e) + dproj_us(m1,:)
         ENDDO
      ENDIF
      DEALLOCATE(dproj_us)
   ENDIF
   !
   CALL stop_clock('dprojdepsilon')
   !
   RETURN
   !
END SUBROUTINE dprojdepsilon_k
!
SUBROUTINE matrix_element_of_dSdepsilon (ik, ipol, jpol, lA, A, lB, B, A_dS_B, lB_s, lB_e, mykey)
   !
   ! This routine computes the matrix element < A | dS/d\epsilon(ipol,jpol) | B >
   ! Written by I. Timrov (2020)
   !
   ! Compute the term <\fi^{at}_{I,m1}|dS/d\epsilon(ipol,jpol)|\psi_{k,v,s}>
   !
   USE kinds,                ONLY : DP
   USE ions_base,            ONLY : nat, ntyp => nsp, ityp
   USE cell_base,            ONLY : tpiba
   USE gvect,                ONLY : g
   USE wvfct,                ONLY : npwx, wg
   USE uspp,                 ONLY : nkb, vkb, qq_at, okvan
   USE uspp_param,           ONLY : nh
   USE wavefunctions,        ONLY : evc
   USE becmod,               ONLY : calbec
   USE klist,                ONLY : xk, igk_k, ngk
   USE force_mod,            ONLY : us_dy, us_dj
   !
   IMPLICIT NONE
   !
   ! Input/Output
   !
   INTEGER, INTENT(IN)      :: ik          ! k point
   INTEGER, INTENT(IN)      :: ipol, jpol  ! Cartesian components
   INTEGER, INTENT(IN)      :: lA, lB, lB_s, lB_e
   ! There is a possibility to parallelize over lB,
   ! where lB_s (start) and lB_e (end)
   COMPLEX(DP), INTENT(IN)  :: A(npwx,lA)
   COMPLEX(DP), INTENT(IN)  :: B(npwx,lB)
   COMPLEX(DP), INTENT(OUT) :: A_dS_B(lA,lB_s:lB_e)
   INTEGER,     INTENT(IN)  :: mykey
   !
   ! Local variables
   !
   INTEGER :: npw, i, nt, na, ih, jh, ig, iA, iB, ijkb0
   REAL(DP) :: gvec
   COMPLEX (DP), ALLOCATABLE :: Adbeta(:,:), Abeta(:,:), &
                                dbetaB(:,:), betaB(:,:), &
                                aux(:,:), qq(:,:)
   REAL (DP) :: q, a1, a2
   REAL (DP), PARAMETER :: eps = 1.0d-8
   REAL (DP), ALLOCATABLE :: gk(:,:), qm1(:)
   !
   A_dS_B(:,:) = (0.0d0, 0.0d0)
   !
   IF (.NOT.okvan .OR. mykey /= 0) RETURN
   !
   npw = ngk(ik)
   !
   ALLOCATE ( qm1(npwx), gk(3,npwx) )
   !
   ! Compute k+G and 1/|k+G|
   DO ig = 1, npw
      DO i = 1, 3
         gk(i,ig) = (xk(i,ik) + g(i,igk_k(ig,ik))) * tpiba
      ENDDO
      q = SQRT(gk(1,ig)**2 + gk(2,ig)**2 + gk(3,ig)**2)
      IF (q.GT.eps) THEN
         qm1(ig)=1.d0/q
      ELSE
         qm1(ig)=0.d0
      ENDIF
   ENDDO
   !
   ijkb0 = 0
   !
   DO nt = 1, ntyp
      !
      ALLOCATE ( Adbeta(lA,nh(nt)) )
      ALLOCATE ( Abeta(lA,nh(nt)) )
      ALLOCATE ( dbetaB(nh(nt),lB) )
      ALLOCATE ( betaB(nh(nt),lB) )
      ALLOCATE ( qq(nh(nt),nh(nt)) )
      !
      DO na = 1, nat
         !
         IF ( ityp(na).EQ.nt ) THEN
            !
            qq(:,:) = CMPLX(qq_at(1:nh(nt),1:nh(nt),na), 0.0d0, kind=DP)
            !
            ! aux is used as a workspace
            ALLOCATE ( aux(npwx,nh(nt)) )
            !
            DO ih = 1, nh(nt)
               ! now we compute the true dbeta function
               DO ig = 1, npw
                  !
                  ! - (k+G)_jpol
                  a1 = -gk(jpol,ig)
                  !
                  ! - (k+G)_ipol * (k+G)_jpol / |k+G|
                  a2 = -gk(ipol,ig)*gk(jpol,ig)*qm1(ig)
                  !
                  aux(ig,ih) = us_dy(ig,ijkb0+ih) * a1 + us_dj(ig,ijkb0+ih) * a2
                  !
                  IF (ipol.EQ.jpol) aux(ig,ih) = aux(ig,ih) - vkb(ig,ijkb0+ih)*0.5d0
                  !
               ENDDO
            ENDDO
            !
            ! Calculate dbetaB = <dbeta|B> 
            CALL calbec(npw, aux, B, dbetaB )
            !
            ! Calculate Adbeta = <A|dbeta>
            CALL calbec(npw, A, aux, Adbeta )
            !
            ! aux is now used as a work space to store vkb
            DO ih = 1, nh(nt)
               DO ig = 1, npw
                  aux(ig,ih) = vkb(ig,ijkb0+ih)
               ENDDO
            ENDDO
            !
            ! Calculate Abeta = <A|beta>
            CALL calbec(npw, A, aux, Abeta )
            !
            ! Calculate betaB = <beta|B>
            CALL calbec( npw, aux, B, betaB )
            !
            DEALLOCATE ( aux )
            !
            ALLOCATE ( aux(nh(nt), lB) )
            ! 
            ! Calculate \sum_jh qq(ih,jh) * dbetaB(jh)
            CALL ZGEMM('N', 'N', nh(nt), lB_e-lB_s+1, nh(nt), (1.0d0,0.0d0), &
                       qq, nh(nt), dbetaB(1,lB_s),    nh(nt), (0.0d0,0.0d0), &
                       aux(1,lB_s), nh(nt))
            dbetaB(:,:) = aux(:,:)
            !
            ! Calculate \sum_jh qq(ih,jh) * betaB(jh)
            CALL ZGEMM('N', 'N', nh(nt), lB_e-lB_s+1, nh(nt), (1.0d0,0.0d0), &
                       qq, nh(nt), betaB(1,lB_s),     nh(nt), (0.0d0,0.0d0), &
                       aux(1,lB_s), nh(nt))
            betaB(:,:) = aux(:,:)
            !
            DEALLOCATE ( aux )
            !
            ijkb0 = ijkb0 + nh(nt)
            !
            ! dproj(iA,iB) = \sum_ih [Adbeta(iA,ih) * betapsi(ih,iB) +
            !                         Abeta(iA,ih)  * dbetaB(ih,iB)] 
            !
            CALL ZGEMM('N', 'N', lA, lB_e-lB_s+1, nh(nt), (1.0d0,0.0d0), &
                       Adbeta, lA, betaB(1,lB_s), nh(nt), (1.0d0,0.0d0), &
                       A_dS_B(1,lB_s), lA)
            CALL ZGEMM('N', 'N', lA, lB_e-lB_s+1, nh(nt), (1.0d0,0.0d0), &
                       Abeta, lA, dbetaB(1,lB_s), nh(nt), (1.0d0,0.0d0), &
                       A_dS_B(1,lB_s), lA)
            !
         ENDIF
         !
      ENDDO
      !
      DEALLOCATE (dbetaB, betaB, Abeta, Adbeta, qq)
      ! 
   ENDDO
   !
   DEALLOCATE ( qm1, gk )
   !
   RETURN
   !
END SUBROUTINE matrix_element_of_dSdepsilon
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
   USE ldaU,                 ONLY : nwfcU, wfcU, is_hubbard, is_hubbard_back,  &
                                    offsetU, offsetU_back, offsetU_back1,      &
                                    oatwfc, oatwfc_back, oatwfc_back1, ldim_u, &
                                    U_projection, Hubbard_l, Hubbard_l_back,   &
                                    backall
   USE lsda_mod,             ONLY : lsda, nspin, isk
   USE wvfct,                ONLY : nbnd, npwx, wg
   USE uspp,                 ONLY : nkb, vkb, qq_at, okvan
   USE uspp_param,           ONLY : upf, nhm, nh
   USE wavefunctions,        ONLY : evc
   USE becmod,               ONLY : becp, calbec
   USE basis,                ONLY : natomwfc
   USE force_mod,            ONLY : at_dy, at_dj, us_dy, us_dj
 
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
   INTEGER :: i, ig, ijkb0, na, ibnd, iwf, nt, ih, jh, npw, &
              offpm, offpmU, m1, ldim_std
   REAL (DP) :: q, a1, a2
   REAL (DP), PARAMETER :: eps=1.0d-8
   COMPLEX (DP), ALLOCATABLE :: &
           dwfc(:,:), dbeta(:,:)
   !       dwfc(npwx,nwfcU),   ! the derivative of the atomic d wfc
   !       dbeta(npwx,nkb),    ! the derivative of the beta function
   REAL (DP), ALLOCATABLE :: &
           betapsi(:,:), dbetapsi(:,:), wfatbeta(:,:), wfatdbeta(:,:), &
           betapsi0(:,:)
   !       betapsi(nhm,nbnd),     ! <beta|evc>
   !       dbetapsi(nhm,nbnd),    ! <dbeta|evc>
   !       wfatbeta(nwfcU,nhm),! <wfc|beta>
   !       wfatdbeta(nwfcU,nhm)! <wfc|dbeta>

   REAL (DP), ALLOCATABLE :: gk(:,:), qm1(:)
   !       gk(3,npwx),
   !       qm1(npwx)
   !
   ! See the implementation in dprojdepsilon_k
   IF (U_projection.EQ."ortho-atomic") CALL errore("dprojdtau_gamma", &
                    " Forces with gamma-only and ortho-atomic are not supported",1)
   !
   ! Number of plane waves at the k point with the index ik
   npw = ngk(ik)
   !
   dproj(:,:) = 0.d0
   !
   ! At first the derivatives of the atomic wfcs: we compute the term
   ! <d\fi^{at}_{I,m1}/d\epsilon(ipol,jpol)|S|\psi_{k,v,s}>
   !
   ALLOCATE ( qm1(npwx), gk(3,npwx) )
   ALLOCATE ( dwfc(npwx,nwfcU) )
   !
   DO ig = 1, npw
      !
      ! Compute k+G and 1/|k+G|
      DO i = 1, 3
         gk(i,ig) = (xk(i,ik) + g(i,igk_k(ig,ik))) * tpiba
      ENDDO
      q = SQRT(gk(1,ig)**2 + gk(2,ig)**2 + gk(3,ig)**2)
      IF (q.GT.eps) THEN
         qm1(ig) = 1.d0/q
      ELSE
         qm1(ig) = 0.d0
      ENDIF
      !
      ! - (k+G)_jpol
      a1 = -gk(jpol,ig)
      !
      ! - (k+G)_ipol * (k+G)_jpol / |k+G|
      a2 = -gk(ipol,ig)*gk(jpol,ig)*qm1(ig)
      !
      DO na = 1, nat
         nt = ityp(na)
         ldim_std = 2*Hubbard_l(nt)+1
         IF (is_hubbard(nt) .OR. is_hubbard_back(nt)) THEN
            DO m1 = 1, ldim_u(nt)
               IF (m1.LE.ldim_std) THEN
                  offpmU = offsetU(na)
                  offpm  = oatwfc(na)
               ELSE
                  offpmU = offsetU_back(na) - ldim_std
                  offpm  = oatwfc_back(na)  - ldim_std
                  IF (backall(nt) .AND. m1.GT.ldim_std+2*Hubbard_l_back(nt)+1) THEN
                     offpmU = offsetU_back1(na) - ldim_std - 2*Hubbard_l_back(nt) - 1
                     offpm  = oatwfc_back1(na)  - ldim_std - 2*Hubbard_l_back(nt) - 1
                  ENDIF
               ENDIF
               dwfc(ig,offpmU+m1) = at_dy(ig,offpm+m1) * a1 + at_dj(ig,offpm+m1) * a2
            ENDDO
         ENDIF
      ENDDO
      !
   ENDDO
   !
   IF (ipol.EQ.jpol) dwfc(1:npw,:) = dwfc(1:npw,:) - wfcU(1:npw,:)*0.5d0
   !
   CALL calbec ( npw, dwfc, spsi, dproj )
   !
   DEALLOCATE (dwfc)
   !
   ! Now the derivatives of the beta functions: we compute the term
   ! <\fi^{at}_{I,m1}|dS/d\epsilon(ipol,jpol)|\psi_{k,v,s}>
   !
   IF (okvan) THEN
    !
    ijkb0 = 0
    !
    DO nt = 1, ntyp
      !
      ALLOCATE (dbeta(npwx,nh(nt)), dbetapsi(nh(nt),nbnd), betapsi(nh(nt),nbnd), &
                wfatbeta(nwfcU,nh(nt)), wfatdbeta(nwfcU,nh(nt)), betapsi0(nh(nt),nbnd) )
      !
      DO na = 1, nat
         !
         IF ( ityp(na).EQ.nt ) THEN
            !
            DO ih = 1, nh(nt)
               ! now we compute the true dbeta function
               DO ig = 1, npw
                  dbeta(ig,ih) = - us_dy(ig,ijkb0+ih)*gk(jpol,ig) - &
                       us_dj(ig,ijkb0+ih) * gk(ipol,ig) * gk(jpol,ig) * qm1(ig)
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
            CALL calbec(npw, dbeta, evc, betapsi0 )
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
                          qq_at(ih,jh,na) * betapsi0(jh,ibnd)
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
      DEALLOCATE (dbeta, dbetapsi, betapsi, wfatbeta, wfatdbeta, betapsi0)
    ENDDO
    ! 
   ENDIF
   !
   DEALLOCATE ( qm1, gk )
   !
   RETURN
   !
END SUBROUTINE dprojdepsilon_gamma
