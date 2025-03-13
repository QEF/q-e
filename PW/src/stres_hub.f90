!
! Copyright (C) 2002-2025 Quantum ESPRESSO group
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
   USE kinds,              ONLY : DP
   USE wavefunctions,      ONLY : evc
   USE ions_base,          ONLY : nat, ityp, ntyp => nsp
   USE cell_base,          ONLY : omega, at, bg
   USE wvfct,              ONLY : nbnd, npwx
   USE ldaU,               ONLY : Hubbard_lmax, Hubbard_l, is_hubbard, &
                                  lda_plus_u_kind, Hubbard_projectors, is_hubbard_back, &
                                  ldim_back, ldmx_b, nsg, v_nsg, max_num_neighbors, &
                                  ldim_u, Hubbard_V, at_sc, neighood, ldmx_tot, &
                                  wfcU, nwfcU, Hubbard_J
   USE becmod,             ONLY : becp, calbec, allocate_bec_type_acc, deallocate_bec_type_acc
   USE lsda_mod,           ONLY : lsda, nspin, current_spin, isk
   USE uspp,               ONLY : nkb, vkb, okvan
   USE klist,              ONLY : nks, xk, ngk, igk_k
   USE basis,              ONLY : natomwfc, wfcatom, swfcatom
   USE io_files,           ONLY : nwordwfc, iunwfc, nwordwfcU
   USE buffers,            ONLY : get_buffer
   USE scf,                ONLY : v, rho
   USE symme,              ONLY : symmatrix
   USE io_global,          ONLY : stdout
   USE mp_pools,           ONLY : inter_pool_comm, me_pool, nproc_pool
   USE mp,                 ONLY : mp_sum
   USE control_flags,      ONLY : gamma_only, offload_type
   USE mp_bands,           ONLY : use_bgrp_in_hpsi, intra_bgrp_comm
   USE noncollin_module,   ONLY : noncolin, npol
   USE force_mod,          ONLY : eigenval, eigenvect, overlap_inv, at_dy, at_dj, &
                                  us_dy, us_dj
   USE uspp_init,          ONLY : init_us_2, gen_us_dj, gen_us_dy
   USE constants,          ONLY : eps16
   !
   IMPLICIT NONE
   !
   REAL(DP), INTENT(OUT) :: sigmah(3,3) 
   !! the Hubbard contribution to stresses
   !
   ! ... local variables
   !
   INTEGER :: ipol, jpol, na, nt, is, is2, m1, m2, na1, nt1, na2, nt2, viz, ik, npw
   INTEGER :: ldim, ldim1, ldim2, ldimb, equiv_na2, nb_s, nb_e, mykey, i
   REAL(DP), ALLOCATABLE :: dns(:,:,:,:), dnsb(:,:,:,:)
   COMPLEX (DP), ALLOCATABLE ::  dns_nc(:,:,:,:)
   REAL(DP) :: xyz(3,3)
   COMPLEX(DP), ALLOCATABLE ::  dnsg(:,:,:,:,:)
   !! the derivative of the atomic occupations
   COMPLEX(DP), ALLOCATABLE :: spsi(:,:)

   LOGICAL :: lhubb
   LOGICAL :: save_flag
   !
   REAL(DP), ALLOCATABLE :: projrd(:,:)
   COMPLEX(DP), ALLOCATABLE :: projkd(:,:)
   !
   CALL start_clock( 'stres_hub' )
   !
   save_flag = use_bgrp_in_hpsi ; use_bgrp_in_hpsi = .false.
   !
   IF (.NOT.((Hubbard_projectors.EQ."atomic") .OR. (Hubbard_projectors.EQ."ortho-atomic"))) &
      CALL errore("stres_hub", &
                   " stress for this Hubbard_projectors type not implemented",1)
   !
   ! IF (noncolin) CALL errore ("stres_hub","Noncollinear case is not supported",1)
   !
   IF (ANY(Hubbard_J(:,:)>eps16)) CALL errore("stres_hub", &
                   " stress in the DFT+U+J scheme is not implemented", 1 ) 
   !
   sigmah(:,:) = 0.d0
   !
   ALLOCATE (spsi(npwx*npol,nbnd))
   ALLOCATE (wfcatom(npwx*npol,natomwfc))
   ALLOCATE (at_dy(npwx*npol,natomwfc), at_dj(npwx*npol,natomwfc))
   IF (okvan) THEN
      ALLOCATE (us_dy(npwx,nkb), us_dj(npwx,nkb))
      !$acc enter data create (us_dy, us_dj)
   END IF
   IF (Hubbard_projectors.EQ."ortho-atomic") THEN
      ALLOCATE (swfcatom(npwx*npol,natomwfc))
      ALLOCATE (eigenval(natomwfc))
      ALLOCATE (eigenvect(natomwfc,natomwfc))
      ALLOCATE (overlap_inv(natomwfc,natomwfc))
      !$acc enter data create(swfcatom,eigenval,eigenvect,overlap_inv)
   ENDIF
   !
   !$acc data create(spsi,wfcatom) present(wfcU)
   !
   IF (gamma_only) THEN
      ALLOCATE( projrd(nwfcU,nbnd))
   ELSE
      ALLOCATE( projkd(nwfcU,nbnd))
   ENDIF
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
      IF (noncolin) then
         ALLOCATE ( dns_nc(ldim, ldim, nspin, nat) )
      ELSE
         ALLOCATE (dns(ldim,ldim,nspin,nat))
      ENDIF
      lhubb = .FALSE.
      DO nt = 1, ntyp
         IF (is_hubbard_back(nt)) lhubb = .TRUE.
      ENDDO
      IF (lhubb) THEN
         IF (noncolin) CALL errore ("stres_hub","Noncollinear and background are not supported",1)         
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
      IF (nks > 1) THEN
        CALL get_buffer (evc, nwordwfc, iunwfc, ik)
        !$acc update device(evc)
      END IF
      !
      CALL init_us_2 (npw, igk_k(1,ik), xk(1,ik), vkb, .TRUE.)
      !
      !$acc update self(vkb)
      !
      ! Compute spsi = S * psi
      CALL allocate_bec_type_acc ( nkb, nbnd, becp)
      !
      CALL calbec( offload_type, npw, vkb, evc, becp )
      CALL s_psi_acc( npwx, npw, nbnd, evc, spsi )
      !
      CALL deallocate_bec_type_acc (becp)
      !
      ! Set up various quantities, in particular wfcU which 
      ! contains Hubbard-U (ortho-)atomic wavefunctions (without ultrasoft S)
      CALL orthoUwfc_k (ik, .TRUE.)
      !
      ! proj=<wfcU|S|evc>
      IF (noncolin) THEN
         !$acc data create(projkd)
         !$acc host_data use_device(spsi,wfcU,projkd)
         CALL MYZGEMM ('C', 'N', nwfcU, nbnd, npwx*npol, (1.0_DP, 0.0_DP), wfcU, &
                    npwx*npol, spsi, npwx*npol, (0.0_DP, 0.0_DP),  projkd, nwfcU)
         CALL mp_sum( projkd( :, 1:nbnd ), intra_bgrp_comm )
         !$acc end host_data
         !$acc update self(projkd)
         !$acc end data
      ELSE
         IF (gamma_only) THEN
            !$acc data create(projrd)
            CALL calbec( offload_type, npw, wfcU, spsi, projrd )
            !$acc update self(projrd)
            !$acc end data
         ELSE
            !$acc data create(projkd)
            CALL calbec( offload_type, npw, wfcU, spsi, projkd )
            !$acc update self(projkd)
            !$acc end data
         ENDIF
      ENDIF
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
               IF (noncolin) THEN          
                  CALL dndepsilon_k_nc (ipol,jpol,ldim,projkd,spsi,ik,nb_s,nb_e,mykey,1,dns_nc )
                  DO na = 1, nat
                     nt = ityp(na)
                     IF ( is_hubbard(nt) ) THEN
                        DO is = 1, npol
                           DO is2 = 1, npol
                              DO m2 = 1, 2*Hubbard_l(nt)+1
                                 DO m1 = 1, 2*Hubbard_l(nt)+1
                                    sigmah(ipol,jpol) = sigmah(ipol,jpol)  -   &
                                           dble( v%ns_nc(m2, m1, npol*(is-1)+is2, na)  *   &
                                           dns_nc(m1, m2, npol*(is2-1)+is, na) )
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDIF
                  ENDDO                    
               ELSE        
                  IF (gamma_only) THEN
                     CALL dndepsilon_gamma(ipol,jpol,ldim,projrd,spsi,ik,nb_s,nb_e,mykey,1,dns)
                  ELSE
                     CALL dndepsilon_k(ipol,jpol,ldim,projkd,spsi,ik,nb_s,nb_e,mykey,1,dns)
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
               ENDIF
               !
               ! The background part
               !
               IF (lhubb) THEN
                  !
                  ! Compute the derivative of the occupation matrix w.r.t epsilon
                  !
                  IF (gamma_only) THEN
                     CALL dndepsilon_gamma(ipol,jpol,ldimb,projrd,spsi,ik,nb_s,nb_e,mykey,2,dnsb)
                  ELSE
                     CALL dndepsilon_k(ipol,jpol,ldimb,projkd,spsi,ik,nb_s,nb_e,mykey,2,dnsb)
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
               IF (noncolin) THEN
                  CALL dngdepsilon_k_nc(ipol,jpol,ldim,projkd,&
                                          spsi,ik,nb_s,nb_e,mykey,dnsg)
                  DO is = 1, npol
                     DO is2 = 1, npol
                        DO na1 = 1, nat
                           nt1 = ityp(na1)
                           IF ( is_hubbard(nt1) ) THEN
                              ldim1 = ldim_u(nt1)
                              DO viz = 1, neighood(na1)%num_neigh
                                 na2 = neighood(na1)%neigh(viz)
                                 equiv_na2 = at_sc(na2)%at
                                 nt2 = ityp (equiv_na2)
                                 ldim2 = ldim_u(nt2)
                                 IF (Hubbard_V(na1,na2,1).NE.0.d0 ) THEN
                                    DO m1 = 1, ldim1
                                       DO m2 = 1, ldim2
                                          sigmah(ipol,jpol) = sigmah(ipol,jpol) - &
                                                DBLE(v_nsg(m2,m1,viz,na1,npol*(is-1)+is2) *&
                                                      dnsg(m2,m1,viz,na1,npol*(is-1)+is2)) 
                                       ENDDO 
                                    ENDDO 
                                 ENDIF
                              ENDDO
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDDO
               ELSE
                  IF (gamma_only) THEN
                     CALL dngdepsilon_gamma(ipol,jpol,ldim,projrd,spsi,ik,nb_s,nb_e,mykey,dnsg)
                  ELSE
                     CALL dngdepsilon_k(ipol,jpol,ldim,projkd,spsi,ik,nb_s,nb_e,mykey,dnsg)
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
                                       sigmah(ipol,jpol) = sigmah(ipol,jpol) - &
                                             DBLE(v_nsg(m2,m1,viz,na1,is) * dnsg(m2,m1,viz,na1,is)) 
                                    ENDDO 
                                 ENDDO 
                              ENDIF
                           ENDDO
                        ENDIF
                     ENDDO
                  ENDDO 
               ENDIF
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
   IF (gamma_only) THEN
      DEALLOCATE( projrd)
   ELSE
      DEALLOCATE( projkd)
   ENDIF
   !
   IF (ALLOCATED(dns))  DEALLOCATE (dns)
   IF (ALLOCATED(dns_nc)) DEALLOCATE(dns_nc)
   IF (ALLOCATED(dnsb)) DEALLOCATE (dnsb)
   IF (ALLOCATED(dnsg)) DEALLOCATE (dnsg)
   !
   !$acc end data
   DEALLOCATE (spsi)
   DEALLOCATE (wfcatom)
   DEALLOCATE (at_dy, at_dj)
   IF (okvan) THEN
      !$acc exit data delete (us_dy, us_dj)
      DEALLOCATE (us_dy, us_dj)
   END IF
   IF (Hubbard_projectors.EQ."ortho-atomic") THEN
      !$acc exit data delete (swfcatom,eigenval,eigenvect,overlap_inv)
      DEALLOCATE (overlap_inv)
      DEALLOCATE (swfcatom)
      DEALLOCATE (eigenval)
      DEALLOCATE (eigenvect)
   ENDIF
   !
   use_bgrp_in_hpsi = save_flag
   !
   CALL stop_clock( 'stres_hub' )
   CALL print_clock('stres_hub')
   CALL print_clock('dprojdepsilon')
   CALL print_clock('dprojdeps1')
   CALL print_clock('dprojdeps2')
   CALL print_clock('dprojdeps3')
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
                                 is_hubbard_back, Hubbard_l2, backall
   !
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
            IF (backall(nt) .AND. m1.GT.2*Hubbard_l2(nt)+1) THEN
               off1 = offsetU_back1(na)
               m11 = m1 - 2*Hubbard_l2(nt) - 1
            ENDIF
            DO m2 = m1, ldim_back(nt)
               off2 = offsetU_back(na)
               m22 = m2
               IF (backall(nt) .AND. m2.GT.2*Hubbard_l2(nt)+1) THEN
                  off2 = offsetU_back1(na)
                  m22 = m2 - 2*Hubbard_l2(nt) - 1
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
SUBROUTINE dndepsilon_k_nc ( ipol,jpol,ldim,proj,spsi,ik,nb_s,nb_e,mykey,lpuk,dns_nc )
   !-----------------------------------------------------------------------------
   !! This routine computes the derivative of the ns_nc atomic occupations with
   !! respect to the strain epsilon(ipol,jpol) used to obtain the Hubbard
   !! contribution to the internal stres tensor in the noncollinear formalism.
   !
   USE kinds,             ONLY : DP
   USE ions_base,         ONLY : nat, ityp
   USE klist,             ONLY : ngk
   USE lsda_mod,          ONLY : nspin, current_spin
   USE wvfct,             ONLY : nbnd, npwx, wg
   USE noncollin_module,  ONLY : npol, noncolin
   USE mp_pools,          ONLY : intra_pool_comm
   USE mp,                ONLY : mp_sum
   USE ldaU,              ONLY : nwfcU, offsetU, Hubbard_l, is_hubbard,  &
                                 ldim_back, offsetU_back, offsetU_back1, &
                                 is_hubbard_back, Hubbard_l2, backall

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
   COMPLEX (DP), INTENT(IN) :: spsi(npwx*npol,nbnd)
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
   COMPLEX(DP), INTENT(OUT) :: dns_nc(ldim,ldim,nspin,nat)
   !! the derivative of the ns atomic occupations
   !
   ! ... local variables
   !
   INTEGER :: ibnd,  & ! count on bands
              is,    & ! count on spins
              npw,   & ! number of plane waves
              na,    & ! atomic index
              nt,    & ! index of the atomic type
              m1, m2, m11, m22,  & ! indices of magnetic quantum numbers
              off1, off2, off22, & ! indices for the offsets
              is1, is2, i, j, ldim1
   COMPLEX(DP), ALLOCATABLE :: dproj(:,:)
   REAL(DP) :: psum
   !
   ALLOCATE ( dproj(nwfcU,nbnd) )
   !
   ! D_Sl for l=1 and l=2 are already initialized, for l=0 D_S0 is 1
   !
   ! Offset of atomic wavefunctions initialized in setup and stored in offsetU
   !
   dns_nc(:,:,:,:) = 0.d0
   !
   npw = ngk(ik)
   !
   ! Calculate the first derivative of proj with respect to epsilon(ipol,jpol)
   !
   CALL dprojdepsilon_k (spsi, ik, ipol, jpol, nb_s, nb_e, mykey, dproj)
   !
   ! Band parallelization. If each band appears more than once
   ! compute its contribution only once (i.e. when mykey=0)
   !
   IF ( mykey /= 0 ) GO TO 10
   !
   DO na = 1, nat
      nt = ityp(na)
      IF ( is_hubbard(nt) ) THEN
         !
         ldim1 = 2*Hubbard_l(nt)+1
         DO is1 = 1, npol
            DO is2 = 1, npol
               i = npol*(is1-1) + is2
               DO m1 = 1, ldim1
                  DO m2 = 1, ldim1
                     DO ibnd = nb_s, nb_e
                        dns_nc(m1,m2,i,na) = dns_nc(m1,m2,i,na) + &
                            wg(ibnd,ik) * dcmplx(CONJG(dproj(offsetU(na)+m2+ldim1*(is2-1),ibnd) )* &
                                                       proj(offsetU(na)+m1+ldim1*(is1-1),ibnd)  + &
                                                 CONJG(proj(offsetU(na)+m2+ldim1*(is2-1),ibnd) )* &
                                                       dproj(offsetU(na)+m1+ldim1*(is1-1),ibnd))
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
       ENDIF
   ENDDO
   !
10 CALL mp_sum(dns_nc, intra_pool_comm)
   !
   !
   ! Impose hermiticity of dns_{m1,m2}
   !
   DO na = 1, nat
      nt = ityp (na)
      IF ( is_hubbard(nt) ) THEN
         DO is1 = 1, npol
           DO is2 = 1, npol
             i = npol*(is1-1) + is2
             j = is1 + npol*(is2-1)
             ldim1 = 2*Hubbard_l(nt)+1
             DO m1 = 1, ldim1
               DO m2 = 1, ldim1
                  psum = ABS( dns_nc(m1,m2,i,na) - CONJG(dns_nc(m2,m1,j,na)) )
                  IF (psum.GT.1.d-10) THEN
                      ! print*, na, m1, m2, is1, is2
                      ! print*, dns_nc(m1,m2,i,na)
                      ! print*, dns_nc(m2,m1,j,na)
                     CALL errore( 'dns_nc', 'non hermitean matrix', 1 )
                  ELSE
                     dns_nc(m2,m1,j,na) = CONJG( dns_nc(m1,m2,i,na) )
                  ENDIF
               ENDDO
             ENDDO
           ENDDO
         ENDDO
      ENDIF
   ENDDO
   !
   DEALLOCATE( dproj )
   !
   RETURN
   !
END SUBROUTINE dndepsilon_k_nc
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
                                 is_hubbard_back, Hubbard_l2, backall
 
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
            IF (backall(nt) .AND. m1.GT.2*Hubbard_l2(nt)+1) THEN
               off1 = offsetU_back1(na)
               m11 = m1 - 2*Hubbard_l2(nt) - 1
            ENDIF
            DO m2 = m1, ldim_back(nt)
               off2 = offsetU_back(na)
               m22 = m2
               IF (backall(nt) .AND. m2.GT.2*Hubbard_l2(nt)+1) THEN
                  off2 = offsetU_back1(na)
                  m22 = m2 - 2*Hubbard_l2(nt)-1
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
   USE ldaU,              ONLY : nwfcU, Hubbard_l, is_hubbard, Hubbard_l2, &
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
                     m1.GT.2*(Hubbard_l(nt1)+Hubbard_l2(nt1)+1) ) &
                     off1 = offsetU_back1(na1) + m1 - &
                            2*(Hubbard_l(nt1)+Hubbard_l2(nt1)+1)
                  DO m2 = 1, ldim2
                     off2 = offsetU(eq_na2) + m2
                     IF (m2.GT.2*Hubbard_l(nt2)+1) & 
                        off2 = offsetU_back(eq_na2)+ m2 - 2*Hubbard_l(nt2) - 1
                     IF (backall(nt2) .AND. &
                        m2.GT.2*(Hubbard_l(nt2)+Hubbard_l2(nt2)+1) ) &
                        off2 = offsetU_back1(eq_na2) + m2 - &
                               2*(Hubbard_l(nt2)+Hubbard_l2(nt2)+1)
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
SUBROUTINE dngdepsilon_k_nc ( ipol,jpol,ldim,proj,spsi,ik,nb_s,nb_e,mykey,dnsg )
   !-----------------------------------------------------------------------
   !! This routine computes the derivative of the nsg atomic occupations with
   !! respect to the strain epsilon(ipol,jpol) used to obtain the noncollinear
   !! generalized Hubbard contribution to the internal stres tensor.
   !
   USE kinds,             ONLY : DP
   USE ions_base,         ONLY : nat, ityp
   USE klist,             ONLY : ngk
   USE lsda_mod,          ONLY : nspin, current_spin
   USE wvfct,             ONLY : nbnd, npwx, wg
   USE mp_pools,          ONLY : intra_pool_comm
   USE mp,                ONLY : mp_sum
   USE ldaU,              ONLY : nwfcU, Hubbard_l, is_hubbard, Hubbard_l2, &
                                 ldim_u, at_sc, neighood, max_num_neighbors,   &
                                 phase_fac, backall, offsetU, offsetU_back,    &
                                 offsetU_back1
   USE noncollin_module,  ONLY : npol 
   USE io_global,          ONLY : stdout
   
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
   COMPLEX (DP), INTENT(IN) :: spsi(npwx*npol,nbnd)
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
              eq_na2, viz,       & ! variables to control neighbors
              is1, is2, i, j       ! for the spin indexes
   COMPLEX(DP), ALLOCATABLE :: dproj(:,:)
   INTEGER, EXTERNAL :: find_viz
   !
   ALLOCATE ( dproj(nwfcU,nbnd) )
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
   CALL dprojdepsilon_k (spsi, ik, ipol, jpol, nb_s, nb_e, mykey, dproj)
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
                     DO is1 = 1, npol
                        DO is2 = 1, npol
                           i = npol*(is2-1) + is1
                           j = npol*(is1-1) + is2
                           dnsg(m2,m1,viz,na1,i) = &
                             CONJG(dnsg(m1,m2,find_viz(na2,na1), na2, j))
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ELSE
               DO is1 = 1, npol
                  DO is2 = 1, npol
                     i = npol*(is2-1) + is1
                     DO m1 = 1, ldim1
                        off1 = offsetU(na1) + m1 +ldim1*(is2-1)
                        DO m2 = 1, ldim2
                           off2 = offsetU(eq_na2) + m2 +ldim2*(is1-1)
                           DO ibnd = nb_s, nb_e
                              dnsg(m2,m1,viz,na1,i) = &
                                 dnsg(m2,m1,viz,na1,i) +  &
                                 wg(ibnd,ik) * CONJG(phase_fac(na2)) * &
                                   dcmplx(proj(off1,ibnd) * CONJG(dproj(off2,ibnd))  + &
                                         dproj(off1,ibnd) *  CONJG(proj(off2,ibnd)) )
                           ENDDO
                        ENDDO 
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
   !
   DEALLOCATE( dproj )
   !
   RETURN
   !
END SUBROUTINE dngdepsilon_k_nc
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
   USE ldaU,              ONLY : nwfcU, Hubbard_l, is_hubbard, Hubbard_l2, &
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
                     m1.GT.2*(Hubbard_l(nt1)+Hubbard_l2(nt1)+1) ) &
                     off1 = offsetU_back1(na1) + m1 - &
                            2*(Hubbard_l(nt1)+Hubbard_l2(nt1)+1)
                  DO m2 = 1, ldim2
                     off2 = offsetU(eq_na2) + m2
                     IF (m2.GT.2*Hubbard_l(nt2)+1) & 
                        off2 = offsetU_back(eq_na2) + m2 - 2*Hubbard_l(nt2) - 1
                     IF (backall(nt2) .AND. &
                        m2.GT.2*(Hubbard_l(nt2)+Hubbard_l2(nt2)+1) ) &
                        off2 = offsetU_back1(eq_na2) + m2 - &
                               2*(Hubbard_l(nt2)+Hubbard_l2(nt2)+1)
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
                                    Hubbard_projectors, Hubbard_l, Hubbard_l2, &
                                    backall
   USE lsda_mod,             ONLY : lsda, nspin, isk
   USE wvfct,                ONLY : nbnd, npwx, wg
   USE uspp,                 ONLY : nkb, vkb, okvan
   USE wavefunctions,        ONLY : evc
   USE becmod,               ONLY : becp, calbec
   USE control_flags,        ONLY : offload_type
   USE basis,                ONLY : natomwfc, wfcatom, swfcatom
   USE force_mod,            ONLY : eigenval, eigenvect, overlap_inv, at_dy, at_dj
   USE mp_bands,             ONLY : intra_bgrp_comm
   USE mp,                   ONLY : mp_sum
   USE noncollin_module,     ONLY : noncolin, npol
   IMPLICIT NONE
   !
   ! I/O variables 
   !
   COMPLEX(DP), INTENT(IN)  :: spsi(npwx*npol,nbnd)
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
   dwfcU(:,:),        & ! the derivative of the (ortho-atomic) wavefunctions
   dwfca(:,:),        & ! the derivative of the atomic wavefunctions
   doverlap(:,:),     & ! derivative of the overlap matrix  
   doverlap_us(:,:),  & ! USPP contribution to doverlap
   doverlap_inv(:,:)    ! derivative of (O^{-1/2})_JI (note the transposition)   
   REAL (DP), ALLOCATABLE :: &
   gk(:,:), & ! k+G
   qm1(:),  & ! 1/|k+G|
   a1_temp(:), a2_temp(:) 
   COMPLEX (DP) :: temp, temp2
   !
   CALL start_clock('dprojdepsilon')
   CALL start_clock('dprojdeps1')
   ! 
   !$acc data present_or_copyin(spsi) present_or_copyout(dproj)
   !
   ! Number of plane waves at the k point with the index ik
   npw = ngk(ik)
   !
   !$acc kernels
   dproj(:,:) = (0.d0, 0.d0)
   !$acc end kernels
   !
   ! At first the derivatives of the atomic wfcs: we compute the term
   ! <d\fi^{at}_{I,m1}/d\epsilon(ipol,jpol)|S|\psi_{k,v,s}>
   !
   ALLOCATE ( qm1(npwx), gk(3,npwx) )
   ALLOCATE ( dwfcU(npwx*npol,nwfcU) )
   ALLOCATE (a1_temp(npw), a2_temp(npw))
   !$acc data create(dwfcU) present(overlap_inv,wfcU)
   !$acc kernels
   dwfcU(:,:) = (0.d0, 0.d0)
   !$acc end kernels
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
      a1_temp(ig) = -gk(jpol,ig)
      !
      ! - (k+G)_ipol * (k+G)_jpol / |k+G|
      a2_temp(ig) = -gk(ipol,ig) * gk(jpol,ig) * qm1(ig)
      !
   ENDDO
   !
   !$acc data copyin(a1_temp, a2_temp, at_dj, at_dy)
   !
   DO na = 1, nat
      nt = ityp(na)
      ldim_std = 2*Hubbard_l(nt)+1
      IF (is_hubbard(nt) .OR. is_hubbard_back(nt)) THEN
         IF (Hubbard_projectors.EQ."atomic") THEN
            DO m1 = 1, ldim_u(nt)
               IF (m1.LE.ldim_std) THEN
                  offpmU = offsetU(na)
                  offpm  = oatwfc(na)
               ELSE
                  offpmU = offsetU_back(na) - ldim_std
                  offpm  = oatwfc_back(na)  - ldim_std
                  IF (backall(nt) .AND. m1.GT.ldim_std+2*Hubbard_l2(nt)+1) THEN
                     offpmU = offsetU_back1(na) - ldim_std - 2*Hubbard_l2(nt) - 1
                     offpm  = oatwfc_back1(na)  - ldim_std - 2*Hubbard_l2(nt) - 1
                  ENDIF
               ENDIF
               IF (m1>ldim_std .and. noncolin) CALL errore("dprojdepsilon_k", &
                             " Stress with background and noncollinear is not supported",1)
               !$acc parallel loop
               DO ig = 1, npw
                  dwfcU(ig,offpmU+m1) = at_dy(ig,offpm+m1) * a1_temp(ig) &
                             + at_dj(ig,offpm+m1) * a2_temp(ig)
                  IF (noncolin) THEN  
                     dwfcU(ig+npwx,offpmU+m1+ldim_std) = &
                                        at_dy(ig+npwx,offpm+m1+ldim_std)*a1_temp(ig) & 
                                      + at_dj(ig+npwx,offpm+m1+ldim_std)*a2_temp(ig)
                  ENDIF
               ENDDO
            ENDDO      
         ELSEIF (Hubbard_projectors.EQ."ortho-atomic") THEN
            DO m1 = 1, ldim_std*npol             
               offpmU = offsetU(na)
               offpm  = oatwfc(na)
               !$acc parallel loop
               DO ig = 1, npw
                  temp = (0.0d0, 0.0d0)
                  temp2 = (0.0d0, 0.0d0)
                  DO m2 = 1, natomwfc
                     temp = temp + overlap_inv(offpm+m1,m2) * &
                         ( at_dy(ig,m2) * a1_temp(ig) + at_dj(ig,m2) * a2_temp(ig) )
                     IF (noncolin) temp2 = temp2 + &
                              overlap_inv(offpm+m1,m2) * &
                              ( at_dy(ig+npwx,m2) * a1_temp(ig) &
                              + at_dj(ig+npwx,m2) * a2_temp(ig) )   
                  ENDDO   
                  dwfcU(ig,offpmU+m1) = dwfcU(ig,offpmU+m1) + temp
                  IF (noncolin) dwfcU(ig+npwx,offpmU+m1) = &
                           dwfcU(ig+npwx,offpmU+m1) + temp2 
               ENDDO
            ENDDO
         ENDIF
      ENDIF
   ENDDO
   !
   ! The diagonal term
   IF (ipol.EQ.jpol) THEN
      !$acc kernels 
      dwfcU(1:npw,:) = dwfcU(1:npw,:) - wfcU(1:npw,:)*0.5d0
      IF (noncolin) dwfcU(1+npwx:npwx+npw,:) = dwfcU(1+npwx:npwx+npw,:) - wfcU(1+npwx:npwx+npw,:)*0.5d0
      !$acc end kernels
   ENDIF   
   CALL stop_clock('dprojdeps1')
   CALL start_clock('dprojdeps2')
   !
   ! 2. Contribution due to the derivative of (O^{-1/2})_JI which
   !    is multiplied by atomic wavefunctions (only for ortho-atomic case)
   !
   IF (Hubbard_projectors.EQ."ortho-atomic") THEN
      !
      ! Compute the derivative dO_IJ/d\epsilon(ipol,jpol)
      !
      ALLOCATE (dwfca(npwx*npol,natomwfc))
      ALLOCATE (doverlap(natomwfc,natomwfc))
      ALLOCATE (doverlap_inv(natomwfc,natomwfc))
      !$acc data create(dwfca,doverlap,doverlap_inv) present_or_copyin(swfcatom, wfcatom)
      !$acc kernels
      dwfca(:,:) = (0.d0, 0.d0)
      doverlap(:,:) = (0.0d0, 0.0d0)
      doverlap_inv(:,:) = (0.0d0, 0.0d0)
      !$acc end kernels
      !
      ! Calculate:
      ! doverlap = < dphi_I/d\epsilon(ipol,jpol) | S | phi_J > 
      !            + < phi_I | S | dphi_J/d\epsilon(ipol,jpol) >
      ! Note that the second term is the hermitian conjugate of the first
      !$acc parallel loop collapse(2)
      DO m1 = 1, natomwfc
         DO ig = 1, npw
            dwfca(ig,m1) = at_dy(ig,m1)*a1_temp(ig) + at_dj(ig,m1)*a2_temp(ig)
            IF (noncolin) THEN
               dwfca(ig+npwx,m1) = at_dy(ig+npwx,m1)*a1_temp(ig) + at_dj(ig+npwx,m1)*a2_temp(ig)
            ENDIF
         ENDDO
      ENDDO
      !
      IF (ipol.EQ.jpol) THEN
         !$acc parallel loop collapse(2)
         DO m1 = 1, natomwfc
            DO ig = 1, npw
               dwfca(ig,m1) = dwfca(ig,m1) - 0.5_dp*wfcatom(ig,m1)
               IF (noncolin) THEN
                  dwfca(ig+npwx,m1) = dwfca(ig+npwx,m1)  - 0.5_dp*wfcatom(ig+npwx,m1)
               ENDIF
            ENDDO
         ENDDO
      END IF
      !
      !$acc host_data use_device(dwfca, swfcatom, doverlap)
      CALL MYZGEMM('C','N', natomwfc, natomwfc, npwx*npol, (1.0_dp,0.0_dp), &
           dwfca, npwx*npol, swfcatom, npwx*npol, (0.0_dp,0.0_dp), &
           doverlap, natomwfc) 
      ! Sum over G vectors
      CALL mp_sum( doverlap, intra_bgrp_comm )
      !$acc end host_data
      !$acc kernels
      doverlap = doverlap + CONJG(TRANSPOSE(doverlap))
      !$acc end kernels
      !
      ! USPP term in dO_IJ/d\epsilon(ipol,jpol)
      !
      IF (okvan) THEN
         ! Calculate doverlap_us = < phi_I | dS/d\epsilon(ipol,jpol) | phi_J >
         ALLOCATE (doverlap_us(natomwfc,natomwfc))
         !$acc data create(doverlap_us) 
         CALL matrix_element_of_dSdepsilon (ik, ipol, jpol, &
              natomwfc, wfcatom, natomwfc, wfcatom, doverlap_us, 1, natomwfc, 0, .false.)
         ! Sum up the "normal" and ultrasoft terms
         !$acc parallel loop collapse(2) 
         DO m2 = 1, natomwfc
            DO m1 = 1, natomwfc
               doverlap(m1,m2) = doverlap(m1,m2) + doverlap_us(m1,m2)
            ENDDO
         ENDDO
         !$acc end data
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
      ! and add it to another term (see above).
      ! Note, doverlap_inv is d(O^{-1/2}) not transposed. The transposition 
      ! of d(O^{-1/2}) is taken into account via a proper usage of the order
      ! of indices in doverlap_inv: 
      ! dwfc(ig,offpmU+m1) = dwfc(ig,offpmU+m1) + wfcatom(ig,m2) * doverlap_inv(m2,offpm+m1)
      ! where m1=1,ldim_u(nt); m2=1,natomwfc; ig=1,npw
      !
      DO na = 1, nat
         nt = ityp(na)
         IF (is_hubbard(nt) .OR. is_hubbard_back(nt)) THEN
            offpmU = offsetU(na)
            offpm  = oatwfc(na)
            !$acc host_data use_device(wfcatom, doverlap_inv, dwfcU)
            CALL MYZGEMM('N','N', npwx*npol, ldim_u(nt)*npol, natomwfc, (1.d0,0.d0), &
                  wfcatom, npwx*npol, doverlap_inv(:,offpm+1:offpm+ldim_u(nt)*npol), &
                  natomwfc, (1.d0,0.d0), dwfcU(:,offpmU+1:offpmU+ldim_u(nt)*npol), npwx*npol)
            !$acc end host_data
         ENDIF
      ENDDO
      !
      !$acc end data
      DEALLOCATE (doverlap_inv)
      DEALLOCATE (doverlap)
      DEALLOCATE (dwfca)
      !
   ENDIF
   !
   ! Compute dproj = <dwfc|S|psi> = <dwfc|spsi>
   IF (noncolin) THEN
      !$acc host_data use_device(dwfcU,spsi,dproj)
      CALL MYZGEMM('C','N', nwfcU, nbnd, npwx*npol, (1.d0,0.d0), &
            dwfcU, npwx*npol, spsi, npwx*npol, (0.d0,0.d0), &
            dproj, nwfcU)   
      CALL mp_sum( dproj, intra_bgrp_comm )
      !$acc end host_data
   ELSE   
      CALL calbec( offload_type, npw, dwfcU, spsi, dproj )
   ENDIF
   !
   !$acc end data
   !$acc end data
   DEALLOCATE ( dwfcU, qm1, gk)
   DEALLOCATE(a1_temp, a2_temp)
   CALL stop_clock('dprojdeps2')
   CALL start_clock('dprojdeps3')
   !
   ! Now the derivatives of the beta functions: we compute the term
   ! <\phi^{at}_{I,m1}|dS/d\epsilon(ipol,jpol)|\psi_{k,v,s}>
   !
   IF (okvan) THEN
      ALLOCATE(dproj_us(nwfcU,nb_s:nb_e))
      !$acc data create(dproj_us) 
      CALL matrix_element_of_dSdepsilon (ik, ipol, jpol, &
                         nwfcU, wfcU, nbnd, evc, dproj_us, nb_s, nb_e, mykey, .true.)
      ! dproj + dproj_us
      !$acc parallel loop 
      DO m1 = 1, nwfcU
         dproj(m1,nb_s:nb_e) = dproj(m1,nb_s:nb_e) + dproj_us(m1,:)
      ENDDO
      !$acc end data
      DEALLOCATE(dproj_us)
   ENDIF
   !
   !$acc end data
   !
   CALL stop_clock('dprojdeps3')
   CALL stop_clock('dprojdepsilon')
   !
   RETURN
   !
END SUBROUTINE dprojdepsilon_k
!
SUBROUTINE matrix_element_of_dSdepsilon (ik, ipol, jpol, lA, A, lB, B, A_dS_B, lB_s, lB_e, mykey, flag)
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
   USE uspp,                 ONLY : nkb, vkb, qq_at, qq_so, okvan
   USE uspp_param,           ONLY : nh, upf
   USE wavefunctions,        ONLY : evc
   USE becmod,               ONLY : calbec
   USE control_flags,        ONLY : offload_type
   USE klist,                ONLY : xk, igk_k, ngk
   USE force_mod,            ONLY : us_dy, us_dj
   USE mp_bands,             ONLY : intra_bgrp_comm
   USE mp,                   ONLY : mp_sum
   USE noncollin_module,     ONLY : noncolin, npol
   USE ldaU,                 ONLY : offsetU, Hubbard_l, is_hubbard
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
   LOGICAL, INTENT(IN)      :: flag  ! noncollinear: controlling whether 
                                     ! calculating <phi|dS|PSI> 
                                     ! or          <phi|dS|PHI> (= .true.)   
   COMPLEX(DP), INTENT(IN)  :: A(npwx*npol,lA)
   COMPLEX(DP), INTENT(IN)  :: B(npwx*npol,lB)
   COMPLEX(DP), INTENT(OUT) :: A_dS_B(lA,lB_s:lB_e)
   INTEGER,     INTENT(IN)  :: mykey
   !
   ! Local variables
   !
   INTEGER :: npw, i, nt, na, nb, ih, jh, ig, iA, iB, ijkb0, mU, mD, ldim, nt1
   REAL(DP) :: gvec
   COMPLEX (DP), ALLOCATABLE :: Adbeta(:,:), Abeta(:,:), &
                                dbetaB(:,:), betaB(:,:), &
                                aux(:,:), qq(:,:)
   REAL (DP) :: q, a1, a2
   REAL (DP), PARAMETER :: eps = 1.0d-8
   REAL (DP), ALLOCATABLE :: gk(:,:), qm1(:)
   INTEGER :: nh_nt
   REAL (DP), ALLOCATABLE :: a1_temp(:), a2_temp(:)
   !
   IF (.NOT.okvan) RETURN
   !
   !$acc data present_or_copyin(A,B) present_or_copyout(A_ds_B)
   !
   !$acc kernels
   A_dS_B(:,:) = (0.0d0, 0.0d0)
   !$acc end kernels
   npw = ngk(ik)
   !
   ALLOCATE ( qm1(npwx), gk(3,npwx) )
   ALLOCATE (a1_temp(npw), a2_temp(npw))
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
      a1_temp(ig) = -gk(jpol,ig)
      a2_temp(ig) = -gk(ipol,ig)*gk(jpol,ig)*qm1(ig)
   ENDDO
   !
   ijkb0 = 0
   !
   !$acc data present(us_dj) copyin(qq_at, a1_temp, a2_temp)
   DO nt = 1, ntyp
      !
      ALLOCATE ( Adbeta(lA,npol*nh(nt)) )
      ALLOCATE ( Abeta(lA,npol*nh(nt)) )
      ALLOCATE ( dbetaB(npol*nh(nt),lB) )
      ALLOCATE ( betaB(npol*nh(nt),lB) )
      ALLOCATE ( qq(npol*nh(nt),npol*nh(nt)) )
      !
      nh_nt = nh(nt)
      !$acc data create(Adbeta,Abeta,dbetaB,betaB,qq)
      !
      DO na = 1, nat
         !
         IF ( ityp(na).EQ.nt ) THEN
            !
            IF (noncolin) THEN
               IF ( upf(nt)%has_so ) THEN
                  !$acc kernels
                  qq(1:nh(nt),1:nh(nt)) = qq_so(:,:,1,nt)
                  qq(1:nh(nt),1+nh(nt):nh(nt)*npol) = qq_so(:,:,2,nt)
                  qq(1+nh(nt):nh(nt)*npol,1:nh(nt)) = qq_so(:,:,3,nt)
                  qq(1+nh(nt):nh(nt)*npol,1+nh(nt):nh(nt)*npol) = qq_so(:,:,4,nt)
                  !$acc end kernels
               ELSE
                  !$acc kernels
                  qq(1:nh(nt),1:nh(nt)) = qq_at(1:nh(nt),1:nh(nt),na)
                  qq(1:nh(nt),1+nh(nt):nh(nt)*npol) = (0.0,0.0)
                  qq(1+nh(nt):nh(nt)*npol,1:nh(nt)) = (0.0,0.0)
                  qq(1+nh(nt):nh(nt)*npol,1+nh(nt):nh(nt)*npol) = qq_at(1:nh(nt),1:nh(nt),na)
                  !$acc end kernels
               ENDIF
            ELSE
               !$acc parallel loop collapse(2) 
               DO jh = 1, nh_nt 
                  DO ih = 1, nh_nt
                     qq(ih,jh) = CMPLX(qq_at(ih,jh,na), 0.0d0, kind=DP)
                  ENDDO
               ENDDO
            ENDIF
            !
            ! aux is used as a workspace
            ALLOCATE ( aux(npwx*npol,nh(nt)*npol) )
            !$acc data create(aux)
            !$acc kernels
            aux=(0.0,0.0)
            !$acc end kernels
            !
            !$acc parallel loop collapse(2) 
            DO ih = 1, nh_nt !nh(nt)
               ! now we compute the true dbeta function
               DO ig = 1, npw
                  aux(ig,ih) = us_dy(ig,ijkb0+ih) * a1_temp(ig) + us_dj(ig,ijkb0+ih) * a2_temp(ig)
                  IF (noncolin) aux(ig+npwx,ih+nh(nt)) = us_dy(ig,ijkb0+ih) * a1_temp(ig) &
                                                + us_dj(ig,ijkb0+ih) * a2_temp(ig)  
               ENDDO
            ENDDO
            !
            IF (ipol.EQ.jpol) THEN
               !$acc parallel loop collapse(2)
               DO ih = 1, nh_nt !nh(nt)
                  DO ig = 1, npw
                     aux(ig,ih) = aux(ig,ih) - vkb(ig,ijkb0+ih)*0.5d0
                     IF (noncolin) aux(ig+npwx,ih+nh(nt)) = aux(ig+npwx,ih+nh(nt)) &
                                   - vkb(ig,ijkb0+ih)*0.5d0
                  ENDDO
               ENDDO   
            ENDIF 
            IF (noncolin) THEN
               ! Calculate betaB = <dbeta|B>
               ! dbetaB(:,1       : nh(nt))      = spin up
               ! dbetaB(:,1+nh(nt): nh(nt)*npol) = spin down
               ! dbetaB=(0.0,0.0)
               !$acc host_data use_device(aux, A, B, Adbeta, dbetaB)
               CALL MYZGEMM ('C', 'N', nh(nt)*npol, lB, npwx*npol, (1.0_DP, 0.0_DP), aux, &
                          npwx*npol, B, npwx*npol, (0.0_DP, 0.0_DP), dbetaB, nh(nt)*npol)
               CALL mp_sum( dbetaB(:, 1:lB) , intra_bgrp_comm )               
               !
               ! Calculate Adbeta = <A|dbeta>
               ! Adbeta(:,1       : nh(nt))      = spin up
               ! Adbeta(:,1+nh(nt): nh(nt)*npol) = spin down      
               !Adbeta=(0.0,0.0)
               CALL MYZGEMM ('C', 'N', lA, nh(nt)*npol, npwx*npol, (1.0_DP, 0.0_DP), A, &
                          npwx*npol, aux, npwx*npol, (0.0_DP, 0.0_DP), Adbeta, lA)
               CALL mp_sum( Adbeta(:, 1:nh(nt)*npol) , intra_bgrp_comm )
               !$acc end host_data 
            ELSE        
               ! Calculate dbetaB = <dbeta|B> 
               CALL calbec(offload_type, npw, aux, B, dbetaB )
               !
               ! Calculate Adbeta = <A|dbeta>
               CALL calbec(offload_type, npw, A, aux, Adbeta )
            ENDIF   
            !
            ! aux is now used as a work space to store vkb
            !$acc parallel loop collapse(2)
            DO ih = 1, nh_nt  !nh(nt)
               DO ig = 1, npw
                  aux(ig,ih) = vkb(ig,ijkb0+ih)
                  IF (noncolin) aux(ig+npwx,ih+nh(nt)) = vkb(ig,ijkb0+ih)
               ENDDO
            ENDDO
            !
            IF (noncolin) THEN
                ! Calculate Abeta = <A|beta>      
                ! (same as Adbeta)
                !
                !Abeta=(0.0,0.0)
                !$acc host_data use_device(aux, A, B, Abeta, betaB)
                CALL MYZGEMM ('C', 'N', lA, nh(nt)*npol, npwx*npol, (1.0_DP, 0.0_DP), A, &
                           npwx*npol, aux, npwx*npol, (0.0_DP, 0.0_DP), Abeta, lA)
                CALL mp_sum( Abeta(:, 1:nh(nt)*npol) , intra_bgrp_comm )
                !
                ! Calculate betaB = <beta|B>
                ! (same as dbetaB)
                !
                !betaB=(0.0,0.0)
                CALL MYZGEMM ('C', 'N', nh(nt)*npol, lB, npwx*npol, (1.0_DP, 0.0_DP), aux, &
                           npwx*npol, B, npwx*npol, (0.0_DP, 0.0_DP), betaB, nh(nt)*npol)
                CALL mp_sum( betaB(:, 1:lB) , intra_bgrp_comm )
               !$acc end host_data 
             ELSE
               ! Calculate Abeta = <A|beta>
               CALL calbec( offload_type, npw, A, aux, Abeta )
               !
               ! Calculate betaB = <beta|B>
               CALL calbec( offload_type, npw, aux, B, betaB )
             ENDIF
            !
            !$acc end data
            DEALLOCATE ( aux )
            !
            ALLOCATE ( aux(nh(nt)*npol, lB) )
            !$acc data create(aux)
            !$acc kernels
            aux(:,:)=(0.0,0.0)            
            !$acc end kernels
            !
            IF (noncolin) THEN
               ! aux(:, 1:nh(nt))             = \sum_jh qq(1,jh)*dbetaB(1,jh) + qq(2,jh)*dbetaB(2,jh)
               ! aux(:, 1+nh(nt):nh(nt)*npol) = \sum_jh qq(3,jh)*dbetaB(1,jh) + qq(4,jh)*dbetaB(2,jh)
               !
               ! spin up
               !$acc host_data use_device(qq,aux,dbetaB)
               CALL MYZGEMM('N', 'N', nh(nt), lB_e-lB_s+1, nh(nt)*npol, (1.0d0,0.0d0), &
                          qq(1, 1), nh(nt)*npol, &
                          dbetaB(1, lB_s), nh(nt)*npol, (0.0d0,0.0d0), &
                          aux(1, lB_s), nh(nt)*npol)
               ! spin down   
               CALL MYZGEMM('N', 'N', nh(nt), lB_e-lB_s+1, nh(nt)*npol, (1.0d0,0.0d0), &
                          qq(1+nh(nt), 1), nh(nt)*npol, &
                          dbetaB(1, lB_s), nh(nt)*npol, (0.0d0,0.0d0), &
                          aux(1+nh(nt), lB_s), nh(nt)*npol)
               !$acc end host_data
            ELSE
               ! Calculate \sum_jh qq_at(ih,jh) * dbetaB(jh)     
               !$acc host_data use_device(qq,dbetaB,aux)
               CALL MYZGEMM('N', 'N', nh(nt), lB_e-lB_s+1, nh(nt), (1.0d0,0.0d0), &
                            qq, nh(nt), dbetaB(1,lB_s),    nh(nt), (0.0d0,0.0d0), &
                            aux(1,lB_s), nh(nt))
               !$acc end host_data
            ENDIF
            !$acc kernels
            dbetaB(:,:) = aux(:,:)
            !$acc end kernels
            !
            IF (noncolin) THEN
               ! (same as dbetaB)     
               !
               ! spin up 
               !$acc host_data use_device(qq,betaB,aux)
               CALL MYZGEMM('N', 'N', nh(nt), lB_e-lB_s+1, nh(nt)*npol, (1.0d0,0.0d0), &
                          qq(1, 1), nh(nt)*npol, &
                          betaB(1, lB_s), nh(nt)*npol, (0.0d0,0.0d0), &
                          aux(1, lB_s), nh(nt)*npol)
               ! spin down   
               CALL MYZGEMM('N', 'N', nh(nt), lB_e-lB_s+1, nh(nt)*npol, (1.0d0,0.0d0), &
                          qq(1+nh(nt), 1), nh(nt)*npol, &
                          betaB(1, lB_s), nh(nt)*npol, (0.0d0,0.0d0), &
                          aux(1+nh(nt), lB_s), nh(nt)*npol)
              !$acc end host_data
            ELSE
               ! Calculate \sum_jh qq_at(ih,jh) * betaB(jh)
               !$acc host_data use_device(qq,betaB,aux)
               CALL MYZGEMM('N', 'N', nh(nt), lB_e-lB_s+1, nh(nt), (1.0d0,0.0d0), &
                            qq, nh(nt), betaB(1,lB_s),     nh(nt), (0.0d0,0.0d0), &
                            aux(1,lB_s), nh(nt))
              !$acc end host_data
            ENDIF
            !
            !$acc kernels
            betaB(:,:) = aux(:,:)
            !$acc end kernels
            !
            !$acc end data
            DEALLOCATE ( aux )
            !
            ijkb0 = ijkb0 + nh(nt)
            !
            ! A_dS_B(iA,iB) = \sum_ih [Adbeta(iA,ih) * betapsi(ih,iB) +
            !                          Abeta(iA,ih)  * dbetaB(ih,iB)] 
            ! Only A_dS_B(:,lB_s:lB_e) are calculated
            !
            IF ( mykey == 0 ) THEN
               IF (noncolin) THEN
                  nt1 = nh(nt) + 1
                  !$acc host_data use_device(Adbeta,betaB,Abeta,dbetaB,A_dS_B)      
                  IF ( flag ) THEN
                     DO nb = 1, nat
                        IF ( is_hubbard(ityp(nb)) ) THEN
                           ldim = 2*hubbard_l(ityp(nb)) + 1
                           mU = offsetU(nb) + 1
                           mD = mU + ldim
                           !
                           ! spin up
                           CALL MYZGEMM('N', 'N', ldim, lB_e-lB_s+1, nh(nt), (1.0d0,0.0d0), &
                                       Adbeta(mU,1), lA, &
                                       betaB(1, lB_s), nh(nt)*npol, (1.0d0,0.0d0), &
                                       A_dS_B(mU, lB_s), lA)
                           CALL MYZGEMM('N', 'N', ldim, lB_e-lB_s+1, nh(nt), (1.0d0,0.0d0), &
                                    Abeta(mU,1), lA, &
                                    dbetaB(1, lB_s), nh(nt)*npol, (1.0d0,0.0d0), &
                                    A_dS_B(mU, lB_s), lA)
                           ! spin down
                           CALL MYZGEMM('N', 'N', ldim, lB_e-lB_s+1, nh(nt), (1.0d0,0.0d0), &
                                       Adbeta(mD, nt1), lA, &
                                       betaB(nt1, lB_s), nh(nt)*npol, (1.0d0,0.0d0), &
                                       A_dS_B(mD,lB_s), lA)
                           CALL MYZGEMM('N', 'N', ldim, lB_e-lB_s+1, nh(nt), (1.0d0,0.0d0), &
                                       Abeta(mD, nt1), lA, &
                                       dbetaB(nt1, lB_s), nh(nt)*npol, (1.0d0,0.0d0), &
                                       A_dS_B(mD,lB_s), lA)
                        ENDIF
                     ENDDO
                  ELSEIF ( .not. flag ) THEN
                     CALL MYZGEMM('N', 'N', lA, lB_e-lB_s+1, nh(nt)*npol, (1.0d0,0.0d0), &
                                  Adbeta, lA, betaB(1,lB_s), nh(nt)*npol, (1.0d0,0.0d0), &
                                  A_dS_B(1,lB_s), lA)
                     CALL MYZGEMM('N', 'N', lA, lB_e-lB_s+1, nh(nt)*npol, (1.0d0,0.0d0), &
                                  Abeta, lA, dbetaB(1,lB_s), nh(nt)*npol, (1.0d0,0.0d0), &
                                  A_dS_B(1,lB_s), lA)
                  ENDIF
                  !$acc end host_data   
               ELSE
                  !$acc host_data use_device(Adbeta,betaB,Abeta,dbetaB,A_dS_B)      
                  CALL MYZGEMM('N', 'N', lA, lB_e-lB_s+1, nh(nt), (1.0d0,0.0d0), &
                             Adbeta, lA, betaB(1,lB_s), nh(nt), (1.0d0,0.0d0), &
                             A_dS_B(1,lB_s), lA)
                  CALL MYZGEMM('N', 'N', lA, lB_e-lB_s+1, nh(nt), (1.0d0,0.0d0), &
                             Abeta, lA, dbetaB(1,lB_s), nh(nt), (1.0d0,0.0d0), &
                             A_dS_B(1,lB_s), lA)
                  !$acc end host_data   
                  !
               ENDIF
            ENDIF
         ENDIF
         !
      ENDDO
      !
      !$acc end data
      DEALLOCATE (dbetaB, betaB, Abeta, Adbeta, qq)
      ! 
   ENDDO
   !$acc end data
   !
   DEALLOCATE (a1_temp, a2_temp)
   DEALLOCATE ( qm1, gk )
   !
   !$acc end data
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
                                    Hubbard_projectors, Hubbard_l, Hubbard_l2, &
                                    backall
   USE lsda_mod,             ONLY : lsda, nspin, isk
   USE wvfct,                ONLY : nbnd, npwx, wg
   USE uspp,                 ONLY : nkb, vkb, qq_at, okvan
   USE uspp_param,           ONLY : nh
   USE wavefunctions,        ONLY : evc
   USE becmod,               ONLY : becp, calbec
   USE control_flags,        ONLY : offload_type
   USE force_mod,            ONLY : at_dy, at_dj, us_dy, us_dj
   !
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
              offpm, offpmU, m1, ldim_std, nh_nt
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

   REAL (DP), ALLOCATABLE :: gk(:,:), qm1(:), a1_temp(:), a2_temp(:)
   !       gk(3,npwx),
   !       qm1(npwx)
   REAL (DP) :: temp
   !
   CALL start_clock('dprojdepsilon')
   !
   ! See the implementation in dprojdepsilon_k
   IF (Hubbard_projectors.EQ."ortho-atomic") CALL errore("dprojdtau_gamma", &
                    " Forces with gamma-only and ortho-atomic are not supported",1)
   !
   !$acc data present_or_copyin(spsi) present_or_copyout(dproj)

   ! Number of plane waves at the k point with the index ik
   npw = ngk(ik)
   !
   !$acc kernels
   dproj(:,:) = 0.d0
   !$acc end kernels
   !
   ! At first the derivatives of the atomic wfcs: we compute the term
   ! <d\fi^{at}_{I,m1}/d\epsilon(ipol,jpol)|S|\psi_{k,v,s}>
   !
   ALLOCATE ( qm1(npwx), gk(3,npwx) )
   ALLOCATE ( dwfc(npwx,nwfcU) )
   ALLOCATE ( a1_temp(npw), a2_temp(npw) )
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
      a1_temp(ig) = -gk(jpol,ig)
      !
      ! - (k+G)_ipol * (k+G)_jpol / |k+G|
      a2_temp(ig) = -gk(ipol,ig)*gk(jpol,ig)*qm1(ig)
      !
   ENDDO
   !
   !$acc data present(us_dy, us_dj, wfcU) copyin(a1_temp, a2_temp, at_dy, at_dj, qq_at)
   !$acc data create(dwfc) 
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
               IF (backall(nt) .AND. m1.GT.ldim_std+2*Hubbard_l2(nt)+1) THEN
                  offpmU = offsetU_back1(na) - ldim_std - 2*Hubbard_l2(nt) - 1
                  offpm  = oatwfc_back1(na)  - ldim_std - 2*Hubbard_l2(nt) - 1
               ENDIF
            ENDIF
            !$acc parallel loop
            DO ig = 1, npw
               dwfc(ig,offpmU+m1) = at_dy(ig,offpm+m1) * a1_temp(ig) + at_dj(ig,offpm+m1) * a2_temp(ig)
            ENDDO
         ENDDO
      ENDIF
   ENDDO
   !
   IF (ipol.EQ.jpol) THEN
      !$acc kernels
      dwfc(1:npw,:) = dwfc(1:npw,:) - wfcU(1:npw,:)*0.5d0
      !$acc end kernels
   ENDIF
   !
   CALL calbec ( offload_type, npw, dwfc, spsi, dproj )
   !
   !$acc end data
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
      !
      DO nt = 1, ntyp
         !
         ALLOCATE (dbeta(npwx,nh(nt)), dbetapsi(nh(nt),nbnd), betapsi(nh(nt),nbnd), &
                   wfatbeta(nwfcU,nh(nt)), wfatdbeta(nwfcU,nh(nt)), betapsi0(nh(nt),nbnd) )
         nh_nt = nh(nt)
         !$acc data create(dbeta, dbetapsi, betapsi, wfatbeta, wfatdbeta, betapsi0)
         !
         DO na = 1, nat
            !
            IF ( ityp(na).EQ.nt ) THEN
               !
               !$acc parallel loop collapse(2)
               DO ih = 1, nh_nt
                  ! now we compute the true dbeta function
                  DO ig = 1, npw
                     dbeta(ig,ih) = us_dy(ig,ijkb0+ih)*a1_temp(ig) + us_dj(ig,ijkb0+ih) * a2_temp(ig)
                  ENDDO
               ENDDO
               !
               IF (ipol.EQ.jpol) THEN
                  !$acc parallel loop collapse(2)     
                  DO ih = 1, nh_nt
                     DO ig = 1, npw
                        dbeta(ig,ih) = dbeta(ig,ih) - vkb(ig,ijkb0+ih)*0.5d0
                     ENDDO
                  ENDDO
               ENDIF   
               !
               CALL calbec(offload_type, npw, dbeta, evc, dbetapsi )
               CALL calbec(offload_type, npw, wfcU, dbeta, wfatdbeta )
               !
               ! dbeta is now used as work space to store vkb
               !$acc parallel loop collapse(2)
               DO ih = 1, nh_nt
                  DO ig = 1, npw
                     dbeta(ig,ih) = vkb(ig,ijkb0+ih)
                  ENDDO
               ENDDO
               !
               CALL calbec(offload_type, npw, wfcU, dbeta, wfatbeta )
               CALL calbec(offload_type, npw, dbeta, evc, betapsi0 )
               !
               ! here starts band parallelization
               ! beta is here used as work space to calculate dbetapsi
               !
               !$acc kernels
               betapsi(:,:) = 0.0_dp
               !$acc end kernels
               !
               !$acc parallel loop collapse(2)
               DO ibnd = nb_s, nb_e
                  DO ih = 1, nh_nt
                     temp = 0.0d0
                     DO jh = 1, nh_nt
                        temp = temp + qq_at(ih,jh,na)  * dbetapsi(jh,ibnd)
                     ENDDO
                     betapsi(ih,ibnd) = betapsi(ih,ibnd) + temp
                  ENDDO
               ENDDO
               !
               !$acc kernels
               dbetapsi(:,:) = betapsi(:,:)
               betapsi(:,:) = 0.0_DP
               !$acc end kernels
               !
               !$acc parallel loop collapse(2)
               DO ibnd = nb_s, nb_e
                  DO ih = 1, nh_nt
                     temp = 0.0d0
                     DO jh = 1, nh_nt
                        temp = temp + qq_at(ih,jh,na)  * betapsi0(jh,ibnd)
                     ENDDO
                     betapsi(ih,ibnd) = betapsi(ih,ibnd) + temp
                  ENDDO
               ENDDO
               !
               ijkb0 = ijkb0 + nh(nt)
               !
               ! dproj(iwf,ibnd) = \sum_ih wfatdbeta(iwf,ih)*betapsi(ih,ibnd) +
               !                           wfatbeta(iwf,ih)*dbetapsi(ih,ibnd) 
               !
               IF ( mykey == 0 .AND. nh(nt) > 0 ) THEN
                  !$acc host_data use_device(wfatdbeta,betapsi,dproj,wfatbeta,dbetapsi)
                  CALL MYDGEMM('N','N',nwfcU, nb_e-nb_s+1, nh(nt), 1.0_dp,  &
                       wfatdbeta, nwfcU, betapsi(1,nb_s), nh(nt), 1.0_dp,&
                       dproj(1,nb_s), nwfcU)
                  CALL MYDGEMM('N','N',nwfcU, nb_e-nb_s+1, nh(nt), 1.0_dp,  &
                       wfatbeta, nwfcU, dbetapsi(1,nb_s), nh(nt), 1.0_dp,&
                       dproj(1,nb_s), nwfcU)
                  !$acc end host_data
               ENDIF
               ! end band parallelization - only dproj(1,nb_s:nb_e) are calculated
            ENDIF
         ENDDO
         !$acc end data
         DEALLOCATE (dbeta, dbetapsi, betapsi, wfatbeta, wfatdbeta, betapsi0)
      ENDDO
      ! 
   ENDIF
   !
   !$acc end data
   DEALLOCATE ( qm1, gk )
   DEALLOCATE ( a1_temp, a2_temp )
   !
   !$acc end data
   CALL stop_clock('dprojdepsilon')
   !
   RETURN
   !
END SUBROUTINE dprojdepsilon_gamma
