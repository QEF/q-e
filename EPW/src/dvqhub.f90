  !
  ! Copyright (C) 2023-2026 EPW-Collaboration
  ! Copyright (C) 2001-2018 Quantum ESPRESSO
  ! This file is distributed under the terms
  ! GNU General Public License. See the file
  ! in the root directory of the present dis
  ! or http://www.gnu.org/copyleft/gpl.txt .
  !
  !
  !-----------------------------------------------------------------------
  MODULE dvqhub
  !-----------------------------------------------------------------------
  !! 
  !! This module contains the routines to compute dV_hub/dtau *|\phi><\phi|\psi> or its
  !! components
  !!
  IMPLICIT NONE
  !
  CONTAINS
    ! 
    !---------------------------------------------------------------------
    SUBROUTINE set_irr_sym_new2 ( t, tmq, npertx )
    !---------------------------------------------------------------------
    !! This subroutine computes the matrices which represent the small
    !! group of q on the pattern basis.
    !
    USE kinds, ONLY : DP
    USE constants, ONLY: tpi
    USE ions_base, ONLY : nat
    USE cell_base, ONLY : at, bg
    USE symm_base, ONLY : s, irt, t_rev
    USE modes,     ONLY : u, nirr, npert
    USE control_flags, ONLY : modenum
    USE qpoint,       ONLY : xq
    USE lr_symm_base, ONLY : nsymq, irotmq, rtau, minus_q
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: npertx
    !! maximum dimension of the irreducible representations
    COMPLEX(DP), INTENT(OUT) :: t(npertx,npertx,48,3*nat)
    !! the symmetry matrices
    COMPLEX(DP), INTENT(OUT) :: tmq (npertx,npertx,3*nat)
    !! the matrices sending \(q \rightarrow -q+G\)
    !
    ! ... the local variables
    !
    INTEGER :: na, imode, jmode, ipert, jpert, kpert, nsymtot, imode0, &
               irr, ipol, jpol, isymq, irot, sna
    ! counters and auxiliary variables
    !
    REAL(DP) :: arg
    !! the argument of the phase
    !
    COMPLEX(DP) :: wrk_u (3, nat), wrk_ru (3, nat), fase, wrk
    ! pattern
    ! rotated pattern
    ! the phase factor
    !
    !   We compute the matrices which represent the symmetry transformation
    !   in the basis of the displacements
    !
    t(:,:,:,:) = (0.d0, 0.d0)
    tmq(:,:,:) = (0.d0, 0.d0)
    IF (minus_q) THEN
       nsymtot = nsymq + 1
    ELSE
       nsymtot = nsymq
    ENDIF
    !
    DO isymq = 1, nsymtot
       IF (isymq.le.nsymq) THEN
          irot = isymq
       ELSE
          irot = irotmq
       ENDIF
       imode0 = 0
       DO irr = 1, nirr
          DO ipert = 1, npert (irr)
             IF (modenum /= 0 .AND. modenum /= irr) CYCLE
             imode = imode0 + ipert
             DO na = 1, nat
                DO ipol = 1, 3
                   jmode = 3 * (na - 1) + ipol
                   wrk_u (ipol, na) = u (jmode, imode)
                ENDDO
             ENDDO
             !
             !     transform this pattern to crystal basis
             !
             DO na = 1, nat
                CALL trnvecc (wrk_u (1, na), at, bg, - 1)
             ENDDO
             !
             !     the patterns are rotated with this symmetry
             !
             wrk_ru(:,:) = (0.d0, 0.d0)
             DO na = 1, nat
                sna = irt (irot, na)
                arg = 0.d0
                DO ipol = 1, 3
                   arg = arg + xq (ipol) * rtau (ipol, irot, na)
                ENDDO
                arg = arg * tpi
                IF ((isymq.eq.nsymtot.and.minus_q).OR.(t_rev(irot)==1)) THEN
                   fase = CMPLX (cos (arg), sin (arg), KIND=DP )
                ELSE
                   fase = CMPLX (cos (arg),-sin (arg), KIND=DP )
                ENDIF
                DO ipol = 1, 3
                   DO jpol = 1, 3
                      wrk_ru (ipol, sna) = wrk_ru (ipol, sna) + s (jpol, ipol, irot) &
                             * wrk_u (jpol, na) * fase
                   ENDDO
                ENDDO
             ENDDO
             !
             !    Transform back the rotated pattern
             !
             DO na = 1, nat
                CALL trnvecc (wrk_ru (1, na), at, bg, 1)
             ENDDO      
             !
             !     Computes the symmetry matrices on the basis of the pattern
             !
             DO jpert = 1, npert (irr)
                imode = imode0 + jpert
                DO na = 1, nat
                   DO ipol = 1, 3
                      jmode = ipol + (na - 1) * 3
                      IF (isymq.eq.nsymtot.and.minus_q) THEN
                         tmq (jpert, ipert, irr) = tmq (jpert, ipert, irr) + CONJG(u ( &
                              jmode, imode) * wrk_ru (ipol, na) )
                      ELSE
                         IF (t_rev(irot)==1) THEN
                            t (jpert,ipert,irot,irr)=t(jpert,ipert,irot,irr) &
                               + CONJG(u (jmode, imode) * wrk_ru (ipol, na))
                         ELSE
                            t (jpert,ipert,irot,irr)=t(jpert,ipert,irot,irr) &
                               + CONJG(u (jmode, imode) ) * wrk_ru (ipol, na)
                         ENDIF
                      ENDIF
                   ENDDO
                ENDDO  
             ENDDO
          ENDDO
          imode0 = imode0 + npert (irr)
          !
          ! If the representations are irreducible, the rotations should be unitary matrices
          ! if this is not the case, the way the representations have been chosen has failed
          ! for some reasons (check set_irr.f90)
          !       
          IF (isymq<=nsymq) THEN
             DO ipert = 1, npert (irr)
                IF (modenum /= 0 .AND. modenum /= irr) CYCLE
                DO jpert = 1, npert (irr)
                   wrk = (0.d0,0.d0)
                   DO kpert = 1, npert (irr)
                      wrk = wrk + t (ipert,kpert,irot,irr) * conjg( t(jpert,kpert,irot,irr))
                   ENDDO
                   IF (jpert.ne.ipert .and. abs(wrk)> 1.d-6 ) &
                         CALL errore('set_irr_sym_new2','wrong representation',100*irr+10*jpert+ipert)
                   IF (jpert.eq.ipert .and. abs(wrk-1.d0)> 1.d-6 ) &
                         CALL errore('set_irr_sym_new2','wrong representation',100*irr+10*jpert+ipert)
                ENDDO
             ENDDO
          ENDIF
       ENDDO
    ENDDO
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE set_irr_sym_new2
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE sym_dns_wrapper2(indsym, ldim, dns_cart)
    !-----------------------------------------------------------------------
    !! This routine symmetrizes dns_cart. This is done in three steps.
    !
    !! Written by I. Timrov using the code by S. de Gironcoli 
    !! and A. Floris (01.10.2018).
    !  
    USE kinds,         ONLY : DP
    USE ep_constants,  ONLY : czero
    USE ions_base,     ONLY : nat
    USE modes,         ONLY : u, nmodes, nirr, npert
    USE lsda_mod,      ONLY : nspin
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: indsym(48)
    !! mapping between the original sym. indices and the new indices
    INTEGER, INTENT(IN) :: ldim
    !! matrix main dimension
    COMPLEX(DP), INTENT(INOUT) :: dns_cart(ldim,ldim,nspin,nat,3,nat)
    !! in/out: dns_cart is the dns matrix in the cartesian coordinates;
    !! on the input it is unsymmetrized, on the output it is symmetrized
    !COMPLEX(DP), INTENT(OUT) :: dns_pattern(ldim,ldim,nspin,nat,3*nat)
    !!! out: dns_pattern is the symmetrized dns matrix in the pattern basis
    !
    ! ... local variables
    !
    INTEGER :: imode, imode0, na, icart, na_icart, irr, npe
    COMPLEX(DP), ALLOCATABLE :: dns_aux(:,:,:,:,:)
    COMPLEX(DP), ALLOCATABLE :: dns_pattern(:,:,:,:,:)
    !
    ALLOCATE(dns_pattern(ldim,ldim,nspin,nat,3*nat))
    !
    dns_pattern = czero
    !
    ! 1 -  transform in the pattern basis
    !
    DO imode = 1, nmodes
       DO na = 1, nat
          DO icart = 1, 3
             na_icart = 3*(na-1) + icart
             dns_pattern(:,:,:,:,imode) = dns_pattern(:,:,:,:,imode) + &
                              dns_cart(:,:,:,:,icart,na) * u(na_icart,imode)
          ENDDO
       ENDDO
    ENDDO
    !
    ! 2 - symmetrize in the pattern basis
    !
    imode0 = 1
    DO irr = 1, nirr
       npe = npert(irr)
       ! allocate
       ALLOCATE (dns_aux(ldim,ldim,nspin,nat,npe))
       ! pack
       dns_aux(:,:,:,:,1:npe) = dns_pattern(:,:,:,:,imode0:imode0-1+npe)
       ! symmetrize
       CALL sym_dns2 (indsym, ldim, npe, irr, dns_aux)
       ! unpack
       dns_pattern(:,:,:,:,imode0:imode0-1+npe) = dns_aux(:,:,:,:,1:npe)
       ! deallocate
       DEALLOCATE (dns_aux)
       ! advance the counter
       imode0 = imode0 + npe
    ENDDO
    !
    ! 3 - back to the cartesian basis
    ! 
    dns_cart = (0.d0, 0.d0)
    !
    DO imode = 1, nmodes
       DO na = 1, nat
          DO icart = 1, 3
             na_icart = 3*(na-1) + icart
             dns_cart(:,:,:,:,icart,na) = dns_cart(:,:,:,:,icart,na) + &
                 dns_pattern(:,:,:,:,imode) * CONJG(u(na_icart,imode))
          ENDDO
       ENDDO
    ENDDO
    !
    RETURN
    !
    CONTAINS
    !-----------------------------------------------------------------------
    SUBROUTINE sym_dns2 (indsym, ldim, npe, irr, dns)
      !-----------------------------------------------------------------------
      !! DFPT+U: This routine symmetrizes the first order variation of 
      !! the occupation matrices dns due to the perturbation caused
      !! by the displacement of atoms.  
      !
      !! This subroutine is from /PHonon/PH/sym_dns.f90
      !
      USE kinds,            ONLY : DP
      USE ep_constants,     ONLY : czero
      USE constants,        ONLY : tpi
      USE ions_base,        ONLY : nat, ityp
      USE ldaU,             ONLY : Hubbard_l, is_hubbard, nwfcU, hubbard_occ
      USE lsda_mod,         ONLY : lsda, nspin
      USE lr_symm_base,     ONLY : nsymq, irgq, minus_q, irotmq, rtau
      USE modes,            ONLY : t, tmq
      USE qpoint,           ONLY : xq
      USE uspp_param,       ONLY : upf
      USE symm_base,        ONLY : d1, d2, d3, nsym, irt, s, invs
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: indsym(48)
      !! mapping between the original sym. indices and the new indices
      INTEGER, INTENT(IN) :: ldim, npe, irr
      COMPLEX(DP), INTENT(INOUT) :: dns(ldim,ldim,nspin,nat,npe)
      !
      ! ... local variables
      !
      INTEGER :: nt, n, counter, l, ip, jp, na, nb, is, m1, m2, &
                 m0, m00, isym, irot
      COMPLEX(DP), ALLOCATABLE :: dnr(:,:,:,:,:), dnraux(:,:,:,:,:)
      COMPLEX(DP) :: phase
      REAL(DP) :: arg
      !
      IF ((nsymq==1) .AND. (.NOT.minus_q)) RETURN
      !
      ! Initialization
      !
      ! D_Sl for l=1, l=2 and l=3 are already initialized, for l=0 D_S0 is 1
      !
      counter = 0  
      DO na = 1, nat  
         nt = ityp(na) 
         IF (.NOT.is_hubbard(nt)) CYCLE
         DO n = 1, upf(nt)%nwfc  
            l = upf(nt)%lchi(n)
            IF (hubbard_occ(nt,1) > 0.d0 .AND. l == Hubbard_l(nt)) &
               counter = counter + 2 * l + 1  
         ENDDO
      ENDDO
      IF (counter.NE.nwfcU) CALL errore ('sym_dns2', 'nwfcU<>counter', 1)
      !
      ALLOCATE (dnraux(ldim,ldim,nspin,nat,npe))
      ALLOCATE (dnr(ldim,ldim,nspin,nat,npe))
      !  
      dnraux = czero
      dnr    = czero
      !
      ! Impose hermiticity of dns_{m1,m2, is,na,ip} for m1<->m2
      ! and put it in dnr zeroing dns
      !
      DO ip = 1, npe
         DO na = 1, nat  
            nt = ityp(na)
            DO is = 1, nspin  
               DO m1 = 1, 2*Hubbard_l(nt)+1
                  DO m2 = 1, 2*Hubbard_l(nt)+1  
                     dnr(m1, m2, is, na,ip) = 0.5d0 * ( dns(m1, m2, is, na, ip) &
                                                      + dns(m2, m1, is, na, ip) )
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      !
      dns = czero
      !
      ! Symmetrize with -q if present (output overwritten on dnr)
      !
      IF (minus_q) THEN
         DO ip = 1, npe
            DO na = 1, nat  
               nt = ityp(na)  
               IF (is_hubbard(nt)) THEN 
                  ! 
                  arg = ( xq(1) * rtau(1,irotmq,na) + &
                          xq(2) * rtau(2,irotmq,na) + &
                          xq(3) * rtau(3,irotmq,na) ) * tpi
                  phase = CMPLX (COS(arg), SIN(arg), kind=DP)
                  !
                  DO is = 1, nspin  
                     DO m1 = 1, 2 * Hubbard_l(nt) + 1  
                        DO m2 = 1, 2 * Hubbard_l(nt) + 1  
                           nb = irt(irotmq, na)  
                           DO m0 = 1, 2 * Hubbard_l(nt) + 1  
                              DO m00 = 1, 2 * Hubbard_l(nt) + 1  
                                 DO jp = 1, npe
                                    IF (Hubbard_l(nt).EQ.0) THEN
                                       dns(m1,m2,is,na,ip) = dns(m1,m2,is,na,ip) + &
                                       dnr(m0,m00,is,nb,jp) * tmq(jp,ip,irr) *     &
                                       phase
                                    ELSE IF (Hubbard_l(nt).EQ.1) THEN
                                       dns(m1,m2,is,na,ip) = dns(m1,m2,is,na,ip) + &
                                       d1(m0 ,m1,irotmq) * d1(m00,m2,irotmq) *     &
                                       dnr(m0,m00,is,nb,jp) * tmq(jp,ip,irr) *     &
                                       phase
                                    ELSE IF (Hubbard_l(nt).EQ.2) THEN
                                       dns(m1,m2,is,na,ip) = dns(m1,m2,is,na,ip) + &
                                       d2(m0 ,m1,irotmq) * d2(m00,m2,irotmq) *     &
                                       dnr(m0,m00,is,nb,jp) * tmq(jp,ip,irr) *     &
                                       phase
                                    ELSE IF (Hubbard_l(nt).EQ.3) THEN
                                       dns(m1,m2,is,na,ip) = dns(m1,m2,is,na,ip) + &
                                       d3(m0 ,m1,irotmq) * d3(m00,m2,irotmq) *     &
                                       dnr(m0,m00,is,nb,jp) * tmq(jp,ip,irr) *     &
                                       phase
                                    ELSE
                                       CALL errore ('sym_dns2', &
                                            'angular momentum not implemented', &
                                            ABS(Hubbard_l(nt)) )
                                    ENDIF
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
         dnr(:,:,:,:,:) = 0.5d0 * ( dnr(:,:,:,:,:) + CONJG(dns(:,:,:,:,:)) )
         dns = czero
      ENDIF
      ! 
      dnraux = czero
      !
      ! Symmetryze dnr -> dns
      !
      DO isym = 1, nsymq
         !
         dns = czero
         !
         irot = isym
         !
         DO ip = 1, npe
            DO na = 1, nat  
               nt = ityp(na)  
               IF (is_hubbard(nt)) THEN  
                  ! 
                  arg = ( xq(1) * rtau(1,irot,na) + &
                          xq(2) * rtau(2,irot,na) + &
                          xq(3) * rtau(3,irot,na) ) * tpi
                  phase = CMPLX (COS(arg), SIN(arg), kind=DP )
                  !
                  DO m1 = 1, 2 * Hubbard_l(nt) + 1  
                     DO m2 = 1, 2 * Hubbard_l(nt) + 1  
                        nb = irt (irot, na)  
                        DO m0 = 1, 2 * Hubbard_l(nt) + 1  
                           DO m00 = 1, 2 * Hubbard_l(nt) + 1  
                              DO jp = 1, npe
                                 IF (Hubbard_l(nt).EQ.0) THEN
                                    dns(m1,m2,:,na,ip) = dns(m1,m2,:,na,ip) +  &
                                    dnr(m0,m00,:,nb,jp) * t(jp,ip,irot,irr) * &
                                    phase 
                                 ELSE IF (Hubbard_l(nt).EQ.1) THEN
                                    dns(m1,m2,:,na,ip) = dns(m1,m2,:,na,ip) + &
                                    d1(m0 ,m1,irot) * d1(m00,m2,irot) *       &
                                    dnr(m0,m00,:,nb,jp) * t(jp,ip,irot,irr) * &
                                    phase 
                                 ELSE IF (Hubbard_l(nt).EQ.2) THEN
                                    dns(m1,m2,:,na,ip) = dns(m1,m2,:,na,ip) + &
                                    d2(m0 ,m1,irot) * d2(m00,m2,irot) *       &
                                    dnr(m0,m00,:,nb,jp) * t(jp,ip,irot,irr) * &
                                    phase 
                                 ELSE IF (Hubbard_l(nt).EQ.3) THEN
                                    dns(m1,m2,:,na,ip) = dns(m1,m2,:,na,ip) + &
                                    d3(m0 ,m1,irot) * d3(m00,m2,irot) *       &
                                    dnr(m0,m00,:,nb,jp) * t(jp,ip,irot,irr) * &
                                    phase 
                                 ELSE
                                    CALL errore ('sym_dns2', &
                                         'angular momentum not implemented', &
                                         ABS(Hubbard_l(nt)) )
                                 ENDIF
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
         dnraux = dnraux + dns / nsymq
      ENDDO
      dns = dnraux
      ! 
      DEALLOCATE (dnr)
      DEALLOCATE (dnraux)
      ! 
      RETURN
      ! 
    END SUBROUTINE sym_dns2
    !---------------------------------------------------------------------
    !-----------------------------------------------------------------------
    END SUBROUTINE sym_dns_wrapper2
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE lr_orthoUwfc2 (xk, xkq, igk, igkq, npw, npwq, wfcatomk, wfcatomkpq, &
        swfcatomk, swfcatomkpq, lgamma, ik)
      !-----------------------------------------------------------------------
      !
      ! This routine computes atomic wavefunctions phi 
      ! at k and k+q. phi(k) is defined in Eq. (A6) in Ref. [1]
      ! [1] Phys. Rev. B 98, 085127 (2018)
      !
      ! This routine will write phi(k) and phi(k+q) to wfcatomk_, and
      ! wfcatomkpq_
      ! Also, write  S(k)*phi(k) and S(k+q)*phi(k+q) to sfwcatomk_, and
      ! swfcatomkpq_
      !
      ! In the norm-conserving case, S(k)=1 and S(k+q)=1.
      ! Note: here the array wfcU is used as a workspace.
      ! Inspired by PW/src/orthoatwfc.f90
      ! Adpated from LR/lr_orthoUwfc
      !
      USE kinds,            ONLY : DP
      USE ep_constants,     ONLY : czero
      USE basis,            ONLY : natomwfc
      USE wvfct,            ONLY : npwx
      USE control_flags,    ONLY : gamma_only
      USE uspp,             ONLY : vkb, nkb, okvan
      USE becmod,           ONLY : allocate_bec_type, deallocate_bec_type, becp
      USE noncollin_module, ONLY : noncolin, npol
      USE ldaU,             ONLY : Hubbard_projectors, wfcU, nwfcU, copy_U_wfc
      USE io_files,         ONLY : iunhub, iunhub_noS, nwordwfcU
      USE buffers,          ONLY : save_buffer
      ! 
      IMPLICIT NONE
      !
      REAL(KIND = DP), INTENT(in) :: xk(3)
      !! k-point coordinate
      REAL(KIND = DP), INTENT(in) :: xkq(3)
      !! k+q point coordinate 
      INTEGER, INTENT(in) :: igk(npw)
      !! k+G mapping
      INTEGER, INTENT(in) :: igkq(npwq)
      !! k+G+q mapping
      INTEGER, INTENT(IN) :: npw
      !! Number of k+G-vectors inside 'ecut sphere'
      INTEGER, INTENT(IN) :: npwq
      !! Number of k+G-vectors inside 'ecut sphere'
      COMPLEX(KIND = DP), INTENT(INOUT) :: wfcatomk(npwx * npol, nwfcU)
      !! \phi (atomic wavefuncton) at k
      COMPLEX(KIND = DP), INTENT(INOUT) :: swfcatomk(npwx * npol, nwfcU)
      !! S * \phi (atomic wavefuncton) at k
      COMPLEX(KIND = DP), INTENT(INOUT) :: wfcatomkpq(npwx * npol, nwfcU)
      !! \phi (atomic wavefuncton) at k+q
      COMPLEX(KIND = DP), INTENT(INOUT) :: swfcatomkpq(npwx * npol, nwfcU)
      !! S * \phi (atomic wavefuncton) at k+q
      LOGICAL, INTENT(IN) :: lgamma
      INTEGER, INTENT(IN), OPTIONAL :: ik
      !
      LOGICAL :: orthogonalize_wfc, normalize_only
      COMPLEX(DP), ALLOCATABLE :: wfcatom(:,:), &  ! atomic wfc
                                  swfcatom(:,:)    ! S * atomic wfc
      !
      IF (Hubbard_projectors=="atomic") THEN
         orthogonalize_wfc = .FALSE.
         normalize_only = .FALSE.
      ELSEIF (Hubbard_projectors=="ortho-atomic") THEN
         orthogonalize_wfc = .TRUE.
         normalize_only = .FALSE.
         IF (gamma_only) CALL errore('lr_orthoUwfc2', &
              'Gamma-only calculation for this case not implemented', 1 )
      ELSEIF (Hubbard_projectors=="norm-atomic") THEN
         orthogonalize_wfc = .TRUE.
         normalize_only = .TRUE.
         IF (gamma_only) CALL errore('lr_orthoUwfc2', &
              'Gamma-only calculation for this case not implemented', 1 )
      ELSE
         CALL errore ("lr_orthoUwfc2"," This Hubbard projectors type is not valid",1)
      ENDIF
      !
      ALLOCATE (wfcatom(npwx*npol,natomwfc))
      ALLOCATE (swfcatom(npwx*npol,natomwfc))
      !
      IF (okvan) CALL allocate_bec_type (nkb,natomwfc,becp)
      !
      wfcatomk(:, :)    = czero
      swfcatomk(:, :)   = czero
      wfcatomkpq(:, :)  = czero
      swfcatomkpq(:, :) = czero
      !
      wfcatom      = czero
      swfcatom     = czero
      !
      ! Determine the atomic orbital at k : phi(k)
      !
      CALL atomic_wfc2 (xk, igk, npw, wfcatom)
      !
      ! Compute S(k)*phi(k) (phi means the atomic orbital)
      !
      CALL s_phi (xk, igk, npw, wfcatom, swfcatom)     
      !
      ! Orthonormalize or normalize the atomic orbitals (if needed) 
      !
      IF (orthogonalize_wfc) THEN
         !$acc data copy(wfcatom, swfcatom)
         CALL ortho_swfc (npw, normalize_only, natomwfc, wfcatom, swfcatom, .TRUE.)
         !$acc end data
      END IF
      !
      ! Copy the result from (orthonormalized) wfcatom 
      ! (which uses the offset oatwfc) to wfcU (which uses the offset offsetU),
      ! and then write wfcU = phi(k) to wfcatomk_.
      !
      IF (ALLOCATED(wfcU)) THEN
        DEALLOCATE (wfcU) 
      ENDIF
      ALLOCATE(wfcU(npwx * npol, nwfcU))
      wfcU(:, :) = czero
      CALL copy_U_wfc (wfcatom, noncolin)
      IF (lgamma) CALL save_buffer (wfcU, nwordwfcU, iunhub_noS, ik)
      wfcatomk(:, :) = wfcU(:, :)
      !
      ! Copy the result from (orthonormalized) swfcatom 
      ! (which uses the offset oatwfc) to wfcU (which uses the offset offsetU),
      ! and then write wfcU = S(k)*phi(k) to swfcatomk_.
      !
      wfcU = czero
      CALL copy_U_wfc (swfcatom, noncolin)
      IF (lgamma) CALL save_buffer (wfcU, nwordwfcU, iunhub, ik)
      swfcatomk = wfcU
      !
      wfcatom  = czero
      swfcatom = czero
      !
      ! Determine the atomic orbital at k+q : phi(k+q)
      !
      CALL atomic_wfc2 (xkq, igkq, npwq, wfcatom)
      !
      ! Compute S(k+q)*phi(k+q) 
      !
      CALL s_phi (xkq, igkq, npwq, wfcatom, swfcatom)
      !
      ! Orthonormalize or normalize the atomic orbitals (if needed)
      !
      IF (orthogonalize_wfc) THEN 
         !$acc data copy(wfcatom, swfcatom)
         CALL ortho_swfc (npwq, normalize_only, natomwfc, wfcatom, swfcatom, .TRUE.)
         !$acc end data
      END IF
      !
      ! If lflag=.TRUE. copy the result from (orthonormalized) wfcatom 
      ! (which uses the offset oatwfc) to wfcU (which uses the offset offsetU),
      ! and then write wfcU = phi(k+q) to wfcatomkpq_.
      !
      wfcU = czero
      CALL copy_U_wfc (wfcatom, noncolin)
      wfcatomkpq = wfcU
      !
      ! Copy the result from (orthonormalized) swfcatom 
      ! (which uses the offset oatwfc) to wfcU (which uses the offset offsetU),
      ! and then write wfcU = S(k+q)*phi(k+q) to swfcatomkpq_.
      !
      wfcU = czero
      CALL copy_U_wfc (swfcatom, noncolin)
      swfcatomkpq = wfcU
      !
      DEALLOCATE (wfcatom)
      DEALLOCATE (swfcatom)
      DEALLOCATE (wfcU)
      ! 
      IF (okvan) CALL deallocate_bec_type (becp)
      !
      RETURN
      !
    CONTAINS
      !
    SUBROUTINE s_phi (xk_, igk_, npw_, wfc, swfc)
      !-----------------------------------------------------------------------
      !
      ! NCPP: swfc = wfc
      ! USPP: swfc = S * wfc
      !
      USE kinds,          ONLY : DP
      USE becmod,         ONLY : calbec, becp
      USE uspp_init,        ONLY : init_us_2
      !
      IMPLICIT NONE
      REAL(KIND = DP), INTENT(IN) :: xk_(3)
      !! k-point coordinate
      INTEGER, INTENT(IN) :: igk_(npw_)
      !! k+G mapping
      INTEGER, INTENT(IN) :: npw_
      !! Number of k+G-vectors inside 'ecut sphere'
      COMPLEX(DP), INTENT(IN) :: wfc( npwx*npol, natomwfc )
      !! Atomic wavefunctions
      COMPLEX(DP), INTENT(OUT) :: swfc( npwx*npol, natomwfc )
      !! S * atomic wavefunctions
      !
      ! NCPP case
      !
      IF ( nkb == 0 .OR. .NOT. okvan ) THEN
         swfc = wfc
         RETURN
      ENDIF
      !
      ! Electron-phonon matrix is not generally Hermitian unless S-matrix = 1.
      ! Therefore, USPP case is commented until problem with Hermicity is
      ! solved. (PRB 101, 184302 (2020))
      !
      !! USPP case
      !!
      !! Compute beta functions vkb at ik_
      !!
      !CALL init_us_2 (npw_, igk_, xk_, vkb, .true.)
      !!$acc update host(vkb)
      !!
      !! Compute the product of beta functions vkb
      !! with the functions wfc : becp = <vkb|wfc>
      !!
      !CALL calbec (npw_, vkb, wfc, becp)
      !!
      !! Calculate S*|wfc> = |wfc> + \sum qq * |vkb> * becp 
      !!  
      !CALL s_psi (npwx, npw_, natomwfc, wfc, swfc)
      !
      RETURN
      !
    END SUBROUTINE s_phi
    ! 
    !-----------------------------------------------------------------------
    SUBROUTINE atomic_wfc2(xk_, igk_, npw_, wfcatom_ )
      !-----------------------------------------------------------------------
      !! This routine computes the superposition of atomic wavefunctions
      !! for k-point "ik" - output in "wfcatom".
      !
      USE kinds,            ONLY : DP
      USE cell_base,        ONLY : omega, tpiba
      USE ions_base,        ONLY : nat, ntyp => nsp, ityp, tau
      USE gvect,            ONLY : mill, eigts1, eigts2, eigts3, g
      USE uspp_param,       ONLY : upf, nwfcm
      USE noncollin_module, ONLY : noncolin, npol 
      USE upf_spinorb,      ONLY : lmaxx
      USE pseudo_types,     ONLY : pseudo_upf ! type
      USE ep_constants,     ONLY : twopi ! parameter
      USE atwfc_mod,        ONLY : interp_atwfc
      !
      IMPLICIT NONE
      REAL(KIND = DP), INTENT(IN) :: xk_(3)
      !! k-point coordinate
      INTEGER, INTENT(in) :: igk_(npw_)
      !! k+G mapping
      INTEGER, INTENT(IN) :: npw_
      !! Number of k+G-vectors inside 'ecut sphere'
      COMPLEX(DP), INTENT(OUT) :: wfcatom_( npwx, npol, natomwfc )
      !! Superposition of atomic wavefunctions
      ! ... local variables
      !
      INTEGER :: n_starting_wfc, lmax_wfc, nt, l, nb, na, m, lm, ig, iig
      REAL(DP) :: xk1, xk2, xk3
      REAL(DP),    ALLOCATABLE :: qg(:), ylm (:,:), chiq (:,:,:), gk (:,:)
      COMPLEX(DP), ALLOCATABLE :: sk (:), aux(:)
      COMPLEX(DP) :: kphase, lphase
      REAL(DP)    :: arg
      !
      ! For noncollinear atomic wavefunctions
      !
      REAL(DP)    :: j
      REAL(DP), ALLOCATABLE :: chiaux(:)
      INTEGER     :: nc, ib
      ! 
      ! calculate max angular momentum required in wavefunctions
      lmax_wfc = 0
      DO nt = 1, ntyp
         lmax_wfc = MAX( lmax_wfc, MAXVAL( upf(nt)%lchi(1:upf(nt)%nwfc) ) )
      END DO
      !
      ALLOCATE( ylm(npw_,(lmax_wfc+1)**2), chiq(npw_,nwfcm,ntyp))
      ALLOCATE( qg(npw_), gk(3,npw_), sk(npw_) )
      ALLOCATE( aux(npw_) )
      !
      xk1 = xk_(1)
      xk2 = xk_(2)
      xk3 = xk_(3)
      !
      DO ig = 1, npw_
         iig = igk_(ig)
         gk (1,ig) = xk1 + g(1,iig)
         gk (2,ig) = xk2 + g(2,iig)
         gk (3,ig) = xk3 + g(3,iig)
         qg(ig) = gk(1,ig)**2 +  gk(2,ig)**2 + gk(3,ig)**2
      END DO
      !
      !  ylm = spherical harmonics
      !
      CALL ylmr2( (lmax_wfc+1)**2, npw_, gk, qg, ylm )
      !
      ! set now q=|k+G| in atomic units
      !
      DO ig = 1, npw_
         qg(ig) = SQRT( qg(ig) )*tpiba
      END DO
      !
      ! chiq = radial fourier transform of atomic orbitals chi
      !
      CALL interp_atwfc ( npw_, qg, nwfcm, chiq )
      !
      wfcatom_(:,:,:) = czero
      n_starting_wfc = 0
      !
      DO na = 1, nat
         arg = (xk1*tau(1,na) + xk2*tau(2,na) + xk3*tau(3,na)) * twopi
         kphase = CMPLX( COS(arg), - SIN(arg) ,KIND=DP)
         !
         !     sk is the structure factor
         !
         DO ig = 1, npw_
            iig = igk_(ig)
            sk (ig) = kphase * eigts1(mill(1,iig),na) * &
                               eigts2(mill(2,iig),na) * &
                               eigts3(mill(3,iig),na)
         END DO
         !
         nt = ityp(na)
         DO nb = 1, upf(nt)%nwfc
            IF ( upf(nt)%oc(nb) >= 0.d0 ) THEN
               l = upf(nt)%lchi(nb)
               lphase = (0.d0,1.d0)**l
               !
               !  the factor i^l MUST BE PRESENT in order to produce
               !  wavefunctions for k=0 that are real in real space
               !
               IF ( noncolin ) THEN
                  !
                  ! noncollinear case
                  IF ( upf(nt)%has_so ) j = upf(nt)%jchi(nb)
                  IF ( upf(nt)%has_so .AND. ABS(j-l+0.5_DP)<1.d-4) CYCLE
                  !
                  ALLOCATE( chiaux(npw_) )
                  !
                  IF ( upf(nt)%has_so) THEN
                     !
                     ! ... Find the index for j=l-1/2
                     !
                     IF (l == 0) THEN
                       chiaux(:) = chiq(:,nb,nt)
                     ELSE
                       DO ib = 1, upf(nt)%nwfc
                          IF ((upf(nt)%lchi(ib) == l) .AND. &
                                       (ABS(upf(nt)%jchi(ib)-l+0.5_DP)<1.d-4)) THEN
                            nc=ib
                            exit
                          ENDIF
                       ENDDO
                       !
                       ! ... Average the two radial functions
                       !
                       chiaux(:) = (chiq(:,nb,nt)*(l+1.0_DP)+chiq(:,nc,nt)*l)/(2.0_DP*l+1.0_DP)
                     ENDIF
                     !
                  ELSE
                     !
                     chiaux(:) = chiq(:,nb,nt)
                     !
                  ENDIF
                  !
                  DO m = 1, 2*l+1
                     lm = l**2 + m
                     n_starting_wfc = n_starting_wfc + 1

                     IF (n_starting_wfc + 2*l+1 > natomwfc) CALL errore &
                           ('atomic_wfc_nc', 'internal error: too many wfcs', 1)
                     DO ig = 1, npw_
                        aux(ig) = sk(ig)*ylm(ig,lm)*chiaux(ig)
                     ENDDO
                     ! 
                     DO ig = 1, npw_
                        !
                        wfcatom_(ig,1,n_starting_wfc) = aux(ig)
                        wfcatom_(ig,2,n_starting_wfc) = 0.d0
                        !
                        wfcatom_(ig,1,n_starting_wfc+2*l+1) = 0.d0
                        wfcatom_(ig,2,n_starting_wfc+2*l+1) = aux(ig)
                        !
                     ENDDO
                  ENDDO
                  !
                  n_starting_wfc = n_starting_wfc + 2*l+1
                  !
                  DEALLOCATE( chiaux )
                  !
               ELSE
                  !
                  ! ... LSDA or nonmagnetic case
                  DO m = 1, 2 * l + 1
                     lm = l**2 + m
                     n_starting_wfc = n_starting_wfc + 1
                     IF ( n_starting_wfc > natomwfc) CALL errore &
                        ('atomic_wfc___', 'internal error: too many wfcs', 1)
                     !
                     DO ig = 1, npw_
                        wfcatom_(ig, 1, n_starting_wfc) = lphase * &
                           sk (ig) * ylm (ig, lm) * chiq (ig, nb, nt)
                     ENDDO
                     !
                  END DO
                  !
               END IF
               !
            END IF
            !
         END DO
         !
      END DO
      !
      IF ( n_starting_wfc /= natomwfc) call errore ('atomic_wfc2', &
         'internal error: some wfcs were lost ', 1 )
      !
      DEALLOCATE( aux, sk, gk, qg, chiq, ylm )
      !
      RETURN
      ! 
    END SUBROUTINE atomic_wfc2
    !-----------------------------------------------------------------------
    END SUBROUTINE lr_orthoUwfc2
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE dvqhub_barepsi_us3(ik, uact, xxkq, sxk, igk, igkq, npw, npwq, evc, &
                                  wfcatomk, wfcatomkpq, swfcatomk, swfcatomkpq, nwfcU)
    !-----------------------------------------------------------------------
    !
    !! DFPT+U 
    !! This routines calculates the BARE derivative of the Hubbard potential
    !! times \(\text{psi}\).
    !
    ! is  = current_spin
    ! isi = opposite of the current_spin 
    !
    !! $$ |\Delta V_{BARE}(k+q,is) \psi(\text{ibnd},k,is)\rangle = 
    !!     + \sum_{I,m1,m2} \text{Hubbard_U}(I) \cdot [0.5\delta(m1,m2) - \text{ns}(m1,m2,is,I)]\cdot
    !!        [ |\text{dqsphi}(imode,I,k+q,m1)\rangle \langle S\phi(I,k,m2)|\psi(\text{ibnd},k,is)\rangle + 
    !!          |S\phi(I,k+q,m1)\rangle\langle\text{dmqsphi}(\text{imode},I,k,m2)|\psi(\text{ibnd},k,is)\rangle ]
    !!     - \sum_{I,m1,m2} \text{Hubbard_U}(I) \cdot \text{dnsbare}(m1,m2,is,I,\text{imode}) \cdot
    !!          |S\phi(I,k+q,m1)\rangle\langle S\phi(I,k,m2)|\psi(\text{ibnd},k,is)\rangle $$
    !
    !! Addition of the J0 terms:
    !
    !! $$ + \sum_{I,m1,m2} \text{Hubbard_J0}(I)\cdot \text{ns}(m1,m2,isi,I)\cdot
    !!       [ |\text{dqsphi}(\text{imode},I,k+q,m1)\rangle\langle S\phi(I,k,m2)|\psi(\text{ibnd},k,is)\rangle + 
    !!         |S\phi(I,k+q,m1)\langle\rangle\text{dmqsphi}(\text{imode},I,k,m2)|\psi(\text{ibnd},k,is)\rangle ]
    !!    + \sum_{I,m1,m2} \text{Hubbard_J0}(I) \cdot \text{dnsbare}(m1,m2,isi,I,\text{imode}) \cdot
    !!                          |S\phi(I,k+q,m1)\rangle\langle S\phi(I,k,m2)|\psi(\text{ibnd},k,is)\rangle $$
    !
    !! Important: in this routine \(\text{vkb}\) is a beta function at k+q, and \(\text{vkb_}\) is beta at k.
    !! This is done so because \(\text{vkb}\) is calculated at k+q in solve_linter (i.e. before calling
    !! this routine), so we keep the same attribution here.
    ! 
    !! Written by A. Floris.  
    !! Modified by I. Timrov (01.10.2018).
    !                 
    USE kinds,            ONLY : DP
    USE io_files,         ONLY : nwordwfcU
    USE ions_base,        ONLY : nat, ityp, ntyp => nsp
    USE ldaU,             ONLY : Hubbard_l, is_hubbard, Hubbard_U, Hubbard_J0, offsetU
    USE ldaU_ph,          ONLY : dnsbare
    USE ldaU_lr,          ONLY : effU
    USE wvfct,            ONLY : npwx
    USE noncollin_module, ONLY : npol
    USE uspp,             ONLY : vkb, nkb, ofsbeta
    USE uspp_param,       ONLY : nh
    USE lsda_mod,         ONLY : lsda, current_spin
    USE eqv,              ONLY : dvpsi
    USE scf,              ONLY : rho
    USE uspp_init,        ONLY : init_us_2
    USE global_var,       ONLY : lower_band, upper_band, nbndep
    USE ep_constants,     ONLY : czero, ryd2ev
    USE input,            ONLY : isk_loc
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ik
    !! Counter on k-point
    INTEGER, INTENT(in) :: npw
    !! Number of k+G-vectors inside 'ecut sphere'
    INTEGER, INTENT(in) :: npwq
    !! Number of k+G-vectors inside 'ecut sphere'
    INTEGER, INTENT(in) :: igk(npw)
    !! k+G mapping
    INTEGER, INTENT(in) :: igkq(npwq)
    !! k+G+q mapping
    REAL(KIND = DP), INTENT(in) :: xxkq(3)
    !! k+q point coordinate 
    REAL(KIND = DP), INTENT(in) :: sxk(3)
    !! Rotated k-point ik
    COMPLEX(KIND = DP), INTENT(in) :: uact(3 * nat)
    !! the pattern of displacements
    COMPLEX(KIND = DP), INTENT(in) :: evc(npwx * npol, nbndep)
    !! Wavefunction at k
    COMPLEX(KIND = DP), INTENT(in) :: wfcatomk(npwx * npol, nwfcU) ! npwx * npol, nwfcU
    !! Atomic wavefunction at k
    COMPLEX(KIND = DP), INTENT(in) :: wfcatomkpq(npwx * npol, nwfcU)
    !! Atomic wavefunction at k+q
    COMPLEX(KIND = DP), INTENT(in) :: swfcatomk(npwx * npol, nwfcU)
    !! S * atomic wavefunction at k
    COMPLEX(KIND = DP), INTENT(in) :: swfcatomkpq(npwx * npol, nwfcU)
    !! S * atomic wavefunction at k+q
    INTEGER, INTENT(in) :: nwfcU
    !! Number of Hubbard corrected atomic wavefunctions
    !
    ! local variables
    !
    INTEGER :: i, j, k, icart, counter, na, nt, l, ih, n, mu, ig, &
               ihubst, ihubst1, ihubst2, nah, m, m1, m2, ibnd, jbnd, op_spin, &
               ibeta
    COMPLEX(DP), ALLOCATABLE :: aux1(:), aux2(:), aux3(:), aux4(:), aux5(:), &
                                dqsphi(:,:), dmqsphi(:,:), dvqi(:,:), dvqhbar(:,:,:,:), &
                                vkbkpq(:,:), vkb_(:,:), dvkb(:,:,:), dvkbkpq(:,:,:), &
                                proj1(:,:), proj2(:,:), &
                                dwfcatom_(:), sdwfcatomk(:,:), sdwfcatomkpq(:,:), dwfcatomkpq(:,:)
    !  
    ALLOCATE (proj1(lower_band:upper_band,nwfcU))
    ALLOCATE (proj2(lower_band:upper_band,nwfcU))
    ALLOCATE (aux1(npwx))
    ALLOCATE (aux2(npwx))
    ALLOCATE (aux3(npwx))
    ALLOCATE (aux4(npwx)) 
    ALLOCATE (aux5(npwx))
    ALLOCATE (dqsphi(npwx,nwfcU))
    ALLOCATE (dmqsphi(npwx,nwfcU))
    ALLOCATE (dvqi(npwx,lower_band:upper_band))
    ALLOCATE (dvqhbar(npwx,lower_band:upper_band,3,nat))
    ALLOCATE (vkbkpq(npwx,nkb))
    ALLOCATE (vkb_(npwx,nkb))
    ALLOCATE (dvkb(npwx,nkb,3))
    ALLOCATE (dvkbkpq(npwx,nkb,3))
    ALLOCATE (dwfcatom_(npwx))
    ALLOCATE (sdwfcatomk(npwx,nwfcU))
    ALLOCATE (sdwfcatomkpq(npwx,nwfcU))
    ALLOCATE (dwfcatomkpq(npwx,nwfcU))
    !
    proj1(:, :)      = czero
    proj2(:, :)      = czero
    aux1(:)          = czero        
    aux2(:)          = czero        
    aux3(:)          = czero        
    aux4(:)          = czero        
    aux5(:)          = czero        
    vkbkpq(:,:)      = czero
    vkb_(:,:)        = czero
    dvkb(:, :, :)    = czero
    dvkbkpq(:, :, :) = czero
    dwfcatom_(:)     = czero
    sdwfcatomk(:,:)  = czero
    sdwfcatomkpq(:,:) = czero
    dwfcatomkpq(:,:) = czero 
    effU             = czero
    ! 
    DO na = 1, nat
      nt = ityp(na)
      !
      effU(nt) = Hubbard_U(nt)
      ! IF Hubbard_J0 /= 0
      IF (Hubbard_J0(nt) .NE. 0.d0) &
        effU(nt) = Hubbard_U(nt) - Hubbard_J0(nt)
    ENDDO
    !
    IF (lsda) THEN 
       current_spin = isk_loc(ik)
       IF (current_spin==1) THEN
          op_spin = 2
       ELSE
          op_spin = 1
       ENDIF
    ELSE        
       op_spin = 1
    ENDIF
    !
    ! Compute the beta function at k and put the result in vkb_
    !
    CALL init_us_2 (npw, igk, sxk, vkb_)
    !
    ! The beta function at k+q. Let us put it in the proper array vkbkpq,
    ! because in the following of this routine the array vkb will be 
    ! used as a workspace in the routine swfc. 
    !
    vkbkpq = vkb
    !
    ! Calculate the derivatives of beta functions 
    ! d^{icart}beta at k and k+q for all the bands and for 
    ! the 3 cartesian directions
    !
    DO icart = 1, 3
       DO na = 1, nat
          nt = ityp(na) 
          DO ih = 1, nh(nt)
             !
             ibeta = ofsbeta(na) + ih
             !
             CALL dwfc2 (npw, igk, sxk, icart, vkb_(:,ibeta), dvkb(:,ibeta,icart))
             CALL dwfc2 (npwq, igkq, xxkq, icart, vkbkpq(:,ibeta), dvkbkpq(:,ibeta,icart))
             !
          ENDDO
       ENDDO
    ENDDO
    !
    dvqhbar = czero
    !
    DO na = 1, nat
       !
       DO icart = 1, 3 
          !  
          dqsphi  = czero   
          dmqsphi = czero
          !
          DO nah = 1, nat
             !
             nt = ityp(nah)
             !
             ! For Hubbard_U - Hubbard_J0
             ! 
             IF (is_hubbard(nt)) THEN
                !
                DO m = 1, 2*Hubbard_l(nt)+1
                   !
                   ihubst = offsetU(nah) + m   ! I m index
                   !
                   IF (nah==na) THEN
                      !
                      ! Calculate | d_icart\phi_(k,I,m)) >
                      !
                      CALL dwfc2 (npw, igk, sxk, icart, wfcatomk(:,ihubst), dwfcatom_) 
                      !
                      ! Apply the S matrix: | S d_^(I,icart)\phi_(k,I,m) >
                      !
                      CALL swfc2 (npw, 1, vkb_, dwfcatom_, sdwfcatomk(:,ihubst))
                      !
                      ! Calculate |d_icart\phi_(k+q,I,m)) >
                      !
                      CALL dwfc2 (npwq, igkq, xxkq, icart, wfcatomkpq(:,ihubst), dwfcatom_)
                      ! 
                      ! Calculate | S d_^(I,icart)\phi_(k+q,I,m) >
                      !
                      CALL swfc2 (npwq, 1, vkbkpq, dwfcatom_, sdwfcatomkpq(:,ihubst))
                      ! 
                   ENDIF
                   !
                   ! Calculate |\Delta_q(S_k \phi_(k,I,m)) >  
                   ! and |\Delta_{-q}(S_{k+q} \phi_(k+q,I,m)) >
                   !
                   CALL delta_sphi2 (npw, npwq, na, icart, nah, ihubst, wfcatomk, wfcatomkpq,  &
                                    sdwfcatomk, sdwfcatomkpq, vkb_, vkbkpq, dvkb(:,:,icart), &
                                    dvkbkpq(:,:,icart), dqsphi, dmqsphi, 1, nwfcU)
                   !
                   ! Calculate:
                   ! proj1 (ihubst,ibnd) = < S_{k}\phi_(k,I,m)| psi(ibnd,k) >
                   ! proj2 (ihubst,ibnd) = < \Delta_{-q}(S_{k+q} \phi_(k+q,I,m)) | psi(ibnd,k) > 
                   !
                   DO ibnd = lower_band, upper_band
                      proj1(ibnd,ihubst) = DOT_PRODUCT(swfcatomk(1:npw,ihubst), evc(1:npw,ibnd))
                      proj2(ibnd,ihubst) = DOT_PRODUCT(dmqsphi(1:npw,ihubst), evc(1:npw,ibnd))
                   ENDDO
                   !
                ENDDO ! m
                !
             ENDIF
             ! 
          ENDDO ! nah
          ! 
          DO nah = 1, nat
             !
             nt = ityp(nah)
             !
             dvqi = czero 
             !
             IF (is_hubbard(nt)) THEN
                !
                DO ibnd = lower_band, upper_band
                   !
                   DO m1 = 1, 2*Hubbard_l(nt)+1
                      !
                      ihubst1 = offsetU(nah) + m1 
                      !
                      DO ig = 1, npwq
                         !
                         aux1(ig) = dqsphi(ig,ihubst1) * proj1(ibnd,ihubst1) 
                         !
                         aux3(ig) = swfcatomkpq(ig,ihubst1) * proj2(ibnd,ihubst1)  
                         !
                         dvqi(ig,ibnd) = dvqi(ig,ibnd) + 0.5d0 * (aux1(ig)+aux3(ig))
                         !
                      ENDDO
                      !
                      DO m2 = 1, 2*Hubbard_l(nt)+1
                         !
                         ihubst2 = offsetU(nah) + m2
                         !
                         DO ig = 1, npwq 
                            !                         
                            aux2(ig) = dqsphi(ig,ihubst1) * rho%ns(m1,m2,current_spin,nah) &
                                       * proj1(ibnd, ihubst2)
                            aux4(ig) = swfcatomkpq(ig,ihubst1) * rho%ns(m1,m2,current_spin,nah) &
                                       * proj2(ibnd, ihubst2)
                            aux5(ig) = swfcatomkpq(ig,ihubst1) &
                                       * dnsbare(m1,m2,current_spin,nah,icart,na) &
                                       * proj1(ibnd, ihubst2)
                            !
                            dvqi(ig, ibnd) = dvqi(ig, ibnd) - aux2(ig) - aux4(ig) - aux5(ig) 
                            ! 
                         ENDDO
                         ! 
                      ENDDO ! m2
                      !
                   ENDDO ! m1
                   !
                ENDDO ! ibnd
                !
                ! effU = Hubbard_U - Hubbard_J0
                !
                dvqi = dvqi * effU(nt)
                ! 
                DO ig = 1, npwq
                   dvqhbar(ig,:,icart,na) = dvqhbar(ig,:,icart,na) + dvqi(ig,:)
                ENDDO
                ! 
             ENDIF
             !
             ! Hubbard_J0 \= 0 case 
             !
             dvqi = (0.d0, 0.d0) 
             ! 
             IF (Hubbard_J0(nt).NE.0.d0) THEN
                !              
                DO ibnd = lower_band, upper_band
                   ! 
                   DO m1 = 1, 2*Hubbard_l(nt)+1
                      ! 
                      ihubst1 = offsetU(nah) + m1 
                      !
                      ! No diagonal term for J0
                      !                          
                      DO m2 = 1, 2*Hubbard_l(nt)+1
                         ! 
                         ihubst2 = offsetU(nah) + m2
                         ! 
                         DO ig = 1, npwq                          
                            aux2(ig) = dqsphi(ig, ihubst1) * rho%ns(m1,m2,op_spin,nah) &
                                       * proj1(ibnd, ihubst2)
                            aux4(ig) = swfcatomkpq(ig,ihubst1) * rho%ns(m1,m2,op_spin,nah) &
                                       * proj2(ibnd, ihubst2)
                            aux5(ig) = swfcatomkpq(ig,ihubst1) &
                                       * dnsbare (m1,m2,op_spin,nah,icart,na) &
                                       * proj1 (ibnd, ihubst2)
                            !
                            ! Note the sign change w.r.t. the case above
                            !
                            dvqi(ig, ibnd) = dvqi(ig, ibnd) + aux2(ig) + aux4(ig) + aux5(ig)
                            ! 
                         ENDDO
                         !
                      ENDDO ! m2
                      !
                   ENDDO ! m1
                   !
                ENDDO ! ibnd
                !
                dvqi = dvqi * Hubbard_J0(nt)
                !   
                DO ig = 1, npwq
                   dvqhbar(ig,:,icart,na) = dvqhbar(ig,:,icart,na) + dvqi(ig,:)
                ENDDO
                !
             ENDIF
             !
          ENDDO ! nah
          !
       ENDDO ! icart
       ! 
    ENDDO ! na
    ! 
    ! Compute the displacement along the pattern uact (for each band).
    ! The result is stored in aux1.
    !
    DO ibnd = lower_band, upper_band    
       !   
       aux1 = (0.d0, 0.d0)
       !
       DO na = 1, nat
          mu = 3 * (na - 1)
          ! Here is the basis transformation from cartesian to pattern uact
          DO icart = 1, 3 
             DO ig = 1, npwq
                aux1(ig) = aux1(ig) + dvqhbar(ig,ibnd,icart,na) * uact(mu+icart) 
             ENDDO
          ENDDO
       ENDDO
       !
       ! Add the result to dvpsi
       !
       DO ig = 1, npwq
          dvpsi(ig,ibnd) = dvpsi(ig,ibnd) + aux1(ig)
       ENDDO
       !
    ENDDO 
    !
    DEALLOCATE (proj1)
    DEALLOCATE (proj2)  
    DEALLOCATE (aux1) 
    DEALLOCATE (aux2)
    DEALLOCATE (aux3) 
    DEALLOCATE (aux4)
    DEALLOCATE (aux5)
    DEALLOCATE (dqsphi)
    DEALLOCATE (dmqsphi)
    DEALLOCATE (dvqi)
    DEALLOCATE (dvqhbar)
    DEALLOCATE (vkb_)
    DEALLOCATE (dwfcatom_) 
    DEALLOCATE (sdwfcatomk)
    DEALLOCATE (sdwfcatomkpq)
    DEALLOCATE (dwfcatomkpq)
    
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE dvqhub_barepsi_us3
    !-----------------------------------------------------------------------
    !
    !-------------------------------------------------------------------
    SUBROUTINE dwfc2 (npw_, igk_, xk_, icart_, func, dfunc)
      !-----------------------------------------------------------------
      !! This routine calculates the first derivative of a wave function 
      !! w.r.t. the position operator r_icart:
      !! $$ d\psi/d(\text{icart}) = \sum_G [i(k+G)_\text{icart}] \psi(G) $$
      !
      USE kinds,      ONLY : DP
      !USE wvfct,      ONLY : npwx
      USE gvect,      ONLY : g
      USE cell_base,  ONLY : tpiba
      USE ep_constants,  ONLY : czero, ci
    
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: npw_
      !! Number of k+G-vectors inside 'ecut sphere'
      INTEGER, INTENT(IN) :: igk_(npw_)
      !! k+G mapping
      INTEGER, INTENT(IN) :: icart_
      !! index of cartesian 
      REAL(KIND = DP), INTENT(IN) :: xk_(3)
      !! k-point coordinate 
      COMPLEX(KIND = DP), INTENT(IN) :: func(npw_)
      !! original 
      COMPLEX(KIND = DP), INTENT(OUT) :: dfunc(npw_)
      !
      ! ... local variables
      !
      INTEGER :: ig
      !
      CALL start_clock( 'dwfc' )
      !
      dfunc = czero
      !
      DO ig = 1, npw_
         dfunc(ig) = -ci * tpiba * (g(icart_, igk_(ig)) + xk_(icart_)) * func(ig)
      ENDDO
      !
      CALL stop_clock( 'dwfc' )
      !
      RETURN 
      ! 
    END SUBROUTINE dwfc2
    !------------------------------------------------------------------------
    !
    !----------------------------------------------------------------
    SUBROUTINE swfc2 (npw_, nbnd_, vkb_, wfc_, swfc_) 
      !---------------------------------------------------------------
      !! This routine applies the S operator to the function \(\text{wfc_}\)
      !! and puts the result in \(\text{swfc_}\), i.e.
      !! \(\text{swfc_} = S * \text{wfc_}\).
      !
      !! Important notice: here, the global array \(\text{vkb}\) is used as 
      !! a workspace, because the routine \(\texttt{s_psi}\) uses 
      !! \(\text{vkb}\) internally.  
      !! \(\text{vkb_}\) can be a beta function at k or k+q, therefore, 
      !! before changing the global array \(\text{vkb}\), we need to be very 
      !! careful: save \(\text{vkb}\) to a temporary array \(\text{vkb_save}\).
      !! Then we copy \(\text{vkb_save}\) back to \(\text{vkb}\), such that
      !! the meaning of \(\text{vkb}\) outside of this routine is restored
      !! (whatever it is).
      !
      !! Written by A. Floris.  
      !! Modified by I. Timrov (01.10.2018).
      !
      USE kinds,  ONLY : DP
      USE becmod, ONLY : calbec
      USE uspp,   ONLY : vkb, nkb
      USE wvfct,  ONLY : npwx
      
      IMPLICIT NONE
      
      INTEGER,     INTENT(IN)  :: npw_
      INTEGER,     INTENT(IN)  :: nbnd_
      COMPLEX(DP), INTENT(IN)  :: vkb_ (npwx, nkb)
      COMPLEX(DP), INTENT(IN)  :: wfc_ (npwx, nbnd_) 
      COMPLEX(DP), INTENT(OUT) :: swfc_(npwx, nbnd_) 
      !
      COMPLEX(DP), ALLOCATABLE :: vkb_save(:, :)
      COMPLEX(DP), ALLOCATABLE :: becp(:, :)
      INTEGER :: ierr
      !! Error status
      !
      CALL start_clock( 'swfc2' )
      !
      swfc_ = (0.d0, 0.d0)
      !
      ALLOCATE(becp(nkb, nbnd_), STAT = ierr)
      IF (ierr /= 0) CALL errore('swfc2', 'Error allocating becp', 1)
      ALLOCATE(vkb_save(npwx,nkb), STAT = ierr)
      IF (ierr /= 0) CALL errore('swfc2', 'Error allocating vkb_save', 1)
      !
      vkb_save = vkb
      vkb = vkb_
      ! 
      CALL calbec (npw_, vkb, wfc_, becp)
      !
      ! s_psi uses vkb
      !
      CALL s_psi (npwx, npw_, nbnd_, wfc_, swfc_)
      !
      vkb = vkb_save
      ! 
      DEALLOCATE(becp)
      IF (ierr /= 0) CALL errore('swfc2', 'Error deallocating becp', 1)
      DEALLOCATE(vkb_save)
      IF (ierr /= 0) CALL errore('swfc2', 'Error deallocating vkb_save', 1)
      !
      CALL stop_clock( 'swfc2' )
      !
      RETURN
      !
    END SUBROUTINE swfc2
    !
    !-----------------------------------------------------------------------------------
    SUBROUTINE delta_sphi2(npw, npwq, na, icart, nah, ihubst, wfcatomk_, wfcatomkpq_,   & 
                           sdwfcatomk_, sdwfcatomkpq_, vkb_, vkbkpq_, dvkb_, dvkbkpq_, &
                           dqsphi, dmqsphi, iflag, nwfcU)  
      !---------------------------------------------------------------------------------
      !
      !! DFPT+U: This routine calculates a vector at k
      !
      ! |\Delta_{-q}(S_{k+q} \phi_(k+q,I,m)) > = S_{k} | \d^{icart}phi_(k,nah,m) > + 
      !       \sum{l1,l2} [ | \dbeta^{icar_}_(k,na_,l1) > * qq_nt(na_, l1 ,l2) * 
      !                           < \beta_(k+q ,na_,l2) | phi_(k+q,nah,m)> + 
      !                      | \beta_(k,na_,l1)> * qq_nt(na_, l1 ,l2) * 
      !                           < \dbeta^{icar_}_(k+q,na_,l2) | phi_(k+q,nah,m) > ]
      !
      !! and also a vector at k+q.  
      !
      ! |\Delta_q(S_{k} \phi_(k,I,m)) > = S_{k+q}| \d^{icart}phi_(k+q,nah,m) > + 
      !       \sum{l1,l2} [ | \dbeta^{icar_}_(k+q,na_,l1) > * qq_nt(na_, l1 ,l2) * 
      !                           < \beta_(k ,na_,l2) | phi_(k,nah,m)> + 
      !                     | \beta_(k+q,na_,l1)> * qq_nt(na_, l1 ,l2) * 
      !                           < \dbeta^{icar_}_(k,na_,l2) | phi_(k ,nah,m) > ]
      !
      !  iflag = 1 : calculate |\Delta_q(S_{k} \phi_(k,I,m)) >  AND 
      !                        |\Delta_{-q}(S_{k+q} \phi_(k+q,I,m)) > 
      !  iflag = 0 : calculate ONLY |\Delta_{-q}(S_{k+q} \phi_(k+q,I,m)) > 
      !
      !! See source comment for details on the implemented formulas.
      !
      !! Written by A. Floris.  
      !! Modified by I. Timrov (01.10.2018).
      ! 
      USE kinds,      ONLY : DP
      USE uspp_param, ONLY : nh, nhm
      USE ions_base,  ONLY : nat, ityp
      USE uspp,       ONLY : nkb, qq_nt, okvan, ofsbeta
      USE wvfct,      ONLY : npwx
      !  
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: npw
      !! number of k+G-vectors inside 'ecut sphere'
      INTEGER, INTENT(IN) :: npwq
      !! number_of k+G-vectors inside 'ecut sphere'
      !
      INTEGER, INTENT(IN) :: na, icart, nah, ihubst 
      ! index of displaced atom
      ! index of cartesian direction of displacement
      ! nah identifies the Hubbard atom I 
      ! index of (I,m) atomic function to which the Dq is applied 
      !
      COMPLEX(DP), INTENT(IN) :: wfcatomk_(npwx,nwfcU),     & 
                                 sdwfcatomk_(npwx,nwfcU),   & 
                                 wfcatomkpq_(npwx,nwfcU),   & 
                                 sdwfcatomkpq_(npwx,nwfcU), &
                                 vkb_(npwx,nkb),            &
                                 vkbkpq_(npwx,nkb),         &
                                 dvkb_(npwx,nkb),           &
                                 dvkbkpq_(npwx,nkb)
      COMPLEX(DP), INTENT(INOUT) :: dqsphi(npwx,nwfcU), &
                                    dmqsphi(npwx,nwfcU)
      INTEGER, INTENT(IN) :: iflag
      INTEGER, INTENT(IN) :: nwfcU
      !
      ! Local variables
      !
      INTEGER :: nt, ih, m3, m4, ig, l
      COMPLEX(DP), ALLOCATABLE :: sc1(:), sc2(:), aux1(:), aux2(:)
      ! 
      CALL start_clock( 'delta_sphi' )
      !
      ALLOCATE (sc1(nhm))
      ALLOCATE (sc2(nhm))
      ALLOCATE (aux1(npwx))
      ALLOCATE (aux2(npwx))
      !
      nt = ityp(na)  
      !  
      ! Calculation of |\Delta_q(S_{k} \phi_(k,I,m)) > 
      ! 
      IF (iflag == 1) THEN
         ! 
         aux1 = (0.d0, 0.d0)
         aux2 = (0.d0, 0.d0)
         !
         ! USPP case 
         !
         IF ( okvan ) THEN
            !
            ! Scalar products in the m3 m4 sum  
            !
            DO ih = 1, nh(nt)
               sc1(ih) = DOT_PRODUCT(vkb_(1:npw,ih+ofsbeta(na)), wfcatomk_(1:npw,ihubst))
               sc2(ih) = DOT_PRODUCT(dvkb_(1:npw,ih+ofsbeta(na)), wfcatomk_(1:npw,ihubst))
            ENDDO
            !
         ENDIF
         ! 
         ! Add to Dq the term |S_{k+q} d_^(na,icart)\phi_(k+q,I,m) > * dkroneker Ina
         ! 
         IF (nah==na) THEN
            DO ig = 1, npwq
               dqsphi(ig,ihubst) = dqsphi(ig,ihubst) + sdwfcatomkpq_(ig,ihubst)
            ENDDO
         ENDIF
         !
         ! USPP case 
         !
         IF ( okvan ) THEN
            DO m3 = 1, nh(nt)
               DO m4 = 1, nh(nt)  
                  DO ig = 1, npwq
                     aux1(ig) = dvkbkpq_(ig,m3+ofsbeta(na)) * qq_nt(m3,m4,nt) * sc1(m4)
                     aux2(ig) =  vkbkpq_(ig,m3+ofsbeta(na)) * qq_nt(m3,m4,nt) * sc2(m4)
                     dqsphi(ig,ihubst) = dqsphi(ig,ihubst) + aux1(ig) + aux2(ig)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
         !
      ENDIF
      !
      ! Calculation of |\Delta_{-q}(S_{k+q} \phi_(k+q,I,m)) >
      !
      IF (iflag == 0 .OR. iflag == 1) THEN
         !  
         aux1 = (0.d0, 0.d0)
         aux2 = (0.d0, 0.d0)  
         !
         ! USPP case
         ! 
         IF ( okvan ) THEN  
            !
            ! Scalar products in the m3 m4 sum
            !
            DO ih = 1, nh(nt)
               sc1(ih) = DOT_PRODUCT(vkbkpq_(1:npwq,ih+ofsbeta(na)), wfcatomkpq_(1:npwq,ihubst))
               sc2(ih) = DOT_PRODUCT(dvkbkpq_(1:npwq,ih+ofsbeta(na)), wfcatomkpq_(1:npwq,ihubst))
            ENDDO
            !
         ENDIF
         !
         ! Add to D-q the term |S_{k} d_^(na,icart)\phi_(k,I,m) > * dkroneker Ina
         !
         IF (nah==na) THEN
            DO ig = 1, npw
               dmqsphi(ig,ihubst) = dmqsphi(ig,ihubst) + sdwfcatomk_(ig,ihubst)
            ENDDO
         ENDIF
         !
         ! USPP case
         !
         IF ( okvan ) THEN
            DO m3 = 1, nh(nt)
               DO m4 = 1, nh(nt)
                  DO ig = 1, npw  
                     aux1(ig) = dvkb_(ig,m3+ofsbeta(na)) * qq_nt(m3,m4,nt) * sc1(m4)
                     aux2(ig) =  vkb_(ig,m3+ofsbeta(na)) * qq_nt(m3,m4,nt) * sc2(m4)
                     dmqsphi(ig,ihubst) = dmqsphi(ig,ihubst) + aux1(ig) + aux2(ig)
                     !
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
         !
      ELSE 
         CALL errore ("delta_sphi"," wrong iflag", 1)
      ENDIF
      !  
      DEALLOCATE (sc1)
      DEALLOCATE (sc2)
      DEALLOCATE (aux1)
      DEALLOCATE (aux2)  
      !
      CALL stop_clock( 'delta_sphi' )
      !
      RETURN
      !
    END SUBROUTINE delta_sphi2
    !
    !--------------------------------------------------------------------------------------
    SUBROUTINE adddvhubscf2(ipert, ik, npw, npwq, evc, swfcatomk, swfcatomkpq, nwfcU)
    !--------------------------------------------------------------------------------------
    ! 
    !! DFPT+U  
    !! This routine calculates the SCF derivative of the Hubbard potential times \(\psi\):
    !! \begin{equation}\notag
    !! \begin{split}
    !!   |\Delta V_{SCF}(k+q,is) \psi(\text{ibnd},k,is)\rangle = 
    !!   - &\sum_{I,m1,m2} \text{Hubbard_U}(I)\cdot\text{dnsscf}(m1,m2,is,I,\text{imode})\cdot \\
    !!     &|S\phi(I,k+q,m1)\rangle\langle S\phi(I,k,m2|\psi(\text{ibnd},k,is)\rangle
    !! \end{split}
    !! \end{equation}
    !
    !! Addition of the \(\text{J0}\) terms:
    !
    !! $$ + \sum_{I,m1,m2} \text{Hubbard_J0}(I)\cdot \text{dnsscf}(m1,m2,\text{isi},I,\text{imode})\cdot 
    !! |S\phi(I,k+q,m1)\rangle\langle S\phi(I,k,m2)|\psi(\text{ibnd},k,is)\rangle $$
    !
    !! Where:  
    !! \(\text{is}\)  = current_spin;  
    !! \(\text{isi}\) = opposite of the current_spin.
    !
    !! Written by A. Floris. Modified by I. Timrov (01.10.2018).
    !
    USE kinds,            ONLY : DP
    USE io_files,         ONLY : nwordwfcU
    USE ions_base,        ONLY : nat, ityp, ntyp => nsp
    USE ldaU,             ONLY : Hubbard_lmax, Hubbard_l, Hubbard_U, offsetU, is_hubbard, &
                                 Hubbard_J0
    USE ldaU_lr,          ONLY : effU, dnsscf
    USE wvfct,            ONLY : npwx 
    USE noncollin_module, ONLY : npol
    USE lsda_mod,         ONLY : lsda, current_spin 
    USE eqv,              ONLY : dvpsi 
    USE global_var,       ONLY : lower_band, upper_band, nbndep
    USE input,            ONLY : isk_loc
    USE ep_constants,     ONLY : czero, ryd2ev
    !  
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: ik
    !! the k point under consideration
    INTEGER, INTENT(IN) :: ipert  
    !! the index of perturbation
    INTEGER, INTENT(in) :: npw
    !! Number of k+G-vectors inside 'ecut sphere'
    INTEGER, INTENT(in) :: npwq
    !! Number of k+G-vectors inside 'ecut sphere'
    COMPLEX(KIND = DP), INTENT(in) :: evc(npwx * npol, nbndep)
    !! Wavefunction at k
    COMPLEX(KIND = DP), INTENT(in) :: swfcatomk(npwx * npol, nwfcU)
    !! S * atomic wavefunction at k
    COMPLEX(KIND = DP), INTENT(in) :: swfcatomkpq(npwx * npol, nwfcU)
    !! S * atomic wavefunction at k+q
    INTEGER, INTENT(in) :: nwfcU
    !! Number of Hubbard corrected atomic wavefunctions
    !  
    ! Local variables
    !
    INTEGER ::  op_spin
    INTEGER ::  i, j, k, nt, l, ih, n, ig, ihubst, ihubst1, ihubst2, &
                nah, m, m1, m2, ibnd, jbnd, ldim
    COMPLEX(DP), ALLOCATABLE :: dvhubscf(:,:), dvqi(:,:), proj1(:, :)
    !
    CALL start_clock ( 'adddvhubscf2' )
    !
    ALLOCATE (proj1(lower_band:upper_band,nwfcU))
    ALLOCATE (dvhubscf(npwx,lower_band:upper_band))
    ALLOCATE (dvqi(npwx,lower_band:upper_band))
    !  
    ! At each iteration dvhubscf is set to zero, 
    ! because it is recomputed through dnsscf 
    !
    dvhubscf = czero
    effU     = czero
    !
    DO nah = 1, nat
      nt = ityp(nah)
      !
      effU(nt) = Hubbard_U(nt)
      ! IF Hubbard_J0 /= 0
      IF (Hubbard_J0(nt) .NE. 0.d0) &
        effU(nt) = Hubbard_U(nt) - Hubbard_J0(nt)
    ENDDO
    !
    IF (lsda) THEN
       current_spin = isk_loc(ik)
       IF (current_spin==1) THEN   
          op_spin = 2
       ELSE
          op_spin = 1
       END IF
    ELSE        
       op_spin = 1
    ENDIF
    !
    DO nah = 1, nat
       !
       nt = ityp(nah)
       !
       IF (is_hubbard(nt)) THEN        
          !   
          DO m = 1, 2*Hubbard_l(nt)+1
             !
             ihubst = offsetU(nah) + m   ! I m index
             !
             ! Calculate proj1(ibnd,ihubst) = < S_{k}\phi_(k,I,m)| psi(inbd,k) >
             !
             DO ibnd = lower_band, upper_band
                proj1(ibnd,ihubst) = DOT_PRODUCT(swfcatomk(1:npw,ihubst), evc(1:npw,ibnd))
             ENDDO
             !
          ENDDO
          !
       ENDIF
       !
    ENDDO
    ! 
    DO nah = 1, nat
       !
       nt = ityp(nah)
       !
       dvqi = (0.d0, 0.d0)
       !
       IF (is_hubbard(nt)) THEN
          !
          DO m1 = 1, 2*Hubbard_l(nt)+1
             !
             ihubst1 = offsetU(nah) + m1
             !
             DO m2 = 1, 2*Hubbard_l(nt)+1
                ! 
                ihubst2 = offsetU(nah) + m2
                ! 
                DO ibnd = lower_band, upper_band
                   DO ig = 1, npwq
                      dvqi(ig,ibnd) = dvqi(ig,ibnd) - effU(nt) * &
                                      dnsscf(m1,m2,current_spin,nah,ipert) * &
                                      swfcatomkpq(ig,ihubst1) * proj1(ibnd,ihubst2)
                   ENDDO
                ENDDO
                !
             ENDDO
             !
          ENDDO
          !
          DO ibnd = lower_band, upper_band
             DO ig = 1, npwq
                dvhubscf(ig,ibnd) = dvhubscf(ig,ibnd) + dvqi(ig,ibnd)
             ENDDO
          ENDDO
          !
       ENDIF
       !
       !  Hubbard_J0
       !
       dvqi = (0.d0, 0.d0)
       !
       IF (Hubbard_J0(nt).NE.0.d0) THEN
          !
          DO m1 = 1, 2*Hubbard_l(nt)+1
             !
             ihubst1 = offsetU(nah) + m1
             !
             DO m2 = 1, 2*Hubbard_l(nt)+1
                !
                ihubst2 = offsetU(nah) + m2
                ! 
                DO ibnd = lower_band, upper_band
                   DO ig = 1, npwq
                      ! Notice that sign change here
                      dvqi(ig,ibnd) = dvqi(ig,ibnd) + Hubbard_J0(nt) * &
                                      dnsscf(m1,m2,op_spin,nah,ipert) * &
                                      swfcatomkpq(ig,ihubst1) * proj1(ibnd,ihubst2)
                   ENDDO
                ENDDO
                !
             ENDDO
             !
          ENDDO
          !
          DO ibnd = lower_band, upper_band
             DO ig = 1, npwq
                dvhubscf(ig,ibnd) = dvhubscf(ig,ibnd) + dvqi(ig,ibnd)
             ENDDO
          ENDDO
          !
       ENDIF
       !
    ENDDO
    !  
    ! Add the result to dvpsi
    !
    DO ibnd = lower_band, upper_band
       DO ig = 1, npwq
          dvpsi(ig,ibnd) = dvpsi(ig,ibnd) + dvhubscf(ig,ibnd)  
       ENDDO
    ENDDO
    !
    DEALLOCATE (proj1)
    DEALLOCATE (dvhubscf)
    DEALLOCATE (dvqi)
    !
    CALL stop_clock ( 'adddvhubscf2' )
    !  
    RETURN
    ! 
    !--------------------------------------------------------------------------------------
    END SUBROUTINE adddvhubscf2
    !--------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------
  END MODULE dvqhub
  !-----------------------------------------------------------------------
