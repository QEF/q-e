!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE hp_symdnsq (dnsq)
  !-----------------------------------------------------------------------
  !
  ! This routine symmetrizes the first order variation of the occupation 
  ! matrices dnsq.
  !
  USE kinds,        ONLY : DP
  USE io_global,    ONLY : stdout
  USE constants,    ONLY : tpi
  USE ions_base,    ONLY : nat, ityp, tau
  USE lsda_mod,     ONLY : nspin
  USE uspp_param,   ONLY : upf
  USE symm_base,    ONLY : d1, d2, d3, nsym, irt, t_rev, s
  USE qpoint,       ONLY : xq
  USE lr_symm_base, ONLY : nsymq, minus_q, irotmq, rtau, gi
  USE ldaU,         ONLY : Hubbard_lmax, Hubbard_l, is_hubbard, d_spin_ldau
  USE ldaU_hp,      ONLY : nah_pert
  USE noncollin_module, ONLY: npol, noncolin, domag
  USE cell_base,        ONLY : bg, at

  IMPLICIT NONE
  ! 
  COMPLEX(DP), INTENT(INOUT) :: dnsq(2*Hubbard_lmax+1, 2*Hubbard_lmax+1, nspin, nat)
  !
  ! Local variables
  !
  INTEGER :: nt, n, l
  INTEGER :: na, nb, is, is2, i, j, m1, m2, m0, m00, ldim, is1, is3, is4, m3, m4
  INTEGER :: isym, irot, isym2
  COMPLEX(DP), ALLOCATABLE :: dnr(:,:,:,:), dnraux(:,:,:,:), &
                              dnr_nc(:,:,:,:,:), dnr1_nc(:,:,:,:,:), dnraux_nc(:,:,:,:,:)
  COMPLEX(DP) :: phase, phase2(0:1)
  REAL(DP) :: arg, arg2, gi_t(3, 48), aq (3), raq (3), wrk (3)
  !
  IF (nsymq == 1 .AND. (.NOT.minus_q)) RETURN
  !
  CALL start_clock('hp_symdnsq')
  !
  ldim = 2 * Hubbard_lmax + 1
  !
  ! Initialization
  IF (noncolin.and.domag) THEN
     IF ( .NOT. ALLOCATED (d_spin_ldau) ) ALLOCATE( d_spin_ldau(2,2,48) )
     CALL comp_dspinldau()
  ENDIF
  !
  ! D_Sl for l=1, l=2 and l=3 are already initialized, for l=0 D_S0 is 1
  !
  ! Allocate auxiliary arrays
  !
  ALLOCATE( dnraux(ldim, ldim, nspin, nat) )  
  dnraux(:,:,:,:) = (0.d0,0.d0)
  !
  IF (noncolin.and.domag) then
    ALLOCATE( dnr_nc(ldim, ldim, npol, npol, nat) )  
    dnr_nc(:,:,:,:,:) = (0.d0,0.d0)  
    ALLOCATE( dnr1_nc(ldim, ldim, npol, npol, nat) )  
    dnr1_nc(:,:,:,:,:) = (0.d0,0.d0) 
    ALLOCATE( dnraux_nc(ldim, ldim, npol, npol, nat) )  
    dnraux_nc(:,:,:,:,:) = (0.d0,0.d0) 
  ELSE
    ALLOCATE( dnr(ldim, ldim, nspin, nat) )  
    dnr(:,:,:,:) = (0.d0,0.d0)
  ENDIF
  !
  ! Impose hermiticity of dnsq_{m1, m2, is, na} for m1<->m2
  ! and put it in dnr zeroing dnsq
  ! IT: Hermiticity is already imposed by construction
  !
  IF (noncolin.and.domag) then
    DO na = 1, nat
       nt = ityp (na)
       IF ( is_hubbard(nt) ) THEN
          DO is = 1, npol
            DO is2 = 1, npol
              i = npol*(is-1) + is2
              j = is + npol*(is2-1)
              ldim = 2*Hubbard_l(nt)+1
              DO m1 = 1, ldim
                DO m2 = 1, ldim
                  !
                  ! Enforce hermiticity: if time-reversal symmetry is broken 
                  ! the occupancy matrix is hermitian (complex), not symmetric (real)  
                   dnr_nc(m1,m2,is,is2,na) =  0.5d0 * ( dnsq(m1, m2, i, na) + CONJG(dnsq(m2, m1, j, na)) )
                ENDDO
              ENDDO
            ENDDO
          ENDDO
       ENDIF
    ENDDO
  ELSE
     DO na = 1, nat  
        nt = ityp(na)
        DO is = 1, nspin  
           DO m1 = 1, 2 * Hubbard_l(nt) + 1
              DO m2 = 1, 2 * Hubbard_l(nt) + 1  
                 dnr(m1, m2, is, na) = 0.5d0 * ( dnsq(m1, m2, is, na) + dnsq(m2, m1, is, na) )
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDIF
  !
  ! Symmetrize with -q if present (output overwritten on dnr)
  !
  IF (minus_q) THEN
     !
     dnsq(:,:,:,:) = (0.d0,0.d0)
     !
     DO na = 1, nat  
        !
        nt = ityp(na)  
        !
        nb = irt(irotmq, na)
        !
        IF (is_hubbard(nt)) THEN 
           !
           ! Compute the phase factor exp(iq*f)
           !
           arg = ( xq(1) * rtau(1, irotmq, na) + &  
                   xq(2) * rtau(2, irotmq, na) + &
                   xq(3) * rtau(3, irotmq, na) ) * tpi
           phase = CMPLX(cos(arg), sin(arg), kind=DP)
           !
           DO is = 1, nspin  
              !
              DO m1 = 1, 2 * Hubbard_l(nt) + 1  
                 DO m2 = 1, 2 * Hubbard_l(nt) + 1  
                    DO m0 = 1, 2 * Hubbard_l(nt) + 1  
                       DO m00 = 1, 2 * Hubbard_l(nt) + 1  
                          !
                          IF (Hubbard_l(nt).EQ.0) THEN
                             dnsq(m1,m2,is,na) = dnsq(m1,m2,is,na) +   &   
                                      dnr(m0,m00,is,nb) * phase
                          ELSEIF (Hubbard_l(nt).EQ.1) THEN
                             dnsq(m1,m2,is,na) = dnsq(m1,m2,is,na)  +          &
                                      d1(m0,m1,irotmq) * d1(m00,m2,irotmq) * &
                                      dnr(m0,m00,is,nb) * phase
                          ELSEIF (Hubbard_l(nt).EQ.2) THEN
                             dnsq(m1,m2,is,na) = dnsq(m1,m2,is,na)  +          &
                                      d2(m0,m1,irotmq) * d2(m00,m2,irotmq) * &
                                      dnr(m0,m00,is,nb) * phase
                          ELSEIF (Hubbard_l(nt).EQ.3) THEN
                             dnsq(m1,m2,is,na) = dnsq(m1,m2,is,na)  +          &
                                      d3(m0,m1,irotmq) * d3(m00,m2,irotmq) * &
                                      dnr(m0,m00,is,nb) * phase
                          ELSE
                             CALL errore ('hp_symdnsq', 'angular momentum not implemented', &
                                         ABS(Hubbard_l(nt)) )
                          ENDIF
                          !
                       ENDDO ! m00
                    ENDDO ! m0
                 ENDDO ! m2
              ENDDO ! m1
              !
           ENDDO ! is
           !
        ENDIF 
        !
     ENDDO ! na
     !
     dnr(:,:,:,:) = 0.5d0 * ( dnr(:,:,:,:) + CONJG(dnsq(:,:,:,:)) )
     !
  ENDIF
  !
  ! Symmetryze dnr -> dnsq
  !
  DO isym = 1, nsymq
    !
    dnsq(:,:,:,:) = (0.d0,0.d0)
    if (noncolin.and.domag) dnr1_nc(:,:,:,:,:) = (0.d0,0.d0)
    !
    irot = isym
    !
    DO na = 1, nat
       !
       nt = ityp(na)  
       !
       nb = irt(irot, na)
       !
       IF (is_hubbard(nt)) THEN  
          !
          ! Compute the phase factor exp(iq*f)
          !
          arg = ( xq(1) * rtau(1, irot, na) + &
                  xq(2) * rtau(2, irot, na) + &
                  xq(3) * rtau(3, irot, na) ) * tpi
          phase = CMPLX (cos(arg), sin(arg), kind=DP)
          !
          ! Compute the phase factor exp(-iG*tau_pert)
          ! where G = Sq - q, and tau_pert is the position
          ! of the perturbed atom 
          !
          IF (noncolin.and.domag) THEN
             gi_t = 0.d0 
             aq   = xq
             call cryst_to_cart (1, aq, at, - 1)
             do isym2 = 1, nsymq
                raq = 0.d0
                do i = 1, 3
                   do j = 1, 3
                      raq (i) = raq (i) + DBLE (s (i, j, isym2) ) * &
                         aq (j)
                   enddo
                enddo
                do i = 1, 3
                   IF (t_rev(isym2)==1) wrk (i) = - aq (i) + raq (i)
                enddo
                call cryst_to_cart (1, wrk, bg, 1)
                gi_t (:, isym2) = wrk (:)
             enddo
          ENDIF
          !
          phase2 = 0.0
          IF (t_rev(irot) == 1) then
            arg2 = ( gi_t(1, irot) * tau(1, nah_pert) + &
                     gi_t(2, irot) * tau(2, nah_pert) + &
                     gi_t(3, irot) * tau(3, nah_pert) ) * tpi
            phase2(1) = CMPLX (cos(arg2), -sin(arg2), kind=DP)
          ELSE
            arg2 = ( gi(1, irot) * tau(1, nah_pert) + &
                     gi(2, irot) * tau(2, nah_pert) + &
                     gi(3, irot) * tau(3, nah_pert) ) * tpi
            phase2(0) = CMPLX (cos(arg2), -sin(arg2), kind=DP)
          ENDIF
          !
          IF (noncolin.and.domag) THEN
            DO m1 = 1, 2 * Hubbard_l(nt) + 1  
               DO m2 = 1, 2 * Hubbard_l(nt) + 1
                  DO is1 = 1, npol
                     DO is2 = 1, npol
                        DO m3 = 1, 2 * Hubbard_l(nt) + 1  
                           DO m4 = 1, 2 * Hubbard_l(nt) + 1
                              DO is3 = 1, npol
                                 DO is4 = 1, npol
                                    IF (Hubbard_l(nt) == 0) THEN
                                       IF (t_rev(irot) == 1) THEN
                                         dnr1_nc(m1,m2,is1,is2,na) = dnr1_nc(m1,m2,is1,is2,na) +  &
                                           CONJG( d_spin_ldau(is1,is3,irot) ) *                 &
                                                  dnr_nc(m4,m3,is4,is3,nb)*(phase * phase2(1)) *    &
                                                  d_spin_ldau(is2,is4,irot)  
                                       ELSE
                                          dnr1_nc(m1,m2,is1,is2,na) = dnr1_nc(m1,m2,is1,is2,na) +  &
                                           CONJG( d_spin_ldau(is1,is3,irot) ) *                 &
                                                  dnr_nc(m3,m4,is3,is4,nb) * phase * phase2(0) *    &
                                                  d_spin_ldau(is2,is4,irot)  
                                       ENDIF
                                    ELSEIF (Hubbard_l(nt) == 1) THEN
                                       IF (t_rev(irot) == 1) THEN
                                          dnr1_nc(m1,m2,is1,is2,na) = dnr1_nc(m1,m2,is1,is2,na) + &
                                            CONJG( d_spin_ldau(is1,is3,irot) )*d1(m1,m3,irot)*  &
                                                   dnr_nc(m4,m3,is4,is3,nb)*(phase * phase2(1)) *    &
                                                   d_spin_ldau(is2,is4,irot)  *d1(m2,m4,irot)
                                        ELSE
                                          dnr1_nc(m1,m2,is1,is2,na) = dnr1_nc(m1,m2,is1,is2,na) + &
                                            CONJG( d_spin_ldau(is1,is3,irot) )*d1(m1,m3,irot)*  &
                                                   dnr_nc(m3,m4,is3,is4,nb)* phase * phase2(0) *    &
                                                   d_spin_ldau(is2,is4,irot)  *d1(m2,m4,irot)
                                        ENDIF
                                      ELSEIF (Hubbard_l(nt) == 2) THEN
                                        IF (t_rev(irot) == 1) THEN
                                          dnr1_nc(m1,m2,is1,is2,na) = dnr1_nc(m1,m2,is1,is2,na) + &
                                            CONJG( d_spin_ldau(is1,is3,irot) )*d2(m1,m3,irot)*  &
                                                   dnr_nc(m4,m3,is4,is3,nb)*(phase * phase2(1)) *    &
                                                   d_spin_ldau(is2,is4,irot)  *d2(m2,m4,irot) 
                                        ELSE
                                          dnr1_nc(m1,m2,is1,is2,na) = dnr1_nc(m1,m2,is1,is2,na) + &
                                            CONJG( d_spin_ldau(is1,is3,irot) )*d2(m1,m3,irot)*  &
                                                   dnr_nc(m3,m4,is3,is4,nb)* phase * phase2(0) *    &
                                                   d_spin_ldau(is2,is4,irot)  *d2(m2,m4,irot)
                                        ENDIF
                                      ELSEIF (Hubbard_l(nt) == 3) THEN
                                        !
                                        IF (t_rev(irot) == 1) THEN
                                          dnr1_nc(m1,m2,is1,is2,na) = dnr1_nc(m1,m2,is1,is2,na) + &
                                            CONJG( d_spin_ldau(is1,is3,irot) )*d3(m1,m3,irot)*  &
                                                   dnr_nc(m4,m3,is4,is3,nb)*(phase * phase2(1)) *  &
                                                   d_spin_ldau(is2,is4,irot)  *d3(m2,m4,irot)
                                        ELSE
                                          dnr1_nc(m1,m2,is1,is2,na) = dnr1_nc(m1,m2,is1,is2,na) + &
                                            CONJG( d_spin_ldau(is1,is3,irot) )*d3(m1,m3,irot)*  &
                                                   dnr_nc(m3,m4,is3,is4,nb)* phase * phase2(0) *    &
                                                   d_spin_ldau(is2,is4,irot)  *d3(m2,m4,irot)
                                        ENDIF
                                        !
                                      ELSE
                                        !
                                        CALL errore( 'hp_symdnsq', &
                                                     'angular momentum not implemented', &
                                                     ABS(Hubbard_l(nt)) )
                                        !
                                      ENDIF
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
          ELSE
            DO is = 1, nspin
               !
               DO m1 = 1, 2 * Hubbard_l(nt) + 1  
                  DO m2 = 1, 2 * Hubbard_l(nt) + 1  
                     DO m0 = 1, 2 * Hubbard_l(nt) + 1  
                        DO m00 = 1, 2 * Hubbard_l(nt) + 1  
                           !
                           IF (Hubbard_l(nt).EQ.0) THEN
                              dnsq(m1,m2,is,na) = dnsq(m1,m2,is,na) + &
                                       dnr(m0,m00,is,nb) * phase * phase2(0)
                           ELSEIF (Hubbard_l(nt).EQ.1) THEN
                              dnsq(m1,m2,is,na) = dnsq(m1,m2,is,na) +       &
                                       d1(m0,m1,irot) * d1(m00,m2,irot) * &
                                       dnr(m0,m00,is,nb) * phase * phase2(0)
                           ELSEIF (Hubbard_l(nt).EQ.2) THEN
                              dnsq(m1,m2,is,na) = dnsq(m1,m2,is,na) +       &
                                       d2(m0,m1,irot) * d2(m00,m2,irot) * &
                                       dnr(m0,m00,is,nb) * phase * phase2(0)
                           ELSEIF (Hubbard_l(nt).EQ.3) THEN
                              dnsq(m1,m2,is,na) = dnsq(m1,m2,is,na) +       &
                                       d3(m0,m1,irot) * d3(m00,m2,irot) * &
                                       dnr(m0,m00,is,nb) * phase * phase2(0)
                           ELSE
                              CALL errore ('hp_symdnsq', 'angular momentum not implemented', &
                                        ABS(Hubbard_l(nt)) )
                           ENDIF
                           !
                        ENDDO ! m00
                     ENDDO ! m0
                  ENDDO ! m2
               ENDDO ! m1
               !
            ENDDO ! is
          ENDIF
          !
       ENDIF 
       !
    ENDDO ! na
    !
    IF (noncolin.and.domag) then
       dnraux_nc = dnraux_nc + dnr1_nc / nsymq
    ELSE
       dnraux = dnraux + dnsq / nsymq
    ENDIF
    !
  ENDDO ! isym
  !
  ! Setup the output matrix ns with combined spin index 
  !
  IF (noncolin.and.domag) THEN
      DO na = 1, nat
         nt = ityp (na)
         IF ( is_hubbard(nt) ) THEN
            DO is1 = 1, npol
               DO is2 = 1, npol
                  i = npol*(is1-1) + is2
                  DO m1 = 1, 2*Hubbard_l(nt)+1
                     DO m2 = 1, 2*Hubbard_l(nt)+1
                        dnsq(m1,m2,i,na) = dnraux_nc(m1,m2,is1,is2,na)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ENDDO
  ELSE
      dnsq = dnraux
  ENDIF
  !
  IF (noncolin.and.domag) THEN
     DEALLOCATE (dnr_nc)
     DEALLOCATE (dnr1_nc)
     DEALLOCATE (dnraux_nc)
     IF ( ALLOCATED (d_spin_ldau) ) DEALLOCATE( d_spin_ldau )
  ELSE
     DEALLOCATE (dnr)
  ENDIF
  !
  DEALLOCATE (dnraux)
  !
  CALL stop_clock('hp_symdnsq')
  !
  RETURN
  !
END SUBROUTINE hp_symdnsq
