!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
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
  USE basis,        ONLY : natomwfc
  USE lsda_mod,     ONLY : nspin
  USE uspp_param,   ONLY : upf
  USE symm_base,    ONLY : d1, d2, d3, nsym, irt
  USE qpoint,       ONLY : xq
  USE lr_symm_base, ONLY : nsymq, minus_q, irotmq, rtau, gi
  USE ldaU,         ONLY : Hubbard_lmax, Hubbard_l, is_hubbard
  USE ldaU_hp,      ONLY : nah_pert

  IMPLICIT NONE
  ! 
  COMPLEX(DP), INTENT(INOUT) :: dnsq(2*Hubbard_lmax+1, 2*Hubbard_lmax+1, nspin, nat)
  !
  ! Local variables
  !
  INTEGER :: nt, n, counter, l
  INTEGER :: na, nb, is, m1, m2, m0, m00, ldim
  INTEGER :: isym, irot
  COMPLEX(DP), ALLOCATABLE :: dnr(:,:,:,:), dnraux(:,:,:,:)
  COMPLEX(DP) :: phase, phase2
  REAL(DP) :: arg, arg2
  !
  IF (nsymq == 1 .AND. (.NOT.minus_q)) RETURN
  !
  CALL start_clock('hp_symdnsq')
  !
  ldim = 2 * Hubbard_lmax + 1
  !
  ! Initialization
  !
  ! D_Sl for l=1, l=2 and l=3 are already initialized, for l=0 D_S0 is 1
  !
  counter = 0  
  DO na = 1, nat  
     nt = ityp (na)  
     DO n = 1, upf(nt)%nwfc  
        IF (upf(nt)%oc(n) >= 0.d0) THEN  
           l = upf(nt)%lchi(n)   
           counter = counter + 2 * l + 1  
        ENDIF
     ENDDO
  ENDDO
  IF (counter.NE.natomwfc) CALL errore ('hp_symdnsq', 'natomwfc<>counter', 1)
  !
  ! Allocate auxiliary arrays
  !
  ALLOCATE( dnraux(ldim, ldim, nspin, nat) )  
  dnraux(:,:,:,:) = (0.d0,0.d0)
  !
  ALLOCATE( dnr(ldim, ldim, nspin, nat) )  
  dnr(:,:,:,:) = (0.d0,0.d0)
  !
  ! Impose hermiticity of dnsq_{m1, m2, is, na} for m1<->m2
  ! and put it in dnr zeroing dnsq
  ! IT: Hermiticity is already imposed by construction
  !
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
          arg2 = ( gi(1, irot) * tau(1, nah_pert) + &
                   gi(2, irot) * tau(2, nah_pert) + &
                   gi(3, irot) * tau(3, nah_pert) ) * tpi
          phase2 = CMPLX (cos(arg2), -sin(arg2), kind=DP)
          !
          DO is = 1, nspin
             !
             DO m1 = 1, 2 * Hubbard_l(nt) + 1  
                DO m2 = 1, 2 * Hubbard_l(nt) + 1  
                   DO m0 = 1, 2 * Hubbard_l(nt) + 1  
                      DO m00 = 1, 2 * Hubbard_l(nt) + 1  
                         !
                         IF (Hubbard_l(nt).EQ.0) THEN
                            dnsq(m1,m2,is,na) = dnsq(m1,m2,is,na) + &
                                     dnr(m0,m00,is,nb) * phase * phase2
                         ELSEIF (Hubbard_l(nt).EQ.1) THEN
                            dnsq(m1,m2,is,na) = dnsq(m1,m2,is,na) +       &
                                     d1(m0,m1,irot) * d1(m00,m2,irot) * &
                                     dnr(m0,m00,is,nb) * phase * phase2
                         ELSEIF (Hubbard_l(nt).EQ.2) THEN
                            dnsq(m1,m2,is,na) = dnsq(m1,m2,is,na) +       &
                                     d2(m0,m1,irot) * d2(m00,m2,irot) * &
                                     dnr(m0,m00,is,nb) * phase * phase2
                         ELSEIF (Hubbard_l(nt).EQ.3) THEN
                            dnsq(m1,m2,is,na) = dnsq(m1,m2,is,na) +       &
                                     d3(m0,m1,irot) * d3(m00,m2,irot) * &
                                     dnr(m0,m00,is,nb) * phase * phase2
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
    dnraux = dnraux + dnsq / nsymq
    !
  ENDDO ! isym
  !
  dnsq = dnraux
  !
  DEALLOCATE (dnr)
  DEALLOCATE (dnraux)
  !
  CALL stop_clock('hp_symdnsq')
  !
  RETURN
  !
END SUBROUTINE hp_symdnsq
