!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE sym_dns (ldim, npe, dns)
  !-----------------------------------------------------------------------
  !! DFPT+U: This routine symmetrizes the first order variation of 
  !! the occupation matrices dns due to the perturbation caused
  !! by the displacement of atoms.  
  !
  !! Written  by S. de Gironcoli and A. Floris.
  !! Modified by I. Timrov (01.10.2018).
  !
  USE kinds,            ONLY : DP
  USE constants,        ONLY : tpi
  USE ions_base,        ONLY : nat, ityp
  USE ldaU,             ONLY : Hubbard_l, is_hubbard
  USE lsda_mod,         ONLY : nspin
  USE lr_symm_base,     ONLY : nsymq, minus_q, irotmq, rtau, upert, upert_mq
  USE qpoint,           ONLY : xq
  USE symm_base,        ONLY : d1, d2, d3, irt
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ldim, npe
  COMPLEX(DP), INTENT(INOUT) :: dns(ldim,ldim,nspin,nat,npe)
  !
  ! ... local variables
  !
  INTEGER :: nt, ip, jp, na, nb, is, m1, m2, &
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
  ALLOCATE (dnraux(ldim,ldim,nspin,nat,npe))
  ALLOCATE (dnr(ldim,ldim,nspin,nat,npe))
  !  
  dnraux = (0.d0, 0.d0)
  dnr    = (0.d0, 0.d0)
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
  dns = (0.d0, 0.d0)
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
                                   dnr(m0,m00,is,nb,jp) * upert_mq(jp,ip) *    &
                                   phase
                                ELSE IF (Hubbard_l(nt).EQ.1) THEN
                                   dns(m1,m2,is,na,ip) = dns(m1,m2,is,na,ip) + &
                                   d1(m0 ,m1,irotmq) * d1(m00,m2,irotmq) *     &
                                   dnr(m0,m00,is,nb,jp) * upert_mq(jp,ip) *    &
                                   phase
                                ELSE IF (Hubbard_l(nt).EQ.2) THEN
                                   dns(m1,m2,is,na,ip) = dns(m1,m2,is,na,ip) + &
                                   d2(m0 ,m1,irotmq) * d2(m00,m2,irotmq) *     &
                                   dnr(m0,m00,is,nb,jp) * upert_mq(jp,ip) *    &
                                   phase
                                ELSE IF (Hubbard_l(nt).EQ.3) THEN
                                   dns(m1,m2,is,na,ip) = dns(m1,m2,is,na,ip) + &
                                   d3(m0 ,m1,irotmq) * d3(m00,m2,irotmq) *     &
                                   dnr(m0,m00,is,nb,jp) * upert_mq(jp,ip) *    &
                                   phase
                                ELSE
                                   CALL errore ('sym_dns', &
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
     dns = (0.d0, 0.d0)
  ENDIF
  ! 
  dnraux = (0.d0, 0.d0)
  !
  ! Symmetryze dnr -> dns
  !
  DO isym = 1, nsymq
     !
     dns = (0.d0, 0.d0)  
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
                          do jp=1, npe
                             IF (Hubbard_l(nt).EQ.0) THEN
                                dns(m1,m2,:,na,ip) = dns(m1,m2,:,na,ip) +  &
                                dnr(m0,m00,:,nb,jp) * upert(jp,ip,irot) * &
                                phase 
                             ELSE IF (Hubbard_l(nt).EQ.1) THEN
                                dns(m1,m2,:,na,ip) = dns(m1,m2,:,na,ip) + &
                                d1(m0 ,m1,irot) * d1(m00,m2,irot) *       &
                                dnr(m0,m00,:,nb,jp) * upert(jp,ip,irot) * &
                                phase 
                             ELSE IF (Hubbard_l(nt).EQ.2) THEN
                                dns(m1,m2,:,na,ip) = dns(m1,m2,:,na,ip) + &
                                d2(m0 ,m1,irot) * d2(m00,m2,irot) *       &
                                dnr(m0,m00,:,nb,jp) * upert(jp,ip,irot) * &
                                phase 
                             ELSE IF (Hubbard_l(nt).EQ.3) THEN
                                dns(m1,m2,:,na,ip) = dns(m1,m2,:,na,ip) + &
                                d3(m0 ,m1,irot) * d3(m00,m2,irot) *       &
                                dnr(m0,m00,:,nb,jp) * upert(jp,ip,irot) * &
                                phase 
                             ELSE
                                CALL errore ('new_ns', &
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
END SUBROUTINE sym_dns
!---------------------------------------------------------------------
