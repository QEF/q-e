!
! Copyright (C) 2001-2018 Quantum_ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE hp_rotate_dnsq (dnr, dns, isym, sxq)
!-----------------------------------------------------------------------
  !
  ! This subroutine rotates the response occupation matrix from the q point 
  ! to the q' point, where q' is from the star of q
  !
  ! Inspired by PH/rotate_and_add_dyn.f90
  !
  USE kinds,        ONLY : DP
  USE io_global,    ONLY : stdout
  USE constants,    ONLY : tpi
  USE ions_base,    ONLY : nat, ityp, tau
  USE lsda_mod,     ONLY : nspin
  USE uspp_param,   ONLY : upf
  USE symm_base,    ONLY : d1, d2, d3, irt, invs
  USE qpoint,       ONLY : xq
  USE lr_symm_base, ONLY : rtau, gi, minus_q, irotmq
  USE ldaU,         ONLY : Hubbard_lmax, Hubbard_l, is_hubbard, nwfcU
  USE ldaU_hp,      ONLY : nah_pert
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN)  :: dnr(2*Hubbard_lmax+1, 2*Hubbard_lmax+1, nspin, nat) 
  COMPLEX(DP), INTENT(OUT) :: dns(2*Hubbard_lmax+1, 2*Hubbard_lmax+1, nspin, nat) 
  INTEGER,  INTENT(IN) :: isym    ! selected symmetry
  REAL(DP), INTENT(IN) :: sxq(3)  ! q point from the star
  !
  ! Local variables
  !
  INTEGER :: counter, n, l, na, sna, nt, ism1, is, m1, m2, m0, m00
  COMPLEX(DP) :: phase, phase2
  REAL(DP) :: arg, arg2
  !
  CALL start_clock('hp_rotate_dnsq')
  !
  dns = (0.0d0, 0.0d0)
  !
  ! D_Sl for l=1, l=2 and l=3 are already initialized, for l=0 D_S0 is 1
  !
  ! Check
  !
  counter = 0
  DO na = 1, nat
     nt = ityp(na)
     IF (.NOT.is_hubbard(nt)) CYCLE
     DO n = 1, upf(nt)%nwfc
        l = upf(nt)%lchi(n)
        IF (upf(nt)%oc(n) >= 0.d0 .AND. l == Hubbard_l(nt)) &
           counter = counter + 2 * l + 1
     ENDDO
  ENDDO
  IF (counter.NE.nwfcU) CALL errore ('hp_rotate_dnsq', 'nwfcU<>counter', 1)
  !
  ! Rotate the response occupation matrix.
  ! Note: here we perform a rotation using the index
  ! of the inverse symmetry operation S^{-1}_isym=S(invs(isym)).
  ! We do so because invs(isym) is also used in star_q.
  !
  ism1 = invs(isym)
  !
  DO na = 1, nat
     !
     nt = ityp(na)
     !
     sna = irt(isym, na)
     !
     IF (is_hubbard(nt)) THEN
        !
        ! Compute the phase factor exp(-iq'*f)
        !
        arg = ( sxq(1) * rtau(1, isym, na) + &
                sxq(2) * rtau(2, isym, na) + &
                sxq(3) * rtau(3, isym, na) ) * tpi
        phase = CMPLX (cos(arg), -sin(arg), kind=DP)
        !
        ! Compute the phase factor exp(-iG*tau_pert)
        ! where G = q' - q, tau_pert is the position
        ! of the perturbed atom, and q' is the point
        ! from the star of q 
        !
        arg2 = ((sxq(1) - xq(1)) * tau(1,nah_pert) + &
                (sxq(2) - xq(2)) * tau(2,nah_pert) + &
                (sxq(3) - xq(3)) * tau(3,nah_pert) ) * tpi
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
                          dns(m1,m2,is,sna) = dns(m1,m2,is,sna) + &
                                   dnr(m0,m00,is,na) * phase * phase2
                       ELSEIF (Hubbard_l(nt).EQ.1) THEN
                          dns(m1,m2,is,sna) = dns(m1,m2,is,sna) +       &
                                   d1(m0,m1,ism1) * d1(m00,m2,ism1) * &
                                   dnr(m0,m00,is,na) * phase * phase2
                       ELSEIF (Hubbard_l(nt).EQ.2) THEN
                          dns(m1,m2,is,sna) = dns(m1,m2,is,sna) +       &
                                   d2(m0,m1,ism1) * d2(m00,m2,ism1) * &
                                   dnr(m0,m00,is,na) * phase * phase2
                       ELSEIF (Hubbard_l(nt).EQ.3) THEN
                          dns(m1,m2,is,sna) = dns(m1,m2,is,sna) +       &
                                   d3(m0,m1,ism1) * d3(m00,m2,ism1) * &
                                   dnr(m0,m00,is,na) * phase * phase2
                       ELSE
                          CALL errore ('hp_rotate_dnsq', 'angular momentum not implemented', &
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
  CALL stop_clock('hp_rotate_dnsq')
  !
  RETURN
  !
END SUBROUTINE hp_rotate_dnsq
