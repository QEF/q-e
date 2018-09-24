!
! Copyright (C) 2001-2018 Quantum_ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE hp_dnstot_sum_q
!-----------------------------------------------------------------------
  !
  ! This subroutine sums up over q points the contributions to
  ! the response occupation matrices (including the phase exp(i q*R)).
  ! If skip_equivalence_q=.FALSE. it also generates and sums over
  ! the q' points from the star of each q point.
  ! Inspired by PH/q2qstar_ph.f90
  ! See Eq. (42) in Ref. [1].
  ! [1] Phys. Rev. B 98, 085127 (2018)
  !
  USE kinds,          ONLY : DP
  USE ions_base,      ONLY : nat, tau, ntyp => nsp, ityp
  USE constants,      ONLY : tpi
  USE qpoint,         ONLY : xq
  USE lsda_mod,       ONLY : nspin
  USE cell_base,      ONLY : at, bg
  USE symm_base,      ONLY : nsym, s, invs, irt, d1, d2, d3
  USE io_global,      ONLY : stdout
  USE control_flags,  ONLY : iverbosity
  USE lr_symm_base,   ONLY : nsymq, invsymq, minus_q, rtau
  USE ldaU,           ONLY : Hubbard_lmax, Hubbard_l, is_hubbard
  USE ldaU_hp,        ONLY : nqsh, Rvect, dnsscf, dns0, dnsscf_tot, dns0_tot, &
                             skip_equivalence_q, nq1, nq2, nq3, x_q, nqs
  !
  IMPLICIT NONE
  !
  INTEGER :: na, nt, m1, m2, is, iq, iq_star, isym
  INTEGER :: nq,      & ! degeneracy of the star of q
             isq(48), & ! index of q in the star of a given sym.op.
             imq,     & ! index of -q in the star of q (0 if not present)
             nqs_all    ! all q points in the grid without symmetry
  REAL(DP) :: sxq(3,48) ! list of vectors in the star of q
  COMPLEX(DP), ALLOCATABLE :: dns0_rot(:,:,:,:),  & 
                              dnsscf_rot(:,:,:,:)
  ! Rotated bare and SCF response occupation matrices
  LOGICAL :: full_q_grid 
  ! If .true. then all q points in the grid (without symmetry) were computed
  ! and hence there is no need to consider the star of q
  !
  CALL start_clock('hp_dnstot_sum_q')
  !
  WRITE( stdout, '(/,5X,"Computing the sum over q of the response occupation matrices...")')
  !
  ALLOCATE (dnsscf_rot(2*Hubbard_lmax+1, 2*Hubbard_lmax+1, nspin, nat))
  ALLOCATE (dns0_rot  (2*Hubbard_lmax+1, 2*Hubbard_lmax+1, nspin, nat))
  !
  nqs_all = nq1*nq2*nq3 
  !
  IF ( nqs == nqs_all ) THEN
     full_q_grid = .TRUE.
  ELSE
     full_q_grid = .FALSE.
  ENDIF
  !
  ! Loop over q points
  !
  DO iq = 1, nqs
     !
     ! Set the coordinates of the current q point
     !
     xq(1:3) = x_q(1:3,iq)
     !
     WRITE( stdout, '(/,5X,"q #",i4," = ",f12.9,2x,f12.9,2x,f12.9)') iq, xq(1), xq(2), xq(3)
     !
     IF ( full_q_grid ) THEN
       !
       nq = 1
       imq = 1
       sxq(:,1) = xq(:)
       !
     ELSE
       !
       ! Set the small group of q
       !
       CALL set_small_group_of_q (nsymq, invsymq, minus_q)
       !
       ! Calculate rtau (the Bravais lattice vector associated to a rotation) 
       ! with the new symmetry order
       !
       CALL sgam_lr (at, bg, nsym, s, irt, tau, rtau, nat)
       !
       ! Since the order of the S matrices is changed (for q\=0) 
       ! we need to re-initialize d1, d2, d3 to rotate the spherical harmonics
       !
       CALL d_matrix( d1, d2, d3 )
       !
       ! Generate the star of q 
       !
       CALL star_q (xq, at, bg, nsym, s, invs, nq, sxq, isq, imq, .TRUE. )
       !
     ENDIF
     !
     ! Loop over q' points in the star of q
     !
     DO iq_star = 1, nq
        !
        IF ( iq_star == 1 ) THEN
           !
           ! There is no need in rotation. Just copy.
           !
           dnsscf_rot(:,:,:,:) = dnsscf(:,:,:,:,iq)
           dns0_rot(:,:,:,:)   = dns0(:,:,:,:,iq)
           ! 
        ELSE
           !
           ! Rotate dnsscf and dns0
           !
           DO isym = 1, nsym
              !
              ! Select the appropriate symmetry which rotates q to q'
              !
              IF (isq(isym)==iq_star) THEN
                 !
                 CALL hp_rotate_dnsq(dnsscf(:,:,:,:,iq), dnsscf_rot(:,:,:,:), isym, sxq(:,iq_star))
                 CALL hp_rotate_dnsq(dns0(:,:,:,:,iq), dns0_rot(:,:,:,:), isym, sxq(:,iq_star))
                 GO TO 100
                 !
              ENDIF
              !
           ENDDO
           !
        ENDIF
        !
100 CONTINUE
        !
        ! Print the response occupation matrices
        !
        IF ( iverbosity > 3 ) CALL print_dns(sxq(:,iq_star), dnsscf_rot)
        !
        ! Add a contribution from this q' point from the star
        !     
        CALL dns_sum_q(sxq(:,iq_star))
        !
        ! Add a contribution from the -q' point if it is not in the star
        ! (i.e. it is in the -q' list)
        !
        IF ( imq == 0 ) THEN
           !
           ! Use the time-reversal symmetry
           !
           CALL dns_conjg()
           !
           IF ( iverbosity > 3 ) THEN
              WRITE(stdout,'(/5x,"Add a contribution from -q which is in a separate list!")')
              CALL print_dns(-sxq(:,iq_star), dnsscf_rot)
           ENDIF
           ! 
           ! Add a contribution from the -q' point 
           !
           CALL dns_sum_q(-sxq(:,iq_star))
           !
        ENDIF
        !
     ENDDO 
     !
  ENDDO
  !
  DEALLOCATE(dns0_rot)
  DEALLOCATE(dnsscf_rot)
  !
  CALL stop_clock('hp_dnstot_sum_q')
  !
  RETURN
  !
CONTAINS
  !
SUBROUTINE dns_sum_q (xq_)
  !
  ! This subroutine sums up contributions from all q points
  ! (including the contributions from the points in the star 
  ! of q if skip_equivalence_q=.FALSE.)
  !
  ! Comment: It may be that the point q' in the star of q 
  ! differs by sign with respect to the corresponding q 
  ! point in the list of the full grid of q points 
  ! (i.e. when skip_equivalence_q=.TRUE.). It seems that
  ! this may happen only when minus_q=.TRUE. However, since
  ! we are interested in the real part of the quantity
  ! exp(iq'*R) dns(q') there is no problem with the sign.
  ! Note that dns(-q') = CONJG(dns(q')), hence
  ! Re[exp(-iq'*R) dns(-q')] = Re[exp(iq'*R) dns(q')].
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: xq_(3) ! q point from the star 
  !
  INTEGER :: icell
  REAL(DP) :: fac, arg
  COMPLEX(DP) :: phase
  !
  ! Normalization factor
  !
  fac = 1.0d0/DBLE(nqsh)
  !
  ! Loop over R points
  !
  DO icell = 1, nqsh
     !
     ! Calculate the phase factor exp(iqR)
     !
     arg = ( xq_(1) * Rvect(1,icell) + &
             xq_(2) * Rvect(2,icell) + &
             xq_(3) * Rvect(3,icell) ) * tpi
     !
     phase = CMPLX (cos(arg), sin(arg), kind=DP)
     !
     ! Accumulate (sum up) over q points in order to obtain 
     ! total occupation matrices for every cell
     !
     DO na = 1, nat
        !
        nt = ityp(na)
        !
        IF (is_hubbard(nt)) THEN
           !
           DO is = 1, nspin
              !
              DO m1 = 1, 2 * Hubbard_l(nt) + 1
                 !
                 DO m2 = m1, 2 * Hubbard_l(nt) + 1
                    !
                    dnsscf_tot(m1,m2,is,na,icell) = dnsscf_tot(m1,m2,is,na,icell) + &
                                                    & fac * phase * dnsscf_rot(m1,m2,is,na)
                    !
                    dns0_tot(m1,m2,is,na,icell)   = dns0_tot(m1,m2,is,na,icell)   + &
                                                    & fac * phase * dns0_rot(m1,m2,is,na)
                    !
                 ENDDO 
                 !
              ENDDO 
              !
           ENDDO 
           ! 
        ENDIF
        !
     ENDDO 
     !
  ENDDO 
  !
  RETURN
  !
END SUBROUTINE dns_sum_q

SUBROUTINE dns_conjg()
  !
  ! Use the time-reversal symmetry: dns = conjg(dns)
  ! 
  IMPLICIT NONE
  !
  DO na = 1, nat
     nt = ityp(na)
     IF (is_hubbard(nt)) THEN
        DO is = 1, nspin
           DO m1 = 1, 2 * Hubbard_l(nt) + 1
              DO m2 = m1, 2 * Hubbard_l(nt) + 1
                 dnsscf_rot(m1,m2,is,na) = CONJG(dnsscf_rot(m1,m2,is,na))
                 dns0_rot(m1,m2,is,na)   = CONJG(dns0_rot(m1,m2,is,na))
              ENDDO
           ENDDO
        ENDDO
     ENDIF
  ENDDO
  !
  RETURN
  !
END SUBROUTINE dns_conjg

SUBROUTINE print_dns(xq_, dns_)
  !
  ! Prints xq_, dns_ (for debugging)
  !
  IMPLICIT NONE
  !
  REAL(DP),    INTENT(IN) :: xq_(3) ! q point from the star
  COMPLEX(DP), INTENT(IN) :: dns_(2*Hubbard_lmax+1, 2*Hubbard_lmax+1, nspin, nat)  
  !
  WRITE( stdout, '(/,5X,"q = ",3(2x,f12.9))') xq_(1), xq_(2), xq_(3)
  ! 
  DO na = 1, nat
     nt = ityp(na)
     IF (is_hubbard(nt)) THEN
        WRITE(stdout,'(5x,"na = ",1x,i3)') na
        DO is = 1, nspin
           WRITE(stdout,'(5x,"is = ",1x,i2)') is
           DO m1 = 1, 2 * Hubbard_l(nt) + 1
              DO m2 = m1, 2 * Hubbard_l(nt) + 1
                 !
                 WRITE(stdout,'(5x,"m1 = ",i1,2x, "m2 = ",i1, 2x, f12.8, 2x, f12.8)') &
                   & m1, m2, DBLE(dns_(m1,m2,is,na)), AIMAG(dns_(m1,m2,is,na))
                 !
              ENDDO
           ENDDO
        ENDDO
     ENDIF
  ENDDO
  !
  RETURN
  !
END SUBROUTINE print_dns

END SUBROUTINE hp_dnstot_sum_q
