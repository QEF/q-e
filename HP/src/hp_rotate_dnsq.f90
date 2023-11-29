!
! Copyright (C) 2001-2023 Quantum_ESPRESSO group
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
  USE symm_base,    ONLY : d1, d2, d3, irt, invs, t_rev
  USE qpoint,       ONLY : xq
  USE lr_symm_base, ONLY : rtau, gi, minus_q, irotmq
  USE ldaU,         ONLY : Hubbard_lmax, Hubbard_l, is_hubbard, nwfcU, d_spin_ldau
  USE ldaU_hp,      ONLY : nah_pert
  USE noncollin_module, ONLY : npol, noncolin, domag
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN)  :: dnr(2*Hubbard_lmax+1, 2*Hubbard_lmax+1, nspin, nat) 
  COMPLEX(DP), INTENT(OUT) :: dns(2*Hubbard_lmax+1, 2*Hubbard_lmax+1, nspin, nat) 
  INTEGER,  INTENT(IN) :: isym    ! selected symmetry
  REAL(DP), INTENT(IN) :: sxq(3)  ! q point from the star
  COMPLEX(DP), ALLOCATABLE :: dnr_nc(:,:,:,:,:), dnr1_nc(:,:,:,:,:)
  !
  ! Local variables
  !
  INTEGER :: n, l, na, sna, nt, ism1, is, m1, m2, m0, m00, ldim, &
             is1, is2, is3, is4, m3, m4, i
  COMPLEX(DP) :: phase, phase2
  REAL(DP) :: arg, arg2
  !
  CALL start_clock('hp_rotate_dnsq')
  !
  dns = (0.0d0, 0.0d0)
  !
  IF (noncolin) THEN
     IF ( .NOT. ALLOCATED (d_spin_ldau) ) ALLOCATE( d_spin_ldau(2,2,48) )
     CALL comp_dspinldau()
     ALLOCATE( dnr_nc(2*Hubbard_lmax+1, 2*Hubbard_lmax+1, npol, npol, nat) )  
     dnr_nc(:,:,:,:,:) = (0.d0,0.d0)  
     ALLOCATE( dnr1_nc(2*Hubbard_lmax+1, 2*Hubbard_lmax+1, npol, npol, nat) )  
     dnr1_nc(:,:,:,:,:) = (0.d0,0.d0) 
     !
     ! D_Sl for l=1, l=2 and l=3 are already initialized, for l=0 D_S0 is 1
     !
     DO na = 1, nat
        nt = ityp (na)
        IF ( is_hubbard(nt) ) THEN
          DO is1 = 1, npol
            DO is2 = 1, npol
              i = npol*(is1-1) + is2
              ldim = 2*Hubbard_l(nt)+1
              DO m1 = 1, ldim
                DO m2 = 1, ldim
                   dnr_nc(m1,m2,is1,is2,na) =  dnr(m1, m2, i, na) 
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF
     ENDDO
  ENDIF
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
        IF (noncolin) THEN
            DO m1 = 1, 2 * Hubbard_l(nt) + 1  
               DO m2 = 1, 2 * Hubbard_l(nt) + 1
                  DO is1 = 1, npol
                     DO is2 = 1, npol
                        DO m3 = 1, 2 * Hubbard_l(nt) + 1  
                           DO m4 = 1, 2 * Hubbard_l(nt) + 1
                              DO is3 = 1, npol
                                 DO is4 = 1, npol
                                    IF (Hubbard_l(nt) == 0) THEN
                                       IF (t_rev(ism1) == 1) THEN
                                         dnr1_nc(m1,m2,is1,is2,sna) = dnr1_nc(m1,m2,is1,is2,sna) +  &
                                           CONJG( d_spin_ldau(is1,is3,ism1) ) *                 &
                                                  dnr_nc(m4,m3,is4,is3,na) * (phase * phase2) *    &
                                                  d_spin_ldau(is2,is4,ism1)  
                                       ELSE
                                          dnr1_nc(m1,m2,is1,is2,sna) = dnr1_nc(m1,m2,is1,is2,sna) +  &
                                           CONJG( d_spin_ldau(is1,is3,ism1) ) *                 &
                                                  dnr_nc(m3,m4,is3,is4,na) * phase * phase2 *    &
                                                  d_spin_ldau(is2,is4,ism1)  
                                       ENDIF
                                    ELSEIF (Hubbard_l(nt) == 1) THEN
                                       IF (t_rev(ism1) == 1) THEN
                                          dnr1_nc(m1,m2,is1,is2,sna) = dnr1_nc(m1,m2,is1,is2,sna) + &
                                            CONJG( d_spin_ldau(is1,is3,ism1) )*d1(m1,m3,ism1)*  &
                                                   dnr_nc(m4,m3,is4,is3,na)* (phase * phase2) *    &
                                                   d_spin_ldau(is2,is4,ism1)  *d1(m2,m4,ism1)
                                        ELSE
                                          dnr1_nc(m1,m2,is1,is2,sna) = dnr1_nc(m1,m2,is1,is2,sna) + &
                                            CONJG( d_spin_ldau(is1,is3,ism1) )*d1(m1,m3,ism1)*  &
                                                   dnr_nc(m3,m4,is3,is4,na)* phase * phase2 *    &
                                                   d_spin_ldau(is2,is4,ism1)  *d1(m2,m4,ism1)
                                        ENDIF
                                      ELSEIF (Hubbard_l(nt) == 2) THEN
                                        IF (t_rev(ism1) == 1) THEN
                                          dnr1_nc(m1,m2,is1,is2,sna) = dnr1_nc(m1,m2,is1,is2,sna) + &
                                          CONJG( d_spin_ldau(is1,is3,ism1) ) * d2(m1,m3,ism1) *    &
                                                   dnr_nc(m4,m3,is4,is3,na)* (phase * phase2) *    &
                                                   d_spin_ldau(is2,is4,ism1)  *d2(m2,m4,ism1)
                                        ELSE
                                          dnr1_nc(m1,m2,is1,is2,sna) = dnr1_nc(m1,m2,is1,is2,sna) + &
                                            CONJG( d_spin_ldau(is1,is3,ism1) )*d2(m1,m3,ism1)*  &
                                                   dnr_nc(m3,m4,is3,is4,na)* phase * phase2 *    &
                                                   d_spin_ldau(is2,is4,ism1)  *d2(m2,m4,ism1)
                                        ENDIF
                                      ELSEIF (Hubbard_l(nt) == 3) THEN
                                        !
                                        IF (t_rev(ism1) == 1) THEN
                                          dnr1_nc(m1,m2,is1,is2,sna) = dnr1_nc(m1,m2,is1,is2,sna) + &
                                            CONJG( d_spin_ldau(is1,is3,ism1) )*d3(m1,m3,ism1)*  &
                                                   dnr_nc(m4,m3,is4,is3,na)* (phase * phase2) *   &
                                                   d_spin_ldau(is2,is4,ism1)  *d3(m2,m4,ism1)
                                        ELSE
                                          dnr1_nc(m1,m2,is1,is2,sna) = dnr1_nc(m1,m2,is1,is2,sna) + &
                                            CONJG( d_spin_ldau(is1,is3,ism1) )*d3(m1,m3,ism1)*  &
                                                   dnr_nc(m3,m4,is3,is4,na)* phase * phase2 *    &
                                                   d_spin_ldau(is2,is4,ism1)  *d3(m2,m4,ism1)
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
               ldim = 2 * Hubbard_l(nt) + 1
               DO m1 = 1, ldim
                  DO m2 = 1, ldim
                     DO m0 = 1, ldim
                        DO m00 = 1, ldim
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
        ENDIF
        !
     ENDIF
     !
  ENDDO ! na
  !
  IF (noncolin) THEN
     DO na = 1, nat
        nt = ityp (na)
        IF ( is_hubbard(nt) ) THEN
           DO is1 = 1, npol
              DO is2 = 1, npol
                 i = npol*(is1-1) + is2
                 DO m1 = 1, 2*Hubbard_l(nt)+1
                    DO m2 = 1, 2*Hubbard_l(nt)+1
                       dns(m1,m2,i,na) = dnr1_nc(m1,m2,is1,is2,na)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDIF
     ENDDO
  ENDIF
  !
  IF (noncolin) THEN
     DEALLOCATE (dnr_nc)
     DEALLOCATE (dnr1_nc)
     IF ( ALLOCATED (d_spin_ldau) ) DEALLOCATE( d_spin_ldau )
  ENDIF
  !
  CALL stop_clock('hp_rotate_dnsq')
  !
  RETURN
  !
END SUBROUTINE hp_rotate_dnsq
