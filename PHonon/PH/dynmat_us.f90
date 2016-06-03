!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE dynmat_us()
  !-----------------------------------------------------------------------
  !
  !  This routine calculates the electronic term: <psi|V"-eS"|psi>
  !  of the dynamical matrix. Eq. B32 of PRB 64, 235118 (2001) is calculated
  !  here. Eqs. B33 and B34 in addusdynmat.
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : tpi
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp, tau
  USE uspp,                 ONLY : nkb, vkb
  USE scf,                  ONLY : rho
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE buffers,              ONLY : get_buffer
  USE gvect,                ONLY : g, ngm, nl, igtongl
  USE wvfct,                ONLY : npwx, nbnd, wg, et
  USE lsda_mod,             ONLY : lsda, current_spin, isk, nspin
  USE vlocal,               ONLY : vloc
  USE klist,                ONLY : xk, ngk, igk_k
  USE wavefunctions_module, ONLY : evc
  USE cell_base,            ONLY : omega, tpiba2
  USE uspp_param,           ONLY : nh, nhm
  USE noncollin_module,     ONLY : noncolin, npol, nspin_lsda
  USE spin_orb,             ONLY : lspinorb
  USE becmod,               ONLY : calbec, bec_type, allocate_bec_type, &
                                   deallocate_bec_type, beccopy
  USE modes,                ONLY : u
  USE dynmat,               ONLY : dyn
  USE phus,                 ONLY : alphap
  USE units_ph,             ONLY : iuwfc, lrwfc
  USE io_global,            ONLY : stdout
  USE mp_pools,             ONLY : my_pool_id, inter_pool_comm
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum

  USE lrus,                 ONLY : becp1
  USE qpoint,               ONLY : nksq, ikks
  USE control_lr,           ONLY : nbnd_occ, lgamma

  IMPLICIT NONE
  INTEGER :: icart, jcart, na_icart, na_jcart, na, ng, nt, ik, &
       ig, is, ibnd, nu_i, nu_j, ijkb0, ikb, jkb, ih, jh, ikk, &
       js,  ijs, npw
  ! counters
  ! ikk: record position of wfc at k

  REAL(DP) :: gtau, fac, wgg
  ! the product G*\tau_s
  ! auxiliary variable
  ! the true weight of a K point

  COMPLEX(DP) :: work, dynwrk (3 * nat, 3 * nat), fact
  ! work space
  TYPE (bec_type) :: gammap(3,3)
  COMPLEX(DP), ALLOCATABLE :: rhog (:), aux1 (:,:), work1 (:), &
               work2 (:), deff_nc(:,:,:,:)
  REAL(DP), ALLOCATABLE :: deff(:,:,:)
  ! fourier transform of rho
  ! the second derivative of the beta
  ! work space

  CALL start_clock ('dynmat_us')
  ALLOCATE (rhog  ( dfftp%nnr))
  ALLOCATE (work1 ( npwx))
  ALLOCATE (work2 ( npwx))
  ALLOCATE (aux1  ( npwx*npol , nbnd))
  IF (noncolin) THEN
     ALLOCATE (deff_nc( nhm, nhm, nat, nspin ))
  ELSE
     ALLOCATE (deff(nhm, nhm, nat ))
  END IF
  DO icart=1,3
     DO jcart=1,3
        CALL allocate_bec_type(nkb,nbnd, gammap(icart,jcart))
     ENDDO
  ENDDO

  dynwrk (:,:) = (0.d0, 0.0d0)
  !
  !   We first compute the part of the dynamical matrix due to the local
  !   potential

  !   ... only the first pool does the calculation (no sum over k needed)
  IF ( my_pool_id /= 0 ) GOTO 100

  !
  rhog (:) = (0.d0, 0.d0)
  DO is = 1, nspin_lsda
     rhog (:) = rhog (:) + CMPLX(rho%of_r(:, is), 0.d0,kind=DP)
  ENDDO

  CALL fwfft ('Dense', rhog, dfftp)
  !
  ! there is a delta ss'
  !
  DO na = 1, nat
     DO icart = 1, 3
        na_icart = 3 * (na - 1) + icart
        DO jcart = 1, 3
           na_jcart = 3 * (na - 1) + jcart
           DO ng = 1, ngm
              gtau = tpi * (g (1, ng) * tau (1, na) + &
                            g (2, ng) * tau (2, na) + &
                            g (3, ng) * tau (3, na) )
              fac = omega * vloc (igtongl (ng), ityp (na) ) * tpiba2 * &
                   ( DBLE (rhog (nl (ng) ) ) * COS (gtau) - &
                    AIMAG (rhog (nl (ng) ) ) * SIN (gtau) )
              dynwrk (na_icart, na_jcart) = dynwrk (na_icart, na_jcart) - &
                   fac * g (icart, ng) * g (jcart, ng)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  CALL mp_sum (dynwrk, intra_bgrp_comm)
  !
  ! each pool contributes to next term
  !
100 CONTINUE
  !
  ! Here we compute  the nonlocal Ultra-soft contribution
  !
  DO ik = 1, nksq
     ikk = ikks(ik)
     IF (lsda) current_spin = isk (ikk)
     npw = ngk(ikk)
     IF (nksq > 1) CALL get_buffer (evc, lrwfc, iuwfc, ikk)
     CALL init_us_2 (npw, igk_k(1,ikk), xk (1, ikk), vkb)
     !
     !    We first prepare the gamma terms, which are the second derivatives
     !    becp terms.
     !
     DO icart = 1, 3
        DO jcart = 1, icart
           aux1=(0.d0,0.d0)
           DO ibnd = 1, nbnd
              DO ig = 1, npw
                 aux1 (ig, ibnd) = - evc (ig, ibnd) * tpiba2 * &
                      (xk (icart, ikk) + g (icart, igk_k(ig,ikk) ) ) * &
                      (xk (jcart, ikk) + g (jcart, igk_k(ig,ikk) ) )
              ENDDO
              IF (noncolin) THEN
                 DO ig = 1, npw
                    aux1 (ig+npwx, ibnd) = - evc (ig+npwx, ibnd) * tpiba2 * &
                      (xk (icart, ikk) + g (icart, igk_k(ig,ikk) ) ) * &
                      (xk (jcart, ikk) + g (jcart, igk_k(ig,ikk) ) )
                 ENDDO
              END IF
           ENDDO

           CALL calbec ( npw, vkb, aux1, gammap(icart,jcart) )
           IF (jcart < icart) &
              CALL beccopy (gammap(icart,jcart),gammap(jcart,icart), nkb, nbnd)
        ENDDO
     ENDDO
     !
     !   And then compute the contribution from the US pseudopotential
     !   which is  similar to the KB one
     !
     DO ibnd = 1, nbnd_occ (ikk)
        wgg = wg (ibnd, ikk)
        IF (noncolin) THEN
           CALL compute_deff_nc(deff_nc,et(ibnd,ikk))
        ELSE
           CALL compute_deff(deff,et(ibnd,ikk))
        ENDIF
        ijkb0 = 0
        DO nt = 1, ntyp
           DO na = 1, nat
              IF (ityp (na) == nt) THEN
                 DO icart = 1, 3
                    na_icart = 3 * (na - 1) + icart
                    DO jcart = 1, 3
                       na_jcart = 3 * (na - 1) + jcart
                       DO ih = 1, nh (nt)
                          ikb = ijkb0 + ih
                          DO jh = 1, nh (nt)
                             jkb = ijkb0 + jh
                             IF (noncolin) THEN
                                ijs=0
                                DO is=1,npol
                                   DO js=1,npol
                                      ijs=ijs+1
                                      dynwrk(na_icart,na_jcart) = &
                                        dynwrk(na_icart,na_jcart) + &
                                             wgg* deff_nc(ih,jh,na,ijs) * &
                                  (CONJG(gammap(icart,jcart)%nc(ikb,is,ibnd))*&
                                     becp1(ik)%nc (jkb, js, ibnd) + &
                                     CONJG(becp1(ik)%nc(ikb, is, ibnd) ) * &
                                     gammap(icart,jcart)%nc(jkb, js, ibnd) + &
                                     CONJG(alphap(icart,ik)%nc(ikb,is,ibnd))* &
                                     alphap(jcart,ik)%nc(jkb, js, ibnd) + &
                                     CONJG(alphap(jcart,ik)%nc(ikb,is,ibnd))*&
                                     alphap(icart,ik)%nc(jkb, js, ibnd) )
                                   END DO
                                END DO
                             ELSE
                                dynwrk(na_icart,na_jcart) = &
                                  dynwrk(na_icart,na_jcart) + &
                                  deff (ih, jh, na)* wgg * &
                                  (CONJG(gammap(icart,jcart)%k(ikb,ibnd)) *&
                                   becp1(ik)%k (jkb, ibnd) + &
                                   CONJG (becp1(ik)%k (ikb, ibnd) ) * &
                                   gammap(icart,jcart)%k(jkb,ibnd) + &
                                   CONJG (alphap(icart,ik)%k(ikb, ibnd) ) * &
                                   alphap(jcart,ik)%k(jkb, ibnd) + &
                                   CONJG (alphap(jcart,ik)%k(ikb, ibnd) ) * &
                                   alphap(icart,ik)%k(jkb, ibnd) )
                             END IF
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDDO
                 ijkb0 = ijkb0 + nh (nt)
              ENDIF
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  !   For true US pseudopotentials there is an additional term in the second
  !   derivative which is due to the change of the self consistent D part
  !   when the atom moves. We compute these terms in an additional routine
  !
  CALL addusdynmat (dynwrk)
  !
  CALL mp_sum ( dynwrk, inter_pool_comm )
  !
  !      do na = 1,nat
  !         do nb = 1,nat
  !           WRITE( stdout, '(2i3)') na,nb
  !            do icart = 1,3
  !              na_icart = 3*(na-1)+icart
  !               WRITE( stdout,'(6f13.8)')  &
  !                     (dynwrk(na_icart,3*(nb-1)+jcart), jcart=1,3)
  !            end do
  !         end do
  !      end do
  !      call stop_ph(.false.)
  !
  !  We rotate the dynamical matrix on the basis of patterns
  !

  CALL rotate_pattern_add(nat, u, dyn, dynwrk)

  IF (noncolin) THEN
     DEALLOCATE (deff_nc)
  ELSE
     DEALLOCATE (deff)
  END IF
  DO icart=1,3
     DO jcart=1,3
        CALL deallocate_bec_type(gammap(icart,jcart))
     ENDDO
  ENDDO
  DEALLOCATE (aux1)
  DEALLOCATE (work2)
  DEALLOCATE (work1)
  DEALLOCATE (rhog)

  CALL stop_clock ('dynmat_us')
  RETURN
END SUBROUTINE dynmat_us
