!
! Copyright (C) 2001-2006 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------
SUBROUTINE dynmat_us()
  !-----------------------------------------------------------------------
  !
  !  This routine calculates the electronic term: <psi|V"|psi>
  !  of the dynamical matrix.
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : tpi
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp, tau
  USE uspp,                 ONLY : deeq, nkb, vkb, qq
  USE scf,                  ONLY : rho
  USE gvect,                ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx
  USE gvect,                ONLY : g, ngm, nl, igtongl
  USE wvfct,                ONLY : npw, npwx, nbnd, igk, wg, et
  USE lsda_mod,             ONLY : nspin, lsda, current_spin, isk
  USE vlocal,               ONLY : vloc
  USE klist,                ONLY : xk
  USE wavefunctions_module, ONLY : evc
  USE cell_base,            ONLY : omega, tpiba2
  USE io_files,             ONLY : iunigk
  USE uspp_param,           ONLY : nh
  USE phcom
  USE mp_global,            ONLY : my_pool_id

  IMPLICIT NONE
  INTEGER :: icart, jcart, na_icart, na_jcart, na, nb, ng, nt, ik, &
       ig, ir, is, ibnd, nu_i, nu_j, ijkb0, ikb, jkb, ih, jh, ikk
  ! counters
  ! ikk: record position of wfc at k

  REAL(DP) :: gtau, fac, wgg
  ! the product G*\tau_s
  ! auxiliary variable
  ! the true weight of a K point

  COMPLEX(DP) :: work, dynwrk (3 * nat, 3 * nat)
  ! work space
  COMPLEX(DP), ALLOCATABLE :: rhog (:), &
       gammap (:,:,:,:), aux1 (:,:), work1 (:), work2 (:)
  ! fourier transform of rho
  ! the second derivative of the beta
  ! work space

  CALL start_clock ('dynmat_us')
  ALLOCATE (rhog  ( nrxx))    
  ALLOCATE (work1 ( npwx))    
  ALLOCATE (work2 ( npwx))    
  ALLOCATE (aux1  ( npwx , nbnd))    
  ALLOCATE (gammap(  nkb, nbnd , 3 , 3))    

  dynwrk (:,:) = (0.d0, 0.0d0)
  !
  !   We first compute the part of the dynamical matrix due to the local
  !   potential

  !   ... only the first pool does the calculation (no sum over k needed)
  IF ( my_pool_id /= 0 ) GOTO 100

  !
  rhog (:) = (0.d0, 0.d0)
  DO is = 1, nspin
     rhog (:) = rhog (:) + CMPLX (rho (:, is), 0.d0)
  ENDDO

  CALL cft3 (rhog, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)
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

  CALL reduce (18 * nat * nat, dynwrk)
  !
  ! each pool contributes to next term
  !
100 CONTINUE
  !      goto 500
  !
  ! Here we compute  the nonlocal Ultra-soft contribution
  !
  IF (nksq > 1) REWIND (unit = iunigk)
  DO ik = 1, nksq
     IF (lgamma) THEN
        ikk = ik
     ELSE
        ikk = 2 * ik - 1
     ENDIF
     IF (lsda) current_spin = isk (ikk)
     IF (nksq > 1) READ (iunigk) npw, igk

     ! npwq and igkq are not actually used
     IF (nksq >1 .AND. .NOT.lgamma) READ (iunigk) npwq, igkq

     IF (nksq > 1) CALL davcio (evc, lrwfc, iuwfc, ikk, - 1)
     CALL init_us_2 (npw, igk, xk (1, ikk), vkb)
     !
     !    We first prepare the gamma terms, which are the second derivatives
     !    becp terms.
     !
     DO icart = 1, 3
        DO jcart = 1, icart
           DO ibnd = 1, nbnd
              DO ig = 1, npw
                 aux1 (ig, ibnd) = - evc (ig, ibnd) * tpiba2 * &
                      (xk (icart, ikk) + g (icart, igk (ig) ) ) * &
                      (xk (jcart, ikk) + g (jcart, igk (ig) ) )
              ENDDO
           ENDDO

           CALL ccalbec(nkb,npwx,npw,nbnd,gammap(1,1,icart,jcart),vkb,aux1)
           IF (jcart < icart) THEN
               CALL ZCOPY (nkb * nbnd, gammap (1, 1, icart, jcart), 1, &
                                       gammap (1, 1, jcart, icart), 1)
            END IF
        ENDDO
     ENDDO
     !
     !   And then compute the contribution from the US pseudopotential
     !   which is  similar to the KB one
     !
     ijkb0 = 0
     DO nt = 1, ntyp
        DO na = 1, nat
           IF (ityp (na) == nt) THEN
              DO icart = 1, 3
                 na_icart = 3 * (na - 1) + icart
                 DO jcart = 1, 3
                    na_jcart = 3 * (na - 1) + jcart
                    DO ibnd = 1, nbnd_occ (ikk)
                       wgg = wg (ibnd, ikk)
                       DO ih = 1, nh (nt)
                          ikb = ijkb0 + ih
                          DO jh = 1, nh (nt)
                             jkb = ijkb0 + jh
                             dynwrk(na_icart,na_jcart) = &
                                  dynwrk(na_icart,na_jcart) + &
                                  (deeq (ih, jh, na, current_spin) - &
                                   et (ibnd, ikk) * qq (ih,jh,nt) ) * wgg * &
                                  (CONJG(gammap(ikb, ibnd, icart, jcart)) *&
                                   becp1 (jkb, ibnd, ik) + &
                                   CONJG (becp1 (ikb, ibnd, ik) ) * &
                                   gammap (jkb, ibnd, icart, jcart) + &
                                   CONJG (alphap (ikb, ibnd, icart, ik) ) * &
                                   alphap (jkb, ibnd, jcart, ik) + &
                                   CONJG (alphap (ikb, ibnd, jcart, ik) ) * &
                                   alphap (jkb, ibnd, icart, ik) )
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
              ijkb0 = ijkb0 + nh (nt)
           ENDIF
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
  CALL poolreduce (18 * nat * nat, dynwrk)
  !
  !      do na = 1,nat
  !         do nb = 1,nat
  !           WRITE( stdout, '(2i3)') na,nb
  !            do icart = 1,3
  !              na_icart = 3*(na-1)+icart
  !               WRITE( stdout,'(6f13.8)')
  !     +               (dynwrk(na_icart,3*(nb-1)+jcart), jcart=1,3)
  !            end do
  !         end do
  !      end do
  !      call stop_ph(.false.)
  !
  !  We rotate the dynamical matrix on the basis of patterns
  !
  DO nu_i = 1, 3 * nat
     DO nu_j = 1, 3 * nat
        work = (0.0d0, 0.0d0)
        DO na_jcart = 1, 3 * nat
           DO na_icart = 1, 3 * nat
              work = work + CONJG (u (na_icart, nu_i) ) * &
                            dynwrk (na_icart, na_jcart) * &
                            u (na_jcart, nu_j)
           ENDDO
        ENDDO
        dyn (nu_i, nu_j) = dyn (nu_i, nu_j) + work
     ENDDO

  ENDDO
  DEALLOCATE (gammap)
  DEALLOCATE (aux1)
  DEALLOCATE (work2)
  DEALLOCATE (work1)
  DEALLOCATE (rhog)

  CALL stop_clock ('dynmat_us')
  RETURN
END SUBROUTINE dynmat_us
