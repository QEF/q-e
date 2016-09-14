!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine incdrhous (drhoscf, weight, ik, dbecsum, evcr, wgg, becq, &
     alpq, mode)
  !-----------------------------------------------------------------------
  !
  !     This routine computes the change of the charge density due
  !     to the displacement of the augmentation charge. Only the
  !     smooth part is computed here.
  !
  USE kinds, only : DP
  USE ions_base, ONLY : ntyp => nsp, nat, ityp
  USE cell_base, ONLY : omega
  USE fft_base,  ONLY : dffts
  USE fft_interfaces, ONLY: invfft
  USE gvecs,   ONLY : nls
  USE noncollin_module, ONLY : npol
  USE uspp,      ONLY : nkb, qq
  USE uspp_param,ONLY : nhm, nh
  USE wvfct,     ONLY : nbnd, npwx
  USE phus,      ONLY : alphap
  USE modes,     ONLY : u
  USE mp_bands,  ONLY : intra_bgrp_comm
  USE mp,        ONLY : mp_sum
  USE becmod,    ONLY : bec_type
  USE klist,     ONLY : ngk, igk_k
  USE lrus,      ONLY : becp1
  USE qpoint,    ONLY : nksq, ikqs, ikks
  USE eqv,       ONLY : evq, dpsi
  USE control_lr, ONLY: nbnd_occ

  implicit none

  integer :: ik, mode
  ! input: the k point
  ! input: the mode which is computed

  real(DP) :: weight, wgg (nbnd, nbnd, nksq)
  ! input: the weight of the k point
  ! input: the weights

  complex(DP) :: evcr (dffts%nnr, nbnd), drhoscf (dffts%nnr), &
       dbecsum(nhm * (nhm + 1) / 2, nat)
  ! input: the wavefunctions at k in real
  ! output: the change of the charge densi
  ! inp/out: the accumulated dbec
  type(bec_type) ::  becq (nksq), &! (nkb, nbnd)
                     alpq (3, nksq)
  ! input: the becp with psi_{k+q}
  ! input: the alphap with psi_{k+q}
  !
  !   here the local variable
  !
  real(DP) :: wgt
  ! the effective weight of the k point

  complex(DP), allocatable :: ps1 (:,:), dpsir (:)
  ! auxiliary space
  ! the change of wavefunctions in real sp

  integer :: ibnd, jbnd, nt, na, mu, ih, jh, ikb, jkb, ijkb0, &
       startb, lastb, ipol, ikk, ir, ig, npwq, ikq
  ! counters

  call start_clock ('incdrhous')
  allocate (dpsir( dffts%nnr))
  allocate (ps1  ( nbnd , nbnd))

  call divide (intra_bgrp_comm, nbnd, startb, lastb)
  ps1 (:,:) = (0.d0, 0.d0)
  ikk=ikks(ik)
  ikq=ikqs(ik)
  npwq = ngk(ikq)
  !
  !   Here we prepare the two terms
  !
  ijkb0 = 0
  do nt = 1, ntyp
     do na = 1, nat
        if (ityp (na) == nt) then
           mu = 3 * (na - 1)
           if (abs(u(mu+1,mode)) + abs(u(mu+2,mode)) &
                                 + abs(u(mu+3,mode)) > 1.0d-12) then
              do ih = 1, nh (nt)
                 ikb = ijkb0 + ih
                 do jh = 1, nh (nt)
                    jkb = ijkb0 + jh
                    do ibnd = 1, nbnd
                       do jbnd = startb, lastb
                          do ipol = 1, 3
                             mu = 3 * (na - 1) + ipol
                             ps1(ibnd,jbnd) = ps1(ibnd,jbnd) - qq(ih,jh,nt) * &
                      ( alphap(ipol,ik)%k(ikb,ibnd) * CONJG(becq(ik)%k(jkb,jbnd)) + &
                        becp1(ik)%k(ikb,ibnd) * CONJG(alpq(ipol,ik)%k(jkb,jbnd)) ) * &
                        wgg (ibnd, jbnd, ik) * u (mu, mode)
                          enddo
                       enddo
                    enddo
                 enddo
              enddo
           endif
           ijkb0 = ijkb0 + nh (nt)
        endif
     enddo
  enddo
#if defined(__MPI)
  call mp_sum (ps1,intra_bgrp_comm)
#endif
  dpsi (:,:) = (0.d0, 0.d0)
  wgt = 2.d0 * weight / omega
  do ibnd = 1, nbnd_occ (ikk)
     do jbnd = 1, nbnd
        call zaxpy (npwq, ps1(ibnd,jbnd), evq(1,jbnd), 1, dpsi(1,ibnd), 1)
     enddo
     dpsir(:) = (0.d0, 0.d0)
     do ig = 1, npwq
        dpsir(nls(igk_k(ig,ikq))) = dpsi (ig, ibnd)
     enddo
     CALL invfft ('Wave', dpsir, dffts)
     do ir = 1, dffts%nnr
        drhoscf(ir) = drhoscf(ir) + wgt * dpsir(ir) * CONJG(evcr(ir,ibnd))
     enddo
  enddo

  call addusdbec (ik, weight, dpsi, dbecsum)
  deallocate (ps1)
  deallocate (dpsir)

  call stop_clock ('incdrhous')
  return
end subroutine incdrhous
