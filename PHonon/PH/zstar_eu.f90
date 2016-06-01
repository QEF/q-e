!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine zstar_eu
  !-----------------------------------------------------------------------
  ! calculate the effective charges Z(E,Us) (E=scf,Us=bare)
  !
  ! epsil =.true. is needed for this calculation to be meaningful
  !
  !
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : bg
  USE ions_base, ONLY : nat, zv, ityp
  USE buffers,   ONLY : get_buffer
  USE klist,     ONLY : wk, xk, ngk, igk_k
  USE symme,     ONLY : symtensor
  USE wvfct,     ONLY : npwx
  USE uspp,      ONLY : okvan, vkb
  use noncollin_module, ONLY : npol
  USE wavefunctions_module,  ONLY: evc

  USE modes,     ONLY : u, nirr, npert
  USE qpoint,    ONLY : npwq, nksq
  USE eqv,       ONLY : dvpsi, dpsi
  USE efield_mod,   ONLY : zstareu0, zstareu
  USE units_ph,  ONLY : iudwf, lrdwf, iuwfc, lrwfc
  USE control_lr,ONLY : nbnd_occ
  USE control_ph,ONLY : done_zeu
  USE ph_restart, ONLY : ph_writefile

  USE mp_pools,              ONLY : inter_pool_comm
  USE mp_bands,              ONLY : intra_bgrp_comm
  USE mp,                    ONLY : mp_sum

  implicit none

  integer :: ibnd, ipol, jpol, icart, na, nu, mu, imode0, irr, &
       imode, nrec, mode, ik, ierr, npw
  ! counters
  real(DP) :: weight
  complex(DP), external :: zdotc
  !  scalar product
  !
  call start_clock ('zstar_eu')

  zstareu0(:,:) = (0.d0,0.d0)
  zstareu (:,:,:) = 0.d0

  do ik = 1, nksq
     npw = ngk(ik)
     npwq = npw
     weight = wk (ik)
     if (nksq > 1) call get_buffer (evc, lrwfc, iuwfc, ik)
     call init_us_2 (npw, igk_k(1,ik), xk (1, ik), vkb)
     imode0 = 0
     do irr = 1, nirr
        do imode = 1, npert (irr)
           mode = imode+imode0
           dvpsi(:,:) = (0.d0, 0.d0)
           !
           ! recalculate  DeltaV*psi(ion) for mode nu
           !
           call dvqpsi_us (ik, u (1, mode), .not.okvan)
           do jpol = 1, 3
              nrec = (jpol - 1) * nksq + ik
              !
              ! read dpsi(scf)/dE for electric field in jpol direction
              !
              call get_buffer (dpsi, lrdwf, iudwf, nrec)
              do ibnd = 1, nbnd_occ(ik)
                 zstareu0(jpol,mode)=zstareu0(jpol, mode)-2.d0*weight*&
                      zdotc(npwx*npol,dpsi(1,ibnd),1,dvpsi(1,ibnd),1)
              enddo
           enddo
        enddo
        imode0 = imode0 + npert (irr)
     enddo
  enddo
  !
  ! Now we add the terms which are due to the USPP
  !
  if (okvan) call zstar_eu_us

  call mp_sum ( zstareu0, intra_bgrp_comm )
  call mp_sum ( zstareu0, inter_pool_comm )
  !
  ! bring the mode index to cartesian coordinates
  ! NOTA BENE: the electric field is in crystal axis
  !
  do jpol = 1, 3
     do mu = 1, 3 * nat
        na = (mu - 1) / 3 + 1
        icart = mu - 3 * (na - 1)
        do nu = 1, 3 * nat
           zstareu (jpol, icart, na) = zstareu (jpol, icart, na) + &
                CONJG(u (mu, nu) ) * ( zstareu0 (1,nu) * bg(jpol,1) + &
                                       zstareu0 (2,nu) * bg(jpol,2) + &
                                       zstareu0 (3,nu) * bg(jpol,3) )
        enddo
     enddo
  enddo
  !
  ! symmetrization
  !
  call symtensor ( nat, zstareu )
  !
  ! add the diagonal part
  !
  do ipol = 1, 3
     do na = 1, nat
        zstareu (ipol, ipol, na) = zstareu (ipol, ipol, na) + zv (ityp ( na) )
     enddo
  enddo

  done_zeu=.TRUE.
  call summarize_zeu()
  CALL ph_writefile('tensors',0,0,ierr)

  call stop_clock ('zstar_eu')
  return
end subroutine zstar_eu
