!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine dielec_test
  !-----------------------------------------------------------------------
  !
  ! Calculates the dielectric tensor using the finite-differences-derivative
  ! of the wavefunctions. This should be used only for testing purposes
  ! while doing a raman calculation
  !
  USE kinds,    ONLY : DP
  USE constants,ONLY : fpi
  USE cell_base,ONLY : omega, at, bg
  USE klist,    ONLY : wk, ngk
  USE symme,    ONLY : symmatrix, crys_to_cart
  USE buffers,  ONLY : get_buffer
  USE wavefunctions_module,  ONLY: evc
  USE efield_mod, ONLY : epsilon
  USE qpoint, ONLY : nksq
  USE eqv, ONLY : dpsi
  USE control_lr, ONLY : nbnd_occ
  USE units_ph, ONLY : lrwfc, iuwfc

  USE ramanm, ONLY : a1j, a2j, lrd2w, iud2w
  USE mp_pools, ONLY : inter_pool_comm
  USE mp_bands, ONLY : intra_bgrp_comm
  USE mp,       ONLY : mp_sum

  implicit none

  INTEGER :: npw
  integer :: ibnd, ipol, jpol, nrec, ik, i1, i2
  real(DP) :: w_, weight, tmp
  complex(DP), external :: zdotc

  epsilon (:,:) = 0.d0
  do ik = 1, nksq
     npw = ngk(ik)
     weight = wk (ik)
     w_ = - fpi * weight / omega
     call get_buffer (evc, lrwfc, iuwfc, ik)
     do ipol = 1, 6
        nrec = (ipol - 1) * nksq + ik
        call davcio (dpsi, lrd2w, iud2w, nrec, -1)
        tmp = 0.d0
        do ibnd = 1, nbnd_occ (ik)
           tmp = tmp + 2.0d0 * w_ *                        &
              real (zdotc (npw, evc (1, ibnd), 1, dpsi (1, ibnd), 1))
        enddo
        i1 = a1j (ipol)
        i2 = a2j (ipol)
        epsilon (i1, i2) = epsilon (i1, i2) + tmp
        if (i1.ne.i2 ) epsilon (i2, i1) = epsilon (i2, i1) + tmp
     enddo
  enddo
  call mp_sum ( epsilon, intra_bgrp_comm )
  call mp_sum ( epsilon, inter_pool_comm )
  !
  !  symmetrize (pass to cartesian axis first)
  !
!  write(6,'(/,10x,''Unsymmetrized in crystal axis '',/)')
!  write(6,'(10x,''('',3f15.5,'' )'')') ((epsilon(ipol,jpol), &
!                                       ipol=1,3),jpol=1,3)

  call crys_to_cart (epsilon)
  call symmatrix(epsilon)
  !
!  write(6,'(/,10x,''Symmetrized in cartesian axis '',/)')
!  write(6,'(10x,''('',3f15.5,'' )'')') ((epsilon(ipol,jpol), &
!                                  ipol=1,3),jpol=1,3)

  !
  ! add the diagonal part
  !
  do ipol = 1, 3
     epsilon (ipol, ipol) = epsilon (ipol, ipol) + 1.d0
  end do
  !
  !  and print the result
  !
  write(6,'(/,10x,''Dielectric constant from finite-differences'',/)')
  write(6,'(10x,''('',3f18.9,'' )'')') ((epsilon(ipol,jpol),       &
                                  ipol=1,3),jpol=1,3)

  return
end subroutine dielec_test

