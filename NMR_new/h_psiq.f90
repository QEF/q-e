!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine h_psiq_nmr (lda, n, m, psi, hpsi, spsi)
  !-----------------------------------------------------------------------
  !
  !     This routine computes the product of the Hamiltonian
  !     and of the S matrix with a m  wavefunctions  contained
  !     in psi. It first computes the bec matrix of these
  !     wavefunctions and then with the routines hus_1psi and
  !     s_psi computes for each band the required products
  !

  use pwcom
  USE wavefunctions_module,  ONLY: psic
  USE becmod, ONLY: becp
  USE kinds, only : DP
  use nmr_module
  implicit none
  !
  !     Here the local variables
  !
  integer :: ibnd
  ! counter on bands

  integer :: lda, n, m
  ! input: the leading dimension of the array psi
  ! input: the real dimension of psi
  ! input: the number of psi to compute
  integer :: j
  ! do loop index

  complex(DP) :: psi (lda, m), hpsi (lda, m), spsi (lda, m)
  ! input: the functions where to apply H and S
  ! output: H times psi
  ! output: S times psi (Us PP's only)


  call start_clock ('h_psiq')
  call start_clock ('init')

  call ccalbec (nkb, npwx, n, m, becp, vkb, psi)
  !
  ! Here we apply the kinetic energy (k+G)^2 psi
  !
  do ibnd = 1, m
     do j = 1, n
        hpsi (j, ibnd) = g2kin (j) * psi (j, ibnd)
     enddo
  enddo
  call stop_clock ('init')
  !
  ! the local potential V_Loc psi. First the psi in real space
  !

  do ibnd = 1, m
     call start_clock ('firstfft')
     psic(:) = (0.d0, 0.d0)
     !do j = 1, n
     !   psic (nls(igk(j))) = psi (j, ibnd)
     !enddo
     psic (nls(igk(1:n))) = psi(1:n, ibnd)
     call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)
     call stop_clock ('firstfft')
     !
     !   and then the product with the potential vrs = (vltot+vr) on the smoo
     !
     !do j = 1, nrxxs
     !   psic (j) = psic (j) * vrs (j, current_spin)
     !enddo
     psic (1:nrxxs) = psic (1:nrxxs) * vrs (1:nrxxs, current_spin)
     !
     !   back to reciprocal space
     !
     call start_clock ('secondfft')
     call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, - 2)
     !
     !   addition to the total product
     !
     !do j = 1, n
     !   hpsi (j, ibnd) = hpsi (j, ibnd) + psic (nls(igk(j)))
     !enddo
     hpsi (1:n, ibnd) = hpsi (1:n, ibnd) + psic (nls(igk(1:n)))
     call stop_clock ('secondfft')
  enddo
  !
  !  Here the product with the non local potential V_NL psi
  !

  call add_vuspsi (lda, n, m, psi, hpsi)

  call s_psi (lda, n, m, psi, spsi)

  call stop_clock ('h_psiq')
  return
end subroutine h_psiq_nmr
