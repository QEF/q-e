!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine h_psi (lda, n, m, psi, hpsi, spsi)
  !-----------------------------------------------------------------------
  !
  !     This routine computes the product of the Hamiltonian
  !     and of the S matrix with m wavefunctions contained in psi
  ! input:
  !     lda   leading dimension of arrays psi, spsi, hpsi
  !     n     true dimension of psi, spsi, hpsi
  !     m     number of states psi
  !     psi
  ! output:
  !     hpsi  H*psi
  !     spsi  S*psi
  !
  use pwcom
  use becmod
  implicit none
  !
  integer :: lda, n, m
  complex(kind=DP) :: psi (lda, m), hpsi (lda, m), spsi (lda, m)
  !
  integer :: ibnd, j
  ! counters

  call start_clock ('h_psi')
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
  ! Here we add the Hubbard potential times psi
  !
  if (lda_plus_u) call vhpsi (lda, n, m, psi, hpsi)
  !
  ! the local potential V_Loc psi. First the psi in real space
  !
  do ibnd = 1, m
     call start_clock ('firstfft')
     psic(1:nrxxs) = (0.d0,0.d0)
     do j = 1, n
        psic(nls(igk(j))) = psi(j,ibnd)
     enddo
     call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)
     call stop_clock ('firstfft')
     !
     !   product with the potential vrs = (vltot+vr) on the smooth grid
     !
     do j = 1, nrxxs
        psic(j) = psic(j) * vrs(j,current_spin)
     enddo
     !
     !   back to reciprocal space
     !
     call start_clock ('secondfft')
     call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -2)
     !
     !   addition to the total product
     !
     do j = 1, n
        hpsi(j,ibnd) = hpsi(j,ibnd) + psic(nls(igk(j)))
     enddo
     call stop_clock ('secondfft')
  enddo
  !
  !  Here the product with the non local potential V_NL psi
  !
  if (nkb.gt.0) call add_vuspsi (lda, n, m, psi, hpsi)
  call s_psi (lda, n, m, psi, spsi)
  call stop_clock ('h_psi')
  return
end subroutine h_psi

