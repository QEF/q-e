!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine h_psi_nc (lda, n, m, psi, hpsi)
  !-----------------------------------------------------------------------
  !
  !     This routine computes the product of the Hamiltonian
  !     matrix with m wavefunctions contained in psi
  ! input:
  !     lda   leading dimension of arrays psi, hpsi
  !     n     true dimension of psi, hpsi
  !     m     number of states psi
  !     psi
  ! output:
  !     hpsi  H*psi
  !
  use uspp, only: vkb, nkb
  use wvfct, only: igk, g2kin
  use gsmooth, only : nls, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, nrxxs
  use ldaU, only : lda_plus_u
  use lsda_mod, only : current_spin
  use scf, only: vrs
  use becmod
  use wavefunctions_module, only: psic_nc
  use noncollin_module, only: noncolin, npol
  implicit none
  !
  integer :: lda, n, m
  complex(DP) :: psi(lda,npol,m), hpsi(lda,npol,m),&
                      sup, sdwn
  !
  integer :: ibnd,j,ipol
  ! counters
  call start_clock ('h_psi')
  call start_clock ('init')
  call ccalbec_nc (nkb, lda, n, npol, m, becp_nc, vkb, psi)
  !
  ! Here we apply the kinetic energy (k+G)^2 psi
  !
  do ibnd = 1, m
     do ipol = 1, npol
        do j = 1, n
           hpsi (j, ipol, ibnd) = g2kin (j) * psi (j, ipol, ibnd)
        enddo
     enddo
  enddo

  call stop_clock ('init')
  !
  ! Here we add the Hubbard potential times psi
  !
  if (lda_plus_u) call vhpsi_nc (lda, n, m, psi(1,1,1), hpsi(1,1,1))
  !
  ! the local potential V_Loc psi. First the psi in real space
  !
  do ibnd = 1, m
     call start_clock ('firstfft')
     psic_nc = (0.d0,0.d0)
     do ipol=1,npol
        do j = 1, n
           psic_nc(nls(igk(j)),ipol) = psi(j,ipol,ibnd)
        enddo
        call cft3s (psic_nc(1,ipol), nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)
     enddo
     call stop_clock ('firstfft')
     !
     !   product with the potential vrs = (vltot+vr) on the smooth grid
     !
     if (noncolin) then
        do j=1, nrxxs
           sup = psic_nc(j,1) * (vrs(j,1)+vrs(j,4)) + &
                 psic_nc(j,2) * (vrs(j,2)-(0.d0,1.d0)*vrs(j,3))
           sdwn = psic_nc(j,2) * (vrs(j,1)-vrs(j,4)) + &
                  psic_nc(j,1) * (vrs(j,2)+(0.d0,1.d0)*vrs(j,3))
           psic_nc(j,1)=sup
           psic_nc(j,2)=sdwn
        end do
     else
        do j = 1, nrxxs
           psic_nc(j,1) = psic_nc(j,1) * vrs(j,current_spin)
        enddo
     endif
     !
     !   back to reciprocal space
     !
     call start_clock ('secondfft')
     do ipol=1,npol
        call cft3s (psic_nc(1,ipol), nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -2)
     enddo
     !
     !   addition to the total product
     !
     do ipol=1,npol
        do j = 1, n
           hpsi(j,ipol,ibnd) = hpsi(j,ipol,ibnd) + psic_nc(nls(igk(j)),ipol)
        enddo
     enddo
     call stop_clock ('secondfft')
  enddo
  !
  !  Here the product with the non local potential V_NL psi
  !
  if (nkb.gt.0) call add_vuspsi_nc (lda, n, m, psi, hpsi(1,1,1))
  call stop_clock ('h_psi')
  return
end subroutine h_psi_nc
