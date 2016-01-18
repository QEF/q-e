!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------

subroutine dpsi_corr (evcq, psidvpsi_x, ik, ikq, nu)
  !-----------------------------------------------------------------------
  ! Used in the metallic case.
  ! If dpsi common variable contains the projection on the conduction
  ! states of the first variation of a wavefunction at a given k-point,
  ! this routine corrects dpsi in such a way that the density matrix
  ! is given by:   Sum_{k,nu} 2 * | dpsi > < psi |
  !
  USE kinds, only : DP
  use pwcom
  use qpoint, ONLY: npwq
  use phcom
  use d3com

  implicit none
  integer :: ik, ikq, nu, ibnd, jbnd
  ! index of the k-point under consideration
  ! index of the corresponding k+q point
  ! mode under consideration
  ! counter on bands
  ! counter on bands

  real (DP) :: wfshift, wgauss, w0gauss, deltae, wg1, wg2, wwg
  ! the shift coefficent for the wave function
  ! function computing the theta function
  ! function computing the derivative of theta
  ! difference of energy
  ! weight for metals
  ! weight for metals
  ! weight for metals

  complex (DP) :: evcq (npwx, nbnd), psidvpsi_x (nbnd, nbnd), &
       psidvpsi
  ! k+q point wavefunction
  ! < psi_{k+q} | V(q) | psi_k >
  !
  ! Multiplies dpsi by the theta function
  !
  do ibnd = 1, nbnd
     wg1 = wgauss ( (ef - et (ibnd, ik) ) / degauss, ngauss)
     call dscal (2 * npwq, wg1, dpsi (1, ibnd), 1)
  enddo
  !
  ! Adds to dpsi the term containing the valence wavefunctions
  !
  do ibnd = 1, nbnd
     do jbnd = 1, nbnd
        deltae = et (ibnd, ik) - et (jbnd, ikq)
        if (abs (deltae) .gt.1.0d-5) then
           wg1 = wgauss ( (ef - et (ibnd, ik) ) / degauss, ngauss)
           wg2 = wgauss ( (ef - et (jbnd, ikq) ) / degauss, ngauss)
           wwg = (wg1 - wg2) / deltae
        else
           wwg = - w0gauss ( (ef - et (ibnd, ik) ) / degauss, ngauss) &
                / degauss
        endif
        psidvpsi = 0.5d0 * wwg * psidvpsi_x (jbnd, ibnd)
        call zaxpy (npwq, psidvpsi, evcq (1, jbnd), 1, dpsi (1, ibnd), &
             1)
     enddo
  enddo
  !
  ! If necessary corrects dpsi with a term depending on FermiEnergy shift
  !
  if (ik.eq.ikq) then
     do ibnd = 1, nbnd_occ (ik)
        wfshift = 0.5d0 * ef_sh (nu) * w0gauss ( (ef - et (ibnd, ik) ) &
             / degauss, ngauss) / degauss
        call daxpy (2 * npw, wfshift, evcq (1, ibnd), 1, dpsi (1, ibnd) &
             , 1)
     enddo

  endif
  return

end subroutine dpsi_corr
