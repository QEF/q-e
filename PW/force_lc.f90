!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine force_lc (nat, tau, ityp, alat, omega, ngm, ngl, &
     igtongl, nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, g, rho, nl, &
     nspin, gstart, gamma_only, vloc, forcelc)
  !----------------------------------------------------------------------
  !
#include "f_defs.h"
  USE kinds
  USE constants, ONLY : tpi
  implicit none
  !
  !   first the dummy variables
  !
  integer :: nat, ngm, nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, nspin, &
       ngl, gstart, igtongl (ngm), nl (ngm), ityp (nat)
  ! input: the number of atoms in the cell
  ! input: the number of G vectors
  ! input: FFT dimensions
  ! input: number of spin polarizations
  ! input: the number of shells
  ! input: correspondence G <-> shell of G
  ! input: the correspondence fft mesh <-> G vec
  ! input: the types of atoms

  logical :: gamma_only

  real(DP) :: tau (3, nat), g (3, ngm), vloc (ngl, * ), &
       rho (nrxx, nspin), alat, omega
  ! input: the coordinates of the atoms
  ! input: the coordinates of G vectors
  ! input: the local potential
  ! input: the valence charge
  ! input: the length measure
  ! input: the volume of the cell

  real(DP) :: forcelc (3, nat)
  ! output: the local forces on atoms

  integer :: ipol, ig, na
  ! counter on polarizations
  ! counter on G vectors
  ! counter on atoms

  real(DP), allocatable :: aux (:,:)
  ! auxiliary space for FFT
  real(DP) :: arg, fact
  !
  ! contribution to the force from the local part of the bare potential
  ! F_loc = Omega \Sum_G n*(G) d V_loc(G)/d R_i
  !
  allocate (aux(2, nrxx))
  aux(1,:) = rho(:,1)
  if (nspin.eq.2) aux(1,:) = aux(1,:) + rho(:,2)
  aux(2,:) = 0.d0
  call cft3 (aux, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)
  !
  !    aux contains now  n(G)
  !
  if (gamma_only) then
     fact = 2.d0
  else
     fact = 1.d0
  end if
  do na = 1, nat
     do ipol = 1, 3
        forcelc (ipol, na) = 0.d0
     enddo
     ! contribution from G=0 is zero
     do ig = gstart, ngm
        arg = (g (1, ig) * tau (1, na) + g (2, ig) * tau (2, na) + &
               g (3, ig) * tau (3, na) ) * tpi
        do ipol = 1, 3
           forcelc (ipol, na) = forcelc (ipol, na) + &
                g (ipol, ig) * vloc (igtongl (ig), ityp (na) ) * &
                (sin (arg) * aux(1,nl(ig)) + cos (arg) * aux(2,nl(ig)) )
        enddo
     enddo
     do ipol = 1, 3
        forcelc (ipol, na) = fact * forcelc (ipol, na) * omega * tpi / alat
     enddo
  enddo
#ifdef __PARA
  call reduce (3 * nat, forcelc)
#endif
  deallocate (aux)
  return
end subroutine force_lc
