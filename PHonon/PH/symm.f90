!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine symm(phi, u, xq, s, isym, rtau, irt, at, bg, nat)
  !-----------------------------------------------------------------------
  !! This routine symmetrizes the matrix of electron-phonon coefficients
  !! written in the basis of the modes.
  !
  USE kinds,     ONLY: DP
  USE constants, ONLY: tpi
  !
  implicit none
  !
  integer, intent (in) :: nat
  !! input: the number of atoms
  integer, intent (in) :: s(3,3,48)
  !! input: the symmetry matrices
  integer, intent (in) :: irt(48,nat)
  !! input: the rotated of each atom
  integer, intent (in) :: isym
  !! input: the small group of q
  real(DP), intent(in) :: xq(3)
  !! input: the coordinates of q
  real(DP), intent(in) :: rtau(3,48,nat)
  !! input: the R associated at each r
  real(DP), intent(in) :: at(3,3)
  !! input: direct lattice vectors
  real(DP), intent(in) :: bg(3,3)
  !! input: reciprocal lattice vectors
  complex(DP), intent(in) :: u(3*nat,3*nat)
  !! input: patterns
  complex(DP), intent(inout) :: phi(3*nat,3*nat)
  !! input: matrix to be symmetrized , output: symmetrized matrix
  !
  ! ... local variables
  !
  integer :: i, j, icart, jcart, na, nb, mu, nu, sna, snb, &
       ipol, jpol, lpol, kpol
  ! counters
  real(DP) :: arg
  !
  complex(DP) :: fase, work, phi0(3*nat,3*nat), phi1(3,3,nat,nat), phi2(3,3,nat,nat)
  ! workspace
  !
  ! First we transform to cartesian coordinates
  !
  phi0 = matmul(u, matmul(phi, conjg(transpose(u))))
  do i = 1, 3 * nat
     na = (i - 1) / 3 + 1
     icart = i - 3 * (na - 1)
     do j = 1, 3 * nat
        nb = (j - 1) / 3 + 1
        jcart = j - 3 * (nb - 1)
        phi1(icart,jcart,na,nb) = phi0(i,j)
     enddo
  enddo
  !
  ! Then we transform to crystal axis
  !
  do na = 1, nat
     do nb = 1, nat
        call trntnsc (phi1(1,1,na,nb), at, bg, - 1)
     enddo
  enddo
  !
  !   And we symmetrize in this basis
  !
  do na = 1, nat
     do nb = 1, nat
        sna = irt (isym, na)
        snb = irt (isym, nb)
        arg = 0.d0
        do ipol = 1, 3
           arg = arg + (xq(ipol)*(rtau(ipol,isym,na) - rtau(ipol,isym,nb)))
        enddo
        arg = arg * tpi
        fase = CMPLX(DCOS (arg), DSIN (arg) ,kind=DP)
        do ipol = 1, 3
           do jpol = 1, 3
              phi2(ipol,jpol,na,nb) = (0.0d0,0.0d0)
              do kpol = 1, 3
                 do lpol = 1, 3
                    phi2(ipol,jpol,na,nb) = phi2(ipol,jpol,na,nb) +  &
                         s(ipol,kpol,isym) * s(jpol,lpol,isym) * &
                         phi1(kpol,lpol,sna,snb) * fase
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  !
  !  Back to cartesian coordinates
  !
  do na = 1, nat
     do nb = 1, nat
        call trntnsc (phi2 (1, 1, na, nb), at, bg, + 1)
     enddo
  enddo
  !
  !  rewrite as an array with dimensions 3nat x 3nat
  !
  do i = 1, 3 * nat
     na = (i - 1) / 3 + 1
     icart = i - 3 * (na - 1)
     do j = 1, 3 * nat
        nb = (j - 1) / 3 + 1
        jcart = j - 3 * (nb - 1)
        phi (i, j) = phi2 (icart, jcart, na, nb)
     enddo
  enddo
  !
  return
end subroutine symm

