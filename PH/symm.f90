!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine symm(dyn1, dyn2, u, xq, s, isym, rtau, irt, at, bg, nat)
  !-----------------------------------------------------------------------
  !
  !    This routine symmetrize the dynamical matrix written in the basis
  !    of the modes
  !
  !
  implicit none
  integer :: nat, s (3,3,48), irt (48, nat), k, isym
  ! input: the number of atoms
  ! input: the symmetry matrices
  ! input: the rotated of each atom
  ! input: the small group of q
  ! input: the inverse of each matrix
  ! input: the order of the small gro
  ! input: the symmetry q -> -q+G

  real(kind(1.d0)) :: xq (3), rtau (3, 48, nat), at (3, 3), bg (3, 3)
  ! input: the coordinates of q
  ! input: the R associated at each r
  ! input: direct lattice vectors
  ! input: reciprocal lattice vectors

  complex(kind((1.d0,1.0d0))) :: dyn1(3*nat,3*nat), dyn2(3*nat,3*nat), &
                                 u(3*nat,3*nat)
  ! inp/out: matrix to symmetrize
  ! input: the patterns

  integer :: i, j, icart, jcart, na, nb, mu, nu
  ! counter on modes
  ! counter on modes
  ! counter on cartesian coordinates
  ! counter on cartesian coordinates
  ! counter on atoms
  ! counter on atoms
  ! counter on modes
  ! counter on modes

  complex(kind((1.d0,1.0d0))) :: work, wrk(3,3), &
            phi1(3,3,nat,nat), phi2(3,3,nat,nat)
  ! auxiliary variable
  ! auxiliary variable
  ! the dynamical matrix
  !
  ! First we transform in the cartesian coordinates
  !

  do i = 1, 3 * nat
     na = (i - 1) / 3 + 1
     icart = i - 3 * (na - 1)
     do j = 1, 3 * nat
        nb = (j - 1) / 3 + 1
        jcart = j - 3 * (nb - 1)
        work = (0.d0, 0.d0)
        do mu = 1, 3 * nat
           do nu = 1, 3 * nat
              work = work + u(i,mu) * dyn1(mu,nu) * conjg(u(j,nu))
           enddo
        enddo
        phi1(icart,jcart,na,nb) = work
     enddo
  enddo
  !
  ! Then we transform to the crystal axis
  !
  do na = 1, nat
     do nb = 1, nat
        call trntnsc (phi1(1,1,na,nb), at, bg, - 1)
     enddo
  enddo
  !
  !   And we symmetrize in this basis
  !
  call symm_q (xq, phi1, phi2, s, isym, rtau, irt, nat) 
  !
  !  Back to cartesian coordinates
  !
  do na = 1, nat
     do nb = 1, nat
        call trntnsc (phi2 (1, 1, na, nb), at, bg, + 1)
     enddo
  enddo
  !
  !  rewrite the dynamical matrix on the array dyn with dimension 3nat x 3
  !
  do i = 1, 3 * nat
     na = (i - 1) / 3 + 1
     icart = i - 3 * (na - 1)
     do j = 1, 3 * nat
        nb = (j - 1) / 3 + 1
        jcart = j - 3 * (nb - 1)
        dyn2 (i, j) = phi2 (icart, jcart, na, nb)
     enddo

  enddo
  return
end subroutine symm

subroutine symm_q (xq, phi1, phi2, s, isym, rtau, irt, nat) 
  !
  !     This routine receives as input an unsymmetrized dynamical
  !     matrix expressed on the crystal axes and imposes the symmetry
  !     of the small group of q. Furthermore it imposes also the symmetry
  !     q -> -q+G if present.
  !
  implicit none
  !
  !    The dummy variables
  !
  integer :: nat, s(3,3,48), irt(48,nat), isym  
  ! input: the number of atoms
  ! input: the symmetry matrices
  ! input: the rotated of each vector
  ! input: the rotation needed in the mapping trom BZ-> IBZ
  real(kind(1.d0)) :: xq(3), rtau(3,48,nat)
  ! input: the q point
  ! input: the R associated at each t

  complex(kind((1.d0,1.0d0))) :: phi1(3,3,nat,nat), phi2(3,3,nat,nat)
  ! inp/out: the matrix to symmetrize
  !
  !     One parameter
  !
  real(kind(1.d0)) :: tpi
  parameter (tpi = 2.0d0 * 3.14159265358979d0)
  !
  !   and the local variables
  !

  integer :: sna, snb, irot, na, nb, ipol, jpol, lpol, kpol
  ! the rotated of the a atom
  ! the rotated of the b atom
  ! counter on rotations
  ! counter on atoms
  ! counter on atoms
  ! counter on polarizations
  ! counter on polarizations
  ! counter on polarizations
  ! counter on polarizations

  real(kind(1.d0)) :: arg
  ! the argument of the phase

  complex(kind((1.d0,1.0d0))) :: fase 
  ! the phase factor
  !
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
           fase = DCMPLX (cos (arg), sin (arg) )
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

  return
end subroutine symm_q

