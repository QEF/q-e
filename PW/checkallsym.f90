!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine checkallsym (nsym, s, nat, tau, ityp, at, bg, nr1, nr2, &
     nr3, irt, ftau, alat, omega)
  !-----------------------------------------------------------------------
  !     given a crystal group this routine checks that the actual
  !     atomic positions and bravais lattice vectors are compatible with
  !     it. Used in relaxation/MD runs to check that atomic motion is
  !     consistent with assumed symmetry.
  !
#include "f_defs.h"
  USE kinds
  implicit none
  !
  !     first the I/O variables
  !
  integer :: nsym, s (3, 3, 48), nat, ityp (nat), nr1, nr2, nr3, &
       irt (48, nat), ftau (3, 48)
  ! input: order of the parent group
  ! input: symmetry operations of parent group
  ! input: number of atoms in the unit cell
  ! input: species of each atom in the unit ce
  !
  ! input: dimensions of the FFT mesh
  !
  ! input: tell which atom is equivalent to
  ! a given atom under a given symmetry operat
  ! input: fractional translation associated t
  ! k-th symmetry operation, referred to the
  ! crystal axis and in units at/nr (0-nr-1)

  real(DP) :: tau (3, nat), at (3, 3), bg (3, 3), alat, omega
  ! input: cartesian coordinates of the atoms
  ! input: basis of the real-space lattice
  ! input:  "   "   "  reciprocal-space lattic
  !
  !    here the local variables
  !
  integer :: na, kpol, isym, i, j, k, l
  ! counter on atoms
  ! counter on polarizations
  ! counter on rotations
  ! counters

  logical :: loksym (48)

  real(DP) :: sx (3, 3), ft (3), ps
  real(DP) , allocatable :: xau(:,:), rau(:,:)
  ! symmetry operation in cartesian axis
  ! atomic coordinate referred to the crystal
  ! crystal coordinates of a rotated atom
  ! actual fractionary translation
  ! scalar product
  real(DP), parameter :: eps = 1.0d-7
  ! a small number

  external checksym

  allocate (xau( 3 , nat))    
  allocate (rau( 3 , nat))    
  !
  !     check that s(i,j, isym) is an orthogonal operation
  !

  do isym = 1, nsym
     do i = 1, 3
        do j = 1, 3
           sx (i, j) = 0.d0
           do l = 1, 3
              do k = 1, 3
                 sx (i, j) = sx (i, j) + bg (i, l) * s (l, k, isym) * at (j, k)
              enddo
           enddo
        enddo

     enddo
     do i = 1, 3
        do j = 1, 3
           ps = sx (1, i) * sx (1, j) + &
                sx (2, i) * sx (2, j) + &
                sx (3, i) * sx (3, j)
           if (i.eq.j) ps = ps - 1.d0
           if (abs (ps) .gt.eps) &
                call errore ('checkallsym', 'not orthogonal operation', isym)
        enddo

     enddo

  enddo
  !
  !     Compute the coordinates of each atom in the basis of the direct la
  !     vectors
  !
  do na = 1, nat
     do kpol = 1, 3
        xau (kpol, na) = bg (1, kpol) * tau (1, na) + &
                         bg (2, kpol) * tau (2, na) + &
                         bg (3, kpol) * tau (3, na)
     enddo

  enddo

  do isym = 1, nsym
     do na = 1, nat
        ! generate the
        do kpol = 1, 3
           ! coordinates
           rau (kpol, na) = s (1, kpol, isym) * xau (1, na) + &
                            s (2, kpol, isym) * xau (2, na) + &
                            s (3, kpol, isym) * xau (3, na)
           ! atom rotated
           ! present symm
           ! operation
        enddo

     enddo

     ft (1) = ftau (1, isym) / DBLE (nr1)
     ft (2) = ftau (2, isym) / DBLE (nr2)
     ft (3) = ftau (3, isym) / DBLE (nr3)

     call checksym (isym, nat, ityp, xau, rau, ft, loksym, irt)
  enddo
  !
  !   deallocate work space
  !
  deallocate(rau)
  deallocate(xau)
  !
  !
  do isym = 1,nsym
     if (.not.loksym (isym) ) call errore ('checkallsym', &
          'the following symmetry operation is not satisfied  ', -isym)
  end do
  if (ANY (.not.loksym (1:nsym) ) ) then
      call symmetrize_at(nsym, s, nat, tau, ityp, at, bg, nr1, nr2, &
                         nr3, irt, ftau, alat, omega)
      call errore ('checkallsym', &
           'some of the original symmetry operations not satisfied ',1)
  end if
  
  !
  return
end subroutine checkallsym

