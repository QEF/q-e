!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine struc_fact (nat, tau, ntyp, ityp, ngm, g, bg, nr1, nr2, &
     nr3, strf, eigts1, eigts2, eigts3)
  !----------------------------------------------------------------------
  !
  !   calculate the structure factors for each type of atoms in the unit
  !   cell
  !
  USE kinds
  USE constants, ONLY : tpi
  USE gvect,     ONLY : mill
  implicit none
  !
  !   Here the dummy variables
  !

  integer :: nat, ntyp, ityp (nat), ngm, nr1, nr2, nr3
  ! input: the number of atom in the unit cel
  ! input: the number of atom types
  ! input: for each atom gives the type
  ! input: the number of G vectors
  ! input: fft dimension along x
  ! input: fft dimension along y
  ! input: fft dimension along z

  real(DP) :: bg (3, 3), tau (3, nat), g (3, ngm)
  ! input: reciprocal crystal basis vectors
  ! input: the positions of the atoms in the c
  ! input: the coordinates of the g vectors

  complex(DP) :: strf (ngm, ntyp),        &
                      eigts1 ( -nr1:nr1, nat), &
                      eigts2 ( -nr2:nr2, nat), &
                      eigts3 ( -nr3:nr3, nat)
  ! output: the structure factor
  !
  ! output: the phases e^{-iG\tau_s}
  !
  !
  !    here the local variables
  !
  integer :: nt, na, ng, n1, n2, n3, ipol
  ! counter over atom type
  ! counter over atoms
  ! counter over G vectors
  ! counter over fft dimension along x
  ! counter over fft dimension along y
  ! counter over fft dimension along z
  ! counter over polarizations

  real(DP) :: arg, bgtau (3)
  ! the argument of the exponent
  ! scalar product of bg and tau

  complex(DP) :: z

!$omp parallel do schedule(static) default(none) &
!$omp private(ipol, arg, bgtau, n1, n2, n3) &
!$omp shared(nr1, nr2, nr3, nat, eigts1, eigts2, eigts3, bg, tau)
  do na = 1, nat
     do ipol = 1, 3
        bgtau (ipol) = bg (1, ipol) * tau (1, na) + &
                       bg (2, ipol) * tau (2, na) + &
                       bg (3, ipol) * tau (3, na)
     enddo
     do n1 = - nr1, nr1
        arg = tpi * n1 * bgtau (1)
        eigts1 (n1, na) = CMPLX(cos (arg), - sin (arg) ,kind=DP)
     enddo
     do n2 = - nr2, nr2
        arg = tpi * n2 * bgtau (2)
        eigts2 (n2, na) = CMPLX(cos (arg), - sin (arg) ,kind=DP)
     enddo
     do n3 = - nr3, nr3
        arg = tpi * n3 * bgtau (3)
        eigts3 (n3, na) = CMPLX(cos (arg), - sin (arg) ,kind=DP)
     enddo
  enddo
!$omp end parallel do

  strf(:,:) = (0.d0,0.d0)
  do nt = 1, ntyp
!$omp parallel do schedule(static) default(none) &
!$omp private(na, z) &
!$omp shared(ngm, nat, nt, ityp, strf, eigts1, eigts2, eigts3, mill)
    do ng = 1, ngm
      z = (0.d0, 0.d0)
      do na = 1, nat
        if (ityp(na).eq.nt) then
          ! compute exp(-i G tau) using product of three exponents
          z = z + eigts1(mill(1, ng), na) * eigts2(mill(2, ng), na) * eigts3(mill(3, ng), na)
        endif
      enddo
      strf (ng, nt) = z
    enddo
!$omp end parallel do
  enddo !nt

  return
end subroutine struc_fact