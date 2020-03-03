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

  integer :: nat, ntyp, ityp (nat), ngm, nr1, nr2, nr3, nmax
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
  integer :: nt, na, ng, n1, n2, n3, ipol, i
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

  complex(DP), allocatable :: eigts(:,:,:)
  ! structure factors for x,y,z components

  complex(DP) :: z

  nmax = max(nr1, max(nr2, nr3))
  allocate(eigts(-nmax:nmax, 3, nat))

  ! compute structure factors exp(i*G_j*B*tau) for j=1,2,3 where
  !  G_j are lattice coordinates of G-vectors,
  !  B - reciprocal lattice vectors and
  !  tau - Cartesian coorinates of atoms
!$omp parallel do schedule(static) default(none) &
!$omp private(ipol,arg,i) &
!$omp shared(nr1,nr2,nr3,nat,eigts,eigts1,eigts2,eigts3,bg,tau,nmax)
  do na = 1, nat
    do ipol = 1, 3
      arg = bg (1, ipol) * tau (1, na) + &
            bg (2, ipol) * tau (2, na) + &
            bg (3, ipol) * tau (3, na)
      do i = -nmax, nmax
        eigts(i, ipol, na) = CMPLX(cos(tpi * arg * i), -sin(tpi * arg * i), kind=DP)
      enddo
    enddo !ipol
    do i = -nr1, nr1
       eigts1(i, na) = eigts(i, 1, na)
    enddo
    do i = -nr2, nr2
       eigts2(i, na) = eigts(i, 2, na)
    enddo
    do i = -nr3, nr3
       eigts3(i, na) = eigts(i, 3, na)
    enddo
  enddo !na
!$omp end parallel do

  strf(:,:) = (0.d0,0.d0)
  do nt = 1, ntyp
!$omp parallel do schedule(static) default(none) &
!$omp private(na, z) &
!$omp shared(ngm, nat, nt, ityp, strf, eigts, mill)
    do ng = 1, ngm
      z = (0.d0, 0.d0)
      do na = 1, nat
        if (ityp(na).eq.nt) then
          ! compute exp(i G tau) using product of three exponents
          z = z + eigts(mill(1, ng), 1, na) * eigts(mill(2, ng), 2, na) * eigts(mill(3, ng), 3, na)
        endif
      enddo
      strf (ng, nt) = z
    enddo
!$omp end parallel do
  enddo !nt

  deallocate(eigts)

  return
end subroutine struc_fact