!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------
subroutine do_voronoy
  !--------------------------------------------------------------------
  !
  !  Calculates charges on atoms by dividing the space into Voronoy polyhe
  !  A Voronoy polyhedron around a given atom is defined as the region
  !  of space that is closer to that atom than to the others
  !
  !  Note that this is a very rough way to associate charges to atoms
  !  and that it is well defined only if all atoms are of the same type!
  !
  !  On input: nr1big, nr2bug, nr3big are FFT grid dimensions larger than
  !  the original ones (the largest, the more accurate the integration)
  !
#include "machine.h"
  use parameters
  use pwcom
  use fft_scalar, only: good_fft_dimension
  implicit none
  integer :: nr1big, nr2big, nr3big, nrx1big
  integer :: n, i, j, ng, na, plot_num
  real(kind=DP) :: total_charge, rhodum
  integer, allocatable :: nlbig (:)
  real(kind=DP), allocatable :: partial_charge (:)
  complex(kind=DP), allocatable :: rhobig (:)
  character(len=80) :: filename
  !
  !
  print '(" Input file > ",$)'
  read (5, '(a)') filename
  !
  ! read file header and allocate objects
  !
  call plot_io (filename, title, nrx1, nrx2, nrx3, nr1, nr2, nr3, &
       nat, ntyp, ibrav, celldm, at, gcutm, dual, ecutwfc, plot_num, &
       atm, ityp, zv, tau, rhodum, 0)
  !
  allocate(tau(3, nat))
  allocate(ityp(nat))
  alat = celldm (1)
  tpiba = 2.d0 * pi / alat
  tpiba2 = tpiba**2
  if (ibrav.gt.0) then
    call latgen (ibrav, celldm, at (1, 1), at (1, 2), at (1, 3), omega )
    at = at / alat !  bring at in units of alat
  end if
  call recips (at (1, 1), at (1, 2), at (1, 3), bg (1, 1), bg (1, 2) &
       , bg (1, 3) )
  call volume (alat, at (1, 1), at (1, 2), at (1, 3), omega)
  call set_fft_dim
  nrxx = nrx1 * nrx2 * nrx3
  nspin = 1
  call allocate_fft
  !
  ! read data from file
  !
  call plot_io (filename, title, nrx1, nrx2, nrx3, nr1, nr2, nr3, &
       nat, ntyp, ibrav, celldm, at, gcutm, dual, ecutwfc, plot_num, &
       atm, ityp, zv, tau, rho, - 1)
  !
  ! calculate g-vectors
  !
  call ggen
  !
  ! interpolate on a larger grid
  !
  print '(" nr1big, nr2big, nr3big (min: ",3i4," ) > ",$)', nr1, &
       nr2, nr3

  read (5, * ) nr1big, nr2big, nr3big

  nrx1big = good_fft_dimension (nr1big)

  allocate (nlbig( ngm))    

  call get_fftindex (g, ngm, nr1big, nr2big, nr3big, nrx1big, &
       nr2big, nr3big, at, nlbig)

  allocate (rhobig(nrx1big * nr2big * nr3big))    

  call rhor_to_rhobig (ngm, nr1, nr2, nr3, nrx1, nl, rho, nr1big, &
       nr2big, nr3big, nrx1big, nlbig, rhobig)

  call calculate_partial_charges (nat, tau, at, bg, nrx1big, nr1big, &
       nr2big, nr3big, rhobig, partial_charge)

  write (6, '(" nr1, nr2, nr3 = ",3i4)') nr1big, nr2big, nr3big
  total_charge = 0.0
  do na = 1, nat
     partial_charge (na) = partial_charge (na) * omega / (nr1big * &
          nr2big * nr3big)
     write (6, '(" Atom # ",i2," tau = ",3f8.4, &
          &       "  charge: ",f8.4) ') na,  (tau (i, na) , i = 1, 3) , &
          partial_charge (na)
     total_charge = total_charge+partial_charge (na)
  enddo

  write (6, '(" Check: total charge = ",f8.4)') total_charge
  stop

end subroutine do_voronoy
!-----------------------------------------------------------------------

subroutine get_fftindex (g, ngm, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
     at, nl)
  !-----------------------------------------------------------------------
  use parameters, only: DP
  implicit none
  integer :: nr1, nr2, nr3, nrx1, nrx2, nrx3, ngm
  real(kind=DP) :: g (3, ngm), at (3, 3)
  integer :: nl (ngm)
  integer :: n1, n2, n3, ng
  !
  do ng = 1, ngm
     n1 = nint (g (1, ng) * at (1, 1) + g (2, ng) * at (2, 1) + g (3, &
          ng) * at (3, 1) ) + 1
     if (n1.lt.1) n1 = n1 + nr1
     n2 = nint (g (1, ng) * at (1, 2) + g (2, ng) * at (2, 2) + g (3, &
          ng) * at (3, 2) ) + 1
     if (n2.lt.1) n2 = n2 + nr2
     n3 = nint (g (1, ng) * at (1, 3) + g (2, ng) * at (2, 3) + g (3, &
          ng) * at (3, 3) ) + 1
     if (n3.lt.1) n3 = n3 + nr3
     if (n1.le.nr1.and.n2.le.nr2.and.n3.le.nr3) then
        nl (ng) = n1 + (n2 - 1) * nrx1 + (n3 - 1) * nrx1 * nrx2
     else
        stop ' Mesh too small '
     endif
  enddo
  !
  return

end subroutine get_fftindex
!-----------------------------------------------------------------------
subroutine rhor_to_rhobig (ngm, nr1, nr2, nr3, nrx1, nl, rho, &
     nr1big, nr2big, nr3big, nrx1big, nlbig, rhobig)
  !-----------------------------------------------------------------------
#include "machine.h"
  use parameters, only: DP
  implicit none
  integer :: ngm, nr1, nr2, nr3, nrx1, nl (ngm), nr1big, nr2big, &
       nr3big, nrx1big, nlbig (ngm)
  real(kind=DP) :: rho (nrx1 * nr2 * nr3)
  complex(kind=DP) :: rhobig (nrx1big * nr2big * nr3big)
  complex(kind=DP), allocatable :: aux (:)
  integer :: nr, ng

  allocate (aux(nrx1 * nr2 * nr3))    
  do nr = 1, nrx1 * nr2 * nr3
     aux (nr) = rho (nr)
  enddo

  call cft3 (aux, nr1, nr2, nr3, nrx1, nr2, nr3, - 1)

  do nr = 1, nrx1big * nr2big * nr3big
     rhobig (nr) = 0.0
  enddo
  do ng = 1, ngm
     rhobig (nlbig (ng) ) = aux (nl (ng) )
  enddo

  deallocate(aux)

  call cft3 (rhobig, nr1big, nr2big, nr3big, nrx1big, nr2big, &
       nr3big, 1)
  return

end subroutine rhor_to_rhobig
!-----------------------------------------------------------------------
subroutine calculate_partial_charges (nat, tau, at, bg, nrx1big, &
     nr1big, nr2big, nr3big, rhobig, partial_charge)
  !-----------------------------------------------------------------------
#include "machine.h"
  use parameters, only: DP
  implicit none
  integer :: nat, nrx1big, nr1big, nr2big, nr3big
  real(kind=DP) :: at (3, 3), bg (3, 3), tau (3, nat), partial_charge ( &
       nat)

  complex(kind=DP) :: rhobig (nrx1big, nr2big, nr3big)
  integer :: atom (nat), equidistant, n1, n2, n3, na, i
  real(kind=DP) :: r (3), dr (3), distance (nat), relative_distance
  !
  !
  do na = 1, nat
     partial_charge (na) = 0.0
  enddo
  do n3 = 1, nr3big
     do n2 = 1, nr2big
        do n1 = 1, nr1big
           ! r is the position of point (n1,n2,n3) on the grid
           do i = 1, 3
              r (i) = (n1 - 1) * at (i, 1) / nr1big + (n2 - 1) * at (i, 2) &
                   / nr2big + (n3 - 1) * at (i, 3) / nr3big
           enddo
           do na = 1, nat
              ! dr is the distance between r and this atom, in crystal axis
              do i = 1, 3
                 dr (i) = (r (1) - tau (1, na) ) * bg (1, i) + (r (2) - tau (2, na) &
                      ) * bg (2, i) + (r (3) - tau (3, na) ) * bg (3, i)
                 ! this brings dr back into the unit cell
                 dr (i) = dr (i) - nint (dr (i) )
              enddo
              ! distance is in cartesioan axis
              distance (na) = (dr (1) * at (1, 1) + dr (2) * at (1, 2) + dr (3) &
                   * at (1, 3) ) **2 + (dr (1) * at (2, 1) + dr (2) * at (2, 2) &
                   + dr (3) * at (2, 3) ) **2 + (dr (1) * at (3, 1) + dr (2) * at (3, &
                   2) + dr (3) * at (3, 3) ) **2
           enddo
           ! initialization needed by hpsort
           atom (1) = 0
           ! sort distances in increasing order
           call hpsort (nat, distance, atom)
           ! find if some atoms are equidistant
           equidistant = 1
           do na = 2, nat
              relative_distance = distance (na) - distance (1)
              if (relative_distance.gt.1.0e-4) goto 10
              equidistant = equidistant + 1
           enddo
10         continue
           ! the charge is assigned to the closest atom or shared among equidistant
           do na = 1, equidistant
              partial_charge (atom (na) ) = partial_charge (atom (na) ) + real ( &
                   rhobig (n1, n2, n3) ) / equidistant
           enddo
        enddo
     enddo
  enddo
  return
end subroutine calculate_partial_charges
