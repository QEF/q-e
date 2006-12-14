!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
subroutine do_shift_ew (alat, nat, ntyp, ityp, zv, delta_zv, at, bg, tau, &
     omega, g, gg, ngm, gcutm, gstart, gamma_only, shift_ion)
  !-----------------------------------------------------------------------
  !
  ! Calculates Ewald energy with both G- and R-space terms.
  ! Determines optimal alpha. Should hopefully work for any structure.
  !
  !
  USE kinds, ONLY : DP
  USE constants, ONLY : tpi, e2
  implicit none
  !
  !   first the dummy variables
  !

  integer :: nat, ntyp, ityp (nat), ngm, gstart
  ! input: number of atoms in the unit cell
  ! input: number of different types of atoms
  ! input: the type of each atom
  ! input: number of plane waves for G sum
  ! input: first non-zero G vector

  logical :: gamma_only

  real(DP) :: tau (3, nat), g (3, ngm), gg (ngm), zv (ntyp), &
       at (3, 3), bg (3, 3), omega, alat, gcutm, delta_zv(ntyp),  &
       shift_ion(nat)
  ! input: the positions of the atoms in the cell
  ! input: the coordinates of G vectors
  ! input: the square moduli of G vectors
  ! input: the charge of each type of atoms
  ! input: the direct lattice vectors
  ! input: the reciprocal lattice vectors
  ! input: the volume of the unit cell
  ! input: lattice parameter
  ! input: cut-off of g vectors
  real(DP) :: ewald
  ! output: the ewald energy
  !
  !    here the local variables
  !
  integer, parameter :: mxr = 50
  ! the maximum number of R vectors included in r
  integer :: ng, nr, na, nb, nt, nrm, ipol
  ! counter over reciprocal G vectors
  ! counter over direct vectors
  ! counter on atoms
  ! counter on atoms
  ! counter on atomic types
  ! number of R vectors included in r sum
  ! counter on polarization

  real(DP) :: charge, tpiba2, ewaldg, ewaldr, dtau (3), alpha, &
       r (3, mxr), r2 (mxr), rmax, rr, upperbound, fact, arg
  ! total ionic charge in the cell
  ! length in reciprocal space
  ! ewald energy computed in reciprocal space
  ! ewald energy computed in real space
  ! the difference tau_s - tau_s'
  ! alpha term in ewald sum
  ! input of the rgen routine ( not used here )
  ! the square modulus of R_j-tau_s-tau_s'
  ! the maximum radius to consider real space sum
  ! buffer variable
  ! used to optimize alpha
  complex(DP), allocatable :: rhon(:)
  real(DP), external :: erfc

  allocate (rhon(ngm))

  shift_ion(:) = 0.d0

  tpiba2 = (tpi / alat) **2
  charge = 0.d0
  do na = 1, nat
     charge = charge+zv (ityp (na) )
  enddo
  alpha = 2.9d0
100 alpha = alpha - 0.1d0
  !
  ! choose alpha in order to have convergence in the sum over G
  ! upperbound is a safe upper bound for the error in the sum over G
  !
  if (alpha.le.0.d0) call errore ('ewald', 'optimal alpha not found', 1)
  upperbound = 2.d0 * charge**2 * sqrt (2.d0 * alpha / tpi) * erfc ( &
       sqrt (tpiba2 * gcutm / 4.d0 / alpha) )
  if (upperbound.gt.1.0d-7) goto 100
  !
  ! G-space sum here.
  ! Determine if this processor contains G=0 and set the constant term
  !
  if (gstart==2) then
     do na =1,nat
        shift_ion(na) = - charge * delta_zv(ityp(na)) /alpha/ 4.0d0
     end do
  endif
  if (gamma_only) then
     fact = 2.d0
  else
     fact = 1.d0
  end if
  do ng = gstart, ngm
     rhon(ng) = (0.d0, 0.d0)
     do na =1, nat
        arg = (g (1, ng) * tau (1, na) + &
               g (2, ng) * tau (2, na) + &
               g (3, ng) * tau (3, na) ) * tpi
        rhon(ng) = rhon(ng) + zv (ityp(na)) * CMPLX (cos (arg), -sin (arg))
     enddo
  end do
  do na=1,nat
     do ng=gstart, ngm
        arg = (g (1, ng) * tau (1, na) + g (2, ng) * tau (2, na) &
             + g (3, ng) * tau (3, na) ) * tpi
        shift_ion(na) = shift_ion(na) + fact * delta_zv(ityp(na)) * &
                        CONJG(rhon(ng)) * CMPLX (cos (arg), -sin (arg)) * &
                        exp ( -gg(ng)*tpiba2/alpha/4.d0) / gg(ng)/tpiba2
     enddo
  enddo
  shift_ion(:) = 2.d0 * tpi / omega * shift_ion(:)
  !
  !  Here add the other constant term
  !
  if (gstart.eq.2) then
     do na = 1, nat
        shift_ion(na) = shift_ion(na) - &
                        zv (ityp (na) ) * delta_zv(ityp(na)) *  &
                        sqrt (8.d0/tpi*alpha)
     enddo
  endif
  !
  ! R-space sum here (only for the processor that contains G=0)
  !
  if (gstart.eq.2) then
     rmax = 4.d0 / sqrt (alpha) / alat
     !
     ! with this choice terms up to ZiZj*erfc(4) are counted (erfc(4)=2x10^-8
     !
     do na = 1, nat
        do nb = 1, nat
           do ipol = 1, 3
              dtau (ipol) = tau (ipol, na) - tau (ipol, nb)
           enddo
           !
           ! generates nearest-neighbors shells
           !
           call rgen (dtau, rmax, mxr, at, bg, r, r2, nrm)
           !
           ! and sum to the real space part
           !
           do nr = 1, nrm
              rr = sqrt (r2 (nr) ) * alat
              shift_ion(na) = shift_ion(na) + &
                       delta_zv(ityp(na)) * zv (ityp (nb) ) *  &
                       erfc ( sqrt (alpha) * rr) / rr
           enddo
        enddo
     enddo
  endif

  shift_ion(:) = e2 * shift_ion(:)

#ifdef __PARA
  call reduce (nat, shift_ion)
#endif
  deallocate (rhon)
  return
end subroutine do_shift_ew

