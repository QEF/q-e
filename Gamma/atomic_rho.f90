!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine atomic_rho (rhoa, nspina)
  !-----------------------------------------------------------------------
  ! This routine calculates rhoa as the superposition of atomic charges.
  !
  ! nspina is the number of spin components to be calculated
  !
  ! if nspina = 1 the total atomic charge density is calculated
  ! if nspina = 2 the spin up and spin down atomic charge densities are
  !               calculated assuming an uniform atomic spin-polarization
  !               equal to starting_magnetization(nt)
  !
  ! NB: nspina may not be equal to nspin because in some cases (as in upda
  ! the total charge only could be needed, even in a LSDA calculation.
  !
  !
#include "machine.h"

  use pwcom
  use gamma
  implicit none
  integer :: nspina
  ! the number of spin polarizations

  real(kind=DP) :: rhoa (nrxx, nspina), rhoneg, rhorea, rhoima, gx
  real(kind=DP), allocatable :: rhocgnt (:), aux (:)
  ! the output atomic charge
  ! negative charge
  ! real charge
  ! imaginary charge
  ! the modulus of G
  ! the value of the integral
  ! the integrand function

  complex(kind=DP), allocatable :: rhocg (:,:)
  ! auxiliary var: charge dens. in G space

  integer :: ir, is, ig, igl, igl0, nt
  ! counter on mesh points
  ! counter on spin polarizations
  ! counter on G vectors
  ! counter on G vectors shells
  ! index of first shell with G != 0
  ! counter on atom types

  !
  ! superposition of atomic charges contained in the array rho_at and
  ! already set in readin-readvan
  !
  allocate (rhocg(  ngm, nspina))    
  allocate (aux( ndm))    
  allocate (rhocgnt( ngl))    

  rhoa(:,:) = 0.d0
  rhocg(:,:) = (0.d0,0.d0)

  do nt = 1, ntyp
     !
     ! Here we compute the G=0 term
     !
     if (gl (1) .lt.1.0d-8) then
        do ir = 1, msh (nt)
           aux (ir) = rho_at (ir, nt)
        enddo
        call simpson (msh (nt), aux, rab (1, nt), rhocgnt (1) )
        igl0 = 2
     else
        igl0 = 1
     endif
     !
     ! Here we compute the G<>0 term
     !
     do igl = igl0, ngl
        gx = sqrt (gl (igl) ) * tpiba
        do ir = 1, msh (nt)
           if (r (ir, nt) .lt.1.0d-8) then
              aux (ir) = rho_at (ir, nt)
           else
              aux (ir) = rho_at (ir, nt) * sin (gx * r (ir, nt) ) / &
                   (r (ir, nt) * gx)
           endif
        enddo
        call simpson (msh (nt), aux, rab (1, nt), rhocgnt (igl) )
     enddo
     !
     ! we compute the 3D atomic charge in reciprocal space
     !
     if (nspina.eq.1) then
        do ig = 1, ngm
           rhocg (ig, 1) = rhocg (ig, 1) + strf (ig, nt) * &
                rhocgnt ( igtongl (ig) ) / omega
        enddo
     else

        do ig = 1, ngm
           rhocg (ig, 1) = rhocg (ig, 1) + 0.5d0 * (1.d0 + &
                starting_magnetization (nt) ) * strf (ig, nt) * &
                rhocgnt ( igtongl (ig) ) / omega

           rhocg (ig, 2) = rhocg (ig, 2) + 0.5d0 * (1.d0 - &
                starting_magnetization (nt) ) * strf (ig, nt) * &
                rhocgnt ( igtongl (ig) ) / omega
        enddo
     endif
  enddo
  deallocate (rhocgnt)

  deallocate (aux)
  do is = 1, nspina
     !
     ! and we return to real space
     !
     psic(:) = (0.d0, 0.d0)
     do ig = 1, ngm
        psic (nl (ig) ) = rhocg (ig, is)
        psic (nlm(ig) ) = conjg( rhocg (ig, is)  )
     enddo
     call cft3 (psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
     !
     ! we check that everything is correct
     !
     rhoneg = 0.d0
     rhoima = 0.d0
     do ir = 1, nrxx
        rhorea = DREAL (psic (ir) )
        rhoneg = rhoneg + min (0.d0, rhorea)
        rhoima = rhoima + abs (DIMAG (psic (ir) ) )
        rhoa (ir, is) = rhorea
     enddo
     rhoneg = rhoneg / (nr1 * nr2 * nr3)
     rhoima = rhoima / (nr1 * nr2 * nr3)
#ifdef __PARA
     call reduce (1, rhoneg)
     call reduce (1, rhoima)
#endif
     if (rhoneg.lt. - 1.0d-4.or.rhoima.gt.1.0d-4) &
          write (6,'(/"  Warning: negative or imaginary starting charge ",&
          &2f12.6,i3)') rhoneg, rhoima, is
  enddo

  deallocate (rhocg)
  return
end subroutine atomic_rho

