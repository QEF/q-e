!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------
subroutine subtract_vxc (raux, spin_component)
  !--------------------------------------------------------------------
  !
  !     This routine computes Exchange and Correlation
  !     potential which corresponds to a given charge density (contained
  !     in the common pwcom) and subtract it from the H+xc potential
  !
  !
#include "machine.h"
  use pwcom
  use funct

  implicit none
  real(kind=DP) :: raux (nrxx)
  ! input: the H+xc_up  potential

  integer :: spin_component
  ! input: the spin to plot in the gga c

  real(kind=DP) :: rhox, arhox, zeta, etxc0, vtxc0, ex, ec, vx (2), vc (2)
  ! the total charge in each point
  ! the absolute value of the charge
  ! the absolute value of the charge
  ! dummy for gga routine
  ! local exchange energy
  ! local correlation energy
  ! local exchange potential
  ! local correlation potential

  integer :: ir, is, negz
  ! counter on mesh points
  ! counter on spin polarizations
  ! number of negative points
  !
  ! exchange-correlation potential and energy
  !
  if (nspin.eq.1) then
     call DCOPY (nrxx, vr (1, nspin), 1, raux, 1)
     !
     ! LDA case
     !
     do ir = 1, nrxx
        rhox = rho (ir, nspin) + rho_core (ir)
        arhox = abs (rhox)
        if (arhox.gt.1.d-30) then
           call xc (arhox, ex, ec, vx, vc)
           raux (ir) = vr (ir, nspin) - e2 * (vx (1) + vc (1) )
        endif
     enddo
  else
     !
     ! LSDA case
     !
     if (spin_component.ne.1.and.spin_component.ne.2) call error ( &
          'subtract vxc', 'spin_component not allowed', 1)
     is = spin_component
     negz = 0
     call DCOPY (nrxx, vr (1, is), 1, raux, 1)
     do ir = 1, nrxx
        rhox = rho (ir, 1) + rho (ir, 2) + rho_core (ir)
        arhox = abs (rhox)
        if (arhox.gt.1.d-30) then
           zeta = (rho (ir, 1) - rho (ir, 2) ) / arhox
           if (abs (zeta) .gt.1.d0) negz = negz + 1
           call xc_spin (arhox, zeta, ex, ec, vx (1), vx (2), vc (1), &
                vc (2) )
           raux (ir) = vr (ir, is) - e2 * (vx (is) + vc (is) )
        endif
     enddo
     if (negz.gt.0) write (6, '(/,4x, " npt with |zeta| > 1:",i8, &
          &                     ", npt tot ",i8, ",",f10.2," % " )') &
          negz, nr1 * nr2 * nr3, float (negz * 100)  / float (nr1 * nr2 * nr3)
  endif
  !
  if (igcx.ne.0.or.igcc.ne.0) then
     call DSCAL (nrxx, - 1.d0, raux, 1)
     call gradcorr (rho, rho_core, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
          nrxx, nl, ngm, g, alat, omega, e2, etxc0, vtxc0, raux, nspin)
     call DSCAL (nrxx, - 1.d0, raux, 1)
  endif
  !
  return
end subroutine subtract_vxc
