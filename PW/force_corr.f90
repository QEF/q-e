!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------

subroutine force_corr (forcescc)
  !-----------------------------------------------------------------------
  !   This routine calculates the force term vanishing at full
  !     self-consistency. It follows the suggestion of Chan-Bohnen-Ho
  !     (PRB 47, 4771 (1993)). The true charge density is approximated
  !     by means of a free atom superposition.
  !     (alessio f.)
  ! Uses superposition of atomic charges contained in the array rho_at
  ! and already set in readin-readvan
  !
#include "machine.h"
  use pwcom
  implicit none
  !
  real(kind=DP) :: forcescc (3, nat)
  !
  real(kind=DP), allocatable :: rhocgnt (:), aux (:)
  ! work space
  real(kind=DP) ::  gx, arg, fact
  ! temp factors
  integer :: ir, isup, isdw, ig, igl0, nt, na, ipol
  ! counters
  !
  ! vnew is V_out - V_in, psic is the temp space
  !
  if (nspin == 1) then
     psic(:) = vnew (:, 1)
  else
     isup = 1
     isdw = 2
     psic(:) = (vnew (:, isup) + vnew (:, isdw)) * 0.5d0
  end if
  !
  allocate ( aux(ndm), rhocgnt(ngl) )

  forcescc(:,:) = 0.d0

  call cft3 (psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)

  if (gamma_only) then
     fact = 2.d0
  else
     fact = 1.d0
  end if

  do nt = 1, ntyp
     !
     ! Here we compute the G.ne.0 term
     !
     do ig = gstart, ngl
        gx = sqrt (gl (ig) ) * tpiba
        do ir = 1, msh (nt)
           if (r (ir, nt) .lt.1.0d-8) then
              aux (ir) = rho_at (ir, nt)
           else
              aux (ir) = rho_at (ir, nt) * sin(gx*r(ir,nt))/(r(ir,nt)*gx)
           endif
        enddo
        call simpson (msh (nt), aux, rab (1, nt), rhocgnt (ig) )
     enddo
     do na = 1, nat
        if (nt.eq.ityp (na) ) then
           do ig = gstart, ngm
              arg = (g (1, ig) * tau (1, na) + g (2, ig) * tau (2, na) &
                   + g (3, ig) * tau (3, na) ) * tpi
              do ipol = 1, 3
                 forcescc (ipol, na) = forcescc (ipol, na) + fact * &
                      rhocgnt (igtongl(ig) ) * DCMPLX(sin(arg),cos(arg)) * &
                      g(ipol,ig) * tpiba * conjg(psic(nl(ig)))
              enddo
           enddo
        endif
     enddo
  enddo
#ifdef __PARA
  call reduce (3 * nat, forcescc)
#endif
  deallocate ( aux, rhocgnt )

  return
end subroutine force_corr

