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
subroutine psymd0rho (nper, irr, dvtosym)
  !-----------------------------------------------------------------------
  !  p-symmetrize the charge density.
  !

#ifdef __PARA
  USE kinds, ONLY : DP
  USE ions_base, ONLY : nat
  USE pwcom
  USE phcom
  USE d3com
  USE mp_global, ONLY : me_pool
  USE fft_base,  ONLY : dfftp

  IMPLICIT NONE

  integer :: nper, irr
  ! the number of perturbations
  ! the representation under consideration

  complex (DP) :: dvtosym (nrxx, nper)
  ! the potential to symmetrize

  ! local variables

  integer :: i, iper, npp0
  complex (DP),pointer  :: ddvtosym (:,:)
  ! the potential to symmetrize

!  if (nsymq.eq.1.and. (.not.minus_q) ) return

  call start_clock ('psymd0rho')

  allocate ( ddvtosym( nrx1 * nrx2 * nrx3, nper))    
  npp0 = 0
  do i = 1, me_pool
     npp0 = npp0 + dfftp%npp (i)
  enddo

  npp0 = npp0 * dfftp%nnp + 1
  do iper = 1, nper
     call cgather_sym (dvtosym (1, iper), ddvtosym (1, iper) )
  enddo

  call symd0rho (max_irr_dim, nper, irr, ddvtosym, s, ftau, nsymg0, irgq, tg0, &
       nat, nr1, nr2, nr3, nrx1, nrx2, nrx3)
  do iper = 1, nper
     call ZCOPY (dfftp%npp (me_pool+1) * dfftp%nnp, ddvtosym (npp0, iper), 1, dvtosym &
          (1, iper), 1)
  enddo
  deallocate(ddvtosym)

  call stop_clock ('psymd0rho')
#endif
  return
end subroutine psymd0rho
