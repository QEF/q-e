!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------

subroutine psyme (dvtosym)
  !-----------------------------------------------------------------------
  !  p-symmetrize the charge density.
  !
#include "f_defs.h"
#ifdef __PARA

  use pwcom
  USE kinds, only : DP
  use phcom
  use para
  implicit none

  complex(kind=DP) :: dvtosym (nrxx, nspin, 3)
  ! the potential to symmetrize
  !-local variable

  integer :: i, is, iper, npp0

  complex(kind=DP), allocatable :: ddvtosym (:,:,:)
  ! the potential to symmet

  allocate (ddvtosym ( nrx1 * nrx2 * nrx3, nspin, 3))    
  npp0 = 0
  do i = 1, me-1
     npp0 = npp0 + npp (i)
  enddo

  npp0 = npp0 * ncplane+1
  do iper = 1, 3
     do is = 1, nspin
        call cgather_sym (dvtosym (1, is, iper), ddvtosym (1, is, iper) )
     enddo

  enddo

  call syme (ddvtosym)
  do iper = 1, 3
     do is = 1, nspin
        call ZCOPY (npp (me) * ncplane, ddvtosym (npp0, is, iper), &
             1, dvtosym (1, is, iper), 1)
     enddo

  enddo

  deallocate (ddvtosym)
#endif
  return
end subroutine psyme
