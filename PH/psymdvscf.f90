!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------

subroutine psymdvscf (nper, irr, dvtosym)
  !-----------------------------------------------------------------------
  !  p-symmetrize the charge density.
  !
#include "machine.h"
#ifdef PARA

  use pwcom
  use parameters, only : DP
  use phcom
  use para
  implicit none

  integer :: nper, irr
  ! the number of perturbations
  ! the representation under consideration


  complex(kind=DP) :: dvtosym (nrxx, nspin, nper)
  ! the potential to symmetrize
  !-local variable

  integer :: i, is, iper, npp0

  complex(kind=DP), allocatable :: ddvtosym (:,:,:)
  ! the potential to symm


  if (nsymq.eq.1.and. (.not.minus_q) ) return
  call start_clock ('psymdvscf')

  allocate (ddvtosym ( nrx1 * nrx2 * nrx3, nspin, nper))    
  npp0 = 1
  do i = 1, me-1
     npp0 = npp0 + npp (i) * ncplane

  enddo
  do iper = 1, nper
     do is = 1, nspin
        call cgather_sym (dvtosym (1, is, iper), ddvtosym (1, is, iper) )
     enddo

  enddo

  call symdvscf (nper, irr, ddvtosym)
  do iper = 1, nper
     do is = 1, nspin
        call ZCOPY (npp (me) * ncplane, ddvtosym (npp0, is, iper), &
             1, dvtosym (1, is, iper), 1)
     enddo

  enddo
  deallocate (ddvtosym)

  call stop_clock ('psymdvscf')
#endif
  return
end subroutine psymdvscf
