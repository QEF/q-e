!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
      subroutine set_dvscf (dvscfs)
  !-----------------------------------------------------------------------
  !
  !   Symmetrizes the variation of the charge
  !   Calculates the local part of the scf potential
  !
#include "f_defs.h"

  use pwcom
  use kinds, only : DP
  use phcom
  implicit none

  complex(kind=DP) :: dvscfs (nrxxs,3)
  complex(kind=DP) , allocatable :: derho (:,:,:)
  integer :: ipl
  !  counter on the polarizations

  allocate (derho ( nrxx, nspin, 3))

  if ( fildrho.eq.' ') call errore ('set_dvscf','unexpected',1)
!
! Symmetrized the variation of the charge; and writes it on the file
!
  do ipl = 1, 3
     call davcio_drho (derho (1, 1, ipl), lrdrho, iudrho, ipl, -1)
  enddo
#ifdef __PARA
  call psyme (derho)
#else
  call syme (derho)
#endif
  do ipl = 1, 3
     call davcio_drho (derho (1, 1, ipl), lrdrho, iudrho, ipl, +1)
  enddo
!
! Calculates the local part of the scf potential
!
  do ipl = 1, 3
     call dv_of_drho (0, derho (1, 1, ipl), .false.)
  enddo
#ifdef __PARA
  call psyme (derho)
#else
  call syme (derho)
#endif

  if (doublegrid) then
     do ipl = 1, 3
        call cinterpolate (derho (1, 1, ipl), dvscfs (1, ipl), -1)
     enddo
  else
     do ipl = 1, 3
        call ZCOPY (nrxx, derho (1, 1, ipl), 1, dvscfs (1, ipl), 1)
     enddo
  endif

  deallocate (derho)

  return
end subroutine set_dvscf
