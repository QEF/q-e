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
  !   Read the variation of the charge and 
  !   calculates the local part of the scf potential
  !
#include "f_defs.h"

  use pwcom
  use kinds, only : DP
  use phcom
  implicit none

  complex(kind=DP) :: dvscfs (nrxxs,3)
  complex(kind=DP) , allocatable :: derho (:,:)
  integer :: ipl
  !  counter on the polarizations

  allocate (derho ( nrxx, nspin))

  if ( fildrho.eq.' ') call errore ('set_dvscf','where is fildrho?',1)
  !
  do ipl = 1, 3
     !
     ! read from file the variation of the charge
     !
     call davcio_drho (derho (1, 1), lrdrho, iudrho, ipl, -1)
     !
     ! Calculates the local part of the scf potential
     !
     call dv_of_drho (0, derho (1, 1), .false.)
     !
     if (doublegrid) then
        call cinterpolate (derho (1, 1), dvscfs (1, ipl), -1)
     else
        call ZCOPY (nrxx, derho (1, 1), 1, dvscfs (1, ipl), 1)
     endif
  end do

  deallocate (derho)

  return
end subroutine set_dvscf
