!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine davcio_drho2 (drho, lrec, iunit, nrec, isw)
  !-----------------------------------------------------------------------
  !
  ! reads/writes variation of the charge with respect to a perturbation
  ! on a file.
  ! isw = +1 : gathers data from the nodes and writes on a single file
  ! isw = -1 : reads data from a single file and distributes them
  !
#include "machine.h"
  use pwcom
  USE kinds, only : DP
  use phcom
#ifdef __PARA
  use para
  USE io_global,     ONLY : ionode_id
  USE mp,            ONLY : mp_bcast  
#endif
  implicit none
#ifdef __PARA
  include 'mpif.h'
#endif
  integer :: iunit, lrec, nrec, isw
  complex(kind=DP) :: drho (nrxx)
#ifdef __PARA
  !
  ! local variables
  !

  integer :: root, errcode, itmp, proc

  complex(kind=DP), allocatable :: ddrho (:)

  allocate (ddrho( nrx1 * nrx2 * nrx3 ))    

  if (isw == 1) then
     !
     ! First task of the pool gathers and writes in the file
     !
     call cgather_sym (drho, ddrho)
     root = 0
     call MPI_barrier (MPI_COMM_POOL, errcode)
     call errore ('davcio_drho2', 'at barrier', errcode)
     if (me.eq.1) call davcio (ddrho, lrec, iunit, nrec, + 1)
  elseif (isw < 0) then
     !
     ! First task of the pool reads ddrho, and broadcasts to all the
     ! processors of the pool
     !
     if (me == 1) call davcio (ddrho, lrec, iunit, nrec, - 1)
     call mp_bcast( ddrho, ionode_id )
     !
     ! Distributes ddrho between between the tasks of the pool
     !
     itmp = 1
     do proc = 1, me-1
        itmp = itmp + ncplane * npp (proc)
     enddo
     drho (:) = (0.d0, 0.d0)
     call ZCOPY (ncplane * npp (me), ddrho (itmp), 1, drho, 1)
  endif

  deallocate(ddrho)
#else
  call davcio (drho, lrec, iunit, nrec, isw)
#endif
  return
end subroutine davcio_drho2
