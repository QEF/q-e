!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------

subroutine davcio_drho (drho, lrec, iunit, nrec, isw)
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
#endif
  implicit none
  integer :: iunit, lrec, nrec, isw
  complex(kind=DP) :: drho (nrxx, nspin)
#ifdef __PARA
  include 'mpif.h'
  !
  ! local variables
  !

  integer :: root, errcode, itmp, proc, is

  complex(kind=DP), allocatable  :: ddrho (:,:)

  allocate (ddrho( nrx1 * nrx2 * nrx3 , nspin))    
  if (isw.eq.1) then
     !
     ! First task of the first pool is the only task allowed to write
     ! the file
     !
     do is = 1, nspin
        call cgather_sym (drho (1, is), ddrho (1, is) )
     enddo
     root = 0
     call MPI_barrier (MPI_COMM_WORLD, errcode)
     call errore ('davcio_drho', 'at barrier', errcode)

     if (me.eq.1) call davcio (ddrho, lrec, iunit, nrec, + 1)
  elseif (isw.lt.0) then
     !
     ! First task of the pool reads ddrho, and broadcasts to all the
     ! processors of the pool
     !
     if (me.eq.1) call davcio (ddrho, lrec, iunit, nrec, - 1)
     call broadcast (2 * nrx1 * nrx2 * nrx3 * nspin, ddrho)
     !
     ! Distributes ddrho between between the tasks of the pool
     !
     itmp = 1
     do proc = 1, me-1
        itmp = itmp + ncplane * npp (proc)
     enddo
     do is = 1, nspin
        call setv (2 * nrxx, 0.d0, drho (1, is), 1)
        call ZCOPY (ncplane * npp (me), ddrho (itmp, is), 1, drho (1, &
             is), 1)
     enddo

  endif
  deallocate(ddrho)
#else
  call davcio (drho, lrec, iunit, nrec, isw)
#endif
  return
end subroutine davcio_drho
