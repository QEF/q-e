!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine io_pot (iop, filename, pot, nc)
  !-----------------------------------------------------------------------
  !
  !     This routine reads (iop=-1) or write (iop=1) the potential
  !     in real space onto a file
  !
#include "machine.h"
  use pwcom
#ifdef __PARA
  use para
#endif
  implicit none
  integer :: iop, nc, ic
  ! option: write if +1, read if -1
  ! number of components and index for them
  character (len=*) :: filename

  real(kind=DP) :: pot (nrxx, nc)
#ifdef __PARA
  real(kind=DP), allocatable :: allv (:,:)
#endif

  logical :: exst
#ifdef __PARA
  if (me.eq.1) allocate( allv(nrx1*nrx2*nrx3, nc) )
  !
  ! On writing: gather the potential on the first node of each pool
  !
  if (iop.eq.1) then
     do ic = 1, nc
        call gather (pot (1, ic), allv (1, ic) )
     enddo
  endif
  if (me.eq.1.and.mypool.eq.1) then
     !
     ! Only the first node of the first pool reads or writes the file
     !
     call seqopn (4, filename, 'unformatted', exst)
     if (iop.eq.1) then
        write (4, err = 10) allv
     else
        read (4, err = 20) allv
     endif
     close (unit = 4)
  endif
  !
  ! On reading: copy the potential on the first node  of all pools
  !             scatter the potential on all nodes of each pool
  !
  if (iop.eq. - 1) then
     if (me.eq.1) call poolbcast (nrx1 * nrx2 * nrx3 * nc, allv)
     do ic = 1, nc
        call scatter (allv (1, ic), pot (1, ic) )
     enddo
  endif
  if (me.eq.1) deallocate(allv)
#else
  call seqopn (4, filename, 'unformatted', exst)
  if (iop.eq.1) then
     write (4, err = 10) pot
  else
     read (4, err = 20) pot
  endif
  close (unit = 4)
#endif
  return
10 call errore ('io_pot', 'error writing '//filename, 1)
  return
20 call errore ('io_pot', 'error reading '//filename, 2)
  return
end subroutine io_pot

