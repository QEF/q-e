!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine macro
  !----------------------------------------------------------------------
  !
#include "machine.h"
  use pwcom
  use cgcom
  !
  implicit none
  integer:: kpoint, ipol
  character(len=7) :: filbar
  logical :: here
  !
  do kpoint=1,nks
     ! NB: this version works only for nks = 1 !
     do ipol=1,3
        write(filbar,'(''filbar'',i1)') ipol
        iubar=ipol
        call seqopn (iubar,filbar,'unformatted',here)
!!!            if (.not.here) then
        ! calculate x * psi  (if not already done)
        call setv(2*nbnd*npwx,0.d0,dvpsi,1)
!!!            else
        ! otherwise restart from x * psi that is present on from file
!!!               read(iubar) dvpsi
!!!            end if
        call dvpsi_e(kpoint,ipol)
        ! write x * psi
        rewind(iubar)
        write(iubar) dvpsi
        close(unit=iubar,status='keep')
     end do
  end do
  !
  return
end subroutine macro
