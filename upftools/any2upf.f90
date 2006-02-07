!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file 'License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
#include "f_defs.h"
!---------------------------------------------------------------------
program mypp2upf  
  !---------------------------------------------------------------------
  !
  !     Convert a pseudopotential written in a user-supplied format 
  !     to unified pseudopotential format - sample program
  !
  implicit none
  character(len=256) filein, fileout
  !
  !
  call get_file ( filein )
  open (unit = 1, file = filein, status = 'old', form = 'formatted')
  call read_mypp(1)
  close (1)

  ! convert variables read from user-supplied format into those needed
  ! by the upf format - add missing quantities

  call convert_mypp

  fileout=trim(filein)//'.UPF'
  print '(''Output PP file in UPF format :  '',a)', fileout

  open(unit=2,file=fileout,status='unknown',form='formatted')
  call write_upf(2)
  close (unit=2)

stop
20 write (6,'("mypp2upf: error reading pseudopotential file name")')
   stop
end program mypp2upf

module mypp
  !
  ! All variables read from user-supplied file format
  ! Must have distinct names from variables in the "upf" module
  !
end module mypp
! 
!     ----------------------------------------------------------
subroutine read_mypp(iunps)
  !     ----------------------------------------------------------
  ! 
  use mypp
  implicit none
  integer :: iunps
  !
  !     ----------------------------------------------------------
  write (6,'(a)') 'Pseudopotential successfully read'
  !     ----------------------------------------------------------
  !
  return
100 write (6,'("read_mypp: error reading pseudopotential file")')
    stop
end subroutine read_mypp

!     ----------------------------------------------------------
subroutine convert_mypp
  !     ----------------------------------------------------------
  !
  use mypp
  use upf
  implicit none
  !     ----------------------------------------------------------
  write (6,'(a)') 'Pseudopotential successfully converted'
  !     ----------------------------------------------------------
  return
end subroutine convert_mypp

