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
  character(len=75) filein, fileout
  logical exst
  integer :: i,ierr,ilen
  integer, external :: iargc
  !
  i = iargc ()  
  if (i.eq.0) then  
5    print '(''  input PP file in your format > '',$)'  
     read (5, '(a)', end = 20, err = 20) filein
     exst=filein.ne.' '
     if (.not. exst) go to 5  
     inquire (file=filein,exist=exst)
     if(.not.exst) go to 5
  elseif (i.eq.1) then  
     call getarg(1, filein)  
  else  
     print '(''   usage: mypp2upf  [input file] '')'  
     stop  
  endif

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
20 call errore ('mypp2upf', 'Reading pseudo file name ', 1)

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
100 call errore ('read_mypp', 'Reading pseudo file', 100 )  

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

