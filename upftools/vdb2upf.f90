!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file 'License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
program vdb2upf
  !---------------------------------------------------------------------
  !
  !     Convert a pseudopotential written in Vanderbilt format
  !     (formatted) to unified pseudopotential format
  !
  implicit none
  character(len=75) filein, fileout
  logical exst
  integer :: i,j
#ifdef ABSOFT
#define iargc  iargc_
#define getarg getarg_
#endif
  integer, external :: iargc  
  !
  i = iargc ()  
  if (i.eq.0) then  
5    print '(''  Input PP file in formatted Vanderbilt format > '',$)'  
     read (5, '(a)', end = 20, err = 20) filein
     exst=filein.ne.' '
     if (.not. exst) go to 5  
     inquire (file=filein,exist=exst)
     if(.not.exst) go to 5
  elseif (i.eq.1) then  
#ifdef T3D
     call pxfgetarg (1, filein, i, j)  
#else
     call getarg (1, filein)  
#endif
  else  
     print '(''   usage: vdb2upf  [input file] '')'  
     stop
  end if

  open(unit=1,file=filein,status='old',form='formatted')
  call read_vdb(1)
  close (unit=1)

  ! convert variables read from Vanderbilt format into those needed
  ! by the upf format - add missing quantities

  call convert_uspp

  fileout=trim(filein)//'.UPF'
  print '(''Output PP file in UPF format :  '',a)', fileout

  open(unit=2,file=fileout,status='unknown',form='formatted')
  call write_upf(2)
  close (unit=2)

  stop
20 call errore ('vdb2upf', 'Reading pseudo file name ', 1)
end program vdb2upf
