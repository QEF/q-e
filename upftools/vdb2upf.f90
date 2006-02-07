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
program vdb2upf
  !---------------------------------------------------------------------
  !
  !     Convert a pseudopotential written in Vanderbilt format
  !     (formatted) to unified pseudopotential format
  !
  implicit none
  character(len=256) filein, fileout
  !
  !
  call get_file ( filein )
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
20 write (6,'("vdb2upf: error reading pseudopotential file name")')
   stop
end program vdb2upf
