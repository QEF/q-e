!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file 'License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
PROGRAM uspp2upf
  !---------------------------------------------------------------------
  !
  !     Convert a pseudopotential written in Vanderbilt format
  !     (unformatted) to unified pseudopotential format
  !
  IMPLICIT NONE
  CHARACTER(len=256) filein, fileout
  !
  !
  IF ( trim(filein) == ' ') &
       CALL errore ('uspp2upf', 'usage: uspp2upf "file-to-be-converted"', 1)
  CALL get_file ( filein )
  OPEN(unit=1,file=filein,status='old',form='unformatted')
  CALL read_uspp(1)
  CLOSE (unit=1)

  ! convert variables read from Vanderbilt format into those needed
  ! by the upf format - add missing quantities

  CALL convert_uspp

  fileout=trim(filein)//'.UPF'
  PRINT '(''Output PP file in UPF format :  '',a)', fileout

  OPEN(unit=2,file=fileout,status='unknown',form='formatted')
  CALL write_upf_v1(2)
  CLOSE (unit=2)

  STOP
20 WRITE (6,'("uspp2upf: error reading pseudopotential file name")')
   STOP
END PROGRAM uspp2upf

