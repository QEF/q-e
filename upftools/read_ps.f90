!
! Copyright (C) 2001-2002 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
PROGRAM read_ps
  !---------------------------------------------------------------------
  !
  !  Read pseudopotentials in the Unified Pseudopotential Format (UPF)
  !
  IMPLICIT NONE
  INTEGER :: is, ios, iunps = 4
  CHARACTER (len=256) :: filein
  !
  is = 0
10 WRITE(*,'("  Input PP file # ",i2," in UPF format > ")',advance="NO") is+1
  READ (5, '(a)', end = 20, err = 20) filein
  OPEN(unit=iunps,file=filein,status='old',form='formatted',iostat=ios)
  IF (ios/=0) STOP
  is = is + 1
  CALL read_pseudo(is, iunps)
  CLOSE (unit=iunps)
  GOTO 10
20 STOP
END PROGRAM read_ps
