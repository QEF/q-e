!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file 'License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
PROGRAM mypp2upf
  !---------------------------------------------------------------------
  !
  !     Convert a pseudopotential written in a user-supplied format
  !     to unified pseudopotential format - sample program
  !
  IMPLICIT NONE
  CHARACTER(len=256) filein, fileout
  !
  !
  CALL get_file ( filein )
  OPEN (unit = 1, file = filein, status = 'old', form = 'formatted')
  CALL read_mypp(1)
  CLOSE (1)

  ! convert variables read from user-supplied format into those needed
  ! by the upf format - add missing quantities

  CALL convert_mypp

  fileout=trim(filein)//'.UPF'
  PRINT '(''Output PP file in UPF format :  '',a)', fileout

  OPEN(unit=2,file=fileout,status='unknown',form='formatted')
  CALL write_upf(2)
  CLOSE (unit=2)

STOP
20 WRITE (6,'("mypp2upf: error reading pseudopotential file name")')
   STOP
END PROGRAM mypp2upf

MODULE mypp
  !
  ! All variables read from user-supplied file format
  ! Must have distinct names from variables in the "upf" module
  !
END MODULE mypp
!
!     ----------------------------------------------------------
SUBROUTINE read_mypp(iunps)
  !     ----------------------------------------------------------
  !
  USE mypp
  IMPLICIT NONE
  INTEGER :: iunps
  !
  !     ----------------------------------------------------------
  WRITE (6,'(a)') 'Pseudopotential successfully read'
  !     ----------------------------------------------------------
  !
  RETURN
100 WRITE (6,'("read_mypp: error reading pseudopotential file")')
    STOP
END SUBROUTINE read_mypp

!     ----------------------------------------------------------
SUBROUTINE convert_mypp
  !     ----------------------------------------------------------
  !
  USE mypp
  USE upf
  IMPLICIT NONE
  !     ----------------------------------------------------------
  WRITE (6,'(a)') 'Pseudopotential successfully converted'
  !     ----------------------------------------------------------
  RETURN
END SUBROUTINE convert_mypp

