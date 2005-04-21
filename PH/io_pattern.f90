!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
SUBROUTINE io_pattern (fildrho,nirr,npert,u,iflag)
!---------------------------------------------------------------------
  !
  USE ions_base, ONLY : nat
  USE io_global,  ONLY : stdout
  USE pwcom
  USE kinds, ONLY : DP
  
  IMPLICIT NONE
!
!   the i/o variables first
!
  INTEGER :: nirr, npert(3*nat), iflag
  COMPLEX(kind=DP) :: u(3*nat,3*nat)
  CHARACTER (len=*)  ::  fildrho   ! name of the file
  CHARACTER (len=256)::  filname   ! complete name of the file
!
!   here the local variables
!
  INTEGER :: i,iunit
  LOGICAL :: exst

  IF (ABS(iflag).NE.1) CALL errore('io_pattern','wrong iflag',1+ABS(iflag))

  iunit = 4
  filname = TRIM(fildrho) //".pat"
  CALL seqopn(iunit,filname,'formatted',exst)

  IF (iflag.GT.0) THEN
     WRITE( stdout,'(5x,"WRITING PATTERNS TO FILE ",a)') TRIM(filname)
     WRITE(iunit,*) nirr
     WRITE(iunit,*) (npert(i),i=1,nirr)
     WRITE(iunit,*) u
  ELSE
     WRITE( stdout,'(5x,"READING PATTERNS FROM FILE ",a)') TRIM(filname)
     READ(iunit,*) nirr
     READ(iunit,*) (npert(i),i=1,nirr)
     READ(iunit,*) u
  END IF

  CLOSE (iunit)

  RETURN
END SUBROUTINE io_pattern
