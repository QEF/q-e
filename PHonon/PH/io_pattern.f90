!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
SUBROUTINE io_pattern (nat,fildrho,nirr,npert,u,xq,directory,iflag)
!---------------------------------------------------------------------
  !
  USE kinds,            ONLY : DP
  USE io_global,        ONLY : stdout
  USE dfile_autoname,   ONLY : dfile_name
  USE io_files,         ONLY : prefix, tmp_dir, seqopn
  USE cell_base,        ONLY : at

  IMPLICIT NONE
!
!   the i/o variables first
!
  INTEGER :: nirr, iflag, nat, npert(3*nat)
  COMPLEX(DP) :: u(3*nat,3*nat)
  REAL(DP) :: xq(3)
  CHARACTER (len=256),INTENT(in)::  directory ! where to read/write the file
  CHARACTER (len=*)  ::  fildrho   ! name of the file
  CHARACTER (len=256)::  filname   ! complete name of the file
!
!   here the local variables
!
  INTEGER :: i,iunit
  INTEGER, EXTERNAL :: find_free_unit
  LOGICAL :: exst

  IF (ABS(iflag).NE.1) CALL errore('io_pattern','wrong iflag',1+ABS(iflag))

  iunit = find_free_unit() 
  filname = TRIM(dfile_name(xq, at, fildrho, TRIM(directory)//prefix, (iflag>0),-1)) //".pat"
  CALL seqopn(iunit,filname,'formatted',exst, directory)

  IF (iflag.GT.0) THEN
     !WRITE( stdout,'(5x,"WRITING PATTERNS TO FILE ",2a)') TRIM(directory), TRIM(filname)
     WRITE(iunit,*) nirr
     WRITE(iunit,*) (npert(i),i=1,nirr)
     WRITE(iunit,*) u
     WRITE(iunit,*) xq
  ELSE
     !WRITE( *,'(5x,"READING PATTERNS FROM FILE ",2a)') TRIM(directory), TRIM(filname)
     READ(iunit,*) nirr
     READ(iunit,*) (npert(i),i=1,nirr)
     READ(iunit,*) u
     READ(iunit,*) xq
  END IF
  !
  CLOSE (iunit)

  RETURN
END SUBROUTINE io_pattern
