!
! Copyright (C) 2008 Simon Binnie
! This file is distributed under the terms of the
! GNU General Public License. See the file 'License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
PROGRAM casino2upf
  !---------------------------------------------------------------------
  !
  !     Convert a pseudopotential written in CASINO tabulated
  !     format to unified pseudopotential format

  USE casino_pp

  IMPLICIT NONE
  CHARACTER(len=256) filein, fileout
  CHARACTER(len=256), ALLOCATABLE:: wavefile(:)
  INTEGER nofiles, i

  PRINT*, 'CASINO2UPF Converter'
  PRINT*, 'Enter CASINO pp.data filename:'

  CALL get_file ( filein )

  PRINT*, 'How many wavefunction *files* are you using?'
  READ(*,*) nofiles
  ALLOCATE(wavefile(nofiles))
  PRINT*, 'Enter wavefunction files, starting with the ground state:'
  DO i=1,nofiles
     CALL get_file ( wavefile(i) )
     OPEN(unit=i,file=wavefile(i),status='old',form='formatted')
  ENDDO
  OPEN(unit=99,file=filein,status='old',form='formatted')

  CALL read_casino(99,nofiles)
  CLOSE (unit=99)
  DO i=1,nofiles
     CLOSE (i)
  ENDDO

  ! convert variables read from CASINO format into those needed
  ! by the upf format - add missing quantities

  CALL convert_casino

  fileout=trim(filein)//'.UPF'
  PRINT '(''Output PP file in US format :  '',a)', fileout

  OPEN(unit=2,file=fileout,status='unknown',form='formatted')
  CALL write_upf(2)
  CLOSE (unit=2)
  STOP
20 CALL errore ('casino2upf', 'Reading pseudo file name ', 1)

END PROGRAM casino2upf
