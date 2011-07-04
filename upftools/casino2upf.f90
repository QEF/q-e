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
  USE upf_module
  
  IMPLICIT NONE
  !
  INTEGER, EXTERNAL :: find_free_unit
  !
  CHARACTER(len=256) :: pp_data = 'pp.data'
  CHARACTER(len=256), ALLOCATABLE:: wavefile(:)
  INTEGER, ALLOCATABLE :: waveunit(:)
  INTEGER nofiles, i, ios, pp_unit
  TYPE(pseudo_upf)      :: upf_out

  NAMELIST / inputpp / &
       pp_data,        &         !CASINO pp filename
       tn_grid,        &         !.true. if Trail and Needs grid is used
       tn_prefac,      &         
       xmin,           &         !xmin for standard QE grid
       dx                        !dx for Trail and Needs and standard QE
                                 !grid

  WRITE(0,*) 'CASINO2UPF Converter'

  READ(*,inputpp,iostat=ios)

  READ(*,*,iostat=ios) nofiles

  ALLOCATE(wavefile(nofiles), waveunit(nofiles))

  !Now read in the awfn file names and open the files

  DO i=1,nofiles
     READ(*,*,iostat=ios) wavefile(:)
     waveunit(i)=find_free_unit()
     OPEN(unit=waveunit(i),file=trim(wavefile(i)),&
          status='old',form='formatted', iostat=ios)
     IF (ios /= 0 ) THEN
        CALL errore ('casino2upf', 'cannot read file', trim(wavefile(i)))
     ENDIF
  ENDDO
  
  pp_unit=find_free_unit()
  OPEN(unit=pp_unit,file=TRIM(pp_data),status='old',form='formatted', iostat=ios)
  IF (ios /= 0 ) THEN
     CALL errore ('casino2upf', 'cannot read file', TRIM(wavefile(i)))
  ENDIF
  
  CALL read_casino(pp_unit,nofiles, waveunit)
  
  CLOSE (unit=pp_unit)
  DO i=1,nofiles
     CLOSE (waveunit(i))
  ENDDO
  
  DEALLOCATE( wavefile, waveunit )

  ! convert variables read from CASINO format into those needed
  ! by the upf format - add missing quantities

  CALL convert_casino(upf_out)

  CALL write_upf(upf=upf_out,unit=6) !filename=fileout)

  STOP

END PROGRAM casino2upf
