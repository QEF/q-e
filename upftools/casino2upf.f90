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
  USE write_upf_v2_module, ONLY :  write_upf_v2
  USE pseudo_types, ONLY : nullify_pseudo_upf, deallocate_pseudo_upf, &
                           pseudo_upf
  IMPLICIT NONE
  !
  INTEGER, EXTERNAL :: find_free_unit
  !
  CHARACTER(len=256) :: pp_data
  CHARACTER(len=256) :: upf_file
  CHARACTER(len=256), ALLOCATABLE:: wavefile(:)
  INTEGER, ALLOCATABLE :: waveunit(:)
  INTEGER nofiles, i, ios, pp_unit
  TYPE(pseudo_upf)      :: upf_out

  NAMELIST / inputpp / &
       pp_data,        &         !CASINO pp filename
       upf_file,        &         !output file
       tn_grid,        &         !.true. if Trail and Needs grid is used
       tn_prefac,      &
       xmin,           &         !xmin for standard QE grid
       dx                        !dx for Trail and Needs and standard QE
                                 !grid
  pp_data= 'pp.data'
  upf_file= 'out.UPF'

  CALL nullify_pseudo_upf( upf_out )

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
  OPEN(unit=pp_unit,file=trim(pp_data),status='old',form='formatted', iostat=ios)
  IF (ios /= 0 ) THEN
     CALL errore ('casino2upf', 'cannot read file', trim(wavefile(i)))
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

  PRINT '(''Output PP file in UPF format :  '',a)', upf_file
  OPEN(unit=2,file=upf_file,status='unknown',form='formatted')
 
  CALL write_upf_v2(u=2,upf=upf_out)

  CLOSE(unit=2,status='keep')
  CALL  deallocate_pseudo_upf( upf_out )

  STOP

END PROGRAM casino2upf
