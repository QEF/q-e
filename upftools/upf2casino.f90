!
! Copyright (C) 2011 Simon Binnie
! This file is distributed under the terms of the
! GNU General Public License. See the file 'License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------

PROGRAM upf2casino
  !---------------------------------------------------------------------
  !
  !     Convert a pseudopotential written in UPF
  !     format to CASINO tabulated format

  USE upf_module
  USE radial_grids, ONLY : radial_grid_type, deallocate_radial_grid, &
       & nullify_radial_grid
  USE pseudo_types, ONLY : pseudo_upf, nullify_pseudo_upf, deallocate_pseudo_upf
  USE casino_pp
  USE upf_module,   ONLY: read_upf 
  USE environment, ONLY: environment_start, environment_end
  USE mp_global, ONLY: mp_startup, mp_global_end
  USE io_global, ONLY: ionode, stdout

  IMPLICIT NONE
  INTEGER :: ierr, ios, prefix_len
  TYPE(pseudo_upf) :: upf_in
  CHARACTER(LEN=256)  :: filein, fileout
  TYPE(radial_grid_type) :: grid

  CALL nullify_pseudo_upf ( upf_in )
  CALL nullify_radial_grid ( grid )
#if defined(__MPI)
  CALL mp_startup()
#endif
  CALL environment_start('UPF2CASINO') 
  IF (ionode) THEN 
      WRITE(0,*) 'UPF2CASINO Converter'
      WRITE(0,*) 'Usage: ./upf2casino  pp.UPF ' 
      WRITE(0,*) 'output printed in pp.out'
      WRITE(0,*) 'All pseudopotential files generated should be &
       &thoroughly checked.'
      WRITE(0,*) 'In paticular make sure the local channel chosen&
       & in the CASINO pp file is what you expected.'

      CALL get_file ( filein ) 
      IF ( INDEX(TRIM(filein),'.UPF' ) > 0) THEN 
         prefix_len = INDEX(TRIM(filein),'.UPF' )
      ELSE IF (INDEX(TRIM(filein),'.upf') > 0 ) THEN
         prefix_len = INDEX(TRIM(filein),'.upf') 
      ELSE 
         prefix_len = LEN(TRIM(filein)) 
      ENDIF
      fileout = filein(1:prefix_len) //'out'

      CALL read_upf( upf_in, IERR = ios, GRID = grid, FILENAME = TRIM(filein)  )  
      IF (upf_in%typ /= 'NC') THEN
         WRITE(0,*) ''
         WRITE(0,*) 'WRONG PSEUDOPOTENTIAL!'
         WRITE(0,*) 'Only norm-conserving pps can be used in CASINO!'
         STOP
      ENDIF

      WRITE(0,*) "Number of grid points: ", grid%mesh
      WRITE(0,*) "Number of KB projectors: ", upf_in%nbeta
      WRITE(0,*) "Channel(s) of KB projectors: ", upf_in%lll
      WRITE(0,*) "Number of channels to be re-constructed: ", upf_in%nbeta+1

      CALL conv_upf2casino(upf_in,grid)
      OPEN ( UNIT=333, FILE=TRIM(fileout), ACTION = 'WRITE') 
      CALL write_casino_tab(upf_in,grid, 333)
      CLOSE (333) 

      DEALLOCATE(vnl)
      CALL deallocate_radial_grid(grid)
      CALL deallocate_pseudo_upf(upf_in)
   END IF
   CALL environment_end('UPF2CASINO')
#if defined(__MPI) 
  CALL mp_global_end()
#endif 

   STOP
END PROGRAM upf2casino
