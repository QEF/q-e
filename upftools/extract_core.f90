!
! Copyright (C) 2011 Simon Binnie
! This file is distributed under the terms of the
! GNU General Public License. See the file 'License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------

PROGRAM extract_core
  !---------------------------------------------------------------------
  !
  ! extract core wavefunctions from GIPAW section of UPF file
  ! and save them to text files

  USE emend_upf_module, ONLY: make_emended_upf_copy
  USE wrappers, ONLY: f_remove
  USE radial_grids, ONLY : radial_grid_type, deallocate_radial_grid, &
       & nullify_radial_grid
  USE pseudo_types, ONLY : pseudo_upf, nullify_pseudo_upf, deallocate_pseudo_upf
  USE upf_module,   ONLY: read_upf 
  USE environment, ONLY: environment_start, environment_end
  USE mp_global, ONLY: mp_startup, mp_global_end
  USE io_global, ONLY: ionode, stdout

  IMPLICIT NONE
  INTEGER :: ierr, ios, prefix_len, i,j
  TYPE(pseudo_upf) :: upf_in
  CHARACTER(LEN=256)  :: filein, fileout
  TYPE(radial_grid_type) :: grid
  LOGICAL :: is_xml

  CALL nullify_pseudo_upf ( upf_in )
  CALL nullify_radial_grid ( grid )
  IF (ionode) THEN 
      WRITE(0,*) 'Usage: ./extract_core.x  pp.UPF ' 

      CALL get_file ( filein ) 
      IF ( INDEX(TRIM(filein),'.UPF' ) > 0) THEN 
         prefix_len = INDEX(TRIM(filein),'.UPF' )
      ELSE IF (INDEX(TRIM(filein),'.upf') > 0 ) THEN
         prefix_len = INDEX(TRIM(filein),'.upf') 
      ELSE 
         prefix_len = LEN(TRIM(filein)) 
      ENDIF
      fileout = filein(1:prefix_len) 

      ios = 0
      CALL read_upf( upf_in, IERR = ios, GRID = grid, FILENAME = TRIM(filein)  )
      !! read UPF v.2
      IF (ios ==-81 ) THEN
         !! Bad characters found: try to remove them and fix the file
         is_xml = make_emended_upf_copy( TRIM(filein), 'tmp.upf' )
         IF (is_xml) THEN
            !! read fixed UPF v.2 file, remove after reading
            CALL read_upf(upf_in, IERR = ios, GRID = grid, FILENAME = 'tmp.upf')
            ios = f_remove('tmp.upf' )
         ELSE
            !! Not an UPF v.2 file: read UPF v.1
            OPEN (UNIT = 999, FILE = filein, STATUS = 'old', FORM = 'formatted')
            CALL  read_upf(upf_in, grid, ios, UNIT = 999 )
            CLOSE (UNIT = 999, STATUS = 'keep')
         END IF
      END IF
 
      IF(upf_in%has_gipaw) THEN 
         DO j =1, upf_in%gipaw_ncore_orbitals
           OPEN(unit=999, file=TRIM(fileout)//TRIM(upf_in%gipaw_core_orbital_el(j))//".out")
           WRITE(stdout,"('writing: ',a)") TRIM(fileout)//TRIM(upf_in%gipaw_core_orbital_el(j))//".out"
           DO i = 1, upf_in%grid%mesh
             WRITE(999,*) upf_in%grid%r(i), upf_in%gipaw_wfs_ae(i,1)
           ENDDO
           CLOSE(999)
         ENDDO
      ELSE
        WRITE(stdout,*) "Error: this pseudopotential does not contains gipaw data"
      ENDIF
 
      CALL deallocate_radial_grid(grid)
      CALL deallocate_pseudo_upf(upf_in)
   END IF

END PROGRAM extract_core
