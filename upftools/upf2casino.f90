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
  !     Convert a pseudopotential written in UFP
  !     format to CASINO tabulated format

  USE upf_module
  USE radial_grids, ONLY : radial_grid_type, deallocate_radial_grid, &
       & nullify_radial_grid
  USE pseudo_types, ONLY : nullify_pseudo_upf
  USE casino_pp

  IMPLICIT NONE
  INTEGER :: ierr
  TYPE(pseudo_upf) :: upf_in
  TYPE(radial_grid_type) :: grid

  CALL nullify_pseudo_upf ( upf_in )
  CALL nullify_radial_grid ( grid )

  WRITE(0,*) 'UPF2CASINO Converter'
  WRITE(0,*) 'Usage: ./upf2casino < pp_in.UPF > pp_out.dat'
  WRITE(0,*) 'All pseudopotential files generated should be &
       &thoroughly checked.'
  WRITE(0,*) 'In paticular make sure the local channel chosen&
       & in the CASINO pp file is what you expected.'

  CALL read_upf(upf_in, grid, ierr, 5)

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
  CALL write_casino_tab(upf_in,grid)

  DEALLOCATE(vnl)
  CALL deallocate_radial_grid(grid)
  CALL deallocate_pseudo_upf(upf_in)

  STOP
END PROGRAM upf2casino
