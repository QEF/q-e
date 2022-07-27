!
! Copyright (C) 2022 Laurent Pizzagalli
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE extffield
  !--------------------------------------------------------------------------
  !
  ! Contains variables and subroutines for external force fields 
  !
  USE kinds,          ONLY : DP
  USE parser,         ONLY : field_count, read_line, get_field
  !
  IMPLICIT NONE
  SAVE
  !  
  INTEGER, PARAMETER :: extff_max=4 
  !! max. number of external force fields
  !
  INTEGER, DIMENSION(extff_max) :: extff_typ = 0
  !! Type of force fields
  !
  REAL(DP), DIMENSION(4,extff_max) :: extff_geo = 0.0_DP
  !! Spatial characteristics of force fields
  ! 
  INTEGER, DIMENSION(extff_max) :: extff_dir = 0
  !! Direction of force fields
  !
  INTEGER, DIMENSION(extff_max) :: extff_axis = 3
  !! Orientation of planar force fields
  !
  REAL(DP), DIMENSION(4,extff_max) :: extff_par = 0.0_DP
  !! Interaction parameters of force fields
  ! 
  CHARACTER(len=15), PARAMETER :: extff_namefile='extffield.dat' 
  !! Name of the file containing extermal force fields specifications
  !
CONTAINS
  !
  !--------------------------------------------------------------------------
  SUBROUTINE read_extffield(nextffield)
    !--------------------------------------------------------------------------
    !! This routine reads external force field parameters 
    !! have been read. If called before reading, allows to read a different
    !! input file without triggering bogus error messages - useful for NEB.
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nextffield
    INTEGER :: i, ios
    INTEGER :: extffunit
    !
    INTEGER, EXTERNAL :: find_free_unit    
    !
    CHARACTER(len=256)         :: input_line
    !
    LOGICAL                    :: tend
    !
    !!  open extffield file
    extffunit = find_free_unit()
    OPEN(unit=extffunit, file=extff_namefile, form='formatted', action='read', iostat=ios)
    IF (ios /= 0) THEN
       CALL errore('read_extffield', 'file ' // extff_namefile // ' not found', ABS(ios))
    END IF
    !
    
  END SUBROUTINE read_extffield
  !
END MODULE extffield
