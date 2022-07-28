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
  !------------------------------------------------------------------------
  !
  ! Contains variables and subroutines for external force fields 
  !
  USE kinds,          ONLY : DP
  USE parser,         ONLY : field_count, read_line, get_field
  USE io_global,      ONLY : ionode, stdout
  USE ions_base,      ONLY : nat, ntyp => nsp, amass, ityp
  USE io_files,       ONLY : tmp_dir, prefix
  !
  IMPLICIT NONE
  SAVE
  !  
  INTEGER, PARAMETER :: extff_max=4 
  !! max. number of external force fields
  !
  INTEGER :: extff_unit 
  !! file unit for printing set by step output information
  !
  INTEGER, DIMENSION(extff_max) :: extff_typ = 0
  !! Type of force fields
  !
  INTEGER, DIMENSION(extff_max) :: extff_atyp = 0
  !! Defines which ion specie feels the external force field (0 = all)
  !
  REAL(DP), DIMENSION(6,extff_max) :: extff_geo = 0.0
  !! Spatial characteristics of force fields
  ! 
  INTEGER, DIMENSION(extff_max) :: extff_dir = 0
  !! Direction of force fields
  !
  INTEGER, DIMENSION(extff_max) :: extff_axis = 3
  !! Orientation of planar force fields (1 = X, 2 = Y, 3 = Z)
  !
  REAL(DP), DIMENSION(4,extff_max) :: extff_par = 0.0
  !! Interaction parameters of force fields
  ! 
  REAL(DP), DIMENSION(extff_max) :: extff_load = 0.0
  !! Computed load for each force field
  ! 
  CHARACTER(len=15), PARAMETER :: extff_namefile='extffield.dat' 
  !! Name of the file containing extermal force fields specifications
  !
CONTAINS
  !
  !--------------------------------------------------------------------------
  SUBROUTINE init_extffield(nextffield)
    !------------------------------------------------------------------------
    !! This routine reads external force field parameters 
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nextffield
    INTEGER :: i, ierr
    INTEGER :: extff_tunit
    !
    INTEGER, EXTERNAL :: find_free_unit    
    !
    CHARACTER(len=256)         :: input_line
    CHARACTER(LEN=256)         :: pprefix, extff_outputfile
    !
    LOGICAL                    :: tend
    !
    !!  open extffield file
    !
    extff_tunit = find_free_unit()
    !
    OPEN(unit=extff_tunit, file=extff_namefile, form='formatted', action='read', iostat=ierr)
    !
    IF (ierr /= 0) THEN
       CALL errore('init_extffield', 'file ' // extff_namefile // ' not found', abs(ierr))
    END IF
    !
    !! read extffield file
    !
    DO i = 1, nextffield
       !
       READ (extff_tunit, *, iostat=ierr ) extff_typ(i)
       CALL errore( 'init_extffield ', 'cannot read external potential type', abs(ierr) )
       !
       SELECT CASE (extff_typ(i))
          !
          ! 1 : Planar indenter (quadratic repulsive force LAMMPS style)
          ! Parameters: axis direction position increment strength
          !
          CASE(1)
             READ (extff_tunit, *, iostat=ierr ) extff_axis(i), & 
                     extff_dir(i), extff_geo(1,i), extff_geo(2,i), extff_par(1,i)
             ! NEED TO MAKE TESTS FOR AXIS AND DIR
          !
          CASE DEFAULT
             CALL errore( 'init_extffield ', 'unknown external potential type', 1 )
          !
       END SELECT
       !
       IF (ierr /= 0) THEN
          CALL errore( 'init_extffield ', 'cannot read external potential parameters', i )
       END IF
       !
    ENDDO
    !
    CLOSE(extff_unit)
    !
    !! print the information in output
    ! 
    IF ( ionode ) THEN
       WRITE( stdout, *) 
       WRITE( stdout, *) '  External potential activation'
       WRITE( stdout, *) '  -----------------------------'
    END IF
    !
    !! open file for printint extffield information
    !
    extff_unit = find_free_unit()
    !
    ! ...  prefix combined with the output path
    pprefix = TRIM( tmp_dir ) // TRIM( prefix )
    extff_outputfile = TRIM(pprefix) // '.extffield'
    !
    OPEN(unit=extff_unit, file=extff_outputfile, form='formatted', action='write', iostat=ierr)
    IF (ierr /= 0) THEN
       CALL errore( 'init_extffield ', 'cannot open external potential output file', i )
    END IF  
      
    IF ( ionode ) THEN
       write(extff_unit, 1010)
       DO i = 1, nextffield
          write(extff_unit, 1020) 'Position', 'Load'
       ENDDO 
       write(extff_unit, 1030)
    END IF

    RETURN
    !
1010  FORMAT('Iterations',2X,$)
1020  FORMAT(A15,2X,A15,2X,$)
1030  FORMAT()
  !
  END SUBROUTINE init_extffield
  !
  !--------------------------------------------------------------------------
  SUBROUTINE apply_extffield(nfi,nextffield,eftau0,effion,efvels)
    !------------------------------------------------------------------------
    !! This routine apply external force field on ions
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nfi,nextffield
    REAL(DP), DIMENSION(3,nat), INTENT(IN):: eftau0
    REAL(DP), DIMENSION(3,nat), INTENT(INOUT):: effion
    REAL(DP), DIMENSION(3,nat), INTENT(IN):: efvels
    !
    INTEGER :: i,ia
    REAL(DP) :: load,for
    ! 
    !
    DO i = 1, nextffield
       !
       load = 0.0
       DO ia = 1, nat
          !
          IF ( extff_dir(i) == 0 .AND. eftau0(extff_axis(i),ia) < extff_geo(1,i) ) THEN 
             for =  extff_par(1,i) * (eftau0(3,ia) - extff_geo(1,i)) ** 2
          ELSE IF ( extff_dir(i) == 1 .AND. eftau0(extff_axis(i),ia) > extff_geo(1,i) ) THEN
             for =  extff_par(1,i) * (eftau0(3,ia) - extff_geo(1,i)) ** 2
          END IF
          !
          load = load + for
          effion(extff_axis(i),ia) = effion(extff_axis(i),ia) + for
          !
       ENDDO
       !
    ENDDO
    !
    CALL print_extffield(nfi,nextffield,extff_geo(1,i),load)
    RETURN
    !
  END SUBROUTINE apply_extffield
  !
  !--------------------------------------------------------------------------
  SUBROUTINE print_extffield(nfi,nextffield,geo,load)
    !------------------------------------------------------------------------
    !! This routine prints information associated with external force field 
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN)   :: nextffield
    INTEGER, INTENT(IN)   :: nfi
    REAL(DP), INTENT(IN)  :: geo    
    REAL(DP), INTENT(IN)  :: load
    !
    IF( ionode ) THEN
       WRITE(extff_unit, 1010) nfi, geo, load
    END IF

1010  FORMAT(I10,2X,6(F15.6,2X))
    
  END SUBROUTINE print_extffield
END MODULE extffield
