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
  USE parameters, ONLY     : ntypx
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
  INTEGER, DIMENSION(ntypx,extff_max) :: extff_atyp = 1
  !! Defines which ion specie feels the external force field (1 = yes, 0 = no)
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
    INTEGER :: i, ierr, j
    INTEGER :: extff_tunit
    INTEGER :: dum = 0 
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
       dum = 0
       DO j = 0, ntyp-1 
          dum = dum + 10 ** j
       ENDDO
       !
       READ (extff_tunit, *, iostat=ierr ) extff_typ(i), dum
       CALL errore( 'init_extffield ', 'cannot read external potential type', abs(ierr) )
       !
       DO j = 1, ntyp
          extff_atyp(j,i) = dum / (10**(ntyp-j))
          dum = dum - extff_atyp(j,i) * (10**(ntyp-j))
       ENDDO       
       !
       SELECT CASE (extff_typ(i))
          !
          CASE(1)
             !
             !! 1 : Planar force field (quadratic repulsive force LAMMPS style)
             !!     Parameters: axis direction position increment strength
             !
             READ (extff_tunit, *, iostat=ierr ) extff_axis(i), & 
                     extff_dir(i), extff_geo(1,i), extff_geo(2,i), extff_par(1,i)
             CALL errore( 'init_extffield ', 'cannot read external potential parameters', ierr )
             IF (extff_axis(i) /= 1 .AND. extff_axis(i) /= 2 .AND. extff_axis(i) /= 3) THEN
                CALL errore( 'init_extffield ', 'incorrect axis for external potential', i )
             END IF
             IF (extff_dir(i) /= 0 .AND. extff_dir(i) /= 1) THEN
                CALL errore( 'init_extffield ', 'incorrect direction for external potential', i )
             END IF             
             !
          CASE(2)
             !
             !! 2 : Viscous drag force field perpendicular to plane
             !!     Parameters: axis direction position increment strength
             !
             READ (extff_tunit, *, iostat=ierr ) extff_axis(i), & 
                     extff_dir(i), extff_geo(1,i), extff_geo(2,i), extff_par(1,i)
             CALL errore( 'init_extffield ', 'cannot read external potential parameters', ierr )
             IF (extff_axis(i) /= 1 .AND. extff_axis(i) /= 2 .AND. extff_axis(i) /= 3) THEN
                CALL errore( 'init_extffield ', 'incorrect axis for external potential', i )
             END IF
             IF (extff_dir(i) /= 0 .AND. extff_dir(i) /= 1) THEN
                CALL errore( 'init_extffield ', 'incorrect direction for external potential', i )
             END IF             
             !
          CASE DEFAULT
             CALL errore( 'init_extffield ', 'unknown external potential type', 1 )
             !
       END SELECT
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
       write(extff_unit, 1000)
       DO i = 1, nextffield
          write(extff_unit, 1010) 'Coordinate', 'Load(X)', 'Load(Y)', 'Load(Z)'
       ENDDO 
       WRITE(extff_unit, *)
    END IF

    RETURN
    !
1000  FORMAT(' Iteration',2X,$)
1010  FORMAT(4(2X,A12),$)
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
    INTEGER :: i,ia,j,k
    REAL(DP), DIMENSION(3,extff_max) :: load = 0.0
    REAL(DP), DIMENSION(3)           :: for = 0.0
    ! 
    !
    DO i = 1, nextffield
       !
       load(1,i) = 0.0
       load(2,i) = 0.0
       load(3,i) = 0.0
       !
       SELECT CASE (extff_typ(i))
       !
          CASE(1)
             !
             !! Planar force field 
             !
             DO ia = 1, nat
                !
                for(extff_axis(i)) = 0.0
                IF ( extff_dir(i) == 0 .AND. eftau0(extff_axis(i),ia) < extff_geo(1,i) ) THEN
                   for(extff_axis(i)) =  extff_atyp(ityp(ia),i) * extff_par(1,i) * (eftau0(3,ia) - extff_geo(1,i)) ** 2
                ELSE IF ( extff_dir(i) == 1 .AND. eftau0(extff_axis(i),ia) > extff_geo(1,i) ) THEN
                   for(extff_axis(i)) =  -extff_atyp(ityp(ia),i) * extff_par(1,i) * (eftau0(3,ia) - extff_geo(1,i)) ** 2
                END IF
                !
                load(extff_axis(i),i) = load(extff_axis(i),i) + for(extff_axis(i))
                effion(extff_axis(i),ia) = effion(extff_axis(i),ia) + for(extff_axis(i))
                !
             ENDDO
             !
          CASE(2)
             !
             !! Viscous drag force field perpendicular to plane
             !
             IF ( extff_axis(i) == 1) THEN
                j = 2
                k = 3
             ELSE IF ( extff_axis(i) == 2) THEN
                j = 1
                k = 3
             ELSE 
                j = 1
                k = 2
             END IF
             DO ia = 1, nat
                !
                for(j) = 0.0
                for(k) = 0.0
                IF ( extff_dir(i) == 0 .AND. eftau0(extff_axis(i),ia) < extff_geo(1,i) ) THEN
                   for(j) =  extff_atyp(ityp(ia),i) * extff_par(1,i) * efvels(j,ia) * amass(ityp(ia))
                   for(k) =  extff_atyp(ityp(ia),i) * extff_par(1,i) * efvels(k,ia) * amass(ityp(ia))
                ELSE IF ( extff_dir(i) == 1 .AND. eftau0(extff_axis(i),ia) > extff_geo(1,i) ) THEN
                   for(j) =  extff_atyp(ityp(ia),i) * extff_par(1,i) * efvels(j,ia) * amass(ityp(ia))
                   for(k) =  extff_atyp(ityp(ia),i) * extff_par(1,i) * efvels(k,ia) * amass(ityp(ia))
                END IF
                !
                load(j,i) = load(j,i) + for(j)
                load(k,i) = load(k,i) + for(k)
                effion(j,ia) = effion(j,ia) + for(j)
                effion(k,ia) = effion(k,ia) + for(k)                
                !
             ENDDO
             !
       END SELECT 
    ENDDO
    !
    IF( ionode ) THEN
       CALL print_extffield(nfi,nextffield,load)
    END IF
    !
    DO i = 1, nextffield
       !
       !! increment extffield position
       !
       extff_geo(1,i) = extff_geo(1,i) + extff_geo(2,i)
       !
    ENDDO
    
    RETURN
    !
  END SUBROUTINE apply_extffield
  !
  !--------------------------------------------------------------------------
  SUBROUTINE print_extffield(nfi,nextffield,load)
    !------------------------------------------------------------------------
    !! This routine prints information associated with external force field 
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN)   :: nextffield
    INTEGER, INTENT(IN)   :: nfi
    REAL(DP), DIMENSION(3,extff_max), INTENT(IN)  :: load
    !
    INTEGER :: i
    !
    WRITE(extff_unit, 2000) nfi
    DO i = 1, nextffield
       WRITE(extff_unit, 2010) extff_geo(1,i),load(1,i),load(2,i),load(3,i)
    ENDDO 
    WRITE(extff_unit, *)

2000  FORMAT(I10,2X,$)
2010  FORMAT(4(2X,F12.6),$)
    
  END SUBROUTINE print_extffield
  !
  !--------------------------------------------------------------------------
  SUBROUTINE close_extffield()
    !------------------------------------------------------------------------
    !! This routine close the output file  
    !
    IMPLICIT NONE
    !
    CLOSE(extff_unit)
    
    END SUBROUTINE close_extffield
    !
END MODULE extffield
