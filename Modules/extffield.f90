!
! Copyright (C) 2002-2022 Quantum ESPRESSO group
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
  ! This module implements external ionic force fields in cp.x and pw.x
  ! Activation is obtained by setting the input parameter 'nextffield' in 
  ! the Namelist &SYSTEM: 
  !         nextffield = INTEGER
  ! nextffield is the number of activated external force fields (default = 0, max = 4)
  ! If nextffield > 0 the file 'extffield.dat' is read
  !  
  ! For each force field, two lines are required. The first one is generic 
  ! to all kinds of force field, with the format:
  ! NTYP  FLAGS
  ! NTYP is an integer defining the type of force fields (see below)
  ! FLAGS is is an integer composed of 1 or 0, as many as ionic species. 
  ! It can be used to restrict the application of force fields to certain species only
  ! For instance, FLAGS = 101 meansthat ionic species 1 and 3 are subject to the force field, but not ionic specie 2 
  ! The second line defines various parameters depending on the force field type.
  !
  ! For a description and an application of the method, see L. Pizzagalli, Phys. Rev. B 102, 094102 (2020)
  !
  ! 
  ! FORCE FIELDS:
  ! -------------
  ! 
  ! -------------------------------------------------------------------------------------------------
  ! NTYP = 1 : Planar quadratic repulsive force field
  ! This force field mimics the command 'FIX INDENT PLANE...' in the LAMMPS code
  ! which is often used in MD studies to model a flat punch indenter
  ! The force expression is $\pm K(z-z_p)^2$
  ! 
  ! The second line in 'extffield.dat' for this potential includes 5 parameters
  ! AXIS  DIR  POS INC STRENGTH
  ! AXIS is an integer defining the axis for the plane (1 = X, 2 = Y, 3 = Z)
  ! DIR is an integer defining the direction of the force (0 is positive, and 1 negative), 
  ! and the selection of ions.
  ! POS is a real defining the position of the plane relative to AXIS (z_p in the formula)
  !   Selected ions are those with a position below POS (DIR = 0) or above POS (DIR = 1)
  ! INC is a real, added to POS at each iteration (dynamic compression)
  ! STRENGTH is a real defining the strength of the repulsion (K in the formula)
  ! Rydberg (Hartree) atomic units are used for pw.x (cp.x). 
  ! 
  ! For instance, with the following two lines
  ! 1   10
  ! 3   0   2.50   0.01   10
  ! one defines a planar repulsive potential acting on the first atomic specie 
  ! (but not on the second one). The potential is applied relatively to a plane normal to Z and of initial position 
  ! 2.50 bohrs along the Z axis, with a positive directionalong the axis Z, with a 
  ! positive direction. All ions of specie 1, 
  ! with an initial z-coordinate below 2.50 will be subject to a positive force along this axis.
  ! The potential strength is 10 a.u.  
  ! The potential threshold is moved up by 0.01 bohr at each ionic iteration
  !
  !
  ! -------------------------------------------------------------------------------------------------
  ! NTYP = 2 : Viscous drag force field perpendicular to a plane
  ! This force field adds a viscous friction for selected atoms, by adding velocity dependent forces in 
  ! the two directions perpendicular to the defined axis 
  ! It can be used in combination with the previous force field, to prevent an excessive rotation 
  ! of the system during the dynamics. The force expression is
  ! $-K*m*v$
  ! The force is proportional and opposed to the ion velocity, and is also proportional to the ion mass 
  ! and the constant K
  ! Selected ions are those with a position below POS (DIR = 0) or above POS (DIR = 1)
  !
  ! The second line in 'extffield.dat' for this potential looks like
  ! AXIS  DIR  POS INC STRENGTH
  ! all parameters have the same meaning than for NTYP = 1
  !
  ! NOTE: NTYP = 2 is only available when using cw.x 
  !
  ! -------------------------------------------------------------------------------------------------
  ! NTYP = 3 : Planar Lennard-Jones potential
  ! This force field allows to impose an interaction of the system of interest with a semi-infinite slab. 
  ! The forces are derived from the well known standard LJ energy formula
  ! 
  ! $V(r) = 4\varepsilon((\sigma/r)^12 - (\sigma/r)^6)$
  ! 
  ! The second line in 'extffield.dat' includes the following 7 parameters:
  ! AXIS  DIR  POS INC \varepsilon \sigma cutoff
  ! The first 4 have the same meaning than for NTYP = 1
  ! \varepsilon and \sigma are the LJ potential parameters (with coherent units for cp.x or pw.x)
  ! cutoff defines the range of the potential
  ! For instance, with cutoff = 5.0, DIR = 0 and AXIS = 3, an ion will feel the potential only if its 
  ! z-coordinate is lower than POS+cutoff
  !
  ! OUTPUT:
  ! -------
  ! 
  ! Information about the defined force fields is written in the standard output file
  ! A file with name 'prefix.extffield' is created, including data per ionic iteration for all defined force fields
  ! For each force field, the position of the plane and the sum of added forces for each axis are written at each step,
  ! For a planar compression, this corresponds to the compression load. 
  ! 
  ! 
  ! RESTRICTIONS:
  ! ------------- 
  ! - No automatic 'restart'. For a follow up calculation, the parameters in 'extffield.dat' have to be set accordingly
  ! - 4 maximum force fields
  ! - Minimal testing done, so use at your own risk
  !
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
  SUBROUTINE init_extffield(prog, nextffield)
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
    CHARACTER(LEN=2)           :: prog
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
    !! print the information in output
    ! 
    IF ( ionode ) THEN
       WRITE( stdout, *) 
       WRITE( stdout, *) '  External force field information'
       WRITE( stdout, *) '  --------------------------------'
       WRITE( stdout, 1020) ntyp 
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
             IF ( ionode) THEN
                WRITE( stdout, *) i,': Repulsive planar (Fix indent lammps style) potential'
                WRITE( stdout, 1030)  extff_axis(i), extff_dir(i), extff_geo(1,i), extff_geo(2,i), extff_par(1,i)
             END IF           
             !
          CASE(2)
             !
             !! 2 : Viscous drag force field perpendicular to plane
             !!     Parameters: axis direction position increment strength
             !
             IF ( prog == 'PW' ) THEN
                CALL errore( 'init_extffield ', 'Viscous force field not available for pw.x', 1 )
             END IF
             READ (extff_tunit, *, iostat=ierr ) extff_axis(i), & 
                     extff_dir(i), extff_geo(1,i), extff_geo(2,i), extff_par(1,i)
             CALL errore( 'init_extffield ', 'cannot read external potential parameters', ierr )
             IF (extff_axis(i) /= 1 .AND. extff_axis(i) /= 2 .AND. extff_axis(i) /= 3) THEN
                CALL errore( 'init_extffield ', 'incorrect axis for external potential', i )
             END IF
             IF (extff_dir(i) /= 0 .AND. extff_dir(i) /= 1) THEN
                CALL errore( 'init_extffield ', 'incorrect direction for external potential', i )
             END IF             
             IF ( ionode) THEN
                WRITE( stdout, *) i,': Viscous drag planar potential'
                WRITE( stdout, 1030)  extff_axis(i), extff_dir(i), extff_geo(1,i), extff_geo(2,i), extff_par(1,i)
             END IF           
             !
          CASE(3)
             !
             !! 3 : Leannard-Jones planar force field
             !!     Parameters: axis direction position increment \varepsilon \sigma cutoff
             !
             READ (extff_tunit, *, iostat=ierr ) extff_axis(i), extff_dir(i), & 
                     extff_geo(1,i), extff_geo(2,i), extff_par(1,i), extff_par(2,i), extff_par(3,i)
             CALL errore( 'init_extffield ', 'cannot read external potential parameters', ierr )
             IF (extff_axis(i) /= 1 .AND. extff_axis(i) /= 2 .AND. extff_axis(i) /= 3) THEN
                CALL errore( 'init_extffield ', 'incorrect axis for external potential', i )
             END IF
             IF (extff_dir(i) /= 0 .AND. extff_dir(i) /= 1) THEN
                CALL errore( 'init_extffield ', 'incorrect direction for external potential', i )
             END IF             
             IF ( ionode) THEN
                WRITE( stdout, *) i,': Lennard-Jones planar potential'
                WRITE( stdout, 1040)  extff_axis(i), extff_dir(i), extff_geo(1,i), extff_geo(2,i), &
                   extff_par(1,i), extff_par(2,i), extff_par(3,i)
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
1020  FORMAT(5X,I1,' external force field(s):')
1030  FORMAT(13X,'axis = ',I1,'  dir = ',I1,' pos = ',F8.4,' inc = ',F8.4,' Strength = ',F10.4)
1040  FORMAT(13X,'axis = ',I1,'  dir = ',I1,' pos = ',F8.4,' inc = ',F8.4,' Eps = ',F8.4,' Sigma = ',F8.4,' Cutoff = ',F10.4)
  !
  END SUBROUTINE init_extffield
  !
  !--------------------------------------------------------------------------
  SUBROUTINE apply_extffield_PW(nfi,nextffield,eftau,effion)
    !------------------------------------------------------------------------
    !! This routine apply external force field on ions  (PW code)
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nfi,nextffield
    REAL(DP), DIMENSION(3,nat), INTENT(IN):: eftau
    REAL(DP), DIMENSION(3,nat), INTENT(INOUT):: effion
    !
    INTEGER :: i,ia,j,k
    REAL(DP), DIMENSION(3,extff_max) :: load = 0.0
    REAL(DP), DIMENSION(3)           :: for = 0.0
    REAL(DP)                         :: alp = 0.0
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
             !! Planar quadratic repulsive force field (Fix indent lammps style)
             !
             DO ia = 1, nat
                !
                for(extff_axis(i)) = 0.0
                IF ( extff_dir(i) == 0 .AND. eftau(extff_axis(i),ia) < extff_geo(1,i) ) THEN
                   for(extff_axis(i)) =  extff_atyp(ityp(ia),i) * extff_par(1,i) * (eftau(extff_axis(i),ia) - extff_geo(1,i)) ** 2
                ELSE IF ( extff_dir(i) == 1 .AND. eftau(extff_axis(i),ia) > extff_geo(1,i) ) THEN
                   for(extff_axis(i)) =  -extff_atyp(ityp(ia),i) * extff_par(1,i) * (eftau(extff_axis(i),ia) - extff_geo(1,i)) ** 2
                END IF
                !
                load(extff_axis(i),i) = load(extff_axis(i),i) + for(extff_axis(i))
                effion(extff_axis(i),ia) = effion(extff_axis(i),ia) + for(extff_axis(i))
                !print *,ia,ityp(ia),load(extff_axis(i),i),eftau(extff_axis(i),ia),extff_geo(1,i)
                !
             ENDDO
             !
!          CASE(2)
!             !
!             !! Viscous drag force field perpendicular to plane
!             !
!             IF ( extff_axis(i) == 1) THEN
!                j = 2
!                k = 3
!             ELSE IF ( extff_axis(i) == 2) THEN
!                j = 1
!                k = 3
!             ELSE 
!                j = 1
!                k = 2
!             END IF
!             DO ia = 1, nat
!                !
!                for(j) = 0.0
!                for(k) = 0.0
!                IF ( extff_dir(i) == 0 .AND. eftau(extff_axis(i),ia) < extff_geo(1,i) ) THEN
!                   for(j) =  -extff_atyp(ityp(ia),i) * extff_par(1,i) * efvels(j,ia) * amass(ityp(ia))
!                   for(k) =  -extff_atyp(ityp(ia),i) * extff_par(1,i) * efvels(k,ia) * amass(ityp(ia))
!                ELSE IF ( extff_dir(i) == 1 .AND. eftau(extff_axis(i),ia) > extff_geo(1,i) ) THEN
!                   for(j) =  -extff_atyp(ityp(ia),i) * extff_par(1,i) * efvels(j,ia) * amass(ityp(ia))
!                   for(k) =  -extff_atyp(ityp(ia),i) * extff_par(1,i) * efvels(k,ia) * amass(ityp(ia))
!                END IF
!                !
!                load(j,i) = load(j,i) + for(j)
!                load(k,i) = load(k,i) + for(k)
!                effion(j,ia) = effion(j,ia) + for(j)
!                effion(k,ia) = effion(k,ia) + for(k)                
!                !
!             ENDDO
             !
          CASE(3)
             !
             !! Planar Lennard-Jones potential
             !
             DO ia = 1, nat
                !
                for(extff_axis(i)) = 0.0
                IF ( extff_dir(i) == 0 .AND. eftau(extff_axis(i),ia) < (extff_geo(1,i) + extff_par(3,i)) ) THEN
                   alp = (extff_par(2,i)/(eftau(extff_axis(i),ia) - extff_geo(1,i)))**6
                   alp = alp - 2.0 * alp**2
                   for(extff_axis(i)) =  -24.0 * extff_par(1,i) / (eftau(extff_axis(i),ia) - extff_geo(1,i)) * alp
                   for(extff_axis(i)) =  extff_atyp(ityp(ia),i) * for(extff_axis(i))
                ELSE IF ( extff_dir(i) == 1 .AND. eftau(extff_axis(i),ia) > ( extff_geo(1,i) - extff_par(3,i)) ) THEN
                   alp = (extff_par(2,i)/(eftau(extff_axis(i),ia) - extff_geo(1,i)))**6
                   alp = alp - 2.0 * alp**2
                   for(extff_axis(i)) =  -24.0 * extff_par(1,i) / (eftau(extff_axis(i),ia) - extff_geo(1,i)) * alp
                   for(extff_axis(i)) =  extff_atyp(ityp(ia),i) * for(extff_axis(i))
                END IF
                !
                load(extff_axis(i),i) = load(extff_axis(i),i) + for(extff_axis(i))
                effion(extff_axis(i),ia) = effion(extff_axis(i),ia) + for(extff_axis(i))
                !print *,ia,ityp(ia),load(extff_axis(i),i),eftau(extff_axis(i),ia),extff_geo(1,i)
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
  END SUBROUTINE apply_extffield_PW
  !
  !--------------------------------------------------------------------------
  SUBROUTINE apply_extffield_CP(nfi,nextffield,eftau,efvels,effion)
    !------------------------------------------------------------------------
    !! This routine apply external force field on ions  (CP code)
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nfi,nextffield
    REAL(DP), DIMENSION(3,nat), INTENT(IN):: eftau
    REAL(DP), DIMENSION(3,nat), INTENT(INOUT):: effion
    REAL(DP), DIMENSION(3,nat), INTENT(IN):: efvels
    !
    INTEGER :: i,ia,j,k
    REAL(DP), DIMENSION(3,extff_max) :: load = 0.0
    REAL(DP), DIMENSION(3)           :: for = 0.0
    REAL(DP)                         :: alp = 0.0
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
             !! Planar quadratic repulsive force field (Fix indent lammps style)
             !
             DO ia = 1, nat
                !
                for(extff_axis(i)) = 0.0
                IF ( extff_dir(i) == 0 .AND. eftau(extff_axis(i),ia) < extff_geo(1,i) ) THEN
                   for(extff_axis(i)) =  extff_atyp(ityp(ia),i) * extff_par(1,i) * (eftau(extff_axis(i),ia) - extff_geo(1,i)) ** 2
                ELSE IF ( extff_dir(i) == 1 .AND. eftau(extff_axis(i),ia) > extff_geo(1,i) ) THEN
                   for(extff_axis(i)) =  -extff_atyp(ityp(ia),i) * extff_par(1,i) * (eftau(extff_axis(i),ia) - extff_geo(1,i)) ** 2
                END IF
                !
                load(extff_axis(i),i) = load(extff_axis(i),i) + for(extff_axis(i))
                effion(extff_axis(i),ia) = effion(extff_axis(i),ia) + for(extff_axis(i))
                !print *,ia,ityp(ia),eftau(extff_axis(i),ia),extff_geo(1,i),for(extff_axis(i)),effion(extff_axis(i),ia)
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
                IF ( extff_dir(i) == 0 .AND. eftau(extff_axis(i),ia) < extff_geo(1,i) ) THEN
                   for(j) =  -extff_atyp(ityp(ia),i) * extff_par(1,i) * efvels(j,ia) * amass(ityp(ia))
                   for(k) =  -extff_atyp(ityp(ia),i) * extff_par(1,i) * efvels(k,ia) * amass(ityp(ia))
                ELSE IF ( extff_dir(i) == 1 .AND. eftau(extff_axis(i),ia) > extff_geo(1,i) ) THEN
                   for(j) =  -extff_atyp(ityp(ia),i) * extff_par(1,i) * efvels(j,ia) * amass(ityp(ia))
                   for(k) =  -extff_atyp(ityp(ia),i) * extff_par(1,i) * efvels(k,ia) * amass(ityp(ia))
                END IF
                !
                load(j,i) = load(j,i) + for(j)
                load(k,i) = load(k,i) + for(k)
                effion(j,ia) = effion(j,ia) + for(j)
                effion(k,ia) = effion(k,ia) + for(k)                
                !
             ENDDO
             !
          CASE(3)
             !
             !! Planar Lennard-Jones potential
             !
             DO ia = 1, nat
                !
                for(extff_axis(i)) = 0.0
                IF ( extff_dir(i) == 0 .AND. eftau(extff_axis(i),ia) < (extff_geo(1,i) + extff_par(3,i)) ) THEN
                   alp = (extff_par(2,i)/(eftau(extff_axis(i),ia) - extff_geo(1,i)))**6
                   alp = alp - 2.0 * alp**2
                   for(extff_axis(i)) =  -24.0 * extff_par(1,i) / (eftau(extff_axis(i),ia) - extff_geo(1,i)) * alp
                   for(extff_axis(i)) =  extff_atyp(ityp(ia),i) * for(extff_axis(i))
                ELSE IF ( extff_dir(i) == 1 .AND. eftau(extff_axis(i),ia) > ( extff_geo(1,i) - extff_par(3,i)) ) THEN
                   alp = (extff_par(2,i)/(eftau(extff_axis(i),ia) - extff_geo(1,i)))**6
                   alp = alp - 2.0 * alp**2
                   for(extff_axis(i)) =  -24.0 * extff_par(1,i) / (eftau(extff_axis(i),ia) - extff_geo(1,i)) * alp
                   for(extff_axis(i)) =  extff_atyp(ityp(ia),i) * for(extff_axis(i))
                END IF
                !
                load(extff_axis(i),i) = load(extff_axis(i),i) + for(extff_axis(i))
                effion(extff_axis(i),ia) = effion(extff_axis(i),ia) + for(extff_axis(i))
                !print *,ia,ityp(ia),load(extff_axis(i),i),eftau(extff_axis(i),ia),extff_geo(1,i)
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
  END SUBROUTINE apply_extffield_CP
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
