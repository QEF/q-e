!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


!  AB INITIO COSTANT PRESSURE MOLECULAR DYNAMICS
!  ----------------------------------------------
!  Car-Parrinello Parallel Program
!  Carlo Cavazzoni - Gerardo Ballabio
!  SISSA, Trieste, Italy - 1997-99
!  Last modified: Sun Nov 14 07:54:47 MET 1999
!  ----------------------------------------------
!  BEGIN manual

!=----------------------------------------------------------------------------=!
      MODULE atoms_type_module
!=----------------------------------------------------------------------------=!

!  this module contains the definitions of several TYPE structures
!  relative to the ionic degrees of freedom,
!  together with their allocation/deallocation routines
!  ----------------------------------------------
!  routines in this module:
!  ----------------------------------------------
!  SUBROUTINE specie_index(isa, na, is, ia)
!  SUBROUTINE allocate_atoms_type(atoms, staur, ismbl, pma, na)
!  SUBROUTINE deallocate_atoms_type(atoms)
!  SUBROUTINE allocate_constrains(constrains, na, cindexes, cvalues, ctype, tol, nc)
!  SUBROUTINE deallocate_constrains(constrains)
!  END manual

        USE kinds
        USE parameters, ONLY: nsx, natx

        IMPLICIT NONE
        SAVE

        PRIVATE

!  BEGIN manual
!  TYPE DEFINITIONS


! ...   title ...
        TYPE atoms_type
          INTEGER :: doft      ! total number of degree_of_freedom
          INTEGER :: nsp       ! number of species
          INTEGER :: nat       ! total number of atoms
          INTEGER :: nax       ! maximum number of atoms per specie
          INTEGER :: dof(nsx)  ! degree_of_freedom for each specie

          CHARACTER(LEN=4) :: label(nsx) !  atomic labels
          INTEGER   :: na(nsx)  !   number of atoms per specie
          INTEGER   :: isa(nsx) !   index of the first atom (in the whole list) of a given specie
          REAL(dbl) :: m(nsx)   !   atomic masses
          REAL(dbl) :: taur(3,natx)  
          REAL(dbl) :: taus(3,natx)  
          ! ... tau: atomic positions, sorted by specie. Atomic positions of specie "is" are
          !          stored in array elements whose index are "isa(is) ... isa(is)+na(is)-1"
          REAL(dbl) :: vels(3,natx)  !  scaled velocities, same layout of "tau"
          REAL(dbl) :: for(3,natx)  !  total force acting on the atom
          INTEGER :: mobile(3,natx) !  atomic freedom, same layout of "tau" ( 1 atom can move )
          INTEGER :: ityp(natx)     !  index of the specie to which the atom belong
          LOGICAL :: tscfor         !  indicate if the force are scaled or real
          REAL(dbl) :: ekin(nsx)    !  kinetic energy per specie
          REAL(dbl) :: ekint        !  total kinetic energy
        END TYPE atoms_type

! .. 4 int + nsx int + 4 char + 2 ( nsx int ) + nsx dbl + 3 ( 3 dbl natx ) + 3 lg natx +
! .. natx int + 3 lg + nsx dbl + dbl 

        TYPE distance_constrain
          INTEGER :: ia1, is1
          INTEGER :: ia2, is2
          REAL(dbl) :: val
        END TYPE distance_constrain

        TYPE constrain_type
          INTEGER :: what
          TYPE (distance_constrain), POINTER :: distance
        END TYPE constrain_type

        TYPE constrains_class
          INTEGER :: nc
          REAL(dbl) :: tolerance
          TYPE (constrain_type), POINTER :: tp(:) 
        END TYPE

!  ----------------------------------------------
!  END manual


        PUBLIC :: atoms_type, constrain_type, constrains_class
        PUBLIC :: distance_constrain
        PUBLIC :: deallocate_atoms_type, allocate_constrains
        PUBLIC :: deallocate_constrains, atoms_type_init

!  end of module-scope declarations
!  ----------------------------------------------

!=----------------------------------------------------------------------------=!
      CONTAINS
!=----------------------------------------------------------------------------=!

!  subroutines

        SUBROUTINE specie_index(isa, na, is, ia)
          INTEGER, INTENT(IN) :: isa, na(:)
          INTEGER, INTENT(OUT) :: is, ia
          INTEGER :: i, nat
          nat = 0
          ia  = 0
          is  = 0
          LOOP: DO i = 1, SIZE( na )
            IF( (nat + na(i) ) >= isa ) THEN
              ia = isa - nat
              is = i
              EXIT LOOP
            ELSE
              nat = nat + na(i)
            END IF
          END DO LOOP
          RETURN
        END SUBROUTINE


        SUBROUTINE atoms_type_init(atoms, staur, ismbl, label, pma, na, nsp, h)
          USE cell_base, ONLY: s_to_r
          TYPE (atoms_type) :: atoms
          REAL(dbl), INTENT(IN) :: staur(:,:)
          LOGICAL, INTENT(IN) :: ismbl(:,:)
          REAL(dbl), INTENT(IN) :: pma(:), h(3,3)
          INTEGER, INTENT(IN) :: na(:), nsp
          CHARACTER(LEN=3), INTENT(IN) :: label(:)
          INTEGER :: nax, nat
          INTEGER :: ierr, is, ia, isa, isatop


          nat = SUM( na( 1 : nsp )  )
          nax = MAXVAL( na( 1 : nsp ) )

          IF( SIZE( na ) < nsp ) &
            CALL errore(' atoms_type_init ', ' wrong na dimensions ', 1)
          IF( SIZE( pma ) < nsp ) &
            CALL errore(' atoms_type_init ', ' wrong pma dimensions ', 1)

          IF( nsp < 1 ) THEN
            CALL errore(' atoms_type_init ', ' nsp less than one ', 3)
          END IF
          IF( nax < 1 ) THEN
            CALL errore(' atoms_type_init ', ' nax less than one ', 4)
          END IF
          IF( nat < 1 ) THEN
            CALL errore(' atoms_type_init ', ' nat less than one ', 5)
          END IF
          IF( ( nat > SIZE(ismbl, 2) ) ) THEN
            CALL errore(' atoms_type_init ', ' invalid nat ', 6)
          END IF
          IF( ( nat > SIZE(staur, 2) ) ) THEN
            CALL errore(' atoms_type_init ', ' invalid nat ', 6)
          END IF

          atoms%nsp = nsp
          atoms%nat = nat
          atoms%nax = nax
          atoms%ekint = 0.0d0

          isa = 1
          atoms%taus = 0.0d0
          atoms%vels = 0.0d0
          atoms%for = 0.0d0
          atoms%mobile = 0
          atoms%ityp = 0
          atoms%tscfor = .FALSE.

          DO is = 1, nsp
            atoms%na(is)    = na(is)
            atoms%m(is)     = pma(is)

            atoms%isa(is)  = isa
            isatop = isa + na(is) - 1

            atoms%label(is)           = TRIM( label(is) )
            atoms%taus(1:3,isa:isatop) = staur(1:3,isa:isatop)
            WHERE( ismbl(1:3,isa:isatop) ) atoms%mobile(1:3,isa:isatop) = 1
            atoms%ityp(isa:isatop)    = is
            atoms%dof(is)             = MAX( COUNT( atoms%mobile(1:3,isa:isatop) == 1 ), 1 )
            atoms%ekin(is)            = 0.0d0

            isa = isa + na(is)
          END DO

          CALL s_to_r( atoms%taus, atoms%taur, atoms%na, atoms%nsp, h )

          atoms%doft = MAX( SUM( atoms%dof(1:nsp) )-3, 1 )

          RETURN
        END SUBROUTINE

!=----------------------------------------------------------------------------=!

        SUBROUTINE deallocate_atoms_type(atoms)
          TYPE (atoms_type) :: atoms
            INTEGER :: is, ierr
          RETURN
        END SUBROUTINE

!=----------------------------------------------------------------------------=!

        SUBROUTINE allocate_constrains(constrains, &
          na, cindexes, cvalues, ctype, tol, nc)
          TYPE (constrains_class) :: constrains
          INTEGER, INTENT(IN) :: na(:)
          INTEGER, INTENT(IN) :: cindexes(:,:)
          REAL(dbl), INTENT(IN) :: cvalues(:)
          INTEGER, INTENT(IN) :: ctype(:)
          REAL(dbl), INTENT(IN) :: tol
          INTEGER, INTENT(IN) :: nc
          INTEGER :: ic, is, ia, isa, ierr, nsp
          
          
          constrains%nc        = nc
          constrains%tolerance = tol
          nsp = SIZE( na )

          IF( nsp < 1 ) THEN
            CALL errore(' allocate_constrains ',' wrong nsp, less than 1 ',7)
          END IF

          IF( nc > 0 ) THEN

            ALLOCATE( constrains%tp( nc ), STAT=ierr )
            IF( ierr /= 0 ) THEN
              CALL errore(' allocate_constrains ',' allocating tp ', ABS(ierr) )
            END IF

            DO ic = 1, nc
              constrains%tp(ic)%what = ctype(ic)
              IF ( ( constrains%tp(ic)%what == 1 ) .OR. &
                   ( constrains%tp(ic)%what == 2 )        ) THEN

                ALLOCATE(constrains%tp(ic)%distance)
                IF( ierr /= 0 ) THEN
                  CALL errore(' allocate_constrains ',' allocating distance ', ABS(ierr) )
                END IF

                CALL specie_index(cindexes(1,ic), na, is, ia)
                constrains%tp(ic)%distance%is1 = is
                constrains%tp(ic)%distance%ia1 = ia

                CALL specie_index(cindexes(2,ic), na, is, ia)
                constrains%tp(ic)%distance%ia2 = ia
                constrains%tp(ic)%distance%is2 = is

                SELECT CASE (constrains%tp(ic)%what)
                  CASE (2)
                    constrains%tp(ic)%distance%val = cvalues(ic)
                  CASE DEFAULT
                    constrains%tp(ic)%distance%val = 0.d0
                END SELECT

              ELSE

                CALL errore(' allocate_constrains ',' unknown constrain type ', 8 )

              END IF

            END DO
          END IF
          RETURN
        END SUBROUTINE allocate_constrains

!=----------------------------------------------------------------------------=!

        SUBROUTINE deallocate_constrains(constrains)
          TYPE (constrains_class) :: constrains
          INTEGER :: ierr, ic
          IF(ASSOCIATED(constrains%tp)) THEN
            DO ic = 1, SIZE(constrains%tp)
              IF(ASSOCIATED(constrains%tp(ic)%distance)) THEN
                DEALLOCATE(constrains%tp(ic)%distance, STAT=ierr)
                IF( ierr /= 0 ) &
                  CALL errore(' deallocate_constrains ', ' deallocating %tp(ic)%distance ', ABS(ierr) )
              END IF
            END DO
            DEALLOCATE(constrains%tp, STAT=ierr)
            IF( ierr /= 0 ) &
              CALL errore(' deallocate_constrains ', ' deallocating %tp ', ABS(ierr) )
          END IF
        END SUBROUTINE deallocate_constrains



!=----------------------------------------------------------------------------=!
      END MODULE 
!=----------------------------------------------------------------------------=!

