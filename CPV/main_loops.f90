!
! Copyright (C) 2002 CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------------=!
!
!  This file contains "sub-main" subroutines that drive the different kind of 
!  "meta", "non-meta" dynamics
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE wf_loop( nloop )
      USE ions_base, ONLY: nat
      USE control_flags, ONLY: tprnfor
      USE input_parameters, ONLY: ion_positions, rd_pos
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nloop
      INTEGER :: iloop
      REAL(kind=8), ALLOCATABLE :: tau( :, : )
      REAL(kind=8), ALLOCATABLE :: fion( :, : )
      REAL(kind=8) :: etot

      IF( nat > 0 ) THEN
        ALLOCATE( tau( 3, nat ), fion( 3, nat ) )
      ELSE
        CALL errore( ' cpr_loop ', ' nat less or equal 0 ', 1 )
      END IF

      ! ... set tprnfor = .true. to get atomic forces
      ! ... even if the atoms do not move

      ! ... set ion_positions = 'from_input'
      ! ... and rd_pos = +your_positions+
      ! ... to force cprmain to compute forces for
      ! ... +your_position+ configuration

      DO iloop = 1, nloop
        call cprmain( tau(1,1), fion(1,1), etot)
        call memstat( 1 )
      END DO

      DEALLOCATE( tau, fion )

      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------------=!

    SUBROUTINE smd_loop( nloop )
      USE ions_base, ONLY: nat
      USE control_flags, ONLY: tprnfor
      USE input_parameters, ONLY: ion_positions, rd_pos, num_of_images
      
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nloop
      INTEGER :: iloop, sm_p 
      REAL(kind=8), ALLOCATABLE :: tau( :, :, : )
      REAL(kind=8), ALLOCATABLE :: fion( :, :, : )
      REAL(kind=8), ALLOCATABLE :: etot( : ) 
      
      sm_p = num_of_images -1
      
      IF( nat > 0 .AND. sm_p > 0 ) THEN
        ALLOCATE( tau( 3, nat, 0:sm_p ), fion( 3, nat, 0:sm_p ), etot( 0:sm_p ) )
      ELSE
        CALL errore( ' smd_loop ', ' nat or sm_p less or equal 0 ', 1 )
      END IF 
      
      ! ... set tprnfor = .true. to get atomic forces
      ! ... even if the atoms do not move
      
      ! ... set ion_positions = 'from_input'
      ! ... and rd_pos = +your_positions+
      ! ... to force cprmain to compute forces for
      ! ... +your_position+ configuration

      DO iloop = 1, nloop
        call smdmain( tau, fion, etot, nat )
        call memstat( 1 )
      END DO

      DEALLOCATE( tau, fion, etot )

      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------------=!

    SUBROUTINE neb_loop( iloop )
      !
      USE kinds
      USE io_global,        ONLY : ionode, stdout
      USE path_variables,   ONLY : conv_path
      USE path_variables,   ONLY : path_deallocation
      USE path_base,        ONLY : initialize_path, search_mep
      USE path_routines,    ONLY : iosys_path
      USE path_io_routines, ONLY : write_output
      USE ions_base,        ONLY : deallocate_ions_base
      !
      IMPLICIT NONE
       !
      INTEGER :: iloop
      !
      !
      ! ... stdout is connected to a file ( specific for each image )
      ! ... via unit 17
      !
      IF( ionode ) THEN
        !
        stdout = 17
        !
      END IF
      !
      CALL iosys_path()
      !
      CALL initialize_path( 'FP' )
      !
      ! ... this routine does all the NEB job
      !
      CALL search_mep()
      !
      ! ... output is written
      !
      CALL write_output()
      !
      CALL deallocate_ions_base( )
      !
      CALL path_deallocation( 'neb' )
      !
      ! ... stdout is reconnected to standard output
      !
      stdout = 6
      !
      RETURN
      !
    END SUBROUTINE neb_loop

!=----------------------------------------------------------------------------=!

    SUBROUTINE cpr_loop( nloop )
      USE ions_base, ONLY: nat
      USE control_flags, ONLY: tprnfor
      USE input_parameters, ONLY: ion_positions, rd_pos
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nloop
      INTEGER :: iloop
      REAL(kind=8), ALLOCATABLE :: tau( :, : )
      REAL(kind=8), ALLOCATABLE :: fion( :, : )
      REAL(kind=8) :: etot

      IF( nat > 0 ) THEN
        ALLOCATE( tau( 3, nat ), fion( 3, nat ) )
      ELSE
        CALL errore( ' cpr_loop ', ' nat less or equal 0 ', 1 )
      END IF

      ! ... set tprnfor = .true. to get atomic forces
      ! ... even if the atoms do not move

      ! ... set ion_positions = 'from_input'
      ! ... and rd_pos = +your_positions+
      ! ... to force cprmain to compute forces for
      ! ... +your_position+ configuration

      DO iloop = 1, nloop
        call cprmain( tau(1,1), fion(1,1), etot)
        call memstat( 1 )
      END DO

      DEALLOCATE( tau, fion )

      RETURN
    END SUBROUTINE


!=----------------------------------------------------------------------------=!


    SUBROUTINE fpmd_loop( iloop )

      USE kinds
      USE main_module,      ONLY: cpmain
      USE input_parameters, ONLY: restart_mode
      USE io_files,         ONLY: outdir
      USE ions_base,        ONLY: nat
      USE io_global,        ONLY: stdout

      IMPLICIT NONE

      INTEGER :: iloop
      CHARACTER(LEN=256), SAVE :: outdir_orig

      REAL(dbl), ALLOCATABLE :: tau( :, : )
      REAL(dbl), ALLOCATABLE :: fion( :, : )
      REAL(dbl) :: etot


      IF( iloop == 1 ) outdir_orig = outdir

      IF( iloop > 1 ) THEN
        restart_mode = 'restart'
      END IF

      SELECT CASE (iloop)
        CASE ( 1 )
          outdir = TRIM( outdir_orig ) // '/' // 'image01'
        CASE ( 2 )
          outdir = TRIM( outdir_orig ) // '/' // 'image02'
      END SELECT

! ... Car-Parrinello main routine

      IF( nat > 0 ) THEN
        ALLOCATE( tau( 3, nat ), fion( 3, nat ) )
      ELSE
        CALL errore( ' cploop ', ' nat less or equal 0 ', 1 )
      END IF

      CALL cpmain( tau, fion, etot )

      ! WRITE(6,*) '  From cploop, etot = ', etot

      RETURN
    END SUBROUTINE fpmd_loop

