!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! ... This file contains "sub-main" subroutines that drive the different
! ... kinds of "meta"/"non-meta" dynamics
!
! ... set tprnfor = .TRUE. to get atomic forces even if the atoms do not move
!
! ... set ion_positions = 'from_input' and rd_pos = +your_positions+
! ... to force cprmain to compute forces for +your_position+ configuration
!
!----------------------------------------------------------------------------
SUBROUTINE smd_loop( nloop )
  !----------------------------------------------------------------------------
  !
  USE kinds,            ONLY : dbl
  USE ions_base,        ONLY : nat
  USE control_flags,    ONLY : tprnfor
  USE input_parameters, ONLY : ion_positions, rd_pos, num_of_images
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nloop
  !
  INTEGER                     :: iloop, sm_p 
  REAL(KIND=dbl), ALLOCATABLE :: tau(:,:,: )
  REAL(KIND=dbl), ALLOCATABLE :: fion(:,:,:)
  REAL(KIND=dbl), ALLOCATABLE :: etot(:) 
  !
  !
  sm_p = num_of_images - 1
  !
  IF ( nat > 0 .AND. sm_p > 0 ) THEN
     !
     ALLOCATE( tau(  3, nat, 0:sm_p ) )
     ALLOCATE( fion( 3, nat, 0:sm_p ) )
     ALLOCATE( etot( 0:sm_p ) )
     !
  ELSE
     !
     CALL errore( ' smd_loop ', ' nat or sm_p less or equal 0 ', 1 )
     !
  END IF
  !
  ! ... initialize g-vectors, fft grids
  !
  CALL init_dimensions()
  !
  DO iloop = 1, nloop
     !
     CALL smdmain( tau, fion, etot, nat )
     !
     CALL memstat( 1 )
     !
  END DO
  !
  DEALLOCATE( tau, fion, etot )
  !
  RETURN
  !
END SUBROUTINE smd_loop
!
!----------------------------------------------------------------------------
SUBROUTINE neb_loop( iloop, program_name )
  !----------------------------------------------------------------------------
  !
  USE path_base,        ONLY : initialize_path, search_mep
  USE path_routines,    ONLY : iosys_path
  USE path_io_routines, ONLY : io_path_start
  !
  IMPLICIT NONE
  !
  INTEGER,          INTENT(IN) :: iloop
  CHARACTER(LEN=*), INTENT(IN) :: program_name
  !
  !
  CALL iosys_path()
  !
  CALL io_path_start()
  !
  CALL initialize_path( 'FP' )
  !
  ! ... this routine does all the NEB job
  !
  CALL search_mep()
  !
  RETURN
  !
END SUBROUTINE neb_loop
!
!----------------------------------------------------------------------------
SUBROUTINE cpr_loop( nloop )
  !----------------------------------------------------------------------------
  !
  USE kinds,            ONLY : dbl
  USE ions_base,        ONLY : nat
  USE control_flags,    ONLY : tprnfor
  USE input_parameters, ONLY : ion_positions, rd_pos
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nloop
  !
  INTEGER                     :: iloop
  REAL(KIND=dbl), ALLOCATABLE :: tau(:,:)
  REAL(KIND=dbl), ALLOCATABLE :: fion(:,:)
  REAL(KIND=dbl)              :: etot
  !
  !
  IF ( nat > 0 ) THEN
     !
     ALLOCATE( tau(  3, nat ) )
     ALLOCATE( fion( 3, nat ) )
     !
  ELSE
     !
     CALL errore( ' cpr_loop ', ' nat less or equal 0 ', 1 )
     !
  END IF
  !
  CALL init_run()
  !
  DO iloop = 1, nloop
     !
     CALL cprmain( tau(1,1), fion(1,1), etot )
     !
     CALL memstat( 1 )
     !
  END DO
  !
  DEALLOCATE( tau, fion )
  !
  RETURN
  !
END SUBROUTINE cpr_loop
!
!----------------------------------------------------------------------------
SUBROUTINE fpmd_loop( iloop )
  !----------------------------------------------------------------------------
  !
  USE kinds,            ONLY : dbl
  USE main_module,      ONLY : cpmain
  USE input_parameters, ONLY : restart_mode
  USE io_files,         ONLY : outdir
  USE ions_base,        ONLY : nat
  USE io_global,        ONLY : stdout
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iloop
  !
  CHARACTER(LEN=256), SAVE :: outdir_orig
  !
  REAL(KIND=dbl), ALLOCATABLE :: tau(:,:)
  REAL(KIND=dbl), ALLOCATABLE :: fion(:,:)
  REAL(KIND=dbl)              :: etot
  !
  !
  IF ( iloop == 1 ) outdir_orig = outdir
  !
  IF ( iloop > 1 ) THEN
     !
     restart_mode = 'restart'
     !
  END IF
  !
  SELECT CASE ( iloop )
  CASE ( 1 )
     !
     outdir = TRIM( outdir_orig ) // '/' // 'image01'
     !
  CASE ( 2 )
     !  
     outdir = TRIM( outdir_orig ) // '/' // 'image02'
     !
  END SELECT
  !
  ! ... Car-Parrinello main routine
  !
  IF ( nat > 0 ) THEN
     !
     ALLOCATE( tau(  3, nat ) )
     ALLOCATE( fion( 3, nat ) )
     !
  ELSE
     !
     CALL errore( ' fpmd_loop ', ' nat less or equal 0 ', 1 )
     !
  END IF
  !
  CALL cpmain( tau, fion, etot )
  !
  RETURN
  !
END SUBROUTINE fpmd_loop
