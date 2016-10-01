!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
FUNCTION read_config_from_file(nat, at_old, omega_old, lmovecell, at, bg, &
     omega, tau) RESULT (ierr)
  !-----------------------------------------------------------------------
  !
  USE kinds,          ONLY : DP
  USE io_global,      ONLY : stdout
  USE io_files,       ONLY : tmp_dir, prefix
!
  USE pw_restart_new,    ONLY  : pw_readschema_file, init_vars_from_schema
  USE qes_types_module,     ONLY :  output_type, input_type, parallel_info_type, general_info_type
  USE qes_libs_module,      ONLY :  qes_reset_output, qes_reset_input, qes_reset_general_info, qes_reset_parallel_info 
!
  USE pw_restart,     ONLY : pw_readfile
!
  !
  IMPLICIT NONE
  !
  LOGICAL,INTENT(in)     :: lmovecell
  INTEGER,INTENT(in)     :: nat
  REAL(DP),INTENT(inout) :: at_old(3,3), omega_old
  REAL(DP),INTENT(inout) :: at(3,3), bg(3,3), omega
  REAL(DP),INTENT(inout) :: tau(3,nat)
  INTEGER :: ierr
!
#if defined(__XSD)
  TYPE ( output_type), ALLOCATABLE   :: output_obj
  TYPE ( input_type ), ALLOCATABLE   :: input_obj 
  TYPE (parallel_info_type),ALLOCATABLE :: parinfo_obj
  TYPE (general_info_type ),ALLOCATABLE :: geninfo_obj 
#endif

  !
  !
  WRITE( stdout, '(/5X,"Atomic positions and unit cell read from directory:", &
                &  /,5X,A)') TRIM( tmp_dir ) // TRIM( prefix ) // ".save/"
  !
  ! ... check if restart file is present, if yes read config parameters
  !
#if defined(__XSD) 
  ALLOCATE (output_obj, input_obj, parinfo_obj, geninfo_obj) 
  CALL pw_readschema_file ( ierr, output_obj, input_obj, parinfo_obj, geninfo_obj)
  IF (ierr == 0 ) THEN 
     CALL init_vars_from_schema ( 'config', ierr, output_obj, input_obj, parinfo_obj, geninfo_obj ) 
     CALL qes_reset_output(output_obj) 
     CALL qes_reset_input (input_obj)
     CALL qes_reset_parallel_info (parinfo_obj) 
     CALL qes_reset_general_info  (geninfo_obj) 
  END IF 
  DEALLOCATE ( output_obj, input_obj, parinfo_obj, geninfo_obj ) 
#else
  CALL pw_readfile( 'config', ierr )
#endif
  !
  IF ( ierr > 0 ) THEN
     !
     WRITE( stdout, '(5X,"Nothing found: ", &
                       & "using input atomic positions and unit cell",/)' )
     RETURN
     !
  END IF
  !
  WRITE( stdout, * )
  !
  IF ( lmovecell ) THEN
     !
     ! ... input value of at and omega (currently stored in xxx_old variables)
     ! ... must be used to initialize G vectors and other things
     ! ... swap xxx and xxx_old variables and scale the atomic position to the
     ! ... input cell shape in order to check the symmetry.
     !
     CALL cryst_to_cart( nat, tau, bg, - 1 )
     !
     CALL dswap( 9, at, 1, at_old,1  )
     CALL dswap( 1, omega, 1, omega_old, 1 )
     !
     CALL cryst_to_cart( nat, tau, at, + 1 )
     !
     CALL recips( at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )
     !
  END IF
  !
  RETURN
  !
END FUNCTION read_config_from_file
