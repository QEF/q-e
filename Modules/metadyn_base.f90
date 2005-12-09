!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!#define __USE_PBC
!
!----------------------------------------------------------------------------
MODULE metadyn_base
  !----------------------------------------------------------------------------
  !
  ! ... this module contains the core methods used to implement meta-dynamics
  !
  ! ... meta-dynamics is implemented following these two references:
  !
  ! ... 1) A. Laio and M. Parrinello; PNAS 99, 12562 (2002);
  ! ... 2) C. Micheletti, A. Laio, and M Parrinello; PRL 92, 17061 (2004).
  !
  ! ... code written by Carlo Sbraccia (2005)
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: set_target, &
            add_gaussians, &
            metadyn_init, &
            evolve_collective_vars
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE metadyn_init( progname, tau )
      !------------------------------------------------------------------------
      !
      USE kinds,              ONLY : DP
      USE input_parameters,   ONLY : restart_mode
      USE constraints_module, ONLY : nconstr, target
      USE control_flags,      ONLY : nstep, ndr
      USE constants,          ONLY : bohr_radius_angs
      USE cell_base,          ONLY : at, alat
      USE metadyn_vars,       ONLY : fe_grad, g_amplitude, max_metadyn_iter, &
                                     metadyn_history, metadyn_fmt, fe_step,  &
                                     first_metadyn_iter, gaussian_pos
      USE metadyn_io,         ONLY : read_metadyn_restart
      USE parser,             ONLY : delete_if_present, int_to_char
      USE io_files,           ONLY : tmp_dir, scradir, prefix, iunaxsf, iunmeta
      USE io_global,          ONLY : stdout, ionode, ionode_id
      USE mp,                 ONLY : mp_bcast
      USE xml_io_base,        ONLY : restart_dir
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)    :: progname
      REAL(DP),         INTENT(INOUT) :: tau(:,:)
      !
      INTEGER            :: idum, i
      REAL(DP)           :: rdum
      LOGICAL            :: file_exists
      CHARACTER(LEN=256) :: dirname
      CHARACTER(LEN=4)   :: c_nconstr
      CHARACTER(LEN=16)  :: fe_step_fmt
      !
      !
      c_nconstr  = int_to_char( nconstr )
      !
      metadyn_fmt = '(I4,' // TRIM( c_nconstr ) // '(2X,F10.5),2X,F14.8,' // &
                  & TRIM( c_nconstr ) // '(2X,F10.5),' // &
                  & TRIM( c_nconstr ) // '(2X,F10.7))'
      !
      IF ( nstep < 1 ) CALL errore( 'metadyn_init', 'nstep < 1', 1 )
      !
      max_metadyn_iter = nstep
      !
      IF ( restart_mode == 'from_scratch' ) THEN
         !
         IF ( ionode ) THEN
            !
            OPEN( UNIT = iunaxsf, &
                  FILE = TRIM( prefix ) // ".axsf", STATUS = 'UNKNOWN' )
            !
            WRITE( UNIT = iunaxsf, &
                   FMT = '(" ANIMSTEPS ",I5)' ) max_metadyn_iter
            !
            WRITE( UNIT = iunaxsf, FMT = '(" CRYSTAL ")' )
            WRITE( UNIT = iunaxsf, FMT = '(" PRIMVEC ")' )
            WRITE( UNIT = iunaxsf, FMT = '(3F14.10)' ) &
                at(1,1) * alat * bohr_radius_angs, &
                at(2,1) * alat * bohr_radius_angs, &
                at(3,1) * alat * bohr_radius_angs
            WRITE( UNIT = iunaxsf, FMT = '(3F14.10)' ) &
                at(1,2) * alat * bohr_radius_angs, &
                at(2,2) * alat * bohr_radius_angs, &
                at(3,2) * alat * bohr_radius_angs
            WRITE( UNIT = iunaxsf, FMT = '(3F14.10)' ) &
                at(1,3) * alat * bohr_radius_angs, &
                at(2,3) * alat * bohr_radius_angs, &
                at(3,3) * alat * bohr_radius_angs
            !
         END IF
         !
         CALL delete_if_present( TRIM( prefix ) // '.metadyn' )
         !
         IF ( ionode ) THEN
            !
            OPEN( UNIT = iunmeta, &
                  FILE = TRIM( prefix ) // '.metadyn', STATUS = 'NEW' )
            !
            WRITE( iunmeta, '(2(2X,I5))' ) nconstr, max_metadyn_iter
            WRITE( iunmeta, '(2(2X,F12.8))' ) g_amplitude
            !
            fe_step_fmt = '(' // TRIM( c_nconstr ) // '(2X,F12.8))'
            !
            WRITE( iunmeta, fe_step_fmt ) fe_step(:)
            !
         END IF
         !
         first_metadyn_iter = 0
         !
      ELSE
         !
         ! ... restarting from file
         !
         IF ( progname == 'PW' ) THEN
            !
            dirname = TRIM( tmp_dir ) // TRIM( prefix ) // '.new-save'
            !
         ELSE IF ( progname == 'CP' ) THEN
            !
            dirname = restart_dir( scradir, ndr )
            !
         ELSE
            !
            CALL errore( 'metadyn_init', &
                         'wrong calling program: ' // TRIM( progname ), 1 )
            !
         END IF
         !
         CALL read_metadyn_restart( dirname, tau, alat )
         !
         IF ( ionode ) THEN
            !
            OPEN( UNIT = iunaxsf, FILE = TRIM( prefix ) // ".axsf", &
                  STATUS = 'UNKNOWN', ACTION = 'WRITE', POSITION = 'APPEND' )
            OPEN( UNIT = iunmeta, FILE = TRIM( prefix ) // '.metadyn', &
                  STATUS = 'UNKNOWN', ACTION = 'WRITE', POSITION = 'APPEND' )
            !
         END IF
         !
      END IF
      !
      IF ( first_metadyn_iter == max_metadyn_iter ) THEN
         !
         WRITE( stdout, '(/,5X,"Simulation already completed",/)' )
         !
         CLOSE( UNIT = iunmeta, STATUS = 'KEEP' )
         !
         CALL stop_run( .FALSE. )
         !
      END IF
      !
      gaussian_pos(:) = target(:)
      !
      RETURN
      !
    END SUBROUTINE metadyn_init
    !
    !------------------------------------------------------------------------
    SUBROUTINE add_gaussians( iter )
      !------------------------------------------------------------------------
      !
      USE constraints_module, ONLY : nconstr
      USE metadyn_vars,       ONLY : metadyn_history, fe_grad, fe_step, &
                                     dfe_acc, g_amplitude
      USE basic_algebra_routines
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: iter
      !
      INTEGER               :: i
      REAL(DP), ALLOCATABLE :: delta(:)
      !
      !
      ! ... history dependent term
      !
      IF ( iter == 1 ) RETURN
      !
      ALLOCATE( delta( nconstr ) )
      !
      dfe_acc = 0.D0
      !
      DO i = 1, iter - 1
         !
         delta = metadyn_history(:,iter) - metadyn_history(:,i)
         !
         dfe_acc(:) = dfe_acc(:) + delta(:) / fe_step(:)**2 * &
                      EXP( - SUM( delta(:)**2 / ( 2.D0 * fe_step(:)**2 ) ) )
         !
      END DO
      !
      fe_grad(:) = fe_grad(:) - g_amplitude * dfe_acc(:)
      !
      DEALLOCATE( delta )
      !
      RETURN
      !
    END SUBROUTINE add_gaussians
    !
    !------------------------------------------------------------------------
    SUBROUTINE evolve_collective_vars( norm_fe_grad )
      !------------------------------------------------------------------------
      !
      USE constraints_module, ONLY : nconstr, constr_type, target, dmax
      USE metadyn_vars,       ONLY : fe_grad, fe_step, new_target, &
                                     to_target, shake_nstep, gaussian_pos
      !
      IMPLICIT NONE
      !
      REAL(DP), INTENT(IN) :: norm_fe_grad
      !
      INTEGER  :: i
      REAL(DP) :: step
      !
      REAL(DP), EXTERNAL :: rndm
      !
      !
      fe_grad(:) = fe_grad(:) / norm_fe_grad
      !
      DO i = 1, nconstr
         !
         gaussian_pos(i) = target(i) - fe_step(i) * fe_grad(i)
         !
         step = ( 1.D0 + 0.5D0*rndm() ) * fe_step(i)
         !
         new_target(i) = target(i) - step * fe_grad(i)
         !
#if defined (__USE_PBC)
         !
         IF ( constr_type(i) == 3 ) THEN
            !
            ! ... for constraints on distances we must use the minimum 
            ! ... image convenction
            !
            new_target(i) = new_target(i) - MAX( 0.D0, new_target(i) - dmax )
            !
         END IF
         !
#endif
         !
      END DO
      !
      to_target(:) = ( new_target(:) - target(:) ) / DBLE( shake_nstep )
      !
      RETURN
      !
    END SUBROUTINE evolve_collective_vars
    !
    !------------------------------------------------------------------------
    SUBROUTINE set_target()
      !------------------------------------------------------------------------
      !
      USE metadyn_vars,       ONLY : to_target, to_new_target, shake_nstep
      USE constraints_module, ONLY : target
      !
      !
      IF ( to_new_target ) THEN
         !
         target(:) = target(:) + to_target(:)
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE set_target
    !
END MODULE metadyn_base
