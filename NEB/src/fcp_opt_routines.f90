!
! Copyright (C) 2003-2017 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE fcp_opt_routines
  !---------------------------------------------------------------------------
  !
  ! ... This module contains all subroutines and functions needed for
  ! ... the optimisation of the FCPs.
  !
  ! ... Written by Carlo Sbraccia ( 2003-2006 )
  !
  ! ... Newton algorithm is implemented by S. Nishihara ( 2016-2017 )
  !
  USE constants,      ONLY : e2
  USE kinds,          ONLY : DP
  USE io_global,      ONLY : meta_ionode, meta_ionode_id
  USE mp,             ONLY : mp_bcast
  USE mp_world,       ONLY : world_comm
  USE path_variables, ONLY : num_of_images
  USE fcp_variables,  ONLY : fcp_mu, lfcp_linmin, lfcp_newton, &
                             fcp_nelec, fcp_ef, fcp_dos, mdiist, &
                             nelec0, force0, firstcall
  USE mdiis,          ONLY : update_by_mdiis
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: fcp_opt_perform, fcp_opt_scale
  !
  CONTAINS
     !
     !----------------------------------------------------------------------
     SUBROUTINE fcp_opt_perform()
        !----------------------------------------------------------------------
        !
        USE fcp_module, ONLY : fcp_check
        !
        IMPLICIT NONE
        !
        REAL(DP) :: step_max
        REAL(DP) :: capacitance
        !
        CALL fcp_check( .TRUE. )
        !
        ! ... evaluate maximum step
        !
        IF ( meta_ionode ) THEN
           !
           CALL fcp_capacitance( capacitance )
           capacitance = e2 * capacitance
           !
           step_max = ABS( capacitance * 0.05_DP ) ! = C * 0.1Ry
           !
        END IF
        !
        CALL mp_bcast( step_max, meta_ionode_id, world_comm )
        !
        IF ( lfcp_linmin ) THEN
           !
           ! ... perform Line-Minimization
           !
           CALL fcp_line_minimization( step_max )
           !
        ELSE IF ( lfcp_newton ) THEN
           !
           ! ... perform Newton-Raphson
           !
           CALL fcp_newton( step_max )
           !
        END IF
        !
     END SUBROUTINE fcp_opt_perform
     !
     !----------------------------------------------------------------------
     SUBROUTINE fcp_line_minimization( step_max )
        !----------------------------------------------------------------------
        !
        USE constants,      ONLY : eps16
        USE path_variables, ONLY : frozen
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: step_max
        !
        INTEGER  :: image
        REAL(DP) :: ef, dos, force, step
        REAL(DP) :: nelec, nelec_new
        !
        IF ( meta_ionode ) THEN
           !
           DO image = 1, num_of_images
              !
              IF ( frozen(image) ) CYCLE
              !
              nelec = fcp_nelec(image)
              ef    = fcp_ef   (image)
              dos   = fcp_dos  (image)
              force = fcp_mu - ef
              !
              IF ( firstcall(image) ) THEN
                 !
                 ! ... initialization
                 !
                 firstcall(image) = .FALSE.
                 !
                 nelec0(image) = nelec
                 force0(image) = force
                 !
              END IF
              !
              IF ( ABS( force0(image) - force ) < eps16 ) THEN
                 !
                 ! ... Steepest-Descent
                 !
                 step = 0.0_DP
                 CALL step_newton( dos, force, step )
                 !
                 nelec_new = nelec + step
                 !
              ELSE
                 !
                 ! ... Line-Minimization
                 !
                 nelec_new = (nelec * force0(image) - nelec0(image) * force) &
                           & / (force0(image) - force)
                 !
              END IF
              !
              ! ... save #electrons and force
              !
              nelec0(image) = nelec
              force0(image) = force
              !
              ! ... update #electrons
              !
              step = nelec_new - nelec
              step = MIN( step, +step_max )
              step = MAX( step, -step_max )
              !
              fcp_nelec(image) = nelec + step
              !
           END DO
           !
        END IF
        !
        CALL mp_bcast( fcp_nelec, meta_ionode_id, world_comm )
        !
        RETURN
        !
     END SUBROUTINE fcp_line_minimization
     !
     !----------------------------------------------------------------------
     SUBROUTINE fcp_newton( step_max )
        !----------------------------------------------------------------------
        !
        USE path_variables, ONLY : frozen
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: step_max
        !
        INTEGER               :: image
        REAL(DP)              :: ef, dos, force, step
        REAL(DP)              :: nelec, nelec_new
        REAL(DP), ALLOCATABLE :: nelec1(:)
        REAL(DP), ALLOCATABLE :: step1(:)
        !
        ALLOCATE(nelec1(num_of_images))
        ALLOCATE(step1(num_of_images))
        !
        IF ( meta_ionode ) THEN
           !
           DO image = 1, num_of_images
              !
              ! ... current #electrons and Newton's steps
              !
              nelec = fcp_nelec(image)
              ef    = fcp_ef   (image)
              dos   = fcp_dos  (image)
              force = fcp_mu - ef
              !
              nelec1(image) = nelec
              CALL step_newton( dos, force, step1(image) )
              !
           END DO
           !
           ! ... apply DIIS
           !
           CALL update_by_mdiis( mdiist, nelec1, step1 )
           !
           DO image = 1, num_of_images
              !
              IF ( frozen(image) ) CYCLE
              !
              ! ... update #electrons
              !
              nelec     = fcp_nelec(image)
              nelec_new = nelec1(image)
              !
              step = nelec_new - nelec
              step = MIN( step, +step_max )
              step = MAX( step, -step_max )
              !
              fcp_nelec(image) = nelec + step
              !
           END DO
           !
        END IF
        !
        CALL mp_bcast( fcp_nelec, meta_ionode_id, world_comm )
        !
        DEALLOCATE(nelec1)
        DEALLOCATE(step1)
        !
        RETURN
        !
     END SUBROUTINE fcp_newton
     !
     !----------------------------------------------------------------------
     SUBROUTINE step_newton( dos, force, step )
       !----------------------------------------------------------------------
       !
       USE constants, ONLY : eps4
       !
       IMPLICIT NONE
       !
       REAL(DP), INTENT(IN)  :: dos
       REAL(DP), INTENT(IN)  :: force
       REAL(DP), INTENT(OUT) :: step
       !
       REAL(DP) :: hess
       REAL(DP) :: capacitance
       !
       hess = dos
       !
       CALL fcp_capacitance( capacitance )
       capacitance = e2 * capacitance
       !
       IF ( capacitance > eps4 ) THEN
          !
          hess = MIN( hess, capacitance )
          !
       END IF
       !
       IF ( hess > eps4 ) THEN
          !
          step = hess * force
          !
       ELSE
          !
          CALL errore( 'step_newton', 'capacitance is not positive', 1 )
          !
          step = 0.0_DP
          !
       END IF
       !
     END SUBROUTINE step_newton
     !
     !----------------------------------------------------------------------
     REAL(DP) FUNCTION fcp_opt_scale()
       !----------------------------------------------------------------------
       !
       ! ... see Modules/bfgs_module.f90
       ! ... this is same as sqrt of the metric.
       !
       USE fcp_variables,  ONLY : fcp_max_volt
       !
       IMPLICIT NONE
       !
       REAL(DP) :: capacitance
       REAL(DP),PARAMETER :: max_step = 0.6_DP
       !
       ! ... set capacitance
       !
       CALL fcp_capacitance( capacitance )
       capacitance = e2 * capacitance
       !
       ! ... set scaling factor
       !
       fcp_opt_scale = max_step / (fcp_max_volt * capacitance)
       !
     END FUNCTION fcp_opt_scale
     !
END MODULE fcp_opt_routines
