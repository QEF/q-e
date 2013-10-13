!
! Copyright (C) 2013 Quantum ESPRESSO groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!==-----------------------------------------------------------------------==!
MODULE qmmm
  !==---------------------------------------------------------------------==!
  USE io_global,        ONLY : ionode, ionode_id, stdout
  USE mp_world,         ONLY : world_comm
  USE mp,               ONLY : mp_bcast, mp_barrier, mp_abort
  USE kinds,            ONLY : DP
  USE parallel_include
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  ! ... MPI communicator to the QM/MM control process, if MPI is used
  INTEGER :: qmmm_comm = MPI_COMM_NULL
  ! ... number of QM/MM steps
  INTEGER :: qmmm_step = -1
  !
  INTEGER :: qmmm_mode  = -1
  ! mode = <0: QM/MM disabled
  ! mode =  0: not properly set up
  ! mode =  1: mechanical coupling
  ! mode =  2: electrostatic coupling
  !
  ! verbosity level
  INTEGER :: qmmm_verb = -1
  !
  ! message tags. keep consistent with MM code
  INTEGER, PARAMETER :: QMMM_TAG_OTHER=0
  INTEGER, PARAMETER :: QMMM_TAG_SIZE=1
  INTEGER, PARAMETER :: QMMM_TAG_COORD=2
  INTEGER, PARAMETER :: QMMM_TAG_FORCE=3
  !
  ! convert forces to LAMMPS "real" units
  REAL(DP), PARAMETER :: QMMM_FORCE_CONV = 592.91102087727177_DP
  !
  ! Number of atoms of the QM/MM systems
  INTEGER :: nat_qm

  ! buffer for converting forces and positions
  REAL(DP), ALLOCATABLE  :: tmp_buf(:,:)
  ! center of mass of the system
  REAL(DP), DIMENSION(3) :: r0 = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
  LOGICAL :: do_init_r0 = .TRUE. 

  PUBLIC :: qmmm_config, qmmm_initialization, qmmm_shutdown
  PUBLIC :: qmmm_update_positions, qmmm_update_forces, qmmm_add_mm_field

CONTAINS

  ! configure the qm/mm interface
  SUBROUTINE qmmm_config( mode, comm, verbose, step )
    IMPLICIT NONE
    INTEGER, OPTIONAL, INTENT(IN) :: mode, comm, verbose, step

    IF (PRESENT(mode)) qmmm_mode = mode
    IF (PRESENT(comm)) qmmm_comm = comm
    IF (PRESENT(verbose)) qmmm_verb = verbose
    IF (PRESENT(step)) qmmm_step = step

  END SUBROUTINE qmmm_config


  SUBROUTINE qmmm_initialization
    USE input_parameters, ONLY : calculation, nstep, nat
    !
    IMPLICIT NONE
    INTEGER :: ierr

    IF (qmmm_mode < 0) RETURN

    ! send global configuration parameters to all ranks
    CALL mp_bcast(qmmm_mode, ionode_id, world_comm)
    CALL mp_bcast(qmmm_step, ionode_id, world_comm)
    nat_qm = nat

    IF (ionode) THEN
        WRITE(stdout,'(/,5X,A)') "QMMM: Initializing QM/MM interface"
        IF (qmmm_comm /= MPI_COMM_NULL) THEN
            WRITE(stdout,'(5X,A)') "QMMM: Using MPI based communication"
        ELSE
            WRITE(stdout,'(5X,A)') "QMMM: Using MS2 daemon based communication"
        END IF

        IF (qmmm_mode == 0) THEN
            WRITE(stdout,'(5X,A)') "QMMM: Running in dummy mode"
        ELSE IF (qmmm_mode == 1) THEN
            WRITE(stdout,'(5X,A)') "QMMM: Using mechanical coupling"
        ELSE IF (qmmm_mode == 2) THEN
            WRITE(stdout,'(5X,A)') "QMMM: Using electrostatic coupling"
        END IF
    END IF

    ! make sure we have sane settings
    IF (TRIM( calculation ) /= 'md' ) THEN
        if (ionode) &
              WRITE(stdout,'(5X,A)') "QMMM Error: 'md' calculation required."
        CALL mp_abort(255,world_comm)
    END IF

    IF (nstep /= qmmm_step) THEN
        IF (ionode) WRITE(stdout,'(5X,A,I6,A,I6)') &
            'QMMM: Adjusting number of steps from', nstep, ' to', qmmm_step
        nstep = qmmm_step
    END IF

    ! only ionode communicates with MM master
    IF (ionode) THEN
        IF (qmmm_comm /= MPI_COMM_NULL) THEN
#if defined(__MPI)
            CALL mpi_send(nat_qm,1,MPI_INTEGER,0,QMMM_TAG_SIZE,qmmm_comm,ierr)
#else
            WRITE(stdout,*) 'Use of QM/MM requires compilation with MPI'
            STOP 200
#endif
        END IF
    END IF
    CALL mp_bcast(nstep, ionode_id, world_comm)

    ! temporary storage
    ALLOCATE( tmp_buf(3,nat_qm) )
    
  END SUBROUTINE qmmm_initialization

  ! private subroutine
  SUBROUTINE qmmm_center_molecule
    USE cell_base, ONLY : alat, at
    USE ions_base, ONLY : nat
    USE ions_base, ONLY : tau
    IMPLICIT NONE
    LOGICAL, SAVE::firstexec = .TRUE.
    INTEGER:: i  
    ! New geometric center
    REAL(DP), DIMENSION(3):: gc = (/0.0d0, 0.0d0, 0.0d0/)
    REAL(DP), DIMENSION(3):: qm_bc = (/0.5d0, 0.5d0, 0.5d0/)

    IF (firstexec) THEN
        ! Take the geometric center during first call
        r0 = SUM(tau, dim = 2) / nat
        WRITE(stdout,'(5X,A,3F10.6)') 'QMMM: r0(old) ', r0
        r0 = MATMUL(at,qm_bc)
        WRITE(stdout,'(5X,A,3F10.6)') 'QMMM: r0(new) ', r0
        firstexec = .FALSE.
    END IF
    ! Recenter the system.
    gc = SUM(tau, dim = 2) / nat
    ! delta = r0 - r1
    DO i = 1, nat
        tau(1,i) = tau(1,i) - gc(1) + r0(1)
        tau(2,i) = tau(2,i) - gc(2) + r0(2)
        tau(3,i) = tau(3,i) - gc(3) + r0(3)
    END DO

  END SUBROUTINE qmmm_center_molecule


  ! update positions of the QM system from MM-master
  SUBROUTINE qmmm_update_positions
    USE constants, ONLY : bohr_radius_angs
    USE cell_base, ONLY : alat
    USE ions_base, ONLY : tau
    IMPLICIT NONE
    INTEGER :: ierr

    IF (qmmm_mode < 0) RETURN
    IF (ionode .and. (qmmm_verb > 0)) &
         WRITE(stdout,'(/,5X,A)') 'QMMM: update positions'

    ! Receive coordinates (from LAMMPS) and broadcast to all processors
    IF (ionode) THEN
#if defined(__MPI)
        CALL mpi_recv(tau(1,1),3*nat_qm,MPI_DOUBLE_PRECISION, &
              0,QMMM_TAG_COORD,qmmm_comm,MPI_STATUS_IGNORE,ierr)
        ! convert from angstrom to alat units
        tau = tau / (alat * bohr_radius_angs)
#else
        WRITE(stdout,*) 'Use of QM/MM requires compilation with MPI support'
        STOP 201
#endif
        CALL qmmm_center_molecule
    END IF

    CALL mp_bcast(tau, ionode_id, world_comm)

  END SUBROUTINE qmmm_update_positions

  ! communicate forces of the QM system to MM-master
  SUBROUTINE qmmm_update_forces(force)
    REAL(DP), INTENT(IN) :: force(:,:)
    INTEGER :: ierr
    IF (qmmm_mode < 0) RETURN
    IF (ionode .and. (qmmm_verb > 0)) &
        WRITE(stdout,'(/,5X,A)') 'QMMM: update forces'

    IF (ionode) THEN
#if defined(__MPI)
        ! convert from atomic to real units
        tmp_buf = force * QMMM_FORCE_CONV
        CALL mpi_send(tmp_buf,3*nat_qm,MPI_DOUBLE_PRECISION, &
              0,QMMM_TAG_FORCE,qmmm_comm,ierr)
#else
        WRITE(stdout,*) 'Use of QM/MM requires compilation with MPI support'
        STOP 201
#endif
    END IF

  END SUBROUTINE qmmm_update_forces

  ! add electrostatic field of MM system to QM system
  SUBROUTINE qmmm_add_mm_field
    IF (qmmm_mode /= 2) RETURN
    IF (ionode .and. (qmmm_verb > 0)) &
        WRITE(stdout,'(/,5X,A)') 'QMMM: add mm field'
  END SUBROUTINE qmmm_add_mm_field

  ! cleanup of QM/MM. free resources
  SUBROUTINE qmmm_shutdown
    IMPLICIT NONE
    !
    IF (qmmm_mode < 0) RETURN
    IF (ionode) THEN
        WRITE(stdout,'(/,5X,A)') "QMMM: Shutting down QM/MM coupling"
    END IF
    deallocate( tmp_buf )
  END SUBROUTINE qmmm_shutdown
END MODULE qmmm

