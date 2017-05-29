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
  USE mp_global,        ONLY : me_pool, intra_pool_comm ! added
  USE mp,               ONLY : mp_bcast, mp_barrier, mp_abort, mp_sum ! added
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
  INTEGER, PARAMETER :: QMMM_TAG_FORCE2=4
  INTEGER, PARAMETER :: QMMM_TAG_CELL=5
  INTEGER, PARAMETER :: QMMM_TAG_RADII=6
  INTEGER, PARAMETER :: QMMM_TAG_CHARGE=7
  INTEGER, PARAMETER :: QMMM_TAG_TYPE=8
  INTEGER, PARAMETER :: QMMM_TAG_MASS=9
  !
  ! convert forces to LAMMPS "real" units
  REAL(DP), PARAMETER :: QMMM_FORCE_CONV = 592.91102087727177_DP
  !
  ! Number of atoms of the QM/MM systems
  ! buffer for converting forces and positions
  REAL(DP), ALLOCATABLE  :: tmp_buf(:,:)
  ! center of mass of the system
  REAL(DP), DIMENSION(3) :: r0 = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
  LOGICAL :: do_init_r0 = .TRUE. 
  !
  REAL(DP), ALLOCATABLE :: charge(:)
  REAL(DP), ALLOCATABLE :: aradii(:)
  REAL(DP), ALLOCATABLE :: tau_mm(:,:)
  REAL(DP), ALLOCATABLE :: force_mm(:,:)
  REAL(DP), ALLOCATABLE :: force_qm(:,:)
  INTEGER, ALLOCATABLE :: tau_mask(:)
  REAL(DP), ALLOCATABLE :: rc_mm(:)
  REAL(DP), ALLOCATABLE :: charge_mm(:)
  REAL(DP), ALLOCATABLE :: mass(:)
  INTEGER, ALLOCATABLE :: types(:)
  REAL(DP) :: cell_data(9) 
  REAL(DP) :: cell_mm(9) 
  INTEGER  :: nat_mm
  INTEGER  :: nat_qm
  INTEGER  :: nat_all
  INTEGER  :: ntypes
  !

  PUBLIC :: qmmm_config, qmmm_initialization, qmmm_shutdown, qmmm_mode
  PUBLIC :: qmmm_update_positions, qmmm_update_forces, qmmm_add_esf, qmmm_force_esf

CONTAINS

  ! configure the qm/mm interface
  SUBROUTINE qmmm_config( mode, comm, verbose, step )
    IMPLICIT NONE
    INTEGER, OPTIONAL, INTENT(IN) :: mode, comm, verbose, step

    IF (PRESENT(mode)) qmmm_mode = mode
    IF (PRESENT(comm)) qmmm_comm = comm
    IF (PRESENT(verbose)) qmmm_verb = verbose
    IF (PRESENT(step)) qmmm_step = step + 1  ! fix step count discrepancy

  END SUBROUTINE qmmm_config

  !---------------------------------------------------------------------!

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
          CALL errore( 'qmmm_initialization', 'Use of QM/MM requires compilation with MPI', 1 )
#endif
       END IF
    END IF
    CALL mp_bcast(nstep, ionode_id, world_comm)
    ! temporary storage
    ALLOCATE( tmp_buf(3,nat_qm) )
    
  END SUBROUTINE qmmm_initialization

  !---------------------------------------------------------------------!

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

    ! Also do the same on tau_mm, if the electrostatic coupling is enabled,
    ! but without keeping the mm atoms in account in the compute of gc/

    IF( qmmm_mode == 2 ) THEN
       DO i = 1, nat_mm
          tau_mm(1, i) = tau_mm(1,i) - gc(1) + r0(1)
          tau_mm(2, i) = tau_mm(2,i) - gc(2) + r0(2)
          tau_mm(3, i) = tau_mm(3,i) - gc(3) + r0(3)
       ENDDO
    ENDIF

  END SUBROUTINE qmmm_center_molecule

  !---------------------------------------------------------------------!

SUBROUTINE qmmm_minimum_image()
  USE constants, ONLY : bohr_radius_angs, eps8
  USE cell_base, ONLY : alat, at
  USE ions_base, ONLY : nat
  USE ions_base, ONLY : tau
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  INTEGER:: i
  REAL(DP), DIMENSION(3):: qm_bc = (/0.5d0, 0.5d0, 0.5d0/)
  REAL(DP), DIMENSION(3):: at_mm
  REAL(DP), DIMENSION(3):: s
  REAL(DP) :: alat_mm
  !only support orthogonal box, xy = xz = yz = 0.d0
  IF( cell_mm(7) .GT. eps8 .OR. cell_mm(7) .LT. -eps8 .OR. &
      cell_mm(8) .GT. eps8 .OR. cell_mm(8) .LT. -eps8 .OR. &
      cell_mm(9) .GT. eps8 .OR. cell_mm(9) .LT. -eps8 ) THEN
     CALL errore("ms2_minimum_image","Only support orthogonal MM box", 1)
  ENDIF
  !
  at_mm(1) = 1.d0
  at_mm(2) = (cell_mm(5) - cell_mm(2)) / (cell_mm(4) - cell_mm(1))
  at_mm(3) = (cell_mm(6) - cell_mm(3)) / (cell_mm(4) - cell_mm(1))
  alat_mm = (cell_mm(4) - cell_mm(1)) / bohr_radius_angs
  !
  DO i = 1,nat_mm 
     s(1) = tau_mm(1,i) - qm_bc(1)
     s(2) = tau_mm(2,i) - qm_bc(2)
     s(3) = tau_mm(3,i) - qm_bc(3)
     !s(:) = matmul(s(:), 1/box)
     s(:) = s(:) / (alat_mm*at_mm(:)/alat)
     s(:) = s(:) - anint(s(:))
     s(:) = s(:) * (alat_mm*at_mm(:)/alat)
     tau_mm(1,i) = s(1) + qm_bc(1)
     tau_mm(2,i) = s(2) + qm_bc(2)
     tau_mm(3,i) = s(3) + qm_bc(3)
  ENDDO
  !
END SUBROUTINE qmmm_minimum_image


  !---------------------------------------------------------------------!
  ! update positions of the QM system from MM-master

  SUBROUTINE qmmm_update_positions
    USE constants, ONLY : bohr_radius_angs
    USE cell_base, ONLY : alat
    USE ions_base, ONLY : tau
    IMPLICIT NONE
    INTEGER :: ierr,i
    INTEGER :: irecv_buf(8)

    IF (qmmm_mode < 0) RETURN

#if defined(__MPI)
   
    IF (ionode .and. (qmmm_verb > 0)) &
        WRITE(stdout,'(/,5X,A)') 'QMMM: update positions'

    IF( ionode ) THEN
        CALL mpi_recv( irecv_buf, 4, MPI_INTEGER, 0, QMMM_TAG_SIZE, qmmm_comm, MPI_STATUS_IGNORE, ierr )
    END IF
    CALL mp_bcast( irecv_buf, ionode_id, world_comm )
    nat_all = irecv_buf(1)
    nat_qm  = irecv_buf(2)
    nat_mm  = irecv_buf(3)
    ntypes  = irecv_buf(4)
    IF (ionode .and. (qmmm_verb > 0 )) THEN
        WRITE(stdout,*) '    QMMM: nat_all = ', nat_all
        WRITE(stdout,*) '    QMMM: nat_qm  = ', nat_qm  ! num_qm in lammps
        WRITE(stdout,*) '    QMMM: nat_mm  = ', nat_mm  ! num_mm in lammps
        WRITE(stdout,*) '    QMMM: ntypes  = ', ntypes  ! num_mm in lammps
    END IF


    IF( .NOT. ALLOCATED( rc_mm ) ) THEN
        ALLOCATE( rc_mm( nat_all ) )
    END IF
    IF( .NOT. ALLOCATED( tau_mm ) ) THEN
        ALLOCATE( tau_mm( 3, nat_mm ) ) 
    END IF
    IF( .NOT. ALLOCATED( tau_mask ) ) THEN
        ALLOCATE( tau_mask( nat_mm ) ) 
    END IF
    IF( .NOT. ALLOCATED( charge_mm ) ) THEN
        ALLOCATE( charge_mm( nat_mm ) ) 
    END IF
    IF( .NOT. ALLOCATED( aradii ) ) THEN
        ALLOCATE( aradii( nat_mm ) ) 
    END IF
    IF( .NOT. ALLOCATED( charge ) ) THEN
        ALLOCATE( charge(nat_qm) )
    END IF
    IF( .NOT. ALLOCATED( force_qm ) ) THEN
        ALLOCATE( force_qm(3,nat_qm) )
    END IF
    IF( .NOT. ALLOCATED( force_mm ) ) THEN
        ALLOCATE( force_mm(3,nat_mm) )
    END IF
    IF( .NOT. ALLOCATED( types ) ) THEN
        ALLOCATE( types( nat_all ) )
    END IF
    IF( .NOT. ALLOCATED( mass ) ) THEN
        ! add 1 to take into account the atom type "0"
        ALLOCATE( mass( ntypes + 1 ) ) 
    END IF

    ! Receive coordinates (from LAMMPS) and broadcast to all processors
    IF (ionode) THEN

        CALL mpi_recv( cell_mm, 9, MPI_DOUBLE_PRECISION, &
              0, QMMM_TAG_CELL, qmmm_comm, MPI_STATUS_IGNORE, ierr )

        CALL mpi_recv(tau(1,1),3*nat_qm,MPI_DOUBLE_PRECISION, &
              0,QMMM_TAG_COORD,qmmm_comm,MPI_STATUS_IGNORE,ierr)
 
        CALL mpi_recv(charge(1),nat_qm,MPI_DOUBLE_PRECISION, &
              0,QMMM_TAG_CHARGE,qmmm_comm,MPI_STATUS_IGNORE,ierr)

        CALL mpi_recv(charge_mm(1),nat_all,MPI_DOUBLE_PRECISION, &
              0,QMMM_TAG_COORD,qmmm_comm,MPI_STATUS_IGNORE,ierr)

        CALL mpi_recv(tau_mm(1,1),3*nat_all,MPI_DOUBLE_PRECISION, &
              0,QMMM_TAG_COORD,qmmm_comm,MPI_STATUS_IGNORE,ierr)

        CALL mpi_recv(tau_mask(1),nat_all,MPI_INTEGER, &
              0,QMMM_TAG_COORD,qmmm_comm,MPI_STATUS_IGNORE,ierr)

        CALL mpi_recv(types(1),nat_all,MPI_INTEGER, &
              0,QMMM_TAG_TYPE,qmmm_comm,MPI_STATUS_IGNORE,ierr)

        CALL mpi_recv(mass(1),ntypes+1,MPI_DOUBLE_PRECISION, &
              0,QMMM_TAG_MASS,qmmm_comm,MPI_STATUS_IGNORE,ierr)

        ! convert from angstrom to alat units
        tau = tau / (alat * bohr_radius_angs)

        tau_mm = tau_mm / (alat * bohr_radius_angs)

        CALL qmmm_center_molecule
        CALL qmmm_minimum_image

        ! set atomic radii
        CALL ec_fill_radii( aradii, nat_mm, mass, types, ntypes, 1 )

    END IF

    CALL mp_bcast(cell_mm, ionode_id, world_comm )
    CALL mp_bcast(aradii, ionode_id, world_comm)    
    CALL mp_bcast(tau, ionode_id, world_comm)
    CALL mp_bcast(charge, ionode_id, world_comm)    
    CALL mp_bcast(charge_mm, ionode_id, world_comm)    
    CALL mp_bcast(tau_mm, ionode_id, world_comm)    
    CALL mp_bcast(tau_mask, ionode_id, world_comm)    
    CALL mp_bcast(types, ionode_id, world_comm)    
    CALL mp_bcast(mass, ionode_id, world_comm)    

    ! clear charge for QM atoms
    DO i = 1, nat_mm
       IF(tau_mask(i) .eq. -1)CYCLE
       charge_mm(i) = 0.0d0
    ENDDO

    rc_mm = aradii
    ! Convert radii to Bohr units
    rc_mm = rc_mm / (alat * bohr_radius_angs)

    IF (ionode) THEN
       WRITE(stdout,*)
       WRITE(stdout,'(5X,A)') 'QMMM: cell_mm'
       WRITE(stdout,'(11X,A,3F6.3)') 'X (lo,hi,len): ',cell_mm(1),cell_mm(4),cell_mm(4)-cell_mm(1)
       WRITE(stdout,'(11X,A,3F6.3)') 'Y (lo,hi,len): ',cell_mm(2),cell_mm(5),cell_mm(5)-cell_mm(2)
       WRITE(stdout,'(11X,A,3F6.3)') 'Z (lo,hi,len): ',cell_mm(3),cell_mm(6),cell_mm(6)-cell_mm(3)
       WRITE(stdout,'(11X,A,3F6.3)') '  (xy,xz,yz) : ',cell_mm(7),cell_mm(8),cell_mm(9)
       WRITE(stdout,*)
       DO i = 1, nat_qm
           WRITE(stdout,'(5X,A,3F10.6,2X,A,F10.6)') &
                'QMMM: tau    ',tau(:,i), ' charge    ',charge(i)
       END DO
       WRITE(stdout,*)
       DO i = 1, nat_all
           WRITE(stdout,'(5X,A,3F10.6,2X,A,F10.6,2X,A,I2)') &
                'QMMM: tau_mm ',tau_mm(:,i),' charge_mm ',charge_mm(i),' QA ',tau_mask(i)
       END DO
    END IF


#else
    CALL errore( 'qmmm_update_positions', 'Use of QM/MM requires compilation with MPI', 1 )
#endif

  END SUBROUTINE qmmm_update_positions

  !---------------------------------------------------------------------!
  ! communicate forces of the QM system to MM-master
  !
  SUBROUTINE qmmm_update_forces( force, rho, nspin, dfftp )
    !
    USE fft_types,          ONLY : fft_type_descriptor
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: force(:,:)
    REAL(DP) :: rho(:,:)
    INTEGER  :: nspin
    TYPE(fft_type_descriptor) :: dfftp
    INTEGER :: ierr, i
    IF (qmmm_mode < 0) RETURN

#if defined(__MPI)

    IF( qmmm_mode == 2 ) THEN
       IF (ionode .and. (qmmm_verb > 0)) &
          WRITE(stdout,'(/,5X,A)') 'QMMM: compute EC forces'
       CALL qmmm_force_esf( rho, nspin, dfftp )
    END IF

    IF (ionode) THEN
        IF (qmmm_verb > 0) WRITE(stdout,'(5X,A)') 'QMMM: update forces'
        ! convert from atomic to real units
        IF( qmmm_mode == 2 ) THEN
           tmp_buf = (force + force_qm) * QMMM_FORCE_CONV
        ELSE
           tmp_buf = force * QMMM_FORCE_CONV
        END IF

        CALL mpi_send(tmp_buf,3*nat_qm,MPI_DOUBLE_PRECISION, 0,QMMM_TAG_FORCE,qmmm_comm,ierr)
        !
        !!!! Note, not used if ec_alg is false. Optimize excluding this send as well
        force_mm = force_mm * QMMM_FORCE_CONV
        CALL mpi_send(force_mm,3*nat_mm,MPI_DOUBLE_PRECISION, 0,QMMM_TAG_FORCE2,qmmm_comm,ierr)
    END IF
#else
    CALL errore( 'qmmm_update_forces', 'Use of QM/MM requires compilation with MPI', 1 )
#endif
  END SUBROUTINE qmmm_update_forces

  !---------------------------------------------------------------------!
  ! add electrostatic field of MM system to QM system

  SUBROUTINE qmmm_add_esf( vltot, dfftp )
    !--------------------------------------------------------------------------
    !
    !   This routine adds an electrostatic field due to MM atoms to the 
    !   local potential.
    !
    USE cell_base,          ONLY : alat, at, omega
    USE ions_base,          ONLY : zv, tau
    USE constants,          ONLY : e2, eps8, bohr_radius_angs
    USE io_global,          ONLY : stdout,ionode
    USE fft_types,          ONLY : fft_type_descriptor
    USE kinds,              ONLY : DP
    !
    USE constraints_module, ONLY : pbc
    !
    IMPLICIT NONE
    !
    REAL(DP) :: vltot(:)
    TYPE(fft_type_descriptor) :: dfftp
    !
    ! local variables
    !
    INTEGER :: index, index0, i, j, k
    INTEGER :: ir
    !
    INTEGER :: i_mm, i_qm, ipol, ii_qm
    ! r_nn is the cutoff for the nearest neighbour
    REAL(DP) :: s(3),r(3), dist, r_nn, fder
    !
    REAL(DP) :: esfcontrib
    REAL(DP),ALLOCATABLE :: esfcontrib_all(:)
    !
    ! if either the MS2 or EC aren't enabled, exit immediately
    IF( qmmm_mode /= 2 ) RETURN
    !
    ! Index for parallel summation
    !
    ALLOCATE(esfcontrib_all(dfftp%nnr))
    esfcontrib_all(:) = 0.D0
    !
    index0 = 0
    !
    DO i = 1, me_pool
       index0 = index0 + dfftp%nr1x*dfftp%nr2x*dfftp%npp(i)
    END DO
    !
    r_nn = 50000.d0 ! cut-off for the nearest neighbour
    !
    r(:) = 0.d0
    ! 
    DO ir = 1, dfftp%nnr
       index = index0 + ir - 1
       k     = index / (dfftp%nr1x*dfftp%nr2x)
       IF ( k .GE. dfftp%nr3 ) CYCLE
       index = index - (dfftp%nr1x*dfftp%nr2x)*k
       j     = index / dfftp%nr1x
       IF ( j .GE. dfftp%nr2 ) CYCLE
       index = index - dfftp%nr1x*j
       i     = index
       IF ( i .GE. dfftp%nr1 ) CYCLE
       !
       s(1) = DBLE(i)/DBLE(dfftp%nr1)
       s(2) = DBLE(j)/DBLE(dfftp%nr2)
       s(3) = DBLE(k)/DBLE(dfftp%nr3)
       !
       r=matmul(at,s)
       ! 
       ! Clear the contribute (it's an accumulator)
       esfcontrib = 0.0D0
       !
       DO i_mm = 1, nat_mm 

          if(tau_mask(i_mm) .ne. -1)cycle ! only MM atoms contribute to ESF

          dist=sqrt((tau_mm(1, i_mm)-r(1))**2 + (tau_mm(2, i_mm)-r(2))**2 + (tau_mm(3, i_mm)-r(3))**2)
          !
          if(dist .LE. r_nn) then
              esfcontrib = esfcontrib - e2*charge_mm(i_mm)*(rc_mm(i_mm)**4 -dist**4)/(rc_mm(i_mm)**5 -dist**5) / alat
          end if
       ENDDO
       !
       ! Add the contribute
       vltot(ir) = vltot(ir) + esfcontrib
       esfcontrib_all(ir) = esfcontrib
       !
    END DO
    !
    r(:) = 0.D0
    force_qm = 0.D0
    !
    !write(stdout, *) "Check QM position"
    !i_qm = 1
    !DO i_mm = 1, nat_mm
    !   IF(tau_mask(i_mm) .eq. -1)CYCLE
    !   write(stdout, '(I5, I5, 3f11.7, " - ", 3f11.7)') &
    !        i_mm, tau_mask(i_mm),tau_mm(1, i_mm),tau_mm(2, i_mm),tau_mm(3, i_mm),tau(1, i_qm),tau(2, i_qm),tau(3, i_qm)
    !   i_qm = i_qm + 1
    !ENDDO
    !write(stdout, *) "All position & charges"
    !DO i_mm = 1, nat_mm
    !   write(stdout, '(I5, I5, 5f11.7)') &
    !        i_mm, tau_mask(i_mm), tau_mm(1, i_mm), tau_mm(2, i_mm), tau_mm(3, i_mm), charge_mm(i_mm), rc_mm(i_mm)
    !ENDDO
    !write(stdout, *) "vltot=", SUM(vltot)
    !write(stdout, *) "esfcontrib_all=", SUM(esfcontrib_all)

    ii_qm = 1
    DO i_qm = 1, nat_mm
       if(tau_mask(i_qm) .eq. -1)cycle
       DO i_mm = 1, nat_mm
          IF(tau_mask(i_mm) .ne. -1)CYCLE
          dist = sqrt((tau_mm(1, i_mm) - tau_mm(1, i_qm))**2 +   &
                      (tau_mm(2, i_mm) - tau_mm(2, i_qm))**2 +   &
                      (tau_mm(3, i_mm) - tau_mm(3, i_qm))**2)

                fder = ( 5.d0*(dist**4)*( rc_mm(i_mm)**4 - dist**4 ) -   &
                         4.d0*(dist**3)*( rc_mm(i_mm)**5 - dist**5 ) ) / & 
                            ( ( rc_mm(i_mm)**5 - dist**5 )**2 )
                DO ipol = 1,3
                      force_qm(ipol,ii_qm) = force_qm(ipol,ii_qm) -  &
                      e2*charge_mm(i_mm)*zv(tau_mask(i_qm)) *       &
                      fder*(tau_mm(ipol, i_qm)-tau_mm(ipol, i_mm))/dist
                ENDDO
       ENDDO
       ii_qm = ii_qm + 1
    ENDDO

    force_qm=force_qm/(alat**2)
    
    !write(stdout, *) "NEW Forces added to QM atoms (Ry / a.u.)"
    !DO i_qm = 1, nat_mm
    !   write(stdout, '(I5, I5, f11.7, f11.7, f11.7, f11.7)') &
    !        i_qm, tau_mask(i_qm),zv(tau_mask(i_qm)),force_qm(IDX1D(1,i_qm): IDX1D(3,i_qm))
    !ENDDO
    !write(stdout, *) "End of NEW forces added to QM atoms (Ry / a.u.)"
    
    DEALLOCATE( esfcontrib_all )
    
    !IF (ionode) THEN
    !PRINT *,"****** END OF ADD_ESF COMPUTATION ******"
    !ENDIF

    RETURN

  END SUBROUTINE qmmm_add_esf

  !---------------------------------------------------------------------!

  SUBROUTINE qmmm_force_esf(rho,nspin,dfftp)
    !
    !   This routine computes the forces on the MM atoms due to the QM part
    !
    
    USE cell_base,          ONLY : alat, at, omega
    USE fft_types,          ONLY : fft_type_descriptor
    USE constants,          ONLY : e2, eps8
    USE io_global,          ONLY : stdout,ionode
    USE ions_base,          ONLY : zv, tau
    USE kinds,              ONLY : DP
    !
    IMPLICIT NONE
    !
    REAL(DP) :: rho(:,:)
    INTEGER  :: nspin
    TYPE(fft_type_descriptor) :: dfftp
    !
    ! local variables
    !
    INTEGER :: index, index0, i, j, k
    INTEGER :: ir
    !
    INTEGER :: i_mm, i_qm, ipol,is
    REAL(DP) :: s(3),r(3), dist, fder, r_nn

    IF( qmmm_mode /= 2 ) RETURN
    !
    ! Index for parallel summation
    !
    index0 = 0
    !
    DO i = 1, me_pool
       index0 = index0 + dfftp%nr1x*dfftp%nr2x*dfftp%npp(i)
    END DO
    !
    r(:) = 0.d0
    r_nn = 5000000.d0 ! cut-off for nearest neighbor
    !
    force_mm = 0.d0
    !
    ! Compute forces on MM atoms due to valence electrons 
    !
    DO i_mm = 1, nat_mm
       !
       if(tau_mask(i_mm) .ne. -1)cycle
       !
       DO is=1,nspin
          DO ir = 1, dfftp%nnr
             !
             ! ... three dimensional indexes
             !
             index = index0 + ir - 1
             k     = index / (dfftp%nr1x*dfftp%nr2x)
             IF ( k .GE. dfftp%nr3 ) CYCLE
             index = index - (dfftp%nr1x*dfftp%nr2x)*k
             j     = index / dfftp%nr1x
             IF ( j .GE. dfftp%nr2 ) CYCLE
             index = index - dfftp%nr1x*j
             i     = index
             IF ( i .GE. dfftp%nr1 ) CYCLE
             !
             s(1) = DBLE(i)/DBLE(dfftp%nr1)
             s(2) = DBLE(j)/DBLE(dfftp%nr2)
             s(3) = DBLE(k)/DBLE(dfftp%nr3)

             !
             r=matmul(at,s)
             dist = sqrt((tau_mm(1, i_mm)-r(1))**2 + (tau_mm(2, i_mm)-r(2))**2 + (tau_mm(3, i_mm)-r(3))**2)
             !
             !  see equation (3) from  j.chem. phys 116 by Laio
             !
             fder = ( 5.d0*(dist**4)*( rc_mm(i_mm)**4 - dist**4 ) -   &
                         4.d0*(dist**3)*( rc_mm(i_mm)**5 - dist**5 ) ) / & 
                            ( ( rc_mm(i_mm)**5 - dist**5 )**2 )
             !
             DO ipol = 1,3
                force_mm(ipol,i_mm) = force_mm(ipol,i_mm) +  &
                        rho(ir,is)*fder*(tau_mm(ipol, i_mm)-r(ipol))/dist
             ENDDO
             !
          END DO
       END DO
       !
       force_mm(1,i_mm) = force_mm(1,i_mm) * charge_mm(i_mm)
       force_mm(2,i_mm) = force_mm(2,i_mm) * charge_mm(i_mm)
       force_mm(3,i_mm) = force_mm(3,i_mm) * charge_mm(i_mm)
    END DO
    ! 
    CALL mp_sum(force_mm, intra_pool_comm)
    !
    force_mm(:,:) = e2*force_mm(:,:)*omega/(dfftp%nr1*dfftp%nr2*dfftp%nr3)

    !write(stdout, *) "RHO = ", SUM(rho)
    !write(stdout, *) "Forces added to MM atoms (Ry / a.u.)"
    !DO i_mm = 1, nat_mm
    !   write(stdout, '(I5, f11.7, f11.7, f11.7, f11.7, f11.7, f11.7)') i_mm, force_mm(1:3,i_mm), tau_mm(1:3, i_mm)
    !ENDDO
    ! 
    DO i_mm = 1, nat_mm
       if(tau_mask(i_mm) .ne. -1)cycle
       DO i_qm = 1, nat_mm
          if(tau_mask(i_qm) .eq. -1)cycle
          dist = sqrt((tau_mm(1, i_mm) - tau_mm(1, i_qm))**2 +   &
                      (tau_mm(2, i_mm) - tau_mm(2, i_qm))**2 +   &
                      (tau_mm(3, i_mm) - tau_mm(3, i_qm))**2)
          fder = ( 5.d0*(dist**4)*( rc_mm(i_mm)**4 - dist**4 ) -   &
                         4.d0*(dist**3)*( rc_mm(i_mm)**5 - dist**5 ) ) / & 
                            ( ( rc_mm(i_mm)**5 - dist**5 )**2 )
          DO ipol = 1,3
             force_mm(ipol,i_mm) = force_mm(ipol,i_mm) -  &
                 e2*charge_mm(i_mm)*zv(tau_mask(i_qm)) *       &
                 fder*(tau_mm(ipol, i_mm)-tau_mm(ipol, i_qm))/dist
          ENDDO
       ENDDO
    ENDDO
    force_mm(:,:)=force_mm(:,:)/(alat**2)
    !
    !write(stdout, *) "Forces added to MM atoms (Ry / a.u.)"
    !DO i_mm = 1, nat_mm
    !   write(stdout, '(I5, f11.7, f11.7, f11.7)') i_mm, force_mm(1:3,i_mm)
    !ENDDO
    !write(stdout, *) "End of forces added to MM atoms (Ry / a.u.)"
    ! convert tau_mm back to angstrom
    ! tau_mm(:) = tau_mm(:)*au2ang
    !
    !
    !IF (ionode) THEN
    !PRINT *, "****** END OF FORCE_ESF COMPUTATION ******"
    !ENDIF
    !
    RETURN
    
  END SUBROUTINE qmmm_force_esf

  !---------------------------------------------------------------------!

  SUBROUTINE qmmm_shutdown
    !
    ! cleanup of QM/MM. free resources
    !
    IMPLICIT NONE
    !
    IF (qmmm_mode < 0) RETURN
    !
    IF (ionode) THEN
        WRITE(stdout,'(/,5X,A)') "QMMM: Shutting down QM/MM coupling"
    END IF
    !
    IF( ALLOCATED( tmp_buf ) ) DEALLOCATE( tmp_buf )
    IF( ALLOCATED( rc_mm ) ) DEALLOCATE( rc_mm )
    IF( ALLOCATED( aradii ) ) DEALLOCATE( aradii )
    IF( ALLOCATED( tau_mm ) ) DEALLOCATE( tau_mm ) 
    IF( ALLOCATED( tau_mask ) ) DEALLOCATE( tau_mask ) 
    IF( ALLOCATED( charge_mm ) ) DEALLOCATE( charge_mm ) 
    IF( ALLOCATED( charge ) ) DEALLOCATE( charge )
    IF( ALLOCATED( force_qm ) ) DEALLOCATE( force_qm )
    IF( ALLOCATED( force_mm ) ) DEALLOCATE( force_mm )
    IF( ALLOCATED( types ) ) DEALLOCATE( types )
    IF( ALLOCATED( mass ) ) DEALLOCATE( mass )

  END SUBROUTINE qmmm_shutdown


END MODULE qmmm
