!
! Copyright (C) 2001-2005 QUANTUM-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------
MODULE vibrations
  !=----------------------------------------------------------------
  ! An implementation of the classical frozen-phonon approach for 
  ! calculation of the Gamma-point dynamical matrix, eigen frequencies
  ! and vectors, and infrared intensities
  !
  ! Programmed by: Silviu Zilberman
  !
  USE kinds,                  ONLY: DP

  IMPLICIT NONE
  SAVE

  PRIVATE
  !
  ! ... private variables
  !
  ! ref_c0      - ground state wave function
  ! ref_tau     - ground state configuration
  ! ref_lambda  - lagrange multipliers at ground state
  ! ref_etot    - total energy of ground state
  ! fion        - forces on the ions
  ! active_atom - frozen/non frozen atom
  !               frozen atoms are not displaced
  !
  COMPLEX (KIND=DP), ALLOCATABLE :: ref_c0(:,:,:,:)

  REAL    (KIND=DP), ALLOCATABLE :: ref_tau(:,:)
  REAL    (KIND=DP), ALLOCATABLE :: ref_lambda(:,:)
  REAL    (KIND=DP)              :: ref_etot
  REAL    (KIND=DP), ALLOCATABLE :: fion(:,:)   

  LOGICAL,           ALLOCATABLE :: active_atom(:)
  !
  ! ... public variables
  !
  PUBLIC :: delta, eigenvals, eigenvecs, born_charge, U, T,              &
       trans_inv_conv_thr, save_freq, nactive_atoms, trans_inv_max_iter, &
       vib_restart_mode, trans_inv_flag, trans_rot_inv_flag, animate
  !
  ! delta              - displacement step
  ! U                  - energy hessian
  ! T                  - diagonal mass matrix
  ! eigenvals          - eigen values
  ! eigenvecs          - eigen vectros
  ! born_charge        - born charge tensor for each atom
  ! trans_inv_thr      - criteria for convergence of the
  !                      repeated application of the symmentrization
  !                      and translational invariance procedure
  ! save_freq          - save restart info every save_freq displacements
  ! nactive_atoms      - # of non-frozen atoms
  ! trans_inv_max_iter - maximal # of iterations to converge
  !                      the translational invariance
  ! vib_restart_mode   - from scrtach = 0,
  !                      restart      = 1,
  !                      auto         = 2
  ! restart_filename   - obvious meaning
  ! trans_inv_flag     - turn on/off translational invariance
  ! trans_rot_inv_flag - turn on/off removal of rigid modes
  ! animate            - creates an xyz animation file for each mode
  !
  REAL(KIND=DP)                  :: delta 
  REAL(KIND=DP), ALLOCATABLE     :: U(:,:)
  REAL(KIND=DP), ALLOCATABLE     :: T(:,:)
  REAL(KIND=DP), ALLOCATABLE     :: eigenvals(:)
  REAL(KIND=DP), ALLOCATABLE     :: eigenvecs(:,:)
  REAL(KIND=DP), ALLOCATABLE     :: born_charge(:,:)
  REAL(KIND=DP)                  :: trans_inv_conv_thr

  INTEGER                        :: save_freq
  INTEGER                        :: nactive_atoms
  INTEGER                        :: trans_inv_max_iter
  INTEGER                        :: vib_restart_mode

  CHARACTER (len=256)            :: restart_filename

  LOGICAL                        :: trans_inv_flag, trans_rot_inv_flag
  LOGICAL                        :: animate
  !
  ! ...Public methods
  !
  PUBLIC :: start_vibrations
  PUBLIC :: calc_hessian
  PUBLIC :: analysis
  PUBLIC :: end_vibrations

  ! ----------------------------------------------------------------
CONTAINS
  ! ----------------------------------------------------------------


  SUBROUTINE start_vibrations (restart_cyc_counter, E_minus,  &
       dip_minus)
    !
    !-------------------------------------------------------------
    !
    USE constants,            ONLY : DIP_DEBYE
    USE cell_base,            ONLY : tpiba2, h
    USE cg_module,            ONLY : tcg
    USE constants,            ONLY : AMU_AU
    USE cp_electronic_mass,   ONLY : emass_precond, emass_cutoff
    USE cp_main_variables,    ONLY : lambda, lambdam, ema0bg, nfi, bec
    USE cp_main_variables,    ONLY : irb, eigrb, rhor, rhog, rhos
    USE cp_main_variables,    ONLY : lambdap, eigr
    USE electrons_base,       ONLY : nbsp, nbspx, nel
    USE electrons_module,     ONLY : cp_eigs
    USE energies,             ONLY : etot, ekin
    USE gvecw,                ONLY : ngw, ggp
    USE io_files,             ONLY : outdir, prefix
    USE io_global,            ONLY : ionode, ionode_id, stdout
    USE ions_base,            ONLY : nsp, nat, iforce, na, pmass
    USE ions_positions,       ONLY : tau0
    USE kinds,                ONLY : DP
    USE mp,                   ONLY : mp_bcast
    USE parameters,           ONLY : natx
    USE printout_base,        ONLY : printout_pos
    USE print_out_module,     ONLY : cp_print_rho
    USE restart_file,         ONLY : writefile
    USE wavefunctions_module, ONLY : c0, cm
    !
    ! ... output variables
    !
    INTEGER,        INTENT(OUT)   :: restart_cyc_counter
    REAL (KIND=DP), INTENT(OUT)   :: dip_minus(3), E_minus
    !
    ! ... local variables
    !
    INTEGER                       :: counter,is,ia,coord,ierr
    INTEGER                       :: dirlen, filep=200, printwfc
    LOGICAL                       :: restart_vib, mass_file_exists
    REAL (KIND=DP)                :: tmp_mass, dipole(3), dipole_moment
    CHARACTER (len=256)           :: mass_file
    !
    !-------------------------------------------------------------
    !
    ! (1) Allocate arrays
    !
    ALLOCATE( ref_c0     ( ngw,   nbspx, 1, 1 ) )
    ALLOCATE( ref_tau    ( 3,     natx        ) )
    ALLOCATE( ref_lambda ( nbsp,  nbsp        ) )
    ALLOCATE( active_atom( nat                ) )
    ALLOCATE( eigenvals  ( 3*nat              ) )
    ALLOCATE( eigenvecs  ( 3*nat, 3*nat       ) )
    ALLOCATE( U          ( 3*nat, 3*nat       ) )
    ALLOCATE( T          ( 3*nat, 3*nat       ) )
    ALLOCATE( born_charge( 3,     3*nat       ) )
    ALLOCATE( fion       ( 3,     natx        ) )
    !
    ! (2) Calculate how many active atoms (non-frozen) are present.
    !     Any atom with ANY frozen coordinate is assume to be totally
    !     frozen.
    !
    nactive_atoms = 0
    DO ia = 1,nat
       IF (iforce(1,ia)+iforce(2,ia)+iforce(3,ia).EQ.3) THEN
          nactive_atoms=nactive_atoms+1
          active_atom(ia)=.TRUE.
       ELSE
          active_atom(ia)=.FALSE.
       END IF
    END DO
    !
    ! (3) initiating variables ...
    !
    U = 0.0
    restart_vib = .FALSE.
    restart_cyc_counter = -1
    !
    ! (4) Setting the T matrix (diagonal masses matrix)
    !
    ! ... set file name
    !
    IF (ionode) THEN
       dirlen    = INDEX(outdir,' ') - 1
       mass_file = TRIM(prefix)//'.vib.isotope'
       mass_file = outdir(1:dirlen) // '/' // mass_file
       !
       ! ... check existance
       !
       INQUIRE (FILE = mass_file, EXIST = mass_file_exists)
       !
       ! ... use it if it's there or creat a new one using masses
       ! ... that were specified in the CP input file
       !
       IF (mass_file_exists) THEN
          !
          ! ... read isotope masses from input.
          ! ... The input file should contain nat lines with the mass 
          ! ... of each individual atom in each line, in AMU units
          !
          WRITE (stdout,*) 'Using existing isotopes file...', mass_file
          OPEN(filep , FILE = mass_file , STATUS = 'old',IOSTAT=ierr)
          IF( ierr /= 0 ) &
               CALL errore(' start_vibrations ', ' opening file '//mass_file, 1 )
          !
          T=0.0
          counter=1
          DO is=1,nat
             READ (filep,*)           tmp_mass
             T(counter,counter)     = tmp_mass * AMU_AU
             T(counter+1,counter+1) = tmp_mass * AMU_AU
             T(counter+2,counter+2) = tmp_mass * AMU_AU
             counter=counter+3
          END DO
          CLOSE (filep)
       ELSE
          !
          ! ... no isotopes file found
          !
          WRITE (stdout,*) 'Creating new isotopes file: ', mass_file
          T=0.
          counter=0
          DO is=1,nsp              
             DO ia=1,na(is)      
                DO coord=1,3      
                   counter=counter+1
                   T(counter,counter)=pmass(is)
                END DO
             END DO
          END DO
          !
          ! write a prefix.vib.isotope file
          !
          !
          OPEN(filep,file=mass_file,status='new',IOSTAT=ierr)
          IF( ierr /= 0 ) &
               CALL errore(' start_vibrations ', ' creating file '//mass_file, 1 )
          !
          DO is=1,nsp            
             DO ia=1,na(is)    
                WRITE (filep,*) pmass(is) / AMU_AU
             END DO
          END DO
          CLOSE(filep)
       END IF
    END IF
    !
    ! (5) restarting from file ? ...
    !
    dirlen           = INDEX(outdir,' ') - 1
    restart_filename = TRIM(prefix)//'.vib.restart'
    restart_filename = outdir(1:dirlen) // '/' // restart_filename
    !
    SELECT CASE ( vib_restart_mode )
    CASE ( 0 )
       restart_vib = .FALSE.
    CASE ( 1 )
       INQUIRE (file=restart_filename, EXIST=restart_vib)
       IF ( .NOT. restart_vib ) &
            CALL errore( 'start_vibrations', &
            'file is absent: '//restart_filename, 1 )
    CASE ( 2 )
       INQUIRE (file=restart_filename, EXIST=restart_vib)
    END SELECT
    !
    IF ( restart_vib ) THEN
       !
       WRITE ( stdout , * ) '------------------------------------'
       WRITE ( stdout , * ) 'CONTINUE A NORMAL MODES CALCULATION '
       WRITE ( stdout , * ) '------------------------------------'
       WRITE ( stdout , * )
       WRITE ( stdout , * ) 'Restart information is read from (and saved to) ', restart_filename
       WRITE ( stdout , * )
       !
       CALL read_restart (E_minus, dip_minus, restart_cyc_counter)
       !
       WRITE ( stdout , * )                                       &
            '... last saved iteration index: ', restart_cyc_counter
       !
    ELSE
       ! Starting vibrational calculation from scratch
       WRITE ( stdout , * ) '---------------------------------------'
       WRITE ( stdout , * ) 'STARTING A NEW NORMAL MODES CALCULATION '
       WRITE ( stdout , * ) '---------------------------------------'
       WRITE ( stdout , * )
       WRITE ( stdout , * ) 'Restart information is saved to ', restart_filename
       WRITE ( stdout , * )
    END IF
    !
    IF (restart_cyc_counter .LT. 3*nat) THEN
       !
       ! ... relaxing wave function and saving to a file
       !
       WRITE ( stdout , * ) '... Initial relaxation of wavefunction...'
       !
       CALL relax_wavefunction (fion)
       !CALL cp_eigs( nfi, bec, c0, irb, eigrb, rhor, &
       !    rhog, rhos, lambdap, lambda, tau0, h )
       !
       !printwfc = 1
       !IF ( printwfc >= 0 ) &
       !     CALL cp_print_rho( nfi, bec, c0, eigr, irb, eigrb, rhor, &
       !     rhog, rhos, lambdap, lambda, tau0, h )
       !
       !
       WRITE ( stdout , * ) 'Done' 
       WRITE ( stdout , * )                                       
       !
       ! (6) Saving refernce ground-state parameters
       !     (coordinates, wavefunctions, energy and polarization)
       !
       ref_c0     = c0
       ref_lambda = lambda
       ref_etot   = etot
       cm         = c0           ! setting zero wavefunction velocity
       lambdam    = lambda
       ref_tau    = tau0
       !
       !     ... printing refernce structure and energy
       !
       CALL calculate_dipole (dipole, dipole_moment,tau0)
       WRITE (stdout,*)
       CALL printout_pos( stdout, ref_tau, nat, 'pos' )
       WRITE (stdout,*)
       WRITE (stdout,110) ref_etot
       WRITE (stdout,111) dipole(1) * DIP_DEBYE, dipole(2) * DIP_DEBYE, &
            dipole(3) * DIP_DEBYE
       WRITE (stdout,112) dipole_moment * DIP_DEBYE
       WRITE (stdout,*)
       !
110    FORMAT(3x,'Ground-state (reference) energy   :',3x,f10.6)
111    FORMAT(3x,'Ground-state dipole vector [debye]:',3x,f10.3,3x,f10.3,3x,f10.3)
112    FORMAT(3x,'Dipole moment [debye]             :',3x,f10.3)
       !
       ! (7) Electronic mass preconditioning, in case damped scf minimization
       !     is used. It assumes a converged ground-state configuration, and 
       !     follows ideas by Pyne and co-workers.
       !     Relevant only of damped dynamics is used for scf
       !
       IF (.NOT. tcg) THEN
          emass_cutoff = ekin/(nel(1)+nel(2))
          CALL emass_precond( ema0bg, ggp, ngw, tpiba2, emass_cutoff )
       END IF
       !
    ELSE
       ref_tau = tau0
    END IF
    !
    RETURN
  END SUBROUTINE start_vibrations
  !
  !
  !  ----------------------------------------------
  !
  !
  SUBROUTINE calc_hessian (restart_cyc_counter, E_minus, dip_minus)
    !
    !-------------------------------------------------------------
    !
    USE constants,            ONLY : DIP_DEBYE, eps4
    USE cp_main_variables,    ONLY : lambda, lambdam, nfi
    USE efield_module,        ONLY : efield_update
    USE energies,             ONLY : etot
    USE from_scratch_module,  ONLY : from_scratch
    USE input_parameters,     ONLY : atom_label
    USE io_global,            ONLY : stdout, ionode
    USE ions_base,            ONLY : nat, ityp
    USE ions_positions,       ONLY : tau0
    USE parameters,           ONLY : natx
    USE phase_factors_module, ONLY : strucf     
    USE wavefunctions_module, ONLY : c0, cm
    !
    ! ... input variables
    !
    INTEGER,        INTENT(IN)    :: restart_cyc_counter
    REAL (KIND=DP), INTENT(INOUT) :: E_minus, dip_minus(3)
    !
    ! ... local variables
    !
    INTEGER                       :: iax, iax2, coord
    INTEGER                       :: index, disp_sign, cyc_counter

    REAL (KIND=DP)                :: E_plus
    REAL (KIND=DP)                :: dipole(3)
    REAL (KIND=DP)                :: tmp1, tmp2, dipole_moment
    CHARACTER                     :: coordinate(3)
    CHARACTER                     :: tmp_string(2)
    !
    cyc_counter = 0
    coordinate(1)='X'
    coordinate(2)='Y'
    coordinate(3)='Z'
    !
    !-----------------------------------------------------------------
    ! CALCULATION OF THE DYNAMICAL MATRIX
    !-----------------------------------------------------------------
    DO iax=1,nat
       IF (active_atom(iax)) THEN
          DO coord=1,3 !x,y,z
             index = (iax-1)*3 + coord
             DO disp_sign=-1,1,2 !displacement sign: +-1
                !
                ! ... Check if displacement was already calculated in previous run
                !
                IF (cyc_counter.LT.restart_cyc_counter) THEN
                   cyc_counter = cyc_counter + 1
                   !
                   ! ... printing information
                   !
                   IF (ionode) THEN
                      IF (disp_sign.EQ.-1) THEN
                         tmp_string(1)='-'
                      ELSE
                         tmp_string(1)='+'
                      END IF
                      tmp_string(2) = coordinate(coord)
                      WRITE (stdout,*)
                      WRITE (stdout,*) 'Skipping previously calculated displacement:'
                      WRITE (stdout,*) 'atom #: ',iax,' Symbol: ',TRIM(atom_label(ityp(iax)))
                      WRITE (stdout,*) 'Displacement in direction: ',tmp_string
                      WRITE (stdout,*)
                   END IF
                ELSE
                   !
                   ! ... setting coordinates of displaced atom and resetting all relevant variables
                   !
                   nfi             = 0
                   tau0            = ref_tau
                   tau0(coord,iax) = tau0(coord,iax) + disp_sign*delta
                   !
                   !
                   IF(disp_sign.EQ.-1) THEN
                      !
                      ! ... use reference wavefunction as initial guess
                      !
                      c0(:,:,1,1) = ref_c0(:,:,1,1)
                      !
                   ELSE
                      !
                      ! ... linearly extrapolate wave function coefficients from their value
                      ! ... at dispalcement at the negative direction
                      !
                      c0(:,:,1,1) = 2*ref_c0(:,:,1,1)-c0(:,:,1,1) 
                      !
                   END IF
                   !
                   ! ... set wavefunction velocity to zero
                   !
                   cm      = c0
                   lambdam = lambda
                   !
                   ! ... printing information
                   !
                   IF (ionode) THEN
                      IF (disp_sign.EQ.-1) THEN
                         tmp_string(1)='-'
                      ELSE
                         tmp_string(1)='+'
                      END IF
                      tmp_string(2) = coordinate(coord)
                      WRITE (stdout,*)
                      WRITE (stdout,*) 'Moving now atom #',iax, &
                           '  Symbol:  ',TRIM(atom_label(ityp(iax)))
                      WRITE (stdout,*) 'Displacement in direction: ',tmp_string
                      WRITE (stdout,*)
                   END IF
                   !
                   !... relax wavefunction in new position
                   !
                   CALL relax_wavefunction (fion)
                   !
                   ! ... calculating the electronic dipole vector
                   !
                   CALL calculate_dipole (dipole, dipole_moment,tau0)
                   !
                   !
                   IF (ionode) THEN
                      !
                      ! ... printing information and various tests on intermediat results
                      !
                      WRITE (stdout,*)
                      WRITE (stdout,*) 'Done scf relaxation for displacement.'
                      WRITE (stdout,113) etot
                      WRITE (stdout,111) dipole        * DIP_DEBYE
                      WRITE (stdout,112) dipole_moment * DIP_DEBYE
                      WRITE (stdout,*)
                      !
                      ! ... verifying that the perturbed structure has higher energy then the
                      !     non perturbed
                      !
                      IF(etot.LT.ref_etot) THEN
                         CALL infomsg('calc_hessian', &
                              'Warning: Reference structure is not converged!!!',-1)
                      END IF
                      !
                      ! ... record observables:
                      !     1. dynamical matrix row elements
                      !     2. Born effective charge (dipole derivative)
                      !
                      IF(disp_sign.EQ.-1) THEN
                         !
                         ! ... disp_sign == -1
                         !
                         E_minus = etot
                         DO iax2 = 1,nat
                            IF (active_atom(iax2)) THEN
                               U((iax2-1)*3+1:(iax2-1)*3+3,index) = fion(1:3,iax2)
                            END IF
                         END DO
                         dip_minus = dipole
                      ELSE
                         !
                         ! ... disp_sign == 1
                         !
                         E_plus = etot
                         DO iax2=1,nat
                            IF (active_atom(iax2)) THEN
                               U((iax2-1)*3+1:(iax2-1)*3+3,index) = U((iax2-1)*3+1:(iax2-1)*3+3,index) &
                                    - fion(1:3,iax2)
                            END IF
                         END DO
                         ! adding negative sign, since electronic charge is positive
                         ! in this CP code.
                         born_charge(:,index)=-(dipole-dip_minus)/(2*delta)
                         !
                         ! A verification on the diagonal element of dynamical matrix:
                         ! ... comparing numerical second derivative of the total energy
                         ! ... to first derivative of the forces
                         !
                         tmp1 = (E_plus+E_minus-2*ref_etot)/(delta*delta)
                         tmp2 = U(index,index)/(2*delta)
                         IF( (ABS(tmp1-tmp2)/(tmp2) > 0.1) .AND. (tmp2 > eps4 ) ) THEN
                            CALL infomsg('calc_hessian','Warning: consistency check',-1)
                            WRITE(stdout,*) '   Numerical second derivative of the total energy, compared to'
                            WRITE(stdout,*) '   first derivative of the forces, for diagonal hessian element,'
                            WRITE(stdout,*) '   deviate by more then 10%:'
                            WRITE(stdout,'(3x,A,f10.4)') '     Energy second derivative : ', tmp1
                            WRITE(stdout,'(3x,A,f10.4)') '     Force  first  derivative : ', tmp2
                         END IF
                      END IF
                   END IF
                   !
                   cyc_counter=cyc_counter+1
                   !
                   ! ... saving restart file
                   !
                   IF(MOD(cyc_counter,save_freq).EQ.0) THEN
                      !
                      ! save restart info
                      !
                      IF (ionode) THEN
                         WRITE (stdout,*)
                         WRITE (stdout,*) 'Saving restart information...'
                         CALL write_restart (E_minus, dip_minus, cyc_counter)
                         WRITE (stdout,*) 'Done.'
                         WRITE (stdout,*)
                      END IF
                      !
                   END IF
                END IF
             END DO
          END DO
       END IF
    END DO
    !
    U=U/(2*delta)
    !
111 FORMAT(3x,'Ground-state dipole vector [debye]:',3x,f10.3,3x,f10.3,3x,f10.3)
112 FORMAT(3x,'Dipole moment [debye]             :',3x,f10.3)
113 FORMAT(3x,'displacement total energy         :',3x,f10.6)
    RETURN
  END SUBROUTINE calc_hessian
  !
  !
  !  ----------------------------------------------
  !
  !
  SUBROUTINE analysis
    !
    USE io_files,           ONLY : outdir, prefix
    USE io_global,          ONLY : ionode, stdout
    USE ions_base,          ONLY : nat
    !
    ! ... local variables
    !
    CHARACTER (len=256)         :: results_file
    INTEGER                     :: filep=200, dirlen, ierr, nmodes, i
    LOGICAL                     :: do_IR_intensity
    REAL (KIND=DP)              :: U_tmp(3*nat,3*nat), T_tmp(3*nat,3*nat)
    REAL (KIND=DP)              :: tot_born_charge(3,3), D(3*nat,3*nat)
    REAL (KIND=DP)              :: intensity(3*nat), mode_mass(3*nat)
    REAL (KIND=DP)              :: mode_force_constant(3*nat)
    REAL (KIND=DP), ALLOCATABLE :: U_internal(:,:)
    REAL (KIND=DP), ALLOCATABLE :: T_internal(:,:)
    REAL (KIND=DP), ALLOCATABLE :: eigenvecs_internal(:,:)
    REAL (KIND=DP), ALLOCATABLE :: eigenvals_internal(:)
    !
    !
    IF (ionode) THEN
       !
       ! ... Setting and opening results file
       !
       dirlen    = INDEX(outdir,' ') - 1
       results_file = TRIM(prefix)//'.vib.analysis'
       results_file = outdir(1:dirlen) // '/' // results_file
       WRITE (stdout,*) 'Results are written to file: ', results_file
       !
       OPEN(filep,file=results_file,status='unknown',IOSTAT=ierr)
       IF( ierr /= 0 ) &
            CALL errore(' analyze_vibrations ', ' opening file '//results_file, 1 )
       !
       ! ... imposing symmetry on the hessian
       !
       CALL symmetrize_matrix(U,3*nat)
       !
       ! ... imposing symmetry on the born charge tensors
       !
       DO i=1,nat
          CALL symmetrize_matrix(born_charge(:,3*(i-1)+1:3*(i-1)+3),3)
       END DO
       !
       U_tmp=U
       T_tmp=T
       !
       ! ------------------------------
       ! *** Analyzing the raw data ***
       ! ------------------------------
       !
       WRITE (filep,*)
       WRITE (filep,*) '**************************'
       WRITE (filep,*) 'Calculations with raw data'
       WRITE (filep,*) '**************************'
       WRITE (filep,*)
       !
       do_IR_intensity = .TRUE.
       CALL analyze_vibrations(U_tmp,T_tmp,3*nat,filep,tot_born_charge,eigenvals, &
            eigenvecs,do_IR_intensity,intensity)
       !
       ! ... computing mode effective mass and spring constant
       !
       CALL get_force_constant &
            (eigenvecs, eigenvals, 3*nat, mode_mass, mode_force_constant)
       !
       ! ... printing harmonic frequencies and eigenmodes
       !
       CALL print_eigenmodes(3*nat,.TRUE.,filep,eigenvals,eigenvecs, &
            mode_mass, mode_force_constant, intensity)
       !
       !------------------------------------------------------------------------------------------
       ! *** REPEATING CALCULATION WITH IMPOSED ACOUSTIC SUM RULE and translational invariance ***
       !------------------------------------------------------------------------------------------
       !
       IF (trans_inv_flag) THEN
          !
          WRITE (filep,*) '************************************************************************'
          WRITE (filep,*) 'Calculations with translational invariance and acoustic sum rule imposed'
          WRITE (filep,*) '************************************************************************'
          !
          U_tmp=U
          T_tmp=T
          !
          ! ... Imposing translational invariance on Hessian
          !
          CALL trans_invariance(U_tmp)
          !
          ! ... Imposing acoustic sum rule on the born effective charges
          !
          CALL apply_asr
          !
          ! ... analyzing
          !
          do_IR_intensity = .TRUE.
          CALL analyze_vibrations(U_tmp,T_tmp,3*nat,filep,tot_born_charge,eigenvals, &
               eigenvecs,do_IR_intensity,intensity)
          !
          ! ... computing mode effective mass and spring constant
          !
          CALL get_force_constant &
               (eigenvecs, eigenvals, 3*nat, mode_mass, mode_force_constant)
          !
          ! ... printing harmonic frequencies and eigenmodes
          !
          CALL print_eigenmodes(3*nat,.TRUE.,filep,eigenvals,eigenvecs, &
               mode_mass, mode_force_constant,intensity)
          !
       END IF
       !
       !------------------------------------------------------------------------------------------
       ! *** REPEATING CALCULATION WITH IMPOSED ACOUSTIC SUM RULE 
       !     translational and rotations are projected out
       !------------------------------------------------------------------------------------------
       !
       IF (trans_rot_inv_flag) THEN
          WRITE (filep,*) '**********************************************************'
          WRITE (filep,*) 'Calculations with rotations and translations projected out'
          WRITE (filep,*) '**********************************************************'
          !
          U_tmp=U
          T_tmp=T
          !
          ! ... Imposing translational invariance on Hessian
          !
          CALL trans_invariance(U_tmp)
          !
          ! ... Imposing acoustic sum rule on the born effective charges
          !
          CALL apply_asr
          !
          ! Internal coordinates representation, rigid translations and rotations are projected out
          !
          CALL internal_hessian(U_tmp,T_tmp,nmodes,ref_tau,D,filep)
          !
          ALLOCATE(U_internal(nmodes,nmodes))
          ALLOCATE(T_internal(nmodes,nmodes))
          ALLOCATE(eigenvecs_internal(nmodes,nmodes))
          ALLOCATE(eigenvals_internal(nmodes))
          !
          U_internal(1:nmodes,1:nmodes)=U_tmp(3*nat-nmodes+1:3*nat,3*nat-nmodes+1:3*nat)
          T_internal(1:nmodes,1:nmodes)=T_tmp(3*nat-nmodes+1:3*nat,3*nat-nmodes+1:3*nat)
          !
          do_IR_intensity = .FALSE.
          CALL analyze_vibrations(U_internal,T_internal,nmodes,filep,tot_born_charge, &
               eigenvals_internal,eigenvecs_internal,do_IR_intensity,intensity)
          !
          ! Calculating the mode effective charge:
          ! 1. First converting back eigenvectors from internal to mass-weighted cartesians
          eigenvecs( : , 1              : 3*nat-nmodes) =        D(: , 1              : 3*nat-nmodes)
          eigenvecs( : , 3*nat-nmodes+1 : 3*nat)        = MATMUL(D(: , 3*nat-nmodes+1 : 3*nat       ),eigenvecs_internal)
          !
          eigenvals(1              : 3*nat-nmodes) = 0.d0
          eigenvals(3*nat-nmodes+1 : 3*nat       ) = eigenvals_internal
          !
          ! 2. converting back the eigen vectors from mass-weighted to cartesian coordinates
          CALL mass_weighted_to_cartesian(3*nat,T,eigenvecs)
          CALL analyze_IR_intensities(tot_born_charge,3*nat,filep,eigenvecs)
          !
          ! ... computing mode effective mass and spring constant
          !
          CALL get_force_constant &
               (eigenvecs, eigenvals, 3*nat, mode_mass, mode_force_constant)

          !
          ! ... printing harmonic frequencies and eigenmodes
          !
          CALL print_eigenmodes(nmodes,.TRUE.,filep,eigenvals_internal,eigenvecs_internal, & 
               mode_mass(3*nat-nmodes+1:3*nat), mode_force_constant(3*nat-nmodes+1:3*nat), &
               intensity(3*nat-nmodes+1:3*nat))
          !
          DEALLOCATE(U_internal          )
          DEALLOCATE(T_internal          )
          DEALLOCATE(eigenvecs_internal  )
          DEALLOCATE(eigenvals_internal  )
          !
       END IF
       CLOSE(filep)
       !
       ! ... creating xyz animation files
       !
       IF (animate) &
            CALL animation (ref_tau,eigenvecs,eigenvals)
       !
    END IF
    !
    RETURN
  END SUBROUTINE analysis
  !
  !
  !  ----------------------------------------------
  !
  !
  SUBROUTINE analyze_vibrations(U_loc,T_loc,dim,filep,tot_born_charge,eigval_loc, &
       eigvec_loc,do_IR_intensity,intensity_loc)
    !
    USE kinds,              ONLY : DP
    !
    ! ... input variables variables
    !
    INTEGER,        INTENT(IN)  :: dim, filep
    REAL (KIND=DP), INTENT(IN)  :: U_loc(dim,dim), T_loc(dim,dim)
    LOGICAL,        INTENT(IN)  :: do_IR_intensity
    !
    ! ... output variables variables
    !
    REAL (KIND=DP), INTENT(OUT) :: tot_born_charge(3,3)
    REAL (KIND=DP), INTENT(OUT) :: eigval_loc(dim), eigvec_loc(dim,dim)
    REAL (KIND=DP), INTENT(OUT), OPTIONAL :: intensity_loc(dim)
    !
    ! ... local variables
    !
    !
    !
    !
    CALL print_matrix(U_loc,dim,dim," Hessian:  ",filep)
    CALL print_matrix(T_loc,dim,dim," Mass:     ",filep)
    !
    ! ... diagonalizing the hessian to obtain eigen frequencies and modes
    !
    CALL rdiaghg(dim,dim,U_loc,T_loc,dim,eigval_loc,eigvec_loc)
    !
    ! ... calculating IR intensities for each mode
    !
    IF (do_IR_intensity) THEN
       CALL analyze_IR_intensities(tot_born_charge,dim,filep, &
            eigvec_loc,intensity_loc)
    ELSE
       intensity_loc = 0.d0
    END IF
    !
    RETURN
  END SUBROUTINE analyze_vibrations
  !
  !
  !  ----------------------------------------------
  !
  !
  SUBROUTINE apply_asr
    !
    ! Impose acoustic sum rule on Born effective charges
    !  ----------------------------------------------
    !
    USE electrons_base,       ONLY : nel, nspin
    USE io_global,            ONLY : stdout
    USE ions_base,            ONLY : na, zv, nsp, nat
    USE kinds,                ONLY : DP
    !
    REAL (KIND=DP) :: tot_born_charge(3,3)
    REAL (KIND=DP) :: tot_charge, tot_system_charge(3,3)
    !
    INTEGER        :: i
    !
    ! ... Calculate the input total system charge
    !
    tot_charge = DOT_PRODUCT(na(1:nsp),zv(1:nsp))-SUM(nel(1:nspin))
    tot_system_charge      = 0.0
    tot_system_charge(1,1) = tot_charge
    tot_system_charge(2,2) = tot_charge
    tot_system_charge(3,3) = tot_charge
    !
    ! ... Test the sum of born charge tensors
    !
    tot_born_charge=0.0
    DO i=1,nat
       tot_born_charge=tot_born_charge+born_charge(:,3*(i-1)+1:3*(i-1)+3)
    END DO
    tot_born_charge=tot_born_charge-tot_system_charge
    !
    ! ... imposing acoustic sum-rule (zero sum) on Born charges
    !
    tot_born_charge=-tot_born_charge/nat
    !
    DO i=1,nat
       born_charge(:,3*(i-1)+1:3*(i-1)+3)=born_charge(:,3*(i-1)+1:3*(i-1)+3)+tot_born_charge
    END DO
    !
    RETURN
  END SUBROUTINE apply_asr
  !
  !
  !  ----------------------------------------------
  !
  !
  SUBROUTINE trans_invariance(U_loc)
    !
    !-------------------------------------------------------------
    !
    USE ions_base,            ONLY : nat
    USE kinds,                ONLY : DP
    !
    ! ... input/output variables
    !
    REAL (KIND=DP), INTENT(INOUT) :: U_loc(nat*3,nat*3)
    !
    ! ... local variables
    !
    INTEGER                       :: na, nb, i, j, iter
    REAL (KIND=DP)                :: sum, max_sum
    !
    iter = 0
    max_sum = 10 !dummy value
    DO WHILE ((max_sum > trans_inv_conv_thr).AND.(iter < trans_inv_max_iter))
       max_sum=0.0
       iter=iter+1
       !
       ! ... imposing symmetry on the dynamical matrix
       !
       CALL symmetrize_matrix(U_loc,3*nat)
       !
       ! ... imposing translational invariance
       !
       DO i=1,3
          DO j=1,3
             DO na=1,nat
                sum=0.d0
                ! ... symmetric implementation
                DO nb=1,nat
                   sum=sum+U_loc((na-1)*3+i,(nb-1)*3+j)
                END DO
                max_sum=MAX(max_sum,sum)
                !
                ! ... distribute the error over all relevant matrix elements
                sum=sum/nat
                DO nb=1,nat
                   U_loc((na-1)*3+i,(nb-1)*3+j)=U_loc((na-1)*3+i,(nb-1)*3+j)-sum
                END DO
                !
             END DO
          END DO
       END DO
       !
    END DO
    !
    RETURN
  END SUBROUTINE trans_invariance
  !
  !
  !  ----------------------------------------------
  !
  !
  SUBROUTINE animation (tau_loc,eigvec_loc,eigval_loc)
    !------------------------------------------------------------------------------------
    ! ANIMATING VIBRATIONS AS XYZ FILES
    !------------------------------------------------------------------------------------
    !
    USE ions_base,            ONLY : nat
    !
    ! ... input variables
    !
    REAL (KIND=DP), INTENT(IN) :: tau_loc(3,nat)
    REAL (KIND=DP), INTENT(IN) :: eigvec_loc(3*nat,3*nat)
    REAL (KIND=DP), INTENT(IN) :: eigval_loc(3*nat)
    !
    ! ... local variables
    !
    INTEGER                    :: i,it,ia, coord
    CHARACTER (len=3)          :: mode_label
    CHARACTER (len=20)         :: free_text
    REAL      (KIND=DP)        :: tau(3,nat)
    !
    DO i=1,3*nat !loop over modes
       !
       WRITE (mode_label,'(i3.3)') i
       WRITE (free_text,'(2x,E7.2,2x,A6)') eigval_loc(i),' cm-1'
       !
       ! ... generate 20 snapshots of along one vibrational period
       !
       tau=0.0
       DO it=0,19 
          DO ia=1,nat
             DO coord=1,3 !x,y,z
                tau(coord,ia)=tau_loc(coord,ia) +        &
                     ( eigvec_loc((ia-1)*3+coord,i) *    &
                       SIN(4.0d0*ASIN(1.0)*it/20)*SQRT(1822.9))
             END DO
          END DO
          CALL write_xyz &
               (tau,free_text,'APPEND','vib_anim_'//mode_label//'.xyz')
       END DO
    END DO
    !
    RETURN
  END SUBROUTINE animation
  !
  !
  !  ----------------------------------------------
  !
  !
  SUBROUTINE internal_hessian(U_loc,T_loc,nmodes,tau,D,filep)
    !
    ! Subroutine that removes the rigid translations and rotation contributions
    ! to the Hessian. The resulting hessian is of dimension 3*nat-6 (or 5), and
    ! reflects only internal modes of motion.
    !
    USE control_flags,        ONLY : iprsta
    USE kinds,                ONLY : DP
    USE ions_base,            ONLY : nat, pmass, ityp, nsp, na, ions_cofmass
    !
    ! ... input variables
    !
    REAL (KIND=DP), INTENT(IN)    :: tau(3,nat)
    REAL (KIND=DP), INTENT(INOUT) :: U_loc(3*nat,3*nat),T_loc(3*nat,3*nat)
    INTEGER,        INTENT(IN)    :: filep
    !
    ! ... output variables
    !
    REAL (KIND=DP), INTENT(OUT)   :: D(3*nat,3*nat) ! transformatiom matrix
    INTEGER,        INTENT(OUT)   :: nmodes
    !
    ! ... local variables
    !
    INTEGER                       :: i, j, coord
    LOGICAL                       :: linear
    !
    !variables for removing transltion and rotations
    REAL (KIND=DP)                :: Rcom(3), tau_com(nat,3), dummy(3,3)
    REAL (KIND=DP)                :: inertia_tensor(3,3), inertia_moments(3), inertia_eigenvecs(3,3)
    REAL (KIND=DP)                :: xx, yy, zz 
    REAL (KIND=DP)                :: P(nat,3)
    INTEGER                       :: index, zero_moment_index
    !
    ! Transforming U_loc to mass weighted coordintaes and T_loc to unity matrix
    !
    DO i=1,3*nat
       DO j=1,3*nat
          IF (i.NE.j) THEN
             T_loc(i,j)=SQRT(T_loc(i,i)*T_loc(j,j))
          END IF
       END DO
    END DO
    U_loc=U_loc/T_loc
    T_loc=0.0
    DO i=1,3*nat
       T_loc(i,i)=1.0
    END DO
    !
    IF(iprsta.GT.4) THEN
       CALL print_matrix(U_loc,3*nat,3*nat, &
            'Mass weighted Hessian with acoustic sum rule:',filep)
       CALL print_matrix(T_loc,3*nat,3*nat, &
            'Unity kinetic tensot:',filep)
    END IF
    !
    ! ... Calculating center of mass and shifting to com coordinates
    !
    IF(iprsta.GT.4) &
         PRINT *,'tau=',tau

    CALL ions_cofmass(tau,pmass,na,nsp,Rcom)
    IF(iprsta.GT.4) &
         PRINT *,'Rcom=',Rcom
    !
    DO i=1,nat
       tau_com(i,:)=tau(:,i)-Rcom(:)
    END DO
    IF(iprsta.GT.4) &
         PRINT *,'tau_com=',tau_com
    !
    ! ... Generating the intertia tensor
    !
    inertia_tensor=0.0
    DO i=1,nat 
       xx=tau_com(i,1)
       yy=tau_com(i,2)
       zz=tau_com(i,3)
       !
       inertia_tensor(1,1)=inertia_tensor(1,1)+pmass(ityp(i))*(yy*yy+zz*zz)
       inertia_tensor(2,2)=inertia_tensor(2,2)+pmass(ityp(i))*(xx*xx+zz*zz)
       inertia_tensor(3,3)=inertia_tensor(3,3)+pmass(ityp(i))*(xx*xx+yy*yy)
       !
       inertia_tensor(1,2)=inertia_tensor(1,2)-pmass(ityp(i))*(xx*yy)
       inertia_tensor(1,3)=inertia_tensor(1,3)-pmass(ityp(i))*(xx*zz)
       inertia_tensor(2,3)=inertia_tensor(2,3)-pmass(ityp(i))*(yy*zz)
    END DO
    !
    inertia_tensor(2,1)=inertia_tensor(1,2)
    inertia_tensor(3,1)=inertia_tensor(1,3)
    inertia_tensor(3,2)=inertia_tensor(2,3)
    !
    IF(iprsta.GT.4) &
         CALL print_matrix(inertia_tensor,3,3,'Inertia tensor:',filep)
    !
    ! ... diagonalizing the inertia tensor
    !
    dummy=0.0
    dummy(1,1)=1.0
    dummy(2,2)=1.0
    dummy(3,3)=1.0
    !CALL rdiaghg(3,3,inertia_tensor,dummy,3,inertia_moments,inertia_eigenvecs)
    CALL rdiagh(3,inertia_tensor,3,inertia_moments,inertia_eigenvecs)
    !
    IF(iprsta.GT.4) THEN
       WRITE (filep,*)
       WRITE (filep,*) 'Inertia moments: '
       WRITE (filep,'(3x,f10.4,3x,f10.4,3x,f10.4)') (inertia_moments(i),i=1,3)
       WRITE (filep,*)
    END IF
    !
    ! ... Checking if the molecule is linear; 
    ! ... A linear molecule have one zero inertia moment
    ! ... and the other two should be identical.
    !
    linear=.FALSE.
    nmodes=3*nat-6
    i=1
    ! looking for a zero moment, arbitrarily chosen to be
    ! smaller than 0.1
    DO WHILE ((i.LE.3).AND.(inertia_moments(i).GT.0.1))
       i=i+1
    END DO
    !
    IF(i.LT.4) THEN
       !
       ! "zero" moment found, testing if the other two are equal
       !
       IF ( inertia_moments(((i+1)/3)*3+MOD(i+1,3))-          &
            inertia_moments(((i+2)/3)*3+MOD(i+2,3) ).LT.1e-1) THEN
          linear=.TRUE.
          zero_moment_index=i
          nmodes=3*nat-5
       END IF
    END IF
    !
    IF (linear) THEN
       WRITE (filep,*)
       WRITE (filep,*) 'linear=',linear,', it is probably a linear molecule.'
       WRITE (filep,*)
    END IF
    !
    IF(iprsta.GT.4) &
         CALL print_matrix(inertia_eigenvecs,3,3,'inertia eigen vectors:',filep)
    !
    ! Generating the coordinates in rotating and translating frame
    ! D is a projection matrix to internal coordinates.
    ! The first 3 vectors are pure translational, by construction.
    ! The next 3 or 2 vectors are pure rotetional, by construction.
    ! The next 3N-6 vectors are constructed to be orthonormal to the
    ! first 5 or 6 ones, by Graham-Schmidt orthonormalization procedure.
    !
    ! Set D initially to contain random number - better results later
    ! when applying the Graham-Schmidt orthonormalization. I don't care
    ! about the seed since these numbers have no meaning other then defining
    ! orthonormal basis.
    !
    CALL RANDOM_NUMBER(D)
    !
    ! ... generating linear translation basis vector
    !
    DO i=1,nat
       D(1:3,(3*(i-1)+1):(3*(i-1)+3))=dummy*SQRT(T(3*(i-1)+1,3*(i-1)+1))
    END DO
    !
    IF(iprsta.GT.4) &
         CALL print_matrix(D,3*nat,3*nat,'D with translations:',filep)
    !
    ! ... generating pure rotation basis vectors
    !
    P=MATMUL(tau_com,inertia_eigenvecs)
    IF(iprsta.GT.4) &
         CALL print_matrix(P,nat,3,'P:',filep)
    !
    IF (.NOT.linear) THEN
       ! non-linear molecule
       DO i=1,nat
          DO coord=1,3
             D(4,3*(i-1)+coord)=(P(i,2)*inertia_eigenvecs(coord,3)-   &
                  P(i,3)*inertia_eigenvecs(coord,2))/SQRT(T(3*(i-1)+coord,3*(i-1)+coord))
             D(5,3*(i-1)+coord)=(P(i,3)*inertia_eigenvecs(coord,1)-   &
                  P(i,1)*inertia_eigenvecs(coord,3))/SQRT(T(3*(i-1)+coord,3*(i-1)+coord))
             D(6,3*(i-1)+coord)=(P(i,1)*inertia_eigenvecs(coord,2)-   &
                  P(i,2)*inertia_eigenvecs(coord,1))/SQRT(T(3*(i-1)+coord,3*(i-1)+coord))
          END DO
       END DO
    ELSE
       ! linear molecule
       DO i=1,nat
          DO coord=1,3
             !
             index=4
             IF (zero_moment_index.NE.1) THEN
                D(index,3*(i-1)+coord)=(P(i,2)*inertia_eigenvecs(coord,3)-   &
                     P(i,3)*inertia_eigenvecs(coord,2))/SQRT(T(3*(i-1)+coord,3*(i-1)+coord))
                index=index+1
             END IF
             !
             IF (zero_moment_index.NE.2) THEN
                D(index,3*(i-1)+coord)=(P(i,3)*inertia_eigenvecs(coord,1)-   &
                     P(i,1)*inertia_eigenvecs(coord,3))/SQRT(T(3*(i-1)+coord,3*(i-1)+coord))
                index=index+1
             END IF
             !
             IF (zero_moment_index.NE.3) THEN
                D(index,3*(i-1)+coord)=(P(i,1)*inertia_eigenvecs(coord,2)-   &
                     P(i,2)*inertia_eigenvecs(coord,1))/SQRT(T(3*(i-1)+coord,3*(i-1)+coord))
                index=index+1
             END IF
             !
          END DO
       END DO
       !
    END IF
    !
    IF(iprsta.GT.4) &
         CALL print_matrix(D,3*nat,3*nat, &
         'D after adding translation and rotations vectors:',filep)
    !
    ! ... Gram-Schmidt orthonormalization.
    !
    D=TRANSPOSE(D)
    CALL orthonormalize(D,3*nat,3*nat)
    !
    IF(iprsta.GT.4) &
         CALL print_matrix(D,3*nat,3*nat,'D orthonormalized:',filep)
    !
    ! ... Transform the Hessian to internal coordinates
    !
    U_loc=MATMUL(TRANSPOSE(D),MATMUL(U_loc,D))
    !
    RETURN
  END SUBROUTINE internal_hessian
  !
  !
  !  ----------------------------------------------
  !
  !
  SUBROUTINE write_restart (E_minus, dip_minus, cyc_counter)
    !
    !  ----------------------------------------------
    !
    USE io_global,         ONLY : ionode
    USE ions_base,         ONLY : nat      
    !
    REAL (KIND=DP), INTENT(IN) :: E_minus, dip_minus(3)
    INTEGER,        INTENT(IN) :: cyc_counter
    INTEGER                    :: i,j
    !
    IF(ionode) THEN
       OPEN  (200,file=restart_filename, status='unknown', form='unformatted')
       WRITE (200) ref_etot
       WRITE (200) E_minus
       WRITE (200) (dip_minus(i),i=1,3)
       WRITE (200) cyc_counter
       WRITE (200) ((born_charge(i,j),i=1,3),j=1,3*nat)
       WRITE (200) ((U(i,j),i=1,3*nat),j=1,3*nat)
       CLOSE (200)
    END IF
    !
    RETURN
  END SUBROUTINE write_restart
  !
  !
  !  ----------------------------------------------
  !
  !
  SUBROUTINE read_restart (E_minus, dip_minus, cyc_counter)
    !
    !  ----------------------------------------------
    !
    USE io_global,           ONLY : ionode, ionode_id
    USE ions_base,           ONLY : nat      
    USE mp,                  ONLY : mp_bcast
    !
    ! ... output variables
    !
    REAL (KIND=DP),  INTENT(OUT) :: E_minus, dip_minus(3)
    INTEGER,         INTENT(OUT) :: cyc_counter
    !
    ! ... local variables
    !
    INTEGER                      :: i,j
    !
    IF(ionode) THEN
       OPEN  (200,file=restart_filename, status='old', form='unformatted')
       READ  (200) ref_etot
       READ  (200) E_minus
       READ  (200) (dip_minus(i),i=1,3)
       READ  (200) cyc_counter
       READ  (200) ((born_charge(i,j),i=1,3),j=1,3*nat)
       READ  (200) ((U(i,j),i=1,3*nat),j=1,3*nat)
       CLOSE (200)
       !
    END IF
    ! BROADCAST TO ALL NODES
    CALL mp_bcast( cyc_counter,            ionode_id )
    !
    RETURN
  END SUBROUTINE read_restart
  !
  !
  !  ----------------------------------------------
  !
  !
  SUBROUTINE end_vibrations
    !
    DEALLOCATE( ref_c0      )
    DEALLOCATE( ref_tau     )
    DEALLOCATE( ref_lambda  )
    DEALLOCATE( active_atom )
    DEALLOCATE( eigenvals   )
    DEALLOCATE( eigenvecs   )
    DEALLOCATE( U           )
    DEALLOCATE( T           )
    DEALLOCATE( born_charge )
    !
    RETURN
  END SUBROUTINE end_vibrations
  !
  !
  !  ----------------------------------------------
  !
  !
  SUBROUTINE mass_weighted_to_cartesian(dim,T,eigenvecs)
    !
    !converting back the eigen vectors from mass-weighted coordinates
    !(the diagonalization procedure converts it to mass-weighted)
    !
    USE kinds,                ONLY : DP
    !
    ! ... input/output variables
    !
    INTEGER,        INTENT(IN)    :: dim
    REAL (KIND=DP), INTENT(IN)    :: T(dim,dim)
    REAL (KIND=DP), INTENT(INOUT) :: eigenvecs(dim,dim)
    !
    ! ... local variables
    !
    INTEGER :: i,j
    !
    DO i=1,dim
       DO j=1,dim
          eigenvecs(j,i)=eigenvecs(j,i)/SQRT(T(j,j))
       END DO
    END DO
  END SUBROUTINE mass_weighted_to_cartesian
  !
  !
  !  ----------------------------------------------
  !
  !
  SUBROUTINE print_eigenmodes(dim,cm_inv_flag,filep,eigval,eigvec, &
       mode_mass, mode_force_constant, intensity)
    !
    USE constants,         ONLY : BOHR_RADIUS_ANGS, AU, AMU_AU
    USE kinds,             ONLY : DP
    !
    ! ... input variables
    !
    INTEGER,        INTENT(IN) :: dim, filep
    LOGICAL,        INTENT(IN) :: cm_inv_flag
    REAL (KIND=DP), INTENT(IN) :: eigval(dim), eigvec(dim,dim)
    REAL (KIND=DP), INTENT(IN) :: mode_mass(dim), mode_force_constant(dim)
    REAL (KIND=DP), INTENT(IN), OPTIONAL   :: intensity(dim)
    !
    ! ... local variables
    !
    INTEGER                    :: i
    REAL (KIND=DP)             :: eigval_loc(dim), intensity_loc(dim)
    !
    !
    eigval_loc = eigval
    IF (PRESENT(intensity)) THEN
       intensity_loc = intensity
    ELSE
       intensity_loc = 0.d0
    END IF
    !
    ! ... converting to frequency
    !
    IF (cm_inv_flag) THEN
       !
       ! ... converting from omega [a.u.] to f in [cm^-1] (omega=2*pi*f)
       !
       DO i=1,dim
          IF(eigval_loc(i).LT.0.0) THEN
             !
             !... imaginary frequency, presented as negative frequency
             !
             eigval_loc(i)=-SQRT(-eigval_loc(i))*2.194746313709741e+05
          ELSE
             eigval_loc(i)= SQRT( eigval_loc(i))*2.194746313709741e+05
          END IF
       END DO
    ELSE
       !
       ! ... omega in a.u.
       !
       DO i=1,dim
          IF(eigval_loc(i).LT.0.0) THEN
             !
             !... imaginary frequency, presented as negative frequency
             !
             eigval_loc(i)=-SQRT(-eigval_loc(i))
          ELSE
             eigval_loc(i)= SQRT( eigval_loc(i))
          END IF
       END DO
    END IF
    !
    !... Writing eigenvalues to file
    !
    WRITE (filep,*)
    IF (cm_inv_flag) THEN
       WRITE (filep,*) 
       WRITE (filep,200) 
       WRITE (filep,201)
       WRITE (filep,202) 
       WRITE (filep,203)
    ELSE
       WRITE (filep,205) 
    END IF
    !
    DO i=1,dim
       WRITE (filep,204,ADVANCE='YES')                                         &
            eigval_loc(i), mode_mass(i)/ AMU_AU,                               &
            mode_force_constant(i) * AU / BOHR_RADIUS_ANGS / BOHR_RADIUS_ANGS, &
            intensity_loc(i)
    END DO
    WRITE (filep,*)
    !
    CALL print_matrix(eigvec,dim,dim, &
         'Eigen vectors (column wise) ordered like the eigen-values:',filep)
    !
200 FORMAT(3x,'Normal modes eigen frequencies: f [cm^-1] (f=omega/2*pi)')
201 FORMAT(3x,'--------------------------------------------------------')
202 FORMAT(3x,'f [cm-1]        mode mass [amu]        effective  force        intensity   ')
203 FORMAT(3x,'                                       constant [eV/A^2]       [hartree au]')
204 FORMAT(3x,f8.2,9x,f14.2,9x,f16.4,9x,E10.3)
205 FORMAT(3x,'Normal modes eigen frequencies: omega [au] (omega=2*pi*f)')
    !
    RETURN
  END SUBROUTINE print_eigenmodes
  !
  !
  !  ----------------------------------------------
  !
  !
  SUBROUTINE analyze_IR_intensities(tot_born_charge,nmodes,filep,eigvec_loc,intensity)
    !
    ! Analysis of the Born effective charges and the infra-red cross-section
    ! Here U is the Hessian, T is the cartesian mass matrix
    !
    USE control_flags,     ONLY  : iprsta
    USE electrons_base,    ONLY  : nel, nspin
    USE input_parameters,  ONLY  : atom_label
    USE ions_base,         ONLY  : nat, ityp, na, nsp, zv 
    !
    ! ... input variables
    !
    INTEGER       , INTENT(IN)  :: filep,nmodes
    REAL (KIND=DP), INTENT(IN)  :: eigvec_loc(nmodes,nmodes)
    !
    ! ... output variables
    !
    REAL (KIND=DP), INTENT(OUT) :: tot_born_charge(3,3)
    REAL (KIND=DP), INTENT(OUT), OPTIONAL :: intensity(nmodes)
    !
    ! ... local variables
    !
    INTEGER                     :: i, ia, coord, index
    REAL (KIND=DP)              :: mode_charge(3,3*nat)
    REAL (KIND=DP)              :: tot_charge, tot_system_charge(3,3)
    !
    !-----------------------------------------------------------------
    ! PRINTING AND TESTING THE BORN CHARGES
    !-----------------------------------------------------------------    
    !
    ! ... Write Born effective charge matrixes
    !
    WRITE (filep,*) ''
    WRITE (filep,*) 'Born effective charges:'
    DO ia=1,nat 
       WRITE (filep,*) ''
       WRITE (filep,*) 'Atom Label: ',TRIM(atom_label(ityp(ia)))
       WRITE (filep,*) ''
       DO coord=1,3             !x,y,z
          index=(ia-1)*3+coord  !index of atom1 coordinate in continueous counting
          WRITE (filep,'(3(f8.3,3X))') born_charge(1,index),born_charge(2,index),born_charge(3,index)
       END DO
    END DO
    !
    ! ... Calculate the input total system charge
    !
    tot_charge=DOT_PRODUCT(na(1:nsp),zv(1:nsp))-SUM(nel(1:nspin))
    tot_system_charge=0.0
    tot_system_charge(1,1)=tot_charge
    tot_system_charge(2,2)=tot_charge
    tot_system_charge(3,3)=tot_charge
    WRITE (filep,*)
    WRITE (filep,'(A20,f10.2)') 'Total system charge:', tot_charge
    WRITE (filep,*)
    !
    ! ... Test the sum of born charge tensors
    !
    tot_born_charge=0.0
    DO i=1,nat
       tot_born_charge=tot_born_charge+born_charge(:,3*(i-1)+1:3*(i-1)+3)
    END DO
    tot_born_charge=tot_born_charge-tot_system_charge
    !
    CALL print_matrix(tot_born_charge,3,3,&
         'Sum of Born charge tensors (minus system charge):',filep)
    !
    !-----------------------------------------------------------------
    ! CALCULATING IR INTENSITIES
    !-----------------------------------------------------------------
    !
    ! ... Calculating the mode effective charge
    !
    mode_charge=MATMUL(born_charge,eigvec_loc)
    !
    IF (iprsta.GT.4) &
         CALL print_matrix(mode_charge,3,3*nat,'Mode effective charge:',filep)
    !
    ! ... IR intensity
    !
    IF (PRESENT(intensity)) THEN
       DO i=3*nat-nmodes+1,3*nat
          intensity(i) = SUM(mode_charge(:,i)*mode_charge(:,i))
          !       WRITE (filep,FMT='(E12.5,3X)',ADVANCE='NO') SUM(mode_charge(:,i)*mode_charge(:,i))
       END DO
    END IF
    !
    RETURN
  END SUBROUTINE analyze_IR_intensities
  !
  !
  !  ----------------------------------------------
  !
  !
  SUBROUTINE get_force_constant &
       (eigvec_loc, eigval_loc, dim, mode_mass_loc, mode_force_constant_loc)
    !
    USE kinds,              ONLY : DP
    !
    ! ... input variables
    !
    INTEGER                     :: dim
    REAL (KIND=DP), INTENT(IN)  :: eigvec_loc(dim,dim), eigval_loc(dim)
    !
    ! ... ouput variables
    !
    REAL (KIND=DP), INTENT(OUT) :: mode_mass_loc(dim)
    REAL (KIND=DP), INTENT(OUT) :: mode_force_constant_loc(dim)    
    !
    ! ... local variables
    !
    INTEGER                     :: i
    REAL (KIND=DP)              :: tmp(dim,dim)
    !
    !
    ! ... compute mode effective mass and force constant
    !
    tmp=MATMUL(TRANSPOSE(eigvec_loc),eigvec_loc)
    DO i=1,dim
       mode_mass_loc(i)=1/tmp(i,i)
       mode_force_constant_loc(i)=eigval_loc(i)*mode_mass_loc(i)
    END DO
    !
    RETURN
  END SUBROUTINE get_force_constant

  !=----------------------------------------------------------------
END MODULE vibrations
!=----------------------------------------------------------------
