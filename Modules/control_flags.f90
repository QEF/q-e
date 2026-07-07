!
! Copyright (C) 2002-2025 Quantum ESPRESSO Fpundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------------=!
MODULE control_flags
  !=--------------------------------------------------------------------------=!
  !! This module contains all basic variables that controls the execution flow.
  !----------------------------------------------
  !
  USE kinds
  USE parameters
  !
  IMPLICIT NONE
  !
  SAVE
  !
  ! ...   declare execution control variables
  !
  LOGICAL :: lforce        = .FALSE. ! compute forces
  LOGICAL :: tv0rd         = .FALSE. ! read ionic velocities from standard input
  LOGICAL :: tstress       = .FALSE. ! compute stress
  LOGICAL :: lecrpa        = .FALSE. ! RPA correlation energy request
  LOGICAL :: tddfpt        = .FALSE. ! use TDDFPT specific tweaks when using the Environ plugin
  LOGICAL :: smallmem      = .FALSE. ! reduce memory by avoiding global sort
  !
  INTEGER :: iprint =10 ! print output every iprint step
  INTEGER :: max_xml_steps =0 ! max number of dynamics included in xml file if all steps are included. 
  !
  ! ... .TRUE. if only gamma point is used
  !
  LOGICAL :: gamma_only = .TRUE.
  !
  ! ... Flag controlling the Nose thermostat for the cell
  !
  LOGICAL :: tnoseh = .FALSE.
  !
  ! ... Flag controlling the Nose thermostat for ions
  !
  LOGICAL  :: tnosep = .FALSE.
  REAL(DP) :: tolp   = 0.0_DP   !  tolerance for temperature variation
  !
  REAL(DP) :: &
       ekin_conv_thr = 0.0_DP, &!  conv. threshold for fictitious e. kinetic energy
       etot_conv_thr = 0.0_DP, &!  conv. threshold for DFT energy
       forc_conv_thr = 0.0_DP   !  conv. threshold for atomic forces
  INTEGER :: &
       ekin_maxiter = 100,   &!  max number of iter. for ekin convergence
       etot_maxiter = 100,   &!  max number of iter. for etot convergence
       forc_maxiter = 100     !  max number of iter. for atomic forces conv.
  !
  ! ... Several variables controlling the run ( used mainly in PW calculations )
  !
  ! ... logical flags controlling the execution
  !
  LOGICAL :: &
    lscf    =.FALSE., &! if .TRUE. the calc. is selfconsistent
    lbfgs   =.FALSE., &! if .TRUE. the calc. is a relaxation based on BFGS
    lmd     =.FALSE., &! if .TRUE. the calc. is a dynamics
    lbands  =.FALSE., &! if .TRUE. the calc. is band structure
    lconstrain=.FALSE.,&! if .TRUE. the calc. is constraint
    llondon =.FALSE., & ! if .TRUE. compute Grimme D2 dispersion corrections
    ldftd3 =.FALSE., & ! if .TRUE. compute Grimme D3 dispersion corrections
    ts_vdw  =.FALSE., & ! as above for Tkatchenko-Scheffler disp.corrections
    mbd_vdw  =.FALSE., &!as above for MBD correction
    lxdm    =.FALSE., & ! if .TRUE. compute XDM dispersion corrections
    lensemb =.FALSE., &! if .TRUE. compute ensemble energies
    restart =.FALSE. ! if .TRUE. restart from results of a preceding run
  !
  ! ... pw self-consistency
  !
  INTEGER :: &
    ngm0,             &! used in mix_rho
    nexxiter,         &! the maximum number of outer iteration (exx)
    niter,            &! the maximum number of iteration
    nmix,             &! the number of iteration kept in the history
    imix               ! the type of mixing (0=plain,1=TF,2=local-TF)
  INTEGER :: &
    n_scf_steps        ! number of scf iterations to reach convergence
  REAL(DP) :: &
    mixing_beta,      &! the mixing parameter
    tr2,              &! the convergence threshold for potential
    scf_error=0.0      ! actual convergence reached

  LOGICAL :: &
    conv_elec          ! if .TRUE. electron convergence has been reached
  ! next 3 variables used for EXX calculations
  LOGICAL :: &
    adapt_thr       ! if .TRUE. an adaptive convergence threshold is used
                       ! for the scf cycle in an EXX calculation.
  REAL(DP)  :: &
    tr2_init,         &! initial value of tr2 for adaptive thresholds
    tr2_multi          ! the dexx multiplier for adaptive thresholds
                       ! tr2 = tr2_multi * dexx after each V_exx update 
  LOGICAL :: scf_must_converge
  !
  ! ... pw diagonalization
  !
  REAL(DP)  :: &
    ethr               ! the convergence threshold for eigenvalues
  INTEGER :: &
    isolve,           &! index selecting Davidson,  CG, ParO or RMM diagonalization
    david,            &! max dimension of subspace in Davidson diagonalization
    max_cg_iter,      &! maximum number of iterations in a CG call
    rmm_ndim,         &! max dimension of subspace in RMM-DIIS diagonalization
    gs_nblock          ! blocking size in Gram-Schmidt orthogonalization
  LOGICAL :: &
    rmm_conv,                     &! if true,  RMM-DIIS is performed up to converge
    rmm_with_davidson  = .TRUE.,  &! if true RMM-DIIS  in alternance with davidson 
    diago_full_acc     = .FALSE.      ! if true,  empty eigenvalues have the same
                                   ! accuracy of the occupied ones
  !
  ! ... ionic dynamics
  !
  INTEGER :: &
    nstep = 1,       &! number of ionic steps
    istep = 0          ! current ionic step
  LOGICAL :: &
    conv_ions          ! if .TRUE. ionic convergence has been reached
  REAL(DP)  :: &
    upscale            ! maximum reduction of convergence threshold
  !
  ! ... system's symmetries
  !
  LOGICAL :: &
    noinv = .FALSE.    ! if .TRUE. q=>-q symmetry not used in k-point generation
  LOGICAL :: symm_by_label = .FALSE. ! use atomic labels to detect symmetry 
  LOGICAL :: use_spinflip = .FALSE.      ! in collinear case add allow rotations + spinflip 
  !
  ! ... phonon calculation
  !
  INTEGER :: &
    modenum            ! for single mode phonon calculation
  !
  ! ... printout control
  !
  INTEGER :: &
    io_level = 1       ! variable controlling the amount of I/O to file
  INTEGER :: & ! variable controlling the amount of I/O to output
    iverbosity = 0     ! -1 minimal, 0 low, 1 medium, 2 high, 3 debug
  !
  ! ... self-interaction correction and scissor operator
  !
  LOGICAL :: sic = .FALSE.
  LOGICAL :: scissor = .FALSE.
  !
  ! ... miscellany
  !
  LOGICAL :: &
    use_para_diag = .FALSE.  ! if .TRUE. a fully distributed memory iteration 
                             ! algorithm and parallel Householder algorithm are used
  !
  LOGICAL :: &
    remove_rigid_rot = .FALSE. ! if .TRUE. the total torque acting on the atoms 
                               ! is removed
  LOGICAL :: &
    do_makov_payne = .FALSE.   ! if .TRUE. makov-payne correction for isolated
                               ! system is used
  LOGICAL :: &
    use_gpu = .FALSE.          ! if .TRUE. selects the accelerated version of the subroutines
                               ! when available (obsolescent, should be removed)
  !
  TYPE(offload_kind_acc) :: offload_acc  ! flag to select CUF/OpenACC offload type
  TYPE(offload_kind_omp) :: offload_omp  ! flag to select OpenMP5 offload type
  TYPE(offload_kind_cpu) :: offload_cpu  ! flag to select no offload type (CPU execution)
#if defined(_OPENACC)
  TYPE(offload_kind_acc) :: offload_type ! flag to point the actual currently used offload type 
#else
  TYPE(offload_kind_cpu) :: offload_type
#endif
  !
  INTEGER :: &
#if defined(__CUDA)
    many_fft = 16              ! the size of FFT batches in vloc_psi and
                               ! sumband. Only use in accelerated subroutines.
#else
    many_fft = 1
#endif
  !
  ! ... Real-space algorithms
  !
  LOGICAL :: tqr=.FALSE. ! if true the Q are in real space
  !
  ! ... Augmentation charge and beta smoothing
  !
  LOGICAL :: tq_smoothing=.FALSE. ! if true the Q are smoothed 
  LOGICAL :: tbeta_smoothing=.FALSE. ! if true the betas are smoothed 
  !
  ! ... External Forces on Ions
  !
  LOGICAL :: textfor = .FALSE.
  !
  LOGICAL :: treinit_gvecs = .FALSE.
  !
END MODULE control_flags

