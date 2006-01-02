!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE read_namelists_module
  !----------------------------------------------------------------------------
  !
  !  ... this module handles the reading of input namelists
  !  ... written by: Carlo Cavazzoni
  !  --------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : factem, kb_au, au_kb, k_boltzman_au, angstrom_au, &
                        amu_au, pi, e2
  USE input_parameters
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  REAL(DP), PARAMETER :: sm_not_set = -20.D0
  !
  PUBLIC :: read_namelists, sm_not_set
  !
  !  ... end of module-scope declarations
  !
  !  ----------------------------------------------
  !
  CONTAINS
     !
     !=-----------------------------------------------------------------------=!
     !
     !  Variables initialization for Namelist CONTROL
     !
     !=-----------------------------------------------------------------------=!
     !
     !----------------------------------------------------------------------
     SUBROUTINE control_defaults( prog )
       !----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       CHARACTER(LEN=2) :: prog   ! ... specify the calling program
       !
       !
       IF ( prog == 'PW' ) THEN
          title = ' ' 
          calculation = 'scf' 
       ELSE
          title = 'MD Simulation' 
          calculation = 'cp' 
       END IF
       verbosity = 'default'
       IF( prog == 'PW' ) restart_mode = 'from_scratch' 
       IF( prog == 'CP' ) restart_mode = 'restart' 
       nstep  = 50
       IF( prog == 'PW' ) iprint = 100000
       IF( prog == 'CP' ) iprint = 10
       IF( prog == 'PW' ) isave = 0
       IF( prog == 'CP' ) isave = 100
       !
       tstress = .FALSE.
       tprnfor = .FALSE.
       !
       IF( prog == 'PW' ) dt  = 20.D0
       IF( prog == 'CP' ) dt  =  1.D0
       !
       ndr = 50
       ndw = 50
       !
       ! ... use the path specified as outdir and the filename prefix to store 
       ! ... the output
       !
       outdir  = './'   
       scradir = './'   
       IF( prog == 'PW' ) prefix = 'pwscf'  
       IF( prog == 'CP' ) prefix = 'cp' 
       !
       ! ... directory containing the pseudopotentials
       !
       pseudo_dir    = './'  
       refg          = 0.05d0
       max_seconds   = 1.D+7
       ekin_conv_thr = 1.D-6
       etot_conv_thr = 1.D-4
       forc_conv_thr = 1.D-3
       disk_io  = 'default'
       dipfield = .FALSE.
       lberry   = .FALSE.
       gdir     = 0
       nppstr   = 0
       wf_collect = .FALSE.
       printwfc = -1
       lelfield = .FALSE.
       nberrycyc  = 1

       ! 
       RETURN
       !
     END SUBROUTINE
     !
     !=----------------------------------------------------------------------=!
     !
     !  Variables initialization for Namelist SYSTEM
     !
     !=----------------------------------------------------------------------=!
     !
     !----------------------------------------------------------------------
     SUBROUTINE system_defaults( prog )
       !---------------------------------------------------------------------- 
       !
       IMPLICIT NONE
       !
       CHARACTER(LEN=2) :: prog   ! ... specify the calling program
       ! 
       !
       ibrav  = -1
       celldm = (/ 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0 /)
       a = 0.D0 
       b = 0.D0
       c = 0.D0
       cosab = 0.D0
       cosac = 0.D0
       cosbc = 0.D0
       nat    = 0
       ntyp   = 0
       nbnd   = 0
       nelec  = 0.D0
       tot_charge = 0.D0
       tot_magnetization = -1
       multiplicity = 0
       ecutwfc = 0.D0
       ecutrho = 0.D0
       nr1  = 0
       nr2  = 0
       nr3  = 0
       nr1s = 0
       nr2s = 0
       nr3s = 0
       nr1b = 3
       nr2b = 3
       nr3b = 3
       occupations = 'fixed'
       smearing = 'gaussian'
       degauss = 0.D0
       nelup = 0.D0
       neldw = 0.D0
       nspin = 1
       nosym = .FALSE.
       ecfixed = 0.D0
       qcutz   = 0.D0
       q2sigma = 0.01D0
       xc_type = 'none'
       input_dft = 'none'
!
! ... set starting_magnetization to an invalid value:
! ... in PW starting_magnetization MUST be set for at least one atomic type
! ... in CP starting_magnetization MUST REMAIN UNSET 
!
       starting_magnetization = sm_not_set

       IF ( prog == 'PW' ) THEN
          !
          starting_ns_eigenvalue = -1.D0
          U_projection_type = 'atomic'
          !
       END IF
       lda_plus_U = .FALSE.
       Hubbard_U = 0.D0
       Hubbard_alpha = 0.D0
       edir = 1
       emaxpos = 0.5D0
       eopreg = 0.1D0
       eamp = 1.0D-3
       !
       !  ... postprocessing of DOS & phonons & el-ph
       la2F = .FALSE.
       !
       ! ... non collinear program variables
       !
       lspinorb = .FALSE.
       noncolin = .FALSE.
       lambda = 1.0
       constrained_magnetization= 'none'
       fixed_magnetization = 0.d0
       B_field = 0.d0
       angle1 = 0.0
       angle2 = 0.0
       report = 1
       ! 
       RETURN
       ! 
     END SUBROUTINE
     !
     !=----------------------------------------------------------------------=!
     !
     !  Variables initialization for Namelist ELECTRONS
     !
     !=----------------------------------------------------------------------=!
     !
     !----------------------------------------------------------------------
     SUBROUTINE electrons_defaults( prog )
       !----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       CHARACTER(LEN=2) :: prog   ! ... specify the calling program
       ! 
       !
       emass = 400.D0
       emass_cutoff = 2.5D0
       orthogonalization = 'ortho'
       ortho_eps = 1.D-8
       ortho_max = 20
       electron_maxstep = 100
       !
       ! ... ( 'sd' | 'cg' | 'damp' | 'verlet' | 'none' | 'diis' )
       !
       electron_dynamics = 'none'  
       electron_damping = 0.1D0
       !
       ! ... ( 'zero' | 'default' )
       !
       electron_velocities = 'default' 
       !
       ! ... ( 'nose' | 'not_controlled' | 'rescaling')
       !
       electron_temperature = 'not_controlled' 
       ekincw = 0.001D0
       fnosee = 1.D0
       ampre  = 0.D0
       grease = 1.D0
       IF ( prog == 'PW' ) THEN
          !
          startingwfc = 'atomic'
          startingpot = 'atomic'          
          !      
       ELSE
          !
          startingwfc = 'random'
          startingpot = ' '
          !
       END IF
       conv_thr = 1.D-6
       empty_states_nbnd = 0
       empty_states_maxstep = 100
       empty_states_delt = 0.D0
       empty_states_emass = 0.D0
       empty_states_ethr = 0.D0
       diis_size = 4
       diis_nreset = 3
       diis_hcut = 1.D0
       diis_wthr = 0.D0
       diis_delt = 0.D0
       diis_maxstep = 100
       diis_rot = .FALSE.
       diis_fthr = 0.D0
       diis_temp = 0.D0
       diis_achmix = 0.D0
       diis_g0chmix = 0.D0
       diis_g1chmix = 0.D0
       diis_nchmix = 3
       diis_nrot = 3
       diis_rothr  = 0.D0
       diis_ethr   = 0.D0
       diis_chguess = .FALSE.
       mixing_mode = 'plain'
       mixing_fixed_ns = 0
       mixing_beta = 0.7D0
       mixing_ndim = 8
       diagonalization = ' '
       diago_thr_init = 0.D0
       diago_cg_maxiter = 20
       diago_david_ndim = 4
       diago_diis_ndim = 3
       !
       sic = 'none' 
       sic_epsilon = 0.D0
       sic_alpha = 0.D0
       force_pairing = .false.
       ! 
       fermi_energy = 0.D0
       n_inner = 2
       rotation_dynamics = "line-minimization"
       occupation_dynamics = "line-minimization"
       rotmass = 0.D0
       occmass = 0.D0
       rotation_damping = 0.D0
       occupation_damping = 0.D0
       !
       tcg     = .FALSE.
       maxiter = 40
       passop  = 0.3D0
       etresh  = 1.D-6
       !
       epol   = 3
       efield = 0.D0
       ! 
       RETURN
       ! 
     END SUBROUTINE
     !
     !=----------------------------------------------------------------------=!
     !
     !  Variables initialization for Namelist IONS
     !
     !=----------------------------------------------------------------------=!
     !
     !----------------------------------------------------------------------
     SUBROUTINE ions_defaults( prog )
       !----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       CHARACTER(LEN=2) :: prog   ! ... specify the calling program
       !
       !
       ! ... ( 'full' | 'coarse-grained' )
       !
       phase_space = 'full'
       !
       ! ... ( 'sd' | 'cg' | 'damp' | 'verlet' | 'none' )
       ! ... ( 'constrained-verlet' | 'bfgs' | 'constrained-damp' | 'beeman' )
       !
       ion_dynamics = 'none'
       ion_radius   = 0.5D0
       ion_damping  = 0.1
       !
       ! ... ( 'default' | 'from_input' )
       !
       ion_positions = 'default'
       !
       ! ... ( 'zero' | 'default' | 'from_input' )
       !
       ion_velocities = 'default' 
       !
       ! ... ( 'nose' | 'not_controlled' | 'rescaling' )
       !
       ion_temperature   = 'not_controlled'
       tempw             = 300.D0
       fnosep            = -1.0D0
       fnosep(1)         = 1.0D0
       nhpcl             = 0
       nhptyp            = 0
       ndega             = 0
       tranp             = .FALSE.
       amprp             = 0.D0
       greasp            = 1.D0
       tolp              = 100.D0
       ion_nstepe        = 1
       ion_maxstep       = 100
       delta_t           = 1.D0
       nraise            = 100
       upscale           = 10
       pot_extrapolation = 'atomic'
       wfc_extrapolation = 'none'
       !
       ! ... defaults for "path" optimisations
       !
       num_of_images   = 0
       first_last_opt  = .FALSE.
       use_masses      = .FALSE.       
       write_save      = .FALSE.
       use_fourier     = .FALSE.
       use_freezing    = .FALSE.
       opt_scheme      = 'quick-min'
       damp            = 1.D0
       temp_req        = 0.D0
       ds              = 1.D0
       path_thr        = 0.05D0
       !
       ! ... NEB specific
       !
       CI_scheme = 'no-CI'
       k_max     = 0.1D0
       k_min     = 0.1D0
       !
       ! ... SMD specific
       !
       init_num_of_images = 3
       use_multistep      = .FALSE.
       fixed_tan          = .FALSE.
       free_energy        = .FALSE.
       !
       ! ... BFGS defaults
       !
       bfgs_ndim        = 1
       trust_radius_max = 0.8D0
       trust_radius_min = 1.D-3
       trust_radius_ini = 0.5D0
       w_1              = 0.01D0
       w_2              = 0.50D0
       !
       sic_rloc = 0.D0
       !
       ! ... SMD defaults (Y.K. 15/04/2004 )
       !
       smd_polm    = .FALSE.
       smd_kwnp    = 2
       smd_linr    = .FALSE.
       smd_stcd    = .FALSE.
       smd_stcd1   = 0
       smd_stcd2   = 0
       smd_stcd3   = 0
       smd_codf    = 50
       smd_forf    = 50
       smd_smwf    = 1
       smd_lmfreq  = 1
       smd_tol     = 1.D-4
       smd_maxlm   = 10
       smd_smcp    = .TRUE.
       smd_smopt   = .FALSE.
       smd_smlm    = .FALSE.
       smd_splc    = .FALSE.
       smd_spal    = 1.D0
       smd_ene_ini = 1.D0
       smd_ene_fin = 1.D0
       !
       ! ... meta-dynamics defaults
       !
       fe_step     = 0.4D0
       fe_nstep    = 100
       shake_nstep = 10
       g_amplitude = 0.005D0
       !
       RETURN
       !
     END SUBROUTINE
     !
     !=----------------------------------------------------------------------=!
     !
     !  Variables initialization for Namelist CELL
     !
     !=----------------------------------------------------------------------=!
     !
     !----------------------------------------------------------------------
     SUBROUTINE cell_defaults( prog )
       !----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       CHARACTER(LEN=2) :: prog   ! ... specify the calling program
       ! 
       !
       cell_parameters = 'default'
       !
       ! ... ( 'sd' | 'pr' | 'none' | 'w' | 'damp-pr' | 'damp-w' )
       !
       cell_dynamics = 'none'
       !
       ! ... ( 'zero' | 'default' )
       !
       cell_velocities = 'default' 
       press = 0.D0
       wmass = 0.D0
       !
       ! ... ( 'nose' | 'not_controlled' | 'rescaling' )
       !
       cell_temperature = 'not_controlled'
       temph = 0.D0
       fnoseh = 1.D0
       greash = 1.D0
       !
       ! ... ('all'* | 'volume' | 'x' | 'y' | 'z' | 'xy' | 'xz' | 'yz' | 'xyz' )
       !
       cell_dofree = 'all' 
       cell_factor = 0.D0
       cell_nstepe = 1
       cell_damping = 0.D0
       !
       RETURN
       !
     END SUBROUTINE
     !
     !=----------------------------------------------------------------------=!
     !
     !  Variables initialization for Namelist PHONON
     !
     !=----------------------------------------------------------------------=!
     !
     !----------------------------------------------------------------------
     SUBROUTINE phonon_defaults( prog )
       !----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       CHARACTER(LEN=2) :: prog   ! ... specify the calling program
       !
       !
       modenum = 0
       xqq = 0.D0
       ! 
       RETURN
       !
     END SUBROUTINE
     !
     !=----------------------------------------------------------------------=!
     !
     !  Variables initialization for Namelist RAMAN
     !
     !----------------------------------------------------------------------
     SUBROUTINE raman_defaults( prog )
       !----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       CHARACTER(LEN=2) :: prog   ! ... specify the calling program
       !
       !
       b_length = 0.d0
       lcart = .false.
       ! 
       RETURN
       !
     END SUBROUTINE
     !
     !=----------------------------------------------------------------------=!
     !
     !  Variables initialization for Namelist WANNIER
     !
     !----------------------------------------------------------------------
     SUBROUTINE wannier_defaults( prog )
       !----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       CHARACTER(LEN=2) :: prog   ! ... specify the calling program
       !
       !
       wf_efield = .FALSE.
       wf_switch = .FALSE.
       !   
       sw_len = 1
       !   
       efx0 = 0.D0
       efy0 = 0.D0
       efz0 = 0.D0
       efx1 = 0.D0
       efy1 = 0.D0
       efz1 = 0.D0
       !   
       wfsd = .FALSE.
       !
       wfdt        = 5.D0
       maxwfdt     = 0.3D0
       wf_q        = 1500.D0
       wf_friction = 0.3D0
       !   
       nit    = 10
       nsd    = 10
       nsteps = 20
       !   
       tolw = 1.D-8
       !
       adapt = .FALSE.
       !   
       calwf  = 3
       nwf    = 0
       wffort = 40
       !
       writev = .FALSE.
       !
       RETURN
       !
     END SUBROUTINE     
     !
     !=----------------------------------------------------------------------=!
     !
     !  Broadcast variables values for Namelist CONTROL
     !
     !=----------------------------------------------------------------------=!
     !
     !----------------------------------------------------------------------
     SUBROUTINE control_bcast()
       !----------------------------------------------------------------------
       !
       USE io_global, ONLY : ionode_id
       USE mp,        ONLY : mp_bcast
       !
       IMPLICIT NONE
       !
       ! 
       CALL mp_bcast( title,         ionode_id )
       CALL mp_bcast( calculation,   ionode_id )
       CALL mp_bcast( verbosity,     ionode_id )
       CALL mp_bcast( restart_mode,  ionode_id )
       CALL mp_bcast( nstep,         ionode_id )
       CALL mp_bcast( iprint,        ionode_id )
       CALL mp_bcast( isave,         ionode_id )
       CALL mp_bcast( tstress,       ionode_id )
       CALL mp_bcast( tprnfor,       ionode_id )
       CALL mp_bcast( dt,            ionode_id )
       CALL mp_bcast( ndr,           ionode_id )
       CALL mp_bcast( ndw,           ionode_id )
       CALL mp_bcast( outdir,        ionode_id )
       CALL mp_bcast( wfcdir,       ionode_id )
       CALL mp_bcast( scradir,       ionode_id )
       CALL mp_bcast( prefix,        ionode_id )
       CALL mp_bcast( max_seconds,   ionode_id )
       CALL mp_bcast( ekin_conv_thr, ionode_id )
       CALL mp_bcast( etot_conv_thr, ionode_id )
       CALL mp_bcast( forc_conv_thr, ionode_id )
       CALL mp_bcast( pseudo_dir,    ionode_id )
       CALL mp_bcast( refg,          ionode_id )
       CALL mp_bcast( disk_io,       ionode_id )
       CALL mp_bcast( tefield,       ionode_id )
       CALL mp_bcast( dipfield,      ionode_id )
       CALL mp_bcast( lberry,        ionode_id )
       CALL mp_bcast( gdir,          ionode_id )
       CALL mp_bcast( nppstr,        ionode_id )
       CALL mp_bcast( wf_collect,    ionode_id )
       CALL mp_bcast( printwfc,      ionode_id )
       CALL mp_bcast( lelfield,      ionode_id )
       CALL mp_bcast( nberrycyc,     ionode_id )

       ! 
       RETURN
       !
     END SUBROUTINE
     !
     !=----------------------------------------------------------------------=!
     !
     !  Broadcast variables values for Namelist SYSTEM
     !
     !=----------------------------------------------------------------------=!
     !
     !----------------------------------------------------------------------
     SUBROUTINE system_bcast()
       !----------------------------------------------------------------------
       !
       USE io_global, ONLY : ionode_id
       USE mp,        ONLY : mp_bcast
       !
       IMPLICIT NONE
       ! 
       !
       CALL mp_bcast( ibrav,                  ionode_id )
       CALL mp_bcast( celldm,                 ionode_id )
       CALL mp_bcast( a,                      ionode_id )
       CALL mp_bcast( b,                      ionode_id )
       CALL mp_bcast( c,                      ionode_id )
       CALL mp_bcast( cosab,                  ionode_id )
       CALL mp_bcast( cosac,                  ionode_id )
       CALL mp_bcast( cosbc,                  ionode_id )
       CALL mp_bcast( nat,                    ionode_id )
       CALL mp_bcast( ntyp,                   ionode_id )
       CALL mp_bcast( nbnd,                   ionode_id )
       CALL mp_bcast( nelec,                  ionode_id )
       CALL mp_bcast( tot_charge,             ionode_id )
       CALL mp_bcast( tot_magnetization,      ionode_id )
       CALL mp_bcast( multiplicity,           ionode_id )
       CALL mp_bcast( ecutwfc,                ionode_id )
       CALL mp_bcast( ecutrho,                ionode_id )
       CALL mp_bcast( nr1,                    ionode_id )            
       CALL mp_bcast( nr2,                    ionode_id )           
       CALL mp_bcast( nr3,                    ionode_id )          
       CALL mp_bcast( nr1s,                   ionode_id )        
       CALL mp_bcast( nr2s,                   ionode_id )       
       CALL mp_bcast( nr3s,                   ionode_id )      
       CALL mp_bcast( nr1b,                   ionode_id )     
       CALL mp_bcast( nr2b,                   ionode_id )    
       CALL mp_bcast( nr3b,                   ionode_id )   
       CALL mp_bcast( occupations,            ionode_id )   
       CALL mp_bcast( smearing,               ionode_id )     
       CALL mp_bcast( degauss,                ionode_id )     
       CALL mp_bcast( nelup,                  ionode_id )
       CALL mp_bcast( neldw,                  ionode_id )
       CALL mp_bcast( nspin,                  ionode_id )
       CALL mp_bcast( nosym,                  ionode_id ) 
       CALL mp_bcast( ecfixed,                ionode_id )
       CALL mp_bcast( qcutz,                  ionode_id )
       CALL mp_bcast( q2sigma,                ionode_id )
       CALL mp_bcast( xc_type,                ionode_id )
       CALL mp_bcast( input_dft,              ionode_id )
#ifdef EXX
       CALL mp_bcast( nqx1,                   ionode_id )
       CALL mp_bcast( nqx2,                   ionode_id )
       CALL mp_bcast( nqx3,                   ionode_id )
#endif
       CALL mp_bcast( starting_magnetization, ionode_id ) 
       CALL mp_bcast( starting_ns_eigenvalue, ionode_id ) 
       CALL mp_bcast( U_projection_type,      ionode_id ) 
       CALL mp_bcast( lda_plus_U,             ionode_id ) 
       CALL mp_bcast( Hubbard_U,              ionode_id ) 
       CALL mp_bcast( Hubbard_alpha,          ionode_id ) 
       CALL mp_bcast( edir,                   ionode_id )
       CALL mp_bcast( emaxpos,                ionode_id )
       CALL mp_bcast( eopreg,                 ionode_id )
       CALL mp_bcast( eamp,                   ionode_id )
       CALL mp_bcast( la2F,                   ionode_id )
       !
       ! ... non collinear broadcast
       !
       CALL mp_bcast( lspinorb,                  ionode_id )
       CALL mp_bcast( noncolin,                  ionode_id )
       CALL mp_bcast( angle1,                    ionode_id )
       CALL mp_bcast( angle2,                    ionode_id )
       CALL mp_bcast( report,                    ionode_id )
       CALL mp_bcast( constrained_magnetization, ionode_id )
       CALL mp_bcast( B_field,                   ionode_id )
       CALL mp_bcast( fixed_magnetization,       ionode_id )
       CALL mp_bcast( lambda,                    ionode_id )       
       ! 
       RETURN
       !
     END SUBROUTINE
     !
     !=----------------------------------------------------------------------=!
     !
     !  Broadcast variables values for Namelist ELECTRONS
     !
     !=----------------------------------------------------------------------=!
     !
     !----------------------------------------------------------------------
     SUBROUTINE electrons_bcast()
       !----------------------------------------------------------------------
       !
       USE io_global, ONLY : ionode_id
       USE mp,        ONLY : mp_bcast
       !
       IMPLICIT NONE
       !
       !
       CALL mp_bcast( emass,                ionode_id )
       CALL mp_bcast( emass_cutoff,         ionode_id )
       CALL mp_bcast( orthogonalization,    ionode_id )
       CALL mp_bcast( electron_maxstep,     ionode_id )
       CALL mp_bcast( ortho_eps,            ionode_id )
       CALL mp_bcast( ortho_max,            ionode_id )
       CALL mp_bcast( electron_dynamics,    ionode_id )
       CALL mp_bcast( electron_damping,     ionode_id )
       CALL mp_bcast( electron_velocities,  ionode_id )
       CALL mp_bcast( electron_temperature, ionode_id )
       CALL mp_bcast( conv_thr,             ionode_id )
       CALL mp_bcast( ekincw,               ionode_id )
       CALL mp_bcast( fnosee,               ionode_id )
       CALL mp_bcast( startingwfc,          ionode_id )
       CALL mp_bcast( ampre,                ionode_id )
       CALL mp_bcast( grease,               ionode_id )
       CALL mp_bcast( startingpot,          ionode_id )
       CALL mp_bcast( empty_states_nbnd,    ionode_id )
       CALL mp_bcast( empty_states_maxstep, ionode_id )
       CALL mp_bcast( empty_states_delt,    ionode_id )
       CALL mp_bcast( empty_states_emass,   ionode_id )
       CALL mp_bcast( empty_states_ethr,    ionode_id )
       CALL mp_bcast( diis_size,            ionode_id )
       CALL mp_bcast( diis_nreset,          ionode_id )
       CALL mp_bcast( diis_hcut,            ionode_id )
       CALL mp_bcast( diis_wthr,            ionode_id )
       CALL mp_bcast( diis_delt,            ionode_id )
       CALL mp_bcast( diis_maxstep,         ionode_id )
       CALL mp_bcast( diis_rot,             ionode_id )
       CALL mp_bcast( diis_fthr,            ionode_id )
       CALL mp_bcast( diis_temp,            ionode_id )
       CALL mp_bcast( diis_achmix,          ionode_id )
       CALL mp_bcast( diis_g0chmix,         ionode_id )
       CALL mp_bcast( diis_g1chmix,         ionode_id )
       CALL mp_bcast( diis_nchmix,          ionode_id )
       CALL mp_bcast( diis_nrot,            ionode_id )
       CALL mp_bcast( diis_rothr,           ionode_id )
       CALL mp_bcast( diis_ethr,            ionode_id )
       CALL mp_bcast( diis_chguess,         ionode_id )
       CALL mp_bcast( mixing_fixed_ns,      ionode_id )
       CALL mp_bcast( mixing_mode,          ionode_id )
       CALL mp_bcast( mixing_beta,          ionode_id )
       CALL mp_bcast( mixing_ndim,          ionode_id )
       CALL mp_bcast( diagonalization,      ionode_id )
       CALL mp_bcast( diago_thr_init,       ionode_id )
       CALL mp_bcast( diago_cg_maxiter,     ionode_id )
       CALL mp_bcast( diago_david_ndim,     ionode_id )
       CALL mp_bcast( diago_diis_ndim,      ionode_id )
       CALL mp_bcast( sic,                  ionode_id )
       CALL mp_bcast( sic_epsilon ,         ionode_id )
       CALL mp_bcast( sic_alpha   ,         ionode_id )
       CALL mp_bcast( force_pairing ,       ionode_id )
       ! 
       ! ... ensemble-DFT
       !
       CALL mp_bcast( fermi_energy,       ionode_id )
       CALL mp_bcast( n_inner,            ionode_id )
       CALL mp_bcast( rotation_dynamics,  ionode_id )
       CALL mp_bcast( occupation_dynamics,ionode_id )
       CALL mp_bcast( rotmass,            ionode_id )
       CALL mp_bcast( occmass,            ionode_id )
       CALL mp_bcast( rotation_damping,   ionode_id )
       CALL mp_bcast( occupation_damping, ionode_id )
       !
       ! ... conjugate gradient
       !
       CALL mp_bcast( tcg,     ionode_id )
       CALL mp_bcast( maxiter, ionode_id )
       CALL mp_bcast( etresh,  ionode_id )
       CALL mp_bcast( passop,  ionode_id )
       !
       ! ... electric field
       !
       CALL mp_bcast( epol,   ionode_id )
       CALL mp_bcast( efield, ionode_id )
       ! 
       RETURN
       !
     END SUBROUTINE
     !
     !=----------------------------------------------------------------------=!
     !
     !  Broadcast variables values for Namelist IONS
     !
     !=----------------------------------------------------------------------=!
     !
     !----------------------------------------------------------------------
     SUBROUTINE ions_bcast()
       !----------------------------------------------------------------------
       !
       USE io_global, ONLY: ionode_id
       USE mp,        ONLY: mp_bcast
       !
       IMPLICIT NONE
       !
       !
       CALL mp_bcast( phase_space,       ionode_id )
       CALL mp_bcast( ion_dynamics,      ionode_id )
       CALL mp_bcast( ion_radius,        ionode_id )
       CALL mp_bcast( ion_damping,       ionode_id )
       CALL mp_bcast( ion_positions,     ionode_id )
       CALL mp_bcast( ion_velocities,    ionode_id )
       CALL mp_bcast( ion_temperature,   ionode_id )
       CALL mp_bcast( tempw,             ionode_id )
       CALL mp_bcast( fnosep,            ionode_id )
       CALL mp_bcast( nhpcl,             ionode_id )
       CALL mp_bcast( nhptyp,            ionode_id )
       CALL mp_bcast( ndega,             ionode_id )
       CALL mp_bcast( tranp,             ionode_id )
       CALL mp_bcast( amprp,             ionode_id )
       CALL mp_bcast( greasp,            ionode_id )
       CALL mp_bcast( tolp,              ionode_id )
       CALL mp_bcast( ion_nstepe,        ionode_id )
       CALL mp_bcast( ion_maxstep,       ionode_id )
       CALL mp_bcast( delta_t,           ionode_id )
       CALL mp_bcast( nraise,            ionode_id )
       CALL mp_bcast( upscale,           ionode_id )
       CALL mp_bcast( pot_extrapolation, ionode_id )
       CALL mp_bcast( wfc_extrapolation, ionode_id )
       !
       ! ... "path" variables broadcast
       !
       CALL mp_bcast( num_of_images,      ionode_id )
       CALL mp_bcast( first_last_opt,     ionode_id )
       CALL mp_bcast( use_masses,         ionode_id )
       CALL mp_bcast( init_num_of_images, ionode_id )
       CALL mp_bcast( use_multistep,      ionode_id )
       CALL mp_bcast( use_fourier,        ionode_id )
       CALL mp_bcast( use_freezing,       ionode_id )
       CALL mp_bcast( fixed_tan,          ionode_id )
       CALL mp_bcast( free_energy,        ionode_id )
       CALL mp_bcast( write_save,         ionode_id )
       CALL mp_bcast( CI_scheme,          ionode_id )
       CALL mp_bcast( opt_scheme,         ionode_id )
       CALL mp_bcast( damp,               ionode_id )
       CALL mp_bcast( temp_req,           ionode_id )
       CALL mp_bcast( ds,                 ionode_id )
       CALL mp_bcast( k_max,              ionode_id )
       CALL mp_bcast( k_min,              ionode_id )
       CALL mp_bcast( path_thr,           ionode_id )
       !
       ! ... BFGS
       !
       CALL mp_bcast( bfgs_ndim,        ionode_id )
       CALL mp_bcast( trust_radius_max, ionode_id )
       CALL mp_bcast( trust_radius_min, ionode_id )
       CALL mp_bcast( trust_radius_ini, ionode_id )
       CALL mp_bcast( w_1,              ionode_id )
       CALL mp_bcast( w_2,              ionode_id )
       !
       CALL mp_bcast( sic_rloc, ionode_id )
       ! 
       ! ... SMD broadcast (Y.K. 15/04/2004)
       !
       CALL mp_bcast( smd_polm,    ionode_id )
       CALL mp_bcast( smd_kwnp,    ionode_id )
       CALL mp_bcast( smd_linr,    ionode_id )
       CALL mp_bcast( smd_stcd,    ionode_id )
       CALL mp_bcast( smd_stcd2,   ionode_id )
       CALL mp_bcast( smd_stcd2,   ionode_id )
       CALL mp_bcast( smd_stcd3,   ionode_id )
       CALL mp_bcast( smd_codf,    ionode_id )
       CALL mp_bcast( smd_forf,    ionode_id )
       CALL mp_bcast( smd_smwf,    ionode_id )
       CALL mp_bcast( smd_lmfreq,  ionode_id )
       CALL mp_bcast( smd_tol,     ionode_id )
       CALL mp_bcast( smd_maxlm,   ionode_id )
       CALL mp_bcast( smd_smcp,    ionode_id )
       CALL mp_bcast( smd_smopt,   ionode_id )
       CALL mp_bcast( smd_smlm,    ionode_id )
       CALL mp_bcast( smd_ene_ini, ionode_id )
       CALL mp_bcast( smd_ene_fin, ionode_id )
       !
       CALL mp_bcast( fe_step,     ionode_id )
       CALL mp_bcast( fe_nstep,    ionode_id )
       CALL mp_bcast( shake_nstep, ionode_id )
       CALL mp_bcast( g_amplitude, ionode_id )
       !
       RETURN
       !
     END SUBROUTINE
     !
     !=----------------------------------------------------------------------=!
     !
     !  Broadcast variables values for Namelist CELL
     !
     !=----------------------------------------------------------------------=!
     !
     !----------------------------------------------------------------------
     SUBROUTINE cell_bcast()
       !----------------------------------------------------------------------
       !
       USE io_global, ONLY: ionode_id
       USE mp, ONLY: mp_bcast
       !
       IMPLICIT NONE
       !
       !
       CALL mp_bcast( cell_parameters,  ionode_id )
       CALL mp_bcast( cell_dynamics,    ionode_id )
       CALL mp_bcast( cell_velocities,  ionode_id )
       CALL mp_bcast( cell_dofree,      ionode_id )
       CALL mp_bcast( press,            ionode_id )
       CALL mp_bcast( wmass,            ionode_id )
       CALL mp_bcast( cell_temperature, ionode_id )
       CALL mp_bcast( temph,            ionode_id )
       CALL mp_bcast( fnoseh,           ionode_id )
       CALL mp_bcast( greash,           ionode_id )
       CALL mp_bcast( cell_factor,      ionode_id )
       CALL mp_bcast( cell_nstepe,      ionode_id )
       CALL mp_bcast( cell_damping,     ionode_id )
       ! 
       RETURN
       !
     END SUBROUTINE
     !
     !=----------------------------------------------------------------------------=!
     !
     !  Broadcast variables values for Namelist PHONON
     !
     !=----------------------------------------------------------------------------=!
     !
     !----------------------------------------------------------------------
     SUBROUTINE phonon_bcast()
       !----------------------------------------------------------------------
       !
       USE io_global, ONLY: ionode_id
       USE mp,        ONLY: mp_bcast
       !
       IMPLICIT NONE
       !
       !
       CALL mp_bcast( modenum, ionode_id )
       CALL mp_bcast( xqq,     ionode_id )
       !
       RETURN
       !
     END SUBROUTINE     
     !
     !=----------------------------------------------------------------------------=!
     !
     !  Broadcast variables values for Namelist RAMAN
     !
     !=----------------------------------------------------------------------------=!
     !
     !----------------------------------------------------------------------
     SUBROUTINE raman_bcast()
       !----------------------------------------------------------------------
       !
       USE io_global, ONLY: ionode_id
       USE mp,        ONLY: mp_bcast
       !
       IMPLICIT NONE
       !
       !
       CALL mp_bcast( b_length, ionode_id )
       CALL mp_bcast( lcart,    ionode_id )
       !
       RETURN
       !
     END SUBROUTINE
     !
     !
     !=----------------------------------------------------------------------------=!
     !
     !  Broadcast variables values for Namelist WANNIER
     !
     !=----------------------------------------------------------------------------=!
     !
     !----------------------------------------------------------------------
     SUBROUTINE wannier_bcast()
       !----------------------------------------------------------------------
       !
       USE io_global, ONLY: ionode_id
       USE mp,        ONLY: mp_bcast
       !
       IMPLICIT NONE
       !
       !
       CALL mp_bcast( wf_efield,   ionode_id )
       CALL mp_bcast( wf_switch,   ionode_id )
       CALL mp_bcast( sw_len,      ionode_id )
       CALL mp_bcast( efx0,        ionode_id )
       CALL mp_bcast( efy0,        ionode_id )
       CALL mp_bcast( efz0,        ionode_id )
       CALL mp_bcast( efx1,        ionode_id )
       CALL mp_bcast( efy1,        ionode_id )
       CALL mp_bcast( efz1,        ionode_id )
       CALL mp_bcast( wfsd,        ionode_id )
       CALL mp_bcast( wfdt,        ionode_id )
       CALL mp_bcast( maxwfdt,     ionode_id )
       CALL mp_bcast( wf_q,        ionode_id )
       CALL mp_bcast( wf_friction, ionode_id )
       CALL mp_bcast( nit,         ionode_id )
       CALL mp_bcast( nsd,         ionode_id )
       CALL mp_bcast( nsteps,      ionode_id )
       CALL mp_bcast( tolw,        ionode_id )
       CALL mp_bcast( adapt,       ionode_id )
       CALL mp_bcast( calwf,       ionode_id )
       CALL mp_bcast( nwf,         ionode_id )
       CALL mp_bcast( wffort,      ionode_id )
       CALL mp_bcast( writev,      ionode_id )
       !
       RETURN
       !
     END SUBROUTINE 
     !
     !=----------------------------------------------------------------------=!
     !
     !  Check input values for Namelist CONTROL
     !
     !=----------------------------------------------------------------------=!
     !
     !----------------------------------------------------------------------
     SUBROUTINE control_checkin( prog )
       !----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       CHARACTER(LEN=2)  :: prog   ! ... specify the calling program
       CHARACTER(LEN=80) :: msg
       CHARACTER(LEN=20) :: sub_name = ' control_checkin '
       INTEGER           :: i
       LOGICAL           :: allowed = .FALSE.
       ! 
       !
       DO i = 1, SIZE( calculation_allowed )
          IF( TRIM(calculation) == calculation_allowed(i) ) allowed = .TRUE.
       END DO
       IF( .NOT. allowed ) &
          CALL errore( sub_name, ' calculation '''// &
                       & TRIM(calculation)//''' not allowed ',1)
       IF( calculation == ' ' ) &
          CALL errore( sub_name,' calculation not specified ',1)
       IF( prog == 'CP' ) THEN
          IF( calculation == 'nscf' .OR. calculation == 'phonon' ) &
             CALL errore( sub_name,' calculation '//TRIM(calculation)// &
                          & ' not implemented ',1)
       END IF
       IF( ndr < 50 ) &
          CALL errore( sub_name,' ndr out of range ', 1 )
       IF( ndw > 0 .AND. ndw < 50 ) &
          CALL errore( sub_name,' ndw out of range ', 1 )
       IF( nstep < 0 ) &
          CALL errore( sub_name,' nstep out of range ', 1 )
       IF( iprint < 1 ) &
          CALL errore( sub_name,' iprint out of range ', 1 )

       IF( prog == 'PW' ) THEN
         IF( isave > 0 ) &
           CALL infomsg( sub_name,' isave not used in PW ', 1 )
       ELSE
         IF( isave < 1 ) &
           CALL errore( sub_name,' isave out of range ', 1 )
       END IF
 
       IF( dt < 0.D0 ) &
          CALL errore( sub_name,' dt out of range ', 1 )
       IF( max_seconds < 0.D0 ) &
          CALL errore( sub_name,' max_seconds out of range ', 1 )

       IF( ekin_conv_thr < 0.D0 ) THEN
          IF( prog == 'PW' ) THEN
            CALL infomsg( sub_name,' ekin_conv_thr not used in PW ', 1 )
          ELSE 
            CALL errore( sub_name,' ekin_conv_thr out of range ', 1 )
          END IF
       END IF

       IF( etot_conv_thr < 0.D0 ) &
          CALL errore( sub_name,' etot_conv_thr out of range ', 1 )
       IF( forc_conv_thr < 0.D0 ) &
          CALL errore( sub_name,' force_conv_thr out of range ', 1 )
       IF( prog == 'CP' ) THEN
          IF( dipfield ) & 
             CALL infomsg( sub_name,' dipfield not implemented yet ', -1)
          IF( lberry ) & 
             CALL infomsg( sub_name,' lberry not implemented yet ', -1)
          IF( disk_io /= 'default' ) &
             CALL infomsg( sub_name,' disk_io not used ', -1)
          IF( gdir /= 0 ) &
             CALL infomsg( sub_name,' gdir not used ', -1)
          IF( nppstr /= 0 ) &
             CALL infomsg( sub_name,' nppstr not used ', -1)
       END IF
       !
       IF( prog == 'PW' .AND. TRIM( restart_mode ) == 'reset_counters' ) THEN
         CALL infomsg ( sub_name, &
                    & ' restart_mode == reset_counters' // &
                    & ' not implemented in PW ', -1 )
       END IF
       !
       IF( refg < 0 ) &
         CALL errore( sub_name, ' wrong table interval refg ', 1 )
       !
       RETURN
       !
     END SUBROUTINE
     !
     !=----------------------------------------------------------------------=!
     !
     !  Check input values for Namelist SYSTEM
     !
     !=----------------------------------------------------------------------=!
     !
     !----------------------------------------------------------------------
     SUBROUTINE system_checkin( prog )
       !---------------------------------------------------------------------- 
       !
       IMPLICIT NONE
       !
       CHARACTER(LEN=2)  :: prog   ! ... specify the calling program
       CHARACTER(LEN=80) :: msg
       CHARACTER(LEN=20) :: sub_name = ' system_checkin '
       !
       !
       IF( ibrav < 0 .OR. ibrav > 14 ) &
          CALL errore( sub_name ,' ibrav out of range ', MAX( 1, ibrav) )
       !
       IF( ( ibrav /= 0 ) .AND. ( celldm(1) == 0.d0 ) .AND. ( a == 0.d0 ) ) &
           CALL errore( ' iosys ', &
                      & ' invalid lattice parameters ( celldm or a )', 1 )
       !
       IF( nat < 1 ) &
          CALL errore( sub_name ,' nat less than one ', MAX( nat, 1) )
       IF( nat > natx ) &
          CALL errore( sub_name , &
                       & ' nat too large, increase NATX ', MAX( nat, 1) )
       !
       IF( ntyp < 1 ) &
          CALL errore( sub_name ,' ntyp less than one ', MAX( ntyp, 1) )
       IF( ntyp < 1 .OR. ntyp > nsx ) &
          CALL errore( sub_name , &
                       & ' ntyp too large, increase NSX ', MAX( ntyp, 1) )
       !
       IF( nspin < 1 .OR. nspin > nspinx ) &
          CALL errore( sub_name ,' nspin out of range ', MAX(nspin, 1 ) )
       !
       IF( ecutwfc <= 0.D0 ) &
          CALL errore( sub_name ,' ecutwfc out of range ',1)
       IF( ecutrho < 0.D0 ) &
          CALL errore( sub_name ,' ecutrho out of range ',1)
       !
       IF( prog == 'CP' ) THEN
          IF( degauss /= 0.D0 ) &
             CALL infomsg( sub_name ,' degauss is not used in CP ', -1)
       END IF
       !
       IF( nelup < 0.d0 .OR. nelup > nelec ) &
          CALL errore( sub_name ,' nelup out of range ',1)
       IF( neldw < 0.d0 .OR. neldw > nelec ) &
          CALL errore( sub_name ,' neldw out of range ',1)
       IF( ecfixed < 0.D0 ) &
          CALL errore( sub_name ,' ecfixed out of range ',1)
       IF( qcutz < 0.D0 ) &
          CALL errore( sub_name ,' qcutz out of range ',1)
       IF( q2sigma < 0.D0 ) &
          CALL errore( sub_name ,' q2sigma out of range ',1)
       IF( prog == 'CP' ) THEN
          IF( ANY(starting_magnetization /= SM_NOT_SET ) ) &
             CALL infomsg( sub_name ,&
                          & ' starting_magnetization is not used in CP ', -1)
          IF( lda_plus_U ) &
             CALL infomsg( sub_name ,' lda_plus_U is not used in CP ', -1)
          IF( la2F ) &
             CALL infomsg( sub_name ,' la2F is not used in CP ', -1)
          IF( ANY(Hubbard_U /= 0.D0) ) &
             CALL infomsg( sub_name ,' Hubbard_U is not used in CP ', -1)
          IF( ANY(Hubbard_alpha /= 0.D0) ) &
             CALL infomsg( sub_name ,' Hubbard_alpha is not used in CP ', -1)
          IF( nosym ) &
             CALL infomsg( sub_name ,' nosym not implemented in CP ', -1)
       END IF
       !
       IF( prog == 'PW' ) THEN
          !
          ! ... stop if starting_magnetization is not set for
          ! ... all atomic types
          !
          IF ( (nspin==2) .AND. ALL(starting_magnetization == sm_not_set) ) &
            CALL errore(sub_name,'some starting_magnetization MUST be set', 1 )
          !
       END IF
       !
       ! ... non collinear check
       !
       IF ( noncolin ) THEN
          !
          IF ( diagonalization == 'cg' ) &
             CALL errore( sub_name ,' cg not allowed with noncolin ', 1 )
          !
          IF ( diagonalization == 'diis' ) &
             CALL errore( sub_name ,' diis not allowed with noncolin ', 1 )
          !
       END IF       
       !
       ! ... control on SIC variables
       !
       IF ( sic /= 'none' ) THEN
          !
          IF (sic_epsilon > 1.D0 )  &
             CALL errore( sub_name, &
                        & ' invalid sic_epsilon, greater than 1.',1 )
          IF (sic_epsilon < 0.d0 )  &
             CALL errore( sub_name, &
                        & ' invalid sic_epsilon, less than 0 ',1 )
          IF (sic_alpha > 1.D0 )  &
             CALL errore( sub_name, &
                        & ' invalid sic_alpha, greater than 1.',1 )
          IF (sic_alpha < 0.d0 )  &
             CALL errore( sub_name, &
                        & ' invalid sic_alpha, less than 0 ',1 )
          !
          IF ( .NOT. force_pairing ) &
             CALL errore( sub_name, &
                        & ' invalid force_pairing with sic activated', 1 )
          IF ( nspin /= 2 ) &
             CALL errore( sub_name, &
                        & ' invalid nspin with sic activated', 1 )
          IF ( ( nelup == 0 ) .AND. ( neldw == 0 ) ) &
             CALL errore( sub_name, &
                      & ' invalid nelup and neldwn spin with sic activated', 1 )
          IF ( nelup /= (neldw + 1) )   &
             CALL errore( sub_name, &
                  & ' invalid nelup /= (neldwn +1) spin with sic activated', 1 )
          !  
       ENDIF
       !
       RETURN
       !
     END SUBROUTINE
     !
     !=----------------------------------------------------------------------=!
     !
     !  Check input values for Namelist ELECTRONS
     !
     !=----------------------------------------------------------------------=!
     !
     !----------------------------------------------------------------------
     SUBROUTINE electrons_checkin( prog )
       !----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       CHARACTER(LEN=2)  :: prog   ! ... specify the calling program
       CHARACTER(LEN=80) :: msg
       CHARACTER(LEN=20) :: sub_name = ' electrons_checkin '
       INTEGER           :: i
       LOGICAL           :: allowed = .FALSE.
       !
       !
       DO i = 1, SIZE(electron_dynamics_allowed)
          IF( TRIM(electron_dynamics) == &
              electron_dynamics_allowed(i) ) allowed = .TRUE.
       END DO
       IF( .NOT. allowed ) &
          CALL errore( sub_name, ' electron_dynamics '''//&
                       & TRIM(electron_dynamics)//''' not allowed ',1)
       IF( emass <= 0.D0 ) &
          CALL errore( sub_name, ' emass less or equal 0 ',1)
       IF( emass_cutoff <= 0.D0 ) &
          CALL errore( sub_name, ' emass_cutoff less or equal 0 ',1)
       IF( ortho_eps <= 0.D0 ) &
          CALL errore( sub_name, ' ortho_eps less or equal 0 ',1)
       IF( ortho_max < 1 ) &
          CALL errore( sub_name, ' ortho_max less than 1 ',1)
       IF( fnosee <= 0.D0 ) &
          CALL errore( sub_name, ' fnosee less or equal 0 ',1)
       IF( ekincw <= 0.D0 ) &
          CALL errore( sub_name, ' ekincw less or equal 0 ',1)
       IF( empty_states_nbnd < 0 ) &
          CALL errore( sub_name, &
                       & ' invalid empty_states_nbnd, less than 0 ',1)
       IF( empty_states_maxstep < 0 ) &
          CALL errore( sub_name,&
                       & ' invalid empty_states_maxstep, less than 0 ',1)
       IF( empty_states_delt < 0.D0 ) &
          CALL errore( sub_name, &
                       & ' invalid empty_states_delt, less than 0 ',1)
       IF( empty_states_emass < 0.D0 ) &
          CALL errore( sub_name, &
                       & ' invalid empty_states_emass, less than 0 ',1)
       IF( empty_states_ethr < 0.D0 ) &
          CALL errore( sub_name, &
                       & ' invalid empty_states_ethr, less than 0 ',1)

!
       RETURN
     END SUBROUTINE
     !
     !=----------------------------------------------------------------------=!
     !
     !  Check input values for Namelist IONS
     !
     !=----------------------------------------------------------------------=!
     !
     !----------------------------------------------------------------------
     SUBROUTINE ions_checkin( prog )
       !----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       CHARACTER(LEN=2)  :: prog   ! ... specify the calling program
       CHARACTER(LEN=80) :: msg
       CHARACTER(LEN=20) :: sub_name = ' ions_checkin '
       INTEGER           :: i
       LOGICAL           :: allowed = .FALSE.
       !
       !
       DO i = 1, SIZE( phase_space_allowed )
          IF( TRIM( phase_space ) == phase_space_allowed(i) ) allowed = .TRUE.
       END DO
       IF ( .NOT. allowed ) &
          CALL errore( sub_name, ' phase_space '''// &
                       & TRIM( phase_space )// ''' not allowed ', 1 )
       !
       DO i = 1, SIZE(ion_dynamics_allowed)
          IF( TRIM(ion_dynamics) == ion_dynamics_allowed(i) ) allowed = .TRUE.
       END DO
       IF( .NOT. allowed ) &
          CALL errore( sub_name, ' ion_dynamics '''// &
                       & TRIM(ion_dynamics)//''' not allowed ',1)
       IF( tempw <= 0.D0 ) &
          CALL errore( sub_name,' tempw out of range ',1)
       IF( fnosep( 1 ) <= 0.D0 ) &
          CALL errore( sub_name,' fnosep out of range ',1)
       IF( nhpcl > nhclm ) &
          CALL infomsg ( sub_name,' nhpcl should be less than nhclm', -1)
       IF( nhpcl < 0 ) &
          CALL infomsg ( sub_name,' nhpcl out of range ', -1)
       IF( ion_nstepe <= 0 ) &
          CALL errore( sub_name,' ion_nstepe out of range ',1)
       IF( ion_maxstep < 0 ) &
          CALL errore( sub_name,' ion_maxstep out of range ',1)
       !
       ! ... general "path" variables checkin
       !
       IF ( ds < 0.D0 ) &
          CALL errore( sub_name,' ds out of range ',1)
       IF ( damp < 0.D0 ) &
          CALL errore( sub_name,' damp out of range ',1)
       IF ( temp_req < 0.D0 ) &
          CALL errore( sub_name,' temp_req out of range ',1)
       !
       DO i = 1, SIZE( opt_scheme_allowed )
          IF ( TRIM( opt_scheme ) == &
               opt_scheme_allowed(i) ) allowed = .TRUE.
       END DO
       IF ( .NOT. allowed ) &
          CALL errore( sub_name, ' opt_scheme '''// &
                     & TRIM( opt_scheme )//''' not allowed ', 1 )
       !
       IF ( calculation == 'neb' .OR. calculation == 'smd' ) THEN
          !
          IF ( phase_space == 'coarse-grained' ) THEN
             !
             full_phs_path_flag = .FALSE.
             cg_phs_path_flag   = .TRUE.
             !
             IF ( calculation /= 'neb' .AND. calculation /= 'smd' ) &
                CALL errore( sub_name, &
                           & ' coarse-grained phase-space is presently' // &
                           & ' allowed only for neb or smd ', 1 )        
             !
          ELSE
             !
             full_phs_path_flag = .TRUE.
             cg_phs_path_flag   = .FALSE.
             !
          END IF
          !
       END IF
       !
       ! ... NEB specific checkin
       !
       IF ( num_of_images > max_num_of_images ) &
          CALL errore( sub_name, &
                     & ' num_of_images too large (increase max_num_of_images)', 1)
       IF ( k_max < 0.D0 ) &
          CALL errore( sub_name, ' k_max out of range', 1 )
       IF ( k_min < 0.D0 ) &
          CALL errore( sub_name, ' k_min out of range', 1 )
       IF ( k_max < k_min ) &
          CALL errore( sub_name, ' k_max < k_min', 1 )
       ! 
       DO i = 1, SIZE( CI_scheme_allowed )
          IF ( TRIM( CI_scheme ) == CI_scheme_allowed(i) ) allowed = .TRUE.
       END DO
       IF ( .NOT. allowed ) &
          CALL errore( sub_name, ' CI_scheme '''// &
                       & TRIM( CI_scheme )//''' not allowed ', 1 )
       !
       ! ... SMD checking ( Y.K. 15/04/2004 )
       !
       IF ( smd_polm .AND. smd_linr ) &
          CALL infomsg( sub_name,' smd_polm & smd_linr  both true ', -1)
       !
       IF ( smd_polm .AND. smd_kwnp < 3 ) &
          CALL infomsg( sub_name,' smd_kwnp < 3 for smd_polm ', -1)
       !
       IF ( smd_stcd .AND. (smd_stcd1==0 .OR. smd_stcd2==0 .OR. smd_stcd3==0) ) & 
          CALL infomsg( sub_name,' smd_stcd not specified ', -1)
       !
       IF( smd_smcp .AND. (smd_smopt .OR. smd_smlm) ) &
          CALL infomsg( sub_name,' smcp ? ', -1)
       !
       IF( smd_smopt .AND. (smd_smcp .OR. smd_smlm) ) &
          CALL infomsg( sub_name,' smopt ? ', -1)
       !
       IF( smd_smlm .AND. (smd_smcp .OR. smd_smopt) ) &
          CALL infomsg( sub_name,' smlm ? ', -1)
       !
       IF (sic /= 'none' .and. sic_rloc == 0.d0) &
          CALL errore( sub_name, ' invalid sic_rloc with sic activated ', 1 )
       !
       RETURN
       !
     END SUBROUTINE
     !
     !=----------------------------------------------------------------------=!
     !
     !  Check input values for Namelist CELL
     !
     !=----------------------------------------------------------------------=!
     !
     !----------------------------------------------------------------------
     SUBROUTINE cell_checkin( prog )
       !----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       CHARACTER(LEN=2)  :: prog   ! ... specify the calling program
       CHARACTER(LEN=80) :: msg
       CHARACTER(LEN=20) :: sub_name = ' cell_checkin '
       INTEGER           :: i
       LOGICAL           :: allowed = .FALSE.
       !
       !
       DO i = 1, SIZE(cell_dynamics_allowed)
          IF( TRIM(cell_dynamics) == &
              cell_dynamics_allowed(i) ) allowed = .TRUE.
       END DO
       IF( .NOT. allowed ) &
          CALL errore( sub_name, ' cell_dynamics '''// &
                       TRIM(cell_dynamics)//''' not allowed ',1)
       IF( wmass < 0.D0 ) &
          CALL errore( sub_name,' wmass out of range ',1)
       IF( prog == 'CP' ) THEN
          IF( cell_factor /= 0.D0 ) &
             CALL infomsg( sub_name,' cell_factor not used in CP ', -1)
       END IF
       IF( cell_nstepe <= 0 ) &
          CALL errore( sub_name,' cell_nstepe out of range ',1)
       !
       RETURN
       !
     END SUBROUTINE
     !
     !=----------------------------------------------------------------------=!
     !
     !  Check input values for Namelist PHONON
     !
     !=----------------------------------------------------------------------=!
     !
     !----------------------------------------------------------------------
     SUBROUTINE phonon_checkin( prog )
       !--------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       CHARACTER(LEN=2) :: prog   ! ... specify the calling program
       !
       !
       RETURN
       !
     END SUBROUTINE
     !
     !=----------------------------------------------------------------------=!
     !
     !  Check input values for Namelist RAMAN
     !
     !=----------------------------------------------------------------------=!
     !
     !----------------------------------------------------------------------
     SUBROUTINE raman_checkin( prog )
       !--------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       CHARACTER(LEN=2) :: prog   ! ... specify the calling program
       !
       !
       RETURN
       !
     END SUBROUTINE
     !
     !=----------------------------------------------------------------------=!
     !
     !  Check input values for Namelist WANNIER
     !
     !=----------------------------------------------------------------------=!
     !
     !----------------------------------------------------------------------
     SUBROUTINE wannier_checkin( prog )
       !--------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       CHARACTER(LEN=2)  :: prog   ! ... specify the calling program
       CHARACTER(LEN=20) :: sub_name = 'wannier_checkin'
       !
       !
       IF ( calwf < 1 .OR. calwf > 5 ) &
          CALL errore( sub_name, ' calwf out of range ', 1 )
       !
       RETURN
       !
     END SUBROUTINE
     !
     !=----------------------------------------------------------------------=!
     !
     !  Set values according to the "calculation" variable
     !
     !=----------------------------------------------------------------------=!
     !
     !----------------------------------------------------------------------
     SUBROUTINE fixval( prog )
       !---------------------------------------------------------------------- 
       !
       IMPLICIT NONE
       !
       CHARACTER(LEN=2)  :: prog   ! ... specify the calling program
       CHARACTER(LEN=80) :: msg
       CHARACTER(LEN=20) :: sub_name = ' fixval '
       !
       !
       IF( prog == 'PW' ) startingpot = 'atomic'
       !       
       SELECT CASE( TRIM( calculation ) )
          CASE ('scf')
             IF( prog == 'CP' ) THEN
                 electron_dynamics = 'damp'
                 ion_dynamics      = 'none'
                 cell_dynamics     = 'none'
             END IF
          CASE ('nscf')
             IF( prog == 'CP' ) &
                CALL errore( sub_name,' calculation '//TRIM(calculation)// &
                             & ' not implemented ',1)
             IF( prog == 'CP' ) occupations = 'bogus'
             IF( prog == 'CP' ) electron_dynamics = 'damp'
             IF( prog == 'PW' ) startingpot = 'file'
          CASE ('phonon')
             IF( prog == 'CP' ) &
                CALL errore( sub_name,' calculation '//TRIM(calculation)// &
                             & ' not implemented ',1)
             IF( prog == 'PW' ) startingpot = 'file'
          CASE ('raman')
             IF( prog == 'CP' ) &
                  CALL errore( sub_name,' calculation '//TRIM(calculation)// &
                  & ' not implemented ',1)
             IF( prog == 'PW' ) startingpot = 'file'
          CASE ( 'cp-wf' )
             IF( prog == 'CP' ) THEN
                electron_dynamics = 'damp'
                ion_dynamics      = 'damp'
             END IF
             IF ( prog == 'PW' ) &
                CALL errore( sub_name, ' calculation ' // &
                           & TRIM( calculation ) // ' not implemented ', 1 )
          CASE ('relax')
             IF( prog == 'CP' ) THEN
                electron_dynamics = 'damp'
                ion_dynamics      = 'damp'
             ELSE IF( prog == 'PW' ) THEN
                ion_dynamics = 'bfgs'
             END IF
          CASE ( 'md', 'cp' )
             IF( prog == 'CP' ) THEN
                electron_dynamics = 'verlet'
                ion_dynamics      = 'verlet'
             ELSE IF( prog == 'PW' ) THEN
                ion_dynamics = 'verlet'
             END IF
          CASE ('vc-relax')
             IF( prog == 'CP' ) THEN
                electron_dynamics = 'damp'
                ion_dynamics      = 'damp'
                cell_dynamics     = 'damp-pr'
             ELSE IF( prog == 'PW' ) THEN
                ion_dynamics = 'damp'
             END IF
          CASE ( 'vc-md', 'vc-cp' )
             IF( prog == 'CP' ) THEN
                electron_dynamics = 'verlet'
                ion_dynamics      = 'verlet'
                cell_dynamics     = 'pr'
             ELSE IF( prog == 'PW' ) THEN
                ion_dynamics = 'beeman'
             END IF
          CASE ( 'neb' )
             !
             ! ... "path" optimizations
             !
             IF( prog == 'CP' ) THEN
                !
                electron_dynamics = 'damp'
                ion_dynamics      = 'none'
                cell_dynamics     = 'none'
                !
             END IF
             !
          CASE ( 'smd' )
             !
             IF( prog == 'CP' ) THEN
                !
                electron_dynamics = 'damp'
                ion_dynamics      = 'damp'
                !
             END IF
             !
          CASE ( 'fpmd' )
             !
             !  Compatibility with old FPMD
             !
             IF ( prog == 'PW' ) &
                CALL errore( sub_name, ' calculation ' // &
                           & TRIM( calculation ) // ' not implemented ', 1 )
             !
             electron_dynamics = 'sd'
             ion_dynamics      = 'none'
             cell_dynamics     = 'none'
             !
          CASE( 'metadyn' )
             !
          CASE DEFAULT
             !
             CALL errore( sub_name,' calculation '// & 
                        & TRIM(calculation)//' not implemented ', 1 )
             !
       END SELECT
       !              
       IF ( prog == 'PW' ) THEN
          !       
          IF ( restart_mode == "from_scratch" ) THEN
             !
             startingwfc = 'atomic'
             startingpot = 'atomic'
             !
          ELSE
             !
             startingwfc = 'file'
             startingpot = 'file'             
             !
          END IF
          !
          IF ( calculation == 'nscf' .OR. &
               calculation == 'raman'.OR. &
               calculation == 'phonon' ) THEN
             !
             startingpot = 'file'
             startingwfc = 'atomic'
             !
          END IF   
          !      
       END IF                   
       !
       IF ( TRIM( sic ) /= 'none' ) THEN 
         IF ( nspin == 2 .AND. nelec > 1 .AND. &
              ( nelup == neldw .OR. nelup == neldw+1 ) ) force_pairing = .TRUE.
       END IF
       !
       IF ( calculation == 'metadyn' .AND. &
            prog == 'CP' ) g_amplitude = g_amplitude / e2
       !
       RETURN
       !
     END SUBROUTINE
     !
     !=----------------------------------------------------------------------=!
     !
     !  Namelist parsing main routine
     !
     !=----------------------------------------------------------------------=!
     !
     !----------------------------------------------------------------------
     SUBROUTINE read_namelists( prog )
       !----------------------------------------------------------------------
       !
       !  this routine reads data from standard input and puts them into
       !  module-scope variables (accessible from other routines by including
       !  this module, or the one that contains them)
       !  ----------------------------------------------
       !
       ! ... declare modules
       !
       USE mp_global, ONLY : mpime, nproc, group
       USE io_global, ONLY : ionode, ionode_id
       USE mp,        ONLY : mp_bcast
       !
       IMPLICIT NONE
       !
       ! ... declare variables
       !
       CHARACTER(LEN=2) :: prog   ! ... specify the calling program
                                  !     prog = 'PW'  pwscf
                                  !     prog = 'CP'  cpr
       !
       ! ... declare other variables
       ! 
       INTEGER :: ios
       !
       ! ... end of declarations
       !
       !  ----------------------------------------------
       !
       !
       IF( prog /= 'PW' .AND. prog /= 'CP' ) &
          CALL errore( ' read_namelists ', ' unknown calling program ', 1 )
       !
       ! ... default settings for all namelists
       !
       CALL control_defaults( prog )       
       CALL system_defaults( prog )
       CALL electrons_defaults( prog )
       CALL ions_defaults( prog )
       CALL cell_defaults( prog )
       CALL phonon_defaults( prog )
       !
       ! ... Here start reading standard input file
       !
       ! ... CONTROL namelist
       !
       ios = 0
       IF( ionode ) THEN
          READ( 5, control, iostat = ios ) 
       END IF
       CALL mp_bcast( ios, ionode_id )
       IF( ios /= 0 ) THEN
          CALL errore( ' read_namelists ', &
                     & ' reading namelist control ', ABS(ios) )
       END IF
       !
       CALL control_bcast( )
       CALL control_checkin( prog )
       !
       ! ... defaults values are changed here according to the CONTROL namelist
       !
       CALL fixval( prog )
       !
       ! ... SYSTEM namelist
       !       
       ios = 0
       IF( ionode ) THEN
          READ( 5, system, iostat = ios ) 
       END IF
       CALL mp_bcast( ios, ionode_id )
       IF( ios /= 0 ) THEN
          CALL errore( ' read_namelists ', &
                     & ' reading namelist system ', ABS(ios) )
       END IF
       !
       CALL system_bcast( )
       CALL system_checkin( prog )
       !
       ! ... ELECTRONS namelist
       !
       ios = 0
       IF( ionode ) THEN
          READ( 5, electrons, iostat = ios ) 
       END IF
       CALL mp_bcast( ios, ionode_id )
       IF( ios /= 0 ) THEN
          CALL errore( ' read_namelists ', &
                     & ' reading namelist electrons ', ABS(ios) )
       END IF
       !
       CALL electrons_bcast( )
       CALL electrons_checkin( prog )
       !
       ! ... IONS namelist   
       !
       ios = 0
       !
       IF ( ionode ) THEN
          !
          IF ( TRIM( calculation ) == 'relax'    .OR. &
               TRIM( calculation ) == 'md'       .OR. &
               TRIM( calculation ) == 'vc-relax' .OR. &
               TRIM( calculation ) == 'vc-md'    .OR. &
               TRIM( calculation ) == 'cp'       .OR. &
               TRIM( calculation ) == 'vc-cp'    .OR. &
               TRIM( calculation ) == 'smd'      .OR. &
               TRIM( calculation ) == 'cp-wf'    .OR. &
               TRIM( calculation ) == 'neb'      .OR. &
               TRIM( calculation ) == 'fpmd'     .OR. &
               TRIM( calculation ) == 'metadyn' ) READ( 5, ions, iostat = ios ) 
          !
       END IF
       CALL mp_bcast( ios, ionode_id )
       IF( ios /= 0 ) THEN
          CALL errore( ' read_namelists ', &
                     & ' reading namelist ions ', ABS(ios) )
       END IF
       !
       CALL ions_bcast( )
       CALL ions_checkin( prog )
       !
       ! ... CELL namelist
       !
       ios = 0
       IF( ionode ) THEN
          IF( TRIM( calculation ) == 'vc-relax' .OR. &
              TRIM( calculation ) == 'vc-cp'    .OR. &
              TRIM( calculation ) == 'vc-md'    .OR. &
              TRIM( calculation ) == 'fpmd'     .OR. &
              TRIM( calculation ) == 'vc-md' ) THEN
             READ( 5, cell, iostat = ios ) 
          END IF
       END IF
       CALL mp_bcast( ios, ionode_id )
       IF( ios /= 0 ) THEN
          CALL errore( ' read_namelists ', &
                     & ' reading namelist cell ', ABS(ios) )
       END IF
       !
       CALL cell_bcast()
       CALL cell_checkin( prog )
       !
       ! ... PHONON namelist
       !
       ios = 0
       IF( ionode ) THEN
          IF( TRIM( calculation ) == 'phonon' ) THEN
             READ( 5, phonon, iostat = ios )
          END IF
       END IF
       CALL mp_bcast( ios, ionode_id )
       IF( ios /= 0 ) THEN
          CALL errore( ' read_namelists ', &
                     & ' reading namelist cell ', ABS(ios) )
       END IF
       !
       CALL phonon_bcast()
       CALL phonon_checkin( prog )
       !
       ! ... RAMAN NAMELIST 
       !
       CALL raman_defaults( prog )
       ios = 0
       IF( ionode ) THEN
          IF( TRIM( calculation ) == 'raman' ) THEN
             READ( 5, raman, iostat = ios )
          END IF
       END IF
       CALL mp_bcast( ios, ionode_id )
       IF( ios /= 0 ) THEN
          CALL errore( ' read_namelists ', &
                     & ' reading namelist raman ', ABS(ios) )
       END IF
       !
       CALL raman_bcast()
       CALL raman_checkin( prog )
       !
       ! ... WANNIER NAMELIST
       !
       CALL wannier_defaults( prog )
       ios = 0
       IF( ionode ) THEN
          IF( TRIM( calculation ) == 'cp-wf' ) THEN
             READ( 5, wannier, iostat = ios )
          END IF
       END IF
       CALL mp_bcast( ios, ionode_id )
       IF( ios /= 0 ) THEN
          CALL errore( ' read_namelists ', &
                     & ' reading namelist wannier ', ABS(ios) )
       END IF
       !
       CALL wannier_bcast()
       CALL wannier_checkin( prog )
       !
       RETURN
       !
     END SUBROUTINE read_namelists
     !
END MODULE read_namelists_module
