!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!==============================================================================
!  Program VIB (programmed by Silviu Zilberman)
!  
!  Calculates the vibrational modes, Born effective charges and infrared
!  cross-section for isolated molecules/clusters in vacuum.
!  The programs uses rither the CP or PW codes as the underlying DFT engine
!  for wave function relaxation (SCF) and the frozen phonon method for
!  construction of the energy hessian (dynamical matrix),
!  i.e. displace each atom in all three cartesian directions to construct
!  the numerical first derivative of the forces.
!  
!  The program uses two mandatory input files
!  1. The regular CP input file, which presumably was used to generate a
!     ground state nuclear configuration (local minimum).
!  2. A file "prefix.vib.inp", where prefix is defined in the CP input
!     file.
!  
!   prefix.vib.inp  contains one name list:
!  
!   &INPUT_VIB
!   ...
!   /
!  
!  Parameters:
!  
!  displacement    real ( default = 0.05 )
!                  displacement step, a uniform value for all atoms
!  
!  save_freq       integer ( default = 1 )
!                  Save restart information every save_freq displacements.
!                  Since it is usually a small file, containing mainly
!                  the intermediate hessian, it is suggested not to
!                  change the default.
!  restart_label   character ( default = 'auto' )
!                  'from_scratch': starts a new calculation
!                  'restart'     : continue a previous calculation
!                                  a file 'prefix.vib.restart' must be present
!                  'auto'        : if a file 'prefix.vib.restart' is present
!                                  the program will pick it up and continue
!                                  the previous calculation, or if absent,
!                                  it will start a new calculation
!  
!  trans_inv_flag  logical ( default = .true. )
!                  impose translational invariance on the energy hessian
!                  to ensure that the modes associated with a rigid
!                  translation would have zero frequency
!  
!  trans_inv_max_iter      integer ( default = 50 )
!                  the maximal number of steps in the iterative procedure
!                  that removes the rigid translation modes
!  
!  trans_inv_conv_thr      real ( default = 1.0D-15 )
!                  the threshhold for convergence of the iterative procedure
!                  for removal of the rigid translation modes
!  
!  trans_rot_inv_flag      logical ( default = .true. )
!                  remove also the rigid rotation modes
!
!  animate         logical ( default = .false. )
!                  generate xyz animation files for all normal modes
!
!  Additional information about the output files:
!  
!  1. After the energy hessian is calculated, the program will write
!     the analysis results to 'prefix.vib.analysis'.
!     The minimal output is:
!          * raw (but symmetrized) hessian
!          * diagonal mass matrix
!          * harmonic frequencies and the associated effective
!            spring constant and effective mode mass
!          * raw Born charge tensors for each atom
!          * the sum over the Born charges (minus total system charge)
!            this quantity should ideally be zero (acoustic sum rule).
!          * infrared intensity of each normal mode.
!  
!     If trans_inv_flag==.true. then the same information as above is repeated
!     but this time with imposed translational invariance on the energy
!     hessian. In addition,  the Born charges are symmetrized and
!     the acoustic sum rule is imposed.
!  
!     if trans_rot_inv_flag==.true. then also rotational modes are projected
!     out.
!  
!  A brief description of the methods:
!  
!  1. The energy hessian is constructed by a simple difference formula for
!     the numerical first derivative of the forces.
!  
!  2. The condition for translational invariance of the hessian is described
!     (for instance) in Gonze and Lee, Phys. Rev. B 55(16), 1997, 10355-10368.
!     However, simple enforcement of translational invariance breaks the
!     symmetry of the hessian. Here symmetrization and translational
!     invariance are repeatedly applied until convergence, i.e., until
!     the largest change in any matrix element of the hessian is smaller
!     than trans_inv_conv_thr parameter.
!     A more sophisticated and general algorithm was suggested by
!     Nicolas Mounet and Nicola Marzari (Ref. ???) and is already implemented
!     in QE.
!  
!  3. The algorithm for removal of rotational rigid modes follows closely the
!     one implemented in Gaussian03 electronic structure package, as described
!     in their white papers (www.gaussian.com).
!
!
!  ==========================================================================
!  ====   units and constants                                            ====
!  ====                                                                  ====
!  ====   1 hartree           = 1 a.u.                                   ====
!  ====   1 bohr radius       = 1 a.u. = 0.529167 Angstrom               ====
!  ====   1 rydberg           = 1/2 a.u.                                 ====
!  ====   1 electron volt     = 1/27.212 a.u.                            ====
!  ====   1 kelvin *k-boltzm. = 1/(27.212*11606) a.u.'='3.2e-6 a.u       ====
!  ====   1 second            = 1/(2.4189)*1.e+17 a.u.                   ====
!  ====   1 proton mass       = 1822.89 a.u.                             ====
!  ====   1 tera              = 1.e+12                                   ====
!  ====   1 pico              = 1.e-12                                   ====
!  ====   1 Volt / meter      = 1/(5.1412*1.e+11) a.u.                   ====
!  ==========================================================================
!
!----------------------------------------------------------------------------
PROGRAM vib
  !----------------------------------------------------------------------------
  !
  USE control_flags,         ONLY : program_name
  USE io_global,             ONLY : stdout
  USE kinds,                 ONLY : DP
  USE vibrations,            ONLY : start_vibrations, calc_hessian, &
       analysis, end_vibrations
  !
  IMPLICIT NONE
  !
  INTEGER                        :: restart_cyc_counter
  REAL    (KIND=DP)              :: E_minus, dip_minus(3)
  !
  ! ... program starts here
  !
#ifdef DFT_CP
  program_name = 'CP90'
  call start_cp
#endif
#ifdef DFT_PW
  program_name = 'PW'
  call start_pw
#endif
  write (stdout,*) 'Using program ',program_name,' as a DFT engine.'
  !
  !
  print *,'reading vib input...'
  CALL read_input_vib()
  print *,'Done...'
  print *,'starting start_vibrations...'
  CALL start_vibrations(restart_cyc_counter, E_minus, dip_minus)
  print *,'Done...'
  CALL calc_hessian    (restart_cyc_counter, E_minus, dip_minus)
  CALL analysis ()
  CALL end_vibrations()
  !
#ifdef DFT_CP
     call end_cp
#endif
#ifdef DFT_PW
     call end_pw
#endif
  !
  STOP
  !
END PROGRAM vib

!----------------------------------------------------------------------------

subroutine start_cp
  !
#ifdef DFT_CP
  !
  USE environment,           ONLY : environment_start
  USE input,                 ONLY : read_input_file, iosys_pseudo, iosys
  USE io_global,             ONLY : io_global_start, io_global_getionode
  USE mp_global,             ONLY : mp_global_start
  USE mp,                    ONLY : mp_end, mp_start, mp_env
  USE check_stop,            ONLY : check_stop_init
  !
  IMPLICIT NONE
  !
  INTEGER                        :: mpime, nproc, gid, ionode_id
  INTEGER,           PARAMETER   :: root = 0
  LOGICAL                        :: ionode
  !
  ! ... initialize MPI (parallel processing handling)
  !
  CALL mp_start()
  CALL mp_env( nproc, mpime, gid )
  CALL mp_global_start( root, mpime, gid, nproc )
  !
  ! ... mpime = processor number, starting from 0
  ! ... nproc = number of processors
  ! ... gid   = group index
  ! ... root  = index of the root processor
  !
  ! ... initialize input output
  !
  CALL io_global_start( mpime, root )
  CALL io_global_getionode( ionode, ionode_id )
  !
  CALL environment_start( )
  !
  ! ... readin the input file
  !
  CALL read_input_file()
  !
  ! ... copy pseudopotential input parameter into internal variables
  ! ... and read in pseudopotentials and wavefunctions files
  !
  CALL iosys_pseudo()
  !
  ! ... copy-in input parameters from input_parameter module
  !
  CALL iosys()
  !
  CALL check_stop_init()
  !
  CALL init_run()
  !     
#endif
  return
end subroutine start_cp
!
!----------------------------------------------------------------------------
!
subroutine end_cp
  IMPLICIT NONE
  !
#ifdef DFT_CP
  !
  CALL terminate_run()
  !
  CALL stop_run( .TRUE. )
  !
#endif
  return
end subroutine end_cp
!
!----------------------------------------------------------------------------
!
subroutine start_pw
  !
#ifdef DFT_PW
  !
  ! ... Plane Wave Self-Consistent Field code
  !
  USE io_global,          ONLY : stdout
  USE parameters,         ONLY : ntypx, npk, lmaxx, nchix, ndmx, nqfx, nbrx
  USE global_version,     ONLY : version_number
  USE wvfct,              ONLY : gamma_only
  USE noncollin_module,   ONLY : noncolin
  USE control_flags,      ONLY : nstep, istep, conv_elec, conv_ions, &
       lpath, lmetadyn
  USE io_files,           ONLY : nd_nmbr, iunpath, tmp_dir
  USE path_variables,     ONLY : conv_path
  USE path_base,          ONLY : initialize_path, search_mep
  USE metadyn_base,       ONLY : metadyn_init
  USE path_io_routines,   ONLY : io_path_start, io_path_stop
  USE io_global,          ONLY : ionode
  USE check_stop,         ONLY : check_stop_init
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  CHARACTER (LEN=9) :: code = 'PWSCF'
  !
  !
  ! ... use ".FALSE." to disable all clocks except the total cpu time clock
  ! ... use ".TRUE."  to enable clocks
  !
  CALL init_clocks( .TRUE. )
  CALL start_clock( code )
  !
  CALL startup( nd_nmbr, code, version_number )
  !
  IF ( ionode ) THEN
     !
     WRITE( UNIT = stdout, &
          FMT = '(/5X,"Ultrasoft (Vanderbilt) Pseudopotentials")')
     !
     WRITE( unit = stdout, FMT = 9010 ) &
          ntypx, npk, lmaxx, nchix, ndmx, nbrx, nqfx
     !
  END IF
  !
  CALL iosys()
  !
  CALL check_stop_init()
  !
  IF ( ionode .AND. noncolin ) &
       WRITE( UNIT = stdout, &
       & FMT = '(/,5X,"non-colinear magnetization allowed",/)' )
  IF ( ionode .AND. gamma_only ) &
       WRITE( UNIT = stdout, &
       & FMT = '(/,5X,"gamma-point specific algorithms are used",/)' )
  !
  !
  CALL init_run()
  !
  istep = 0
  !
  !
  !
9010 FORMAT( /,5X,'Current dimensions of program pwscf are:', /, &
       & /,5X,'ntypx = ',I2,'   npk = ',I5,'  lmax = ',I2   &
       & /,5X,'nchix = ',I2,'  ndmx = ',I5,'  nbrx = ',I2,'  nqfx = ',I2 )
  !
#endif
  return
end subroutine start_pw
!
!----------------------------------------------------------------------------
!
subroutine end_pw
  !
#ifdef DFT_PW
  USE control_flags,      ONLY : conv_ions
  !
  IMPLICIT NONE
  !
  CALL stop_run( conv_ions )
  !
#endif
  return
end subroutine end_pw
