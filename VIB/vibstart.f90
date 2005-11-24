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
!  Program CPVIB (programmed by Silviu Zilberman)
!  
!  Calculates the vibrational modes, Born effective charges and infrared
!  cross-section for isolated molecules/clusters in vacuum.
!  The programs uses the CP code as the underlying DFT engine
!  for wave function relaxation (SCF) and the frozen phonon method for
!  construction of the energy hessian (dynamical matrix),
!  i.e. displace each atom in all three cartesian directions to construct
!  the numerical first derivative of the forces.
!  
!  The program uses two mandatory input files (optional files are discussed
!  further below):
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
!  delta           real ( default = 0.05 )
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
!  Additional information about the input/output files:
!  
!  1. The program will look for the presence of 'prefix.vib.isotope' file.
!     It is a text file containing NAT lines, where NAT is the number of atoms.
!     Each line contains a single number, the mass of this particular atom
!     in the same order as in the CP input file. The masses are in AMU, i.e.
!     the hydrogen is ~1, carbon is 12 etc. If this file is present, these masses
!     are used in the calculation. If not, the default masses are used,
!     as specified in the CP input, and an isotopes file is created.
!     This is useful for testing the isotope shifts, without having to
!     redo the full calculation.
!  
!  2. After the energy hessian is calculated, the program will write
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
PROGRAM cpvib
  !----------------------------------------------------------------------------
  !
  USE control_flags,         ONLY : program_name
  USE environment,           ONLY : environment_start
  USE input,                 ONLY : read_input_file, iosys_pseudo, iosys
  USE io_global,             ONLY : io_global_start, io_global_getionode
  USE kinds,                 ONLY : DP
  USE mp_global,             ONLY : mp_global_start
  USE mp,                    ONLY : mp_end, mp_start, mp_env
  USE vibrations,            ONLY : start_vibrations, calc_hessian, &
       analysis, end_vibrations
  !
  IMPLICIT NONE
  !
  INTEGER                        :: mpime, nproc, gid, ionode_id, &
       restart_cyc_counter
  INTEGER,           PARAMETER   :: root = 0
  LOGICAL                        :: ionode
  REAL    (KIND=DP)              :: E_minus, dip_minus(3)
  !
  ! ... program starts here
  !
  program_name = 'CP90'
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
  CALL init_run()
  !     
  !  CALL memstat( 1 )
  !
  CALL read_input_vib()
  CALL start_vibrations(restart_cyc_counter, E_minus, dip_minus)
  CALL calc_hessian    (restart_cyc_counter, E_minus, dip_minus)
  CALL analysis ()
  CALL end_vibrations()
  !
  CALL terminate_run()
  !
  CALL stop_run( .TRUE. )
  !
  STOP
  !
END PROGRAM cpvib
