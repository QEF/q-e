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
!***  Molecular Dynamics using Density-Functional Theory                   ****
!***  this is the main routine driver for Car-Parrinello simulations       ****
!******************************************************************************
!***  based on
!***  version 11 of cpv code including ggapw 07/04/99
!***  copyright Alfredo Pasquarello 10/04/1996
!***  parallelized and converted to f90 by Paolo Giannozzi (2000),
!***  PBE added by Michele Lazzeri (2000)
!***  variable-cell dynamics by Andrea Trave (1998-2000)
!***  Makov Payne Correction for charged systems by Filippo De Angelis
!******************************************************************************
!***  also based on
!***  FPMD code written by Carlo Cavazzoni, Paolo Focher, Guido Chiarotti, 
!***  Sandro Scandolo
!******************************************************************************
!***  appropriate citation for use of this code:
!***  Car-Parrinello method    R. Car and M. Parrinello, PRL 55, 2471 (1985) 
!***  current implementation   A. Pasquarello, K. Laasonen, R. Car, 
!***                           C. Lee, and D. Vanderbilt, PRL 69, 1982 (1992);
!***                           K. Laasonen, A. Pasquarello, R. Car, 
!***                           C. Lee, and D. Vanderbilt, PRB 47, 10142 (1993).
!***  implementation gga       A. Dal Corso, A. Pasquarello, A. Baldereschi,
!***                           and R. Car, PRB 53, 1180 (1996).
!***  implementation Wannier   M. Sharma, Y. Wu and R. Car, Int.J.Quantum.Chem.
!***  function dynamics        95, 821, (2003).
!***
!***  implementation           M. Sharma and R.Car, ???
!***  Electric Field
!***  ensemble-DFT
!***  cf. "Ensemble Density-Functional Theory for Ab Initio Molecular Dynamics
!***       of Metals and Finite-Temperature Insulators"  PRL v.79,nbsp.7 (1997)
!***       nbsp. Marzari, D. Vanderbilt and M.C. Payne
!***  string methods           Yosuke Kanai et al.  J. Chem. Phys., 
!***                           2004, 121, 3359-3367
!***
!***  parallelization          C. Cavazzoni and G. L. Chiarotti "A Parallel and
!***                           Modular Car-Parrinello Code", Computer Physics 
!***                           Communication 123, 56-76, (1999).
!***                           P. Giannozzi, F. De Angelis and R. Car ...
!***  variable cell            P. Focher et al., Europhys. Lett. 36, 345-351
!***                           (1994)
!******************************************************************************
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
PROGRAM main
  !----------------------------------------------------------------------------
  !
  USE input,         ONLY : read_input_file, iosys_pseudo, iosys
  USE io_global,     ONLY : io_global_start, io_global_getmeta
  USE mp_global,     ONLY : mp_global_start, init_pool
  USE mp,            ONLY : mp_end, mp_start, mp_env, mp_bcast
  USE control_flags, ONLY : lneb, lsmd, lmetadyn, program_name
  USE environment,   ONLY : environment_start
  USE check_stop,    ONLY : check_stop_init
  !
  IMPLICIT NONE
  !
  INTEGER            :: mpime, nproc, world, meta_ionode_id
  INTEGER            :: nimage, ntask_groups
  LOGICAL            :: meta_ionode
  INTEGER, PARAMETER :: root = 0
  !
  ! ... program starts here
  !
  program_name = 'CP90'
  !
  ! ... initialize MPI (parallel processing handling)
  !
  CALL mp_start()
  !
  ! ... get from communication sub-sistem basic parameters
  ! ... to handle processors
  !
  CALL mp_env( nproc, mpime, world )
  !
  ! ... now initialize module holding processors and groups
  ! ... variables
  !
  CALL mp_global_start( root, mpime, world, nproc )
  !
  ! ... mpime = processor number, starting from 0
  ! ... nproc = number of processors
  ! ... world = group index of all processors
  ! ... root  = index of the root processor
  !
  ! ... initialize input output
  !
  CALL io_global_start( mpime, root )
  !
  ! ... get the "meta" io node
  !
  CALL io_global_getmeta( meta_ionode, meta_ionode_id )
  !
  IF ( meta_ionode ) THEN
     !
     ! ... check for command line arguments
     !
     CALL get_arg_nimage( nimage )
     !
     nimage = MAX( nimage, 1 )
     nimage = MIN( nimage, nproc )
     !
     CALL get_arg_ntg( ntask_groups )
     !
  END IF
  !
  CALL mp_bcast( nimage,       meta_ionode_id, world )
  CALL mp_bcast( ntask_groups, meta_ionode_id, world )
  !
  ! ... start the environment
  !
  CALL environment_start( )
  !
  ! ... readin the input file
  !
  CALL read_input_file()
  !
  ! ... read in pseudopotentials files and then
  ! ... copy pseudopotential parameters into internal variables
  !
  CALL iosys_pseudo()
  !
  ! ... copy-in input parameters from input_parameter module
  !
  CALL iosys()
  !
  CALL check_stop_init()
  !
  ! ... here reorganize processors in groups
  !
  CALL init_pool( nimage, ntask_groups )
  !
  IF ( lneb ) THEN
     !
     CALL neb_loop( )
     !
  ELSE IF ( lsmd ) THEN
     !
     CALL smd_loop( 1 )
     !
  ELSE
     !
     CALL cpr_loop( 1 )
     !
  END IF
  !
  CALL stop_run( .TRUE. )
  !
  STOP
  !
END PROGRAM main
