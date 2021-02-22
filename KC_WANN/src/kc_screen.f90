!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-------------------------------------------------------------------
PROGRAM kc_screen
  !-----------------------------------------------------------------
  !
  !!  This is the main driver of the kc_screen.x code. It computes the 
  !!  screaning parameters alpha for KOOPMANS. It reads the ouput of 
  !!  a PWscf calculation and the U matrices from W90 and computes
  !!  the orbital dependent screening coefficients as described in 
  !!  N. Colonna et al. JCTC 14, 2549 (2018) 
  !!  https://pubs.acs.org/doi/10.1021/acs.jctc.7b01116
  !!
  !!  Code written by Nicola Colonna (EPFL April 2019) 
  !
  USE mp_global,             ONLY : mp_startup,mp_global_end 
  USE environment,           ONLY : environment_start, environment_end
  USE check_stop,            ONLY : check_stop_init
  USE klist,                 ONLY : nkstot
  USE lsda_mod,              ONLY : nspin
  USE control_kc_wann,       ONLY : calculation
  USE mp_global,             ONLY : mp_startup
  USE mp_world,              ONLY : world_comm
  USE mp_pools,              ONLY : intra_pool_comm
  USE mp_bands,              ONLY : intra_bgrp_comm, inter_bgrp_comm
  USE command_line_options,  ONLY : ndiag_
  !
  IMPLICIT NONE
  !
  include 'laxlib.fh'
  !
  CHARACTER(LEN=9) :: code='KC_WANN'
  !
  ! 1) Initialize MPI, clocks, print initial messages
  CALL mp_startup ( )
  CALL laxlib_start ( ndiag_, world_comm, intra_bgrp_comm, &
       do_distr_diag_inside_bgrp_ = .true. )
  CALL set_mpi_comm_4_solvers( intra_pool_comm, intra_bgrp_comm, &
       inter_bgrp_comm )
  CALL environment_start ( code )
  !
  CALL check_stop_init ( )
  !
  calculation = 'screen'
  !
  ! 2) Read the input file and the PW output
  CALL kc_readin( ) 
  !
  ! 3) Set up for the KC calculation. 
  CALL kc_setup_screen( )
  !
  ! 4) Compute the screening coefficient via Linear Respoonse
  CALL screen_coeff ( )
  ! 
  CALL clean_pw( .TRUE. )
  CALL close_kc ( ) 
  !
  !
  IF (nkstot/nspin .gt. 1) CALL print_clock_pw ( )
  CALL print_clock_kc ( )
  !
  ! 5) Clean and Close 
  CALL mp_global_end()
  CALL environment_end( code )
  !
END PROGRAM kc_screen
