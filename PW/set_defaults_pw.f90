!
! Copyright (C) 2008 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE set_defaults_pw
  !-----------------------------------------------------------------------------
  !
  ! ...  this subroutine sets the default values for the variables
  ! ...  read from input by pw which are not saved into the xml file. 
  ! ...  It has to be called by all the programs that
  ! ...  run pw. It is possible to change 
  ! ...  the default values initialized here. The variables which
  ! ...  are in input_parameters are initialized by that module
  ! ...  and do not need to be initialized here. 
  ! ...  Actually this routine should not be needed. All the variables
  ! ...  which are read from input should have a default value and should
  ! ...  be in a module that initializes them. The pw code should then
  ! ...  use those variables. Unfortunately the routine input now 
  ! ...  initializes several variables that are in pw modules and
  ! ...  have just a different name from the variable contained in 
  ! ...  input parameters, or are calculated starting from the variable
  ! ...  contained in input_parameters.  
  ! ...  Moreover many variables contained in control flags are not initialized
  ! ...  and needs to be initialized here ... 
  !
  !
  USE kinds,         ONLY : DP
  USE bp,            ONLY : lberry,   &
                            lelfield
  !
  USE basis,         ONLY : startingwfc, &
                            startingpot
  !
  USE char,          ONLY : crystal
  !
  USE cellmd,        ONLY : calc, lmovecell
  !
  USE force_mod,     ONLY : lforce, lstres
  !
  USE gvect,         ONLY : ecfixed, qcutz, q2sigma
  !
  USE klist,         ONLY : lxkcry, tot_charge, &
                            tot_magnetization, &
                            multiplicity

  USE relax,         ONLY : starting_scf_threshold
  !
  USE control_flags, ONLY : isolve, max_cg_iter, tr2, imix, &
                            nmix, iverbosity, niter, pot_order, wfc_order, &
                            assume_isolated, &
                            diago_full_acc, &
                            mixing_beta, &
                            upscale, &
                            nstep, &
                            iprint, &
                            nosym, &
                            io_level, lscf, lbfgs, lmd, lpath, lneb,   &
                            lsmd, lphonon, ldamped, lbands, lmetadyn, llang, &
                            lconstrain, lcoarsegrained, restart, &
                            use_para_diag

  USE bfgs_module,   ONLY : bfgs_ndim, &
                            trust_radius_max, &
                            trust_radius_min, &
                            trust_radius_ini, &
                            w_1, &
                            w_2
  USE us, ONLY : spline_ps
  USE a2F, ONLY : la2F

  !
  IMPLICIT NONE
  !
  iprint = 100000
  lberry   = .FALSE.
  lelfield = .FALSE.
  lxkcry=.FALSE.
  tot_charge = 0.0_DP
  tot_magnetization = -1
  multiplicity = 0
  nosym = .FALSE.
  ecfixed = 0.0_DP
  qcutz   = 0.0_DP
  q2sigma = 0.01_DP
  !
  !  ... postprocessing of DOS & phonons & el-ph
  la2F = .FALSE.
  !
  ! ... non collinear program variables
  !
  assume_isolated = .FALSE.
  !
  spline_ps = .FALSE.
  !
  diago_full_acc = .FALSE.
  !
  upscale           = 10.0_DP
  mixing_beta       = 0.7
  !
  ! ... BFGS defaults
  !
  bfgs_ndim        = 1
  trust_radius_max = 0.8_DP   ! bohr
  trust_radius_min = 1.E-4_DP ! bohr
  trust_radius_ini = 0.5_DP   ! bohr
  w_1              = 0.01_DP
  w_2              = 0.50_DP
  !
  startingpot = 'file'
  startingwfc = 'atomic'
  !
  restart        = .FALSE.
  !
  io_level = 1
  !
  ! ... various initializations of control variables
  !
  lscf      = .FALSE.
  lmd       = .FALSE.
  lmetadyn  = .FALSE.
  lpath     = .FALSE.
  lneb      = .FALSE.
  lsmd      = .FALSE.
  lmovecell = .FALSE.
  lphonon   = .FALSE.
  lbands    = .FALSE.
  lbfgs     = .FALSE.
  ldamped   = .FALSE.
  lforce    = .FALSE.
  lstres    = .FALSE.
  !
  nstep = 1
  !
  isolve = 1
  max_cg_iter = 100
  use_para_diag = .FALSE.
  !
  niter = 1000
  !
  pot_order = 0
  wfc_order = 0
  !
  tr2=1.D-6
  starting_scf_threshold = tr2
  imix = 0
  nmix = 0
  !
  iverbosity = 0
  !
  crystal     = ' '
  calc      = ' '
  !
  RETURN
  !
END SUBROUTINE set_defaults_pw
