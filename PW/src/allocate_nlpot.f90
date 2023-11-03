!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE allocate_nlpot
  !-----------------------------------------------------------------------
  !! This routine allocates arrays containing the non-local part of the
  !! pseudopotential for each atom or atomic species.
  !
  !! Requires in input:  
  !! * dimensions: nhm, nsp, nat, lmaxkb, nbetam, nwfcm, nspin
  !! * parameters: ecutrho, qnorm, dq, ecutwfc, cell_factor
  !! * options: tqr, noncolin, lspinorb
  !
  !! Computes the following global quantities:  
  !! * nqx: number of points of the interpolation table
  !! * nqxq: as above, for q-function interpolation table
  !
  USE control_flags,    ONLY : tqr, use_gpu
  USE ions_base,        ONLY : nat
  USE cellmd,           ONLY : cell_factor
  USE klist,            ONLY : qnorm
  USE lsda_mod,         ONLY : nspin
  USE noncollin_module, ONLY : noncolin, lspinorb
  USE gvect,            ONLY : ecutrho
  USE gvecw,            ONLY : ecutwfc
  USE uspp_data,        ONLY : dq, nqx, nqxq, allocate_uspp_data
  USE uspp,             ONLY : allocate_uspp
  USE uspp_param,       ONLY : upf, lmaxq, lmaxkb, nh, nhm, nsp, nbetam, nwfcm
  IMPLICIT NONE
  !
  ! Note: computation of the number of beta functions for
  ! each atomic type and the maximum number of beta functions
  ! and the number of beta functions of the solid has been
  ! moved to init_run.f90 : pre_init()
  !
  call allocate_uspp(use_gpu,noncolin,lspinorb,tqr,nhm,nsp,nat,nspin)
  !
  ! This routine is called also by the phonon code, in which case it should
  ! allocate an array that includes q+G vectors up to |q+G|_max <= |Gmax|+|q|
  !
  nqxq = INT( ( (SQRT(ecutrho) + qnorm) / dq + 4) * cell_factor )
  !
  ! Calculate dimensions for array tab (including a possible factor
  ! coming from cell contraction during variable cell relaxation/MD)
  !
  nqx = INT( (SQRT(ecutwfc) / dq + 4) * cell_factor )
  !
  ! uspp_data  actual allocation
  !
  call allocate_uspp_data(use_gpu,nqxq,nqx,nbetam,nwfcm,lmaxq,nsp)
  !
  return
  !
END SUBROUTINE allocate_nlpot
