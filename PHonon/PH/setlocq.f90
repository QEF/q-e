!
! Copyright (C) 2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE init_vlocq ( )
  !----------------------------------------------------------------------
  !
  USE kinds,                ONLY : dp
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE gvect,                ONLY : g, ngm, ecutrho
  USE cell_base,            ONLY : omega, tpiba2
  USE ions_base,            ONLY : nsp
  USE klist,                ONLY : qnorm
  USE qpoint,               ONLY : xq
  USE vloc_mod,             ONLY : vloc_of_g, init_tab_vloc
  USE Coul_cut_2D,          ONLY : do_cutoff_2D     
  USE Coul_cut_2D_ph,       ONLY : cutoff_lr_Vlocq , cutoff_fact_qg
  USE eqv,                  ONLY : vlocq
  !
  IMPLICIT NONE
  INTEGER :: ig, nt, ierr
  REAL(dp) :: qmax
  REAL(dp), ALLOCATABLE :: qg(:)
  !
  qmax = (sqrt(ecutrho)+qnorm)
  CALL init_tab_vloc (qmax, do_cutoff_2d, omega, intra_bgrp_comm, ierr )
  !
  ALLOCATE ( qg(ngm) )
  DO ig = 1 , ngm
     qg(ig) = (xq(1)+G(1,ig))**2 + (xq(2)+G(2,ig))**2 + (xq(3)+G(3,ig))**2
  END DO
  DO nt = 1, nsp
     CALL vloc_of_g( nt, ngm, qg, tpiba2, do_cutoff_2d, omega, vlocq(:,nt) )
  END DO
  DEALLOCATE ( qg )
  !
  ! for 2d calculations, we need to initialize the fact for the q+G 
  ! component of the cutoff of the COulomb interaction
  !
  IF (do_cutoff_2D) call cutoff_fact_qg()
  !
  !  in 2D calculations the long range part of vlocq(g) (erf/r part)
  ! was not re-added in g-space because everything is calculated in
  ! radial coordinates, which is not compatible with 2D cutoff. 
  ! It will be re-added each time vlocq(g) is used in the code. 
  ! Here, this cutoff long-range part of vlocq(g) is computed only once
  ! by the routine below and stored
  !
  IF (do_cutoff_2D) call cutoff_lr_Vlocq() 
  !
  RETURN
  !
END SUBROUTINE init_vlocq
