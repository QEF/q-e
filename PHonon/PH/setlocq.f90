!
! Copyright (C) 2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE init_vlocq ( xq )
  !----------------------------------------------------------------------
  !
  USE kinds,                ONLY : dp
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE gvect,                ONLY : g, ngm, ecutrho
  USE cell_base,            ONLY : omega, tpiba2
  USE ions_base,            ONLY : nsp
  USE vloc_mod,             ONLY : vloc_of_g, init_tab_vloc
  USE Coul_cut_2D,          ONLY : do_cutoff_2D     
  USE Coul_cut_2D_ph,       ONLY : cutoff_lr_Vlocq , cutoff_fact_qg, lr_Vlocq
  USE eqv,                  ONLY : vlocq
  !
  IMPLICIT NONE
  INTEGER :: ig, nt, ierr
  REAL(dp) :: qmax
  REAL(dp), INTENT(IN) :: xq(3)
  REAL(dp), ALLOCATABLE :: qg(:)
  !
  qmax = sqrt ( tpiba2 * ( xq(1)**2+xq(2)**2+xq(3)**2 ) ) + sqrt(ecutrho)
  CALL init_tab_vloc (qmax, do_cutoff_2d, omega, intra_bgrp_comm, ierr )
  IF ( ierr == 1 ) THEN
     CALL errore('init_vloc','Coulomb or GTH PPs incompatible with 2D cutoff &
             & (see upflib/vloc_mod.f90)',ierr)
  ELSE IF ( ierr == -1 ) THEN
     CALL infomsg('init_vloc','Interpolation table for Vloc re-allocated')
  END IF
  !
  ALLOCATE ( qg(ngm) )
  DO ig = 1 , ngm
     qg(ig) = (xq(1)+g(1,ig))**2 + (xq(2)+g(2,ig))**2 + (xq(3)+g(3,ig))**2
  END DO
  DO nt = 1, nsp
     CALL vloc_of_g( nt, ngm, qg, tpiba2, do_cutoff_2d, omega, vlocq(:,nt) )
  END DO
  DEALLOCATE ( qg )
  !
  ! For 2d calculations, we need to initialize the fact for the q+G
  ! component of the cutoff of the Coulomb interaction
  !
  IF (do_cutoff_2D) call cutoff_fact_qg()
  !
  ! In 2D calculations the long range part of vlocq(g) (erf/r part)
  ! was not re-added in g-space because everything is calculated in
  ! radial coordinates, which is not compatible with 2D cutoff. 
  ! Here, this cutoff long-range part of vlocq(g) is computed only once
  ! by the routine below, added to vlocq, and stored for possible separate use.
  !
  IF (do_cutoff_2D) THEN
     CALL cutoff_lr_Vlocq()
     vlocq = vlocq + lr_Vlocq
  ENDIF
  !
  RETURN
  !
END SUBROUTINE init_vlocq
