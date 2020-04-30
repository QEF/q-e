!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE init_vloc()
  !----------------------------------------------------------------------
  !! This routine computes the fourier coefficient of the local
  !! potential vloc(ig,it) for each type of atom.
  !
  USE atom,           ONLY : msh, rgrid
  USE m_gth,          ONLY : vloc_gth
  USE kinds,          ONLY : DP
  USE uspp_param,     ONLY : upf
  USE ions_base,      ONLY : ntyp => nsp
  USE cell_base,      ONLY : omega, tpiba2
  USE vlocal,         ONLY : vloc
  USE gvect,          ONLY : ngl, gl
  USE Coul_cut_2D,    ONLY : do_cutoff_2D, cutoff_lr_Vloc
  !
  IMPLICIT NONE
  !
  INTEGER :: nt
  ! counter on atomic types
  !
  CALL start_clock( 'init_vloc' )
  !
  vloc(:,:) = 0._DP
  !
  DO nt = 1, ntyp
     !
     ! compute V_loc(G) for a given type of atom
     !
     IF ( upf(nt)%is_gth ) THEN
        !
        ! special case: GTH pseudopotential
        !
        CALL vloc_gth( nt, upf(nt)%zp, tpiba2, ngl, gl, omega, vloc(1,nt) )
        !
     ELSE IF ( upf(nt)%tcoulombp ) THEN
        !
        ! special case: pseudopotential is coulomb 1/r potential
        !
        CALL vloc_coul( upf(nt)%zp, tpiba2, ngl, gl, omega, vloc(1,nt) )
        !
     ELSE
        !
        ! normal case
        !
        CALL vloc_of_g( rgrid(nt)%mesh, msh(nt), rgrid(nt)%rab, rgrid(nt)%r, &
                        upf(nt)%vloc(1), upf(nt)%zp, tpiba2, ngl, gl, omega, &
                        vloc(1,nt) )
        !
     ENDIF
     !
  ENDDO
  !
  ! in 2D calculations the long range part of vloc(g) (erf/r part)
  ! was not re-added in g-space because everything is caclulated in
  ! radial coordinates, which is not compatible with 2D cutoff. 
  ! It will be re-added each time vloc(g) is used in the code. 
  ! Here, this cutoff long-range part of vloc(g) is computed only once
  ! by the routine below and stored
  !
  IF (do_cutoff_2D) CALL cutoff_lr_Vloc() 
  !
  CALL stop_clock( 'init_vloc' )
  !
  RETURN
  !
END SUBROUTINE init_vloc

