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
  USE kinds,          ONLY : dp
  USE vloc_mod,       ONLY : vloc_of_g
  USE uspp_param,     ONLY : upf
  USE ions_base,      ONLY : ntyp => nsp
  USE cell_base,      ONLY : omega, tpiba2
  USE gvect,          ONLY : ngl, gl
  USE atom,           ONLY : msh, rgrid
  USE Coul_cut_2D,    ONLY : do_cutoff_2D, cutoff_lr_Vloc, lz
  USE esm,            ONLY : do_comp_esm, esm_bc
  USE vlocal,         ONLY : vloc
  !
  IMPLICIT NONE
  !
  INTEGER :: nt
  !! counter on atomic types
  INTEGER :: ierr
  !! counter on atomic types
  LOGICAL :: modified_coulomb
  !
  CALL start_clock( 'init_vloc' )
  !
  vloc(:,:) = 0._dp
  !
  modified_coulomb = do_cutoff_2D .OR. (do_comp_esm .and. ( esm_bc .ne. 'pbc' ))
  !
  DO nt = 1, ntyp
     IF (do_cutoff_2D .AND. rgrid(nt)%r(msh(nt)) > lz) THEN 
        CALL errore('init_vloc','2D cutoff smaller than pseudo cutoff radius: &
           & increase interlayer distance (or see Modules/read_pseudo.f90)',nt)
     END IF
     !
     ! compute V_loc(G) for a given type of atom
     !
     CALL vloc_of_g( nt, ngl, gl, tpiba2, modified_coulomb, omega, &
             vloc(1,nt), ierr )
     CALL errore('init_vloc','Coulomb or GTH PPs incompatible with 2D cutoff &
             & or ESM (see upflib/vloc_mod.f90)',ierr)
     !
  ENDDO
  !
  ! in 2D calculations the long range part of vloc(g) (erf/r part)
  ! was not re-added in g-space because everything is calculated in
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

