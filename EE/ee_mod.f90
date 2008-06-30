!
! Copyright (C) 2007-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! original version by I. Dabo and N. Marzari (MIT)
!
! contributions by E. Lamas and S. de Gironcoli (SISSA/DEMOCRITOS)
!
!----------------------------------------------------------------------------
MODULE ee_mod
  !----------------------------------------------------------------------------
  !
  ! ... this modules contains the variables needed for the 
  ! ... ELECTROSTATIC EMBEDDING method
  !
  !       Ismaila Dabo and Nicola Marzari MIT
  !
  USE kinds, ONLY :  DP
  SAVE
  !
  LOGICAL ::                        &
       do_comp,                     &
       do_coarse,                   &
       do_mltgrid                            
  REAL (KIND=DP), ALLOCATABLE ::    &
       vcomp(:),                    &
       vloccoul(:),                 &
       rhoion(:),                   &
       vcoul(:),                    &
       vsolvation(:),               &
       atomicspread(:)
  INTEGER ::                        &
       n_self_interaction,          &
       n_charge_compensation,       &
       n_cycle,                     &
       ncomp,                       &
       poisson_maxiter,             &
       mr1, mr2, mr3,               & ! MG grid dimensions
       ncompx,                      &
       ncompy,                      &
       ncompz,                      &
       nlev,                        &
       itmax,                       &
       n_smoothing,                 &
       icomp,                       &
       whichbc(3)                       
  REAL(KIND=DP) ::                  &
       ecomp,                       &
       poisson_thr,                 &
       mixing_charge_compensation,  &
       ecutcoarse,                  &
       errtol,                      &
       vloc_of_g_zero,              &
       madelung,                    &
       rhocut,                      &
       rhocution,                   &
       cellmin(3),                  &
       cellmax(3),                  &
       omegafact,                   &
       comp_thr,                    &
       rhothr,                      &
       rhothrbis,                   &
       depsthr,                     &
       epsthr,                      &
       deltapot,                    &
       smoothspr( 3 ) ,             &
       rhoionmax,                   &
       tbeta,                       &
       epsinfty

  CHARACTER (LEN=256) ::            &
       which_compensation,          &
       which_smoothing
  !
END MODULE ee_mod
!
