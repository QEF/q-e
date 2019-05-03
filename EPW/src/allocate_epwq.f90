  !                                                                        
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino     
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                        
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! Code adapted from PH/allocate_phq - Quantum-ESPRESSO group
  ! 09/2009 There is a lot of excess in this file.  
  !
  !----------------------------------------------------------------------- 
  SUBROUTINE allocate_epwq
  !-----------------------------------------------------------------------
  !!
  !! Dynamical allocation of arrays: quantities needed for the linear
  !! response problem
  !!
  !! RM - Nov/Dec - 2014 - Imported the noncolinear case implemented by xlzhang
  !! SP - 2016 - Updated for QE 5
  !! RM - Jan 2019 - Updated based on QE 6.3
  !!
  USE ions_base,    ONLY : nat, ntyp => nsp
  USE pwcom,        ONLY : nbnd, nspin
  USE gvect,        ONLY : ngm
  USE noncollin_module, ONLY : noncolin, nspin_mag
  USE spin_orb,     ONLY : lspinorb
  USE phus,         ONLY : int1, int1_nc, int2, int2_so, &
                           int4, int4_nc, int5, int5_so, & 
                           alphap
  USE lrus,         ONLY : becp1
  USE elph2,        ONLY : elph, el_ph_mat
  USE becmod,       ONLY : becp, allocate_bec_type
  USE uspp_param,   ONLY : nhm
  USE uspp,         ONLY : okvan, nkb
  USE klist,        ONLY : nks
  USE fft_base,     ONLY : dfftp
  USE epwcom,       ONLY : nstemp  
  ! 
  IMPLICIT NONE
  !  
  INTEGER :: ik
  !! k-point
  INTEGER :: ipol
  !! Polarization index
  !
  !  ALLOCATE space for the quantities needed in EPW
  !
  ! SP: nrxx is not used in QE 5 ==> tg_nnr is the maximum among nnr
  !     This should have the same dim as nrxx had.
  !     ALLOCATE (dmuxc ( nrxx, nspin, nspin))  
  ! SP: Again a new change in QE (03/08/2016)  
  !     ALLOCATE (dmuxc ( dffts%tg_nnr, nspin, nspin))    
  ! SP: Following new FFT restructuration from Aug. 2017 (SdG)
  !     nnr = local number of FFT grid elements  ( ~nr1*nr2*nr3/nproc )
  !     nnr_tg = local number of grid elements for task group FFT ( ~nr1*nr2*nr3/proc3 )  
  !           --> tg = task group    
  !     ALLOCATE (dmuxc ( dffts%nnr, nspin, nspin))    
  !
  ! 
  RETURN
  ! 
  END SUBROUTINE allocate_epwq
