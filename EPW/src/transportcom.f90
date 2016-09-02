  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Roxana Margine
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !--------------------------------------------------------------------------
  MODULE transport_common
  !-------------------------------------------------------------------------- 
  USE kinds, ONLY :  DP
  !
  SAVE
  !
  ! Global variables for transport
  !
  INTEGER :: lower_bnd
  !! lower bound for the k-depend index among the mpi pools

  REAL(DP) :: mobilityh_save, mobilityel_save

  REAL(DP), ALLOCATABLE :: tdf_sigma(:,:,:), transp_temp(:), & 
                           SigmaS(:,:), SigmaS2(:,:), Seebeck(:,:), & 
                           Kappael(:,:), Kappa(:,:)
  !
  ! tdf_sigma(9,nbnd,nkf) : transport distribution function
  ! transp_temp(nstemp) : temperature array
  ! SigmaS(9,nstemp) :
  ! SigmaS2(9,nstemp) : power factor
  ! Seebeck(9,nstemp) : Seebeck coefficient
  ! Kappael(9,nstemp) : elec contribution to thermal conductivity
  ! Kappa(9,nstemp) : thermal conductivity
  !
  END MODULE transport_common
  !
  !--------------------------------------------------------------------------
  MODULE transportcom
  !-------------------------------------------------------------------------- 
  !
  USE transport_common
  !
  END MODULE transportcom
