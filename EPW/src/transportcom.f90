  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
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
  INTEGER :: upper_bnd
  !! lower bound for the k-depend index among the mpi pools
  INTEGER, ALLOCATABLE :: ixkqf_tr(:,:)
  !! Mapping matrix from k+q (where q is full BZ) to IBZ
  INTEGER, ALLOCATABLE :: s_BZtoIBZ_full(:,:,:,:)
  !! Rotation that brink that k-point from BZ to IBZ
  !
  REAL(kind=DP), ALLOCATABLE :: transp_temp(:), & 
                           SigmaS(:,:), SigmaS2(:,:), Seebeck(:,:), & 
                           Kappael(:,:), Kappa(:,:),  mobilityh_save(:), &
                           mobilityel_save(:)
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
