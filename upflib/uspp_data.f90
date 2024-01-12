!
! Copyright (C) 2001-2023 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE uspp_data
  !
  !! These parameters are needed with the US pseudopotentials.
  !
  USE upf_kinds,      ONLY : DP
  !
  SAVE
  PRIVATE
  !
  PUBLIC :: nqxq, nqx, dq
  PUBLIC :: tab_at
  !
  PUBLIC :: allocate_uspp_data
  PUBLIC :: deallocate_uspp_data
  PUBLIC :: scale_uspp_data
  !
  INTEGER :: nqxq
  !! size of interpolation table
  INTEGER :: nqx
  !! number of interpolation points
  REAL(DP), PARAMETER:: dq = 0.01D0
  !! space between interpolation points 
  REAL(DP), ALLOCATABLE :: tab_at(:,:,:)
  !! interpolation table for atomic wfc
  !
contains
  !
  subroutine allocate_uspp_data(use_gpu,nqxq_,nqx_,nbetam,nwfcm,lmaxq,nsp)
     implicit none
     logical, intent(in) :: use_gpu
     integer, intent(in) :: nqxq_,nqx_,nbetam,nwfcm,lmaxq,nsp
     !
     if (nqxq_/=nqxq) call upf_error("allocate_uspp_data","invalid nqxq_",1)
     if (nqx_/=nqx)   call upf_error("allocate_uspp_data","invalid nqx_",1)
     !
     allocate(tab_at(nqx_,nwfcm,nsp))
     !$acc enter data create(tab_at)
     !
  end subroutine allocate_uspp_data
  !
  subroutine deallocate_uspp_data()
     implicit none
     !$acc exit data delete(tab_at)
     if( allocated( tab_at ) )    deallocate( tab_at )
     !
  end subroutine
  !
  subroutine scale_uspp_data( vol_ratio_m1 )
     ! vol_ratio_m1 = omega_old / omega
     implicit none
     real(DP), intent(in) :: vol_ratio_m1
     !
     tab_at(:,:,:)   = tab_at(:,:,:) * SQRT(vol_ratio_m1)
!$acc update device (tab_at)
  end subroutine scale_uspp_data
  !
END MODULE uspp_data

