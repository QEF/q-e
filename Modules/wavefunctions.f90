!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!=----------------------------------------------------------------------------=!
   MODULE wavefunctions_module
!=----------------------------------------------------------------------------=!
     USE kinds
     USE kinds, ONLY :  DP

     IMPLICIT NONE
     SAVE

     !
     COMPLEX(KIND=DP), ALLOCATABLE, TARGET :: &
       evc(:,:)     ! wavefunctions in the PW basis

     !
     COMPLEX(KIND=DP) , ALLOCATABLE, TARGET :: &
       psic(:)      ! additional memory for FFT


     ! three index wave functions, for non colinear code
     !
     COMPLEX(KIND=DP), ALLOCATABLE, TARGET :: &
       evc_nc(:,:,:)     !  wavefunctions in the PW basis
     !
     COMPLEX(KIND=DP), ALLOCATABLE, TARGET :: &
       psic_nc(:,:)      !  additional memory for FFT
     !     

     ! electronic wave functions, FPMD code
     !
     COMPLEX(dbl), ALLOCATABLE :: c0(:,:,:,:)  ! wave functions at time t
     COMPLEX(dbl), ALLOCATABLE :: cm(:,:,:,:)  ! wave functions at time t-delta t
     COMPLEX(dbl), ALLOCATABLE :: cp(:,:,:,:)  ! wave functions at time t+delta t
     COMPLEX(dbl), ALLOCATABLE :: ce(:,:,:,:)  ! empty states wave func. at time t

!=----------------------------------------------------------------------------=!
   END MODULE wavefunctions_module
!=----------------------------------------------------------------------------=!
