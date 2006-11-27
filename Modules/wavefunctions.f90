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
     USE kinds, ONLY :  DP

     IMPLICIT NONE
     SAVE

     !
     COMPLEX(DP), ALLOCATABLE, TARGET :: &
       evc(:,:)     ! wavefunctions in the PW basis set
                    ! noncolinear case: first index
                    ! is a combined PW + spin index
     !
     COMPLEX(DP) , ALLOCATABLE, TARGET :: &
       psic(:), &      ! additional memory for FFT
       psic_nc(:,:)    ! as above for the noncolinear case
     !
     ! electronic wave functions, FPMD code
     !
     COMPLEX(DP), ALLOCATABLE :: c0(:,:)  ! wave functions at time t
     COMPLEX(DP), ALLOCATABLE :: cm(:,:)  ! wave functions at time t-delta t
     COMPLEX(DP), ALLOCATABLE :: cp(:,:)  ! wave functions at time t+delta t

   CONTAINS

     SUBROUTINE deallocate_wavefunctions
       IF( ALLOCATED( c0 ) ) DEALLOCATE( c0 )
       IF( ALLOCATED( cm ) ) DEALLOCATE( cm )
       IF( ALLOCATED( cp ) ) DEALLOCATE( cp )
       IF( ALLOCATED( psic_nc ) ) DEALLOCATE( psic_nc )
       IF( ALLOCATED( psic ) ) DEALLOCATE( psic )
       IF( ALLOCATED( evc ) ) DEALLOCATE( evc )
     END SUBROUTINE deallocate_wavefunctions

!=----------------------------------------------------------------------------=!
   END MODULE wavefunctions_module
!=----------------------------------------------------------------------------=!
