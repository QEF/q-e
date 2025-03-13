!
! Copyright (C) 2002-2024 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!=----------------------------------------------------------------------------=!
   MODULE wavefunctions
!=----------------------------------------------------------------------------=!
     !! Wavefunction arrays.
     !
     USE kinds, ONLY :  DP
     USE control_flags, ONLY : use_gpu
#if defined (__CUDA)
     USE cudafor
#endif

     IMPLICIT NONE
     SAVE

     !
#if defined(__CUDA)
     COMPLEX(DP), POINTER, CONTIGUOUS :: evc(:,:)
#else
     COMPLEX(DP), ALLOCATABLE, TARGET :: evc(:,:)
#endif
       !! wavefunctions in the PW basis set.  
       !! noncolinear case: first index is a combined PW + spin index
       !
     COMPLEX(DP) , ALLOCATABLE, TARGET :: psic(:)
     !! additional memory for FFT
     COMPLEX(DP) , ALLOCATABLE, TARGET :: psic_nc(:,:)
     !! additional memory for FFT for the noncolinear case
     !
   CONTAINS

!----------------------------------------------------------------------------
SUBROUTINE allocate_wfc(npwx, npol, nbnd)
  !----------------------------------------------------------------------------
  !! Dynamical allocation of wavefunctions.  
  !! Requires dimensions: \(\text{npwx}\), \(\text{nbnd}\), \(\text{npol}\)
  !
#if defined (__CUDA)
  use, intrinsic :: iso_c_binding
#endif
  INTEGER, INTENT(IN) :: npwx, npol, nbnd
  !
  INTEGER :: istat
  !
  !
  ALLOCATE( evc(npwx*npol,nbnd) )
!civn: PIN evc memory here
#if defined(__CUDA)
  IF(use_gpu) istat = cudaHostRegister(C_LOC(evc(1,1)), sizeof(evc), cudaHostRegisterMapped)
  !$acc enter data create(evc)
#endif
  !
END SUBROUTINE allocate_wfc
!
!----------------------------------------------------------------------------
SUBROUTINE deallocate_wfc()
  !----------------------------------------------------------------------------
#if defined (__CUDA)
  use, intrinsic :: iso_c_binding
#endif
  !
  IMPLICIT NONE
  INTEGER :: istat
  !
#if defined(__CUDA)
  !$acc exit data delete(evc)
  IF(use_gpu) istat = cudaHostUnregister(C_LOC(evc(1,1)))
#endif
  IF( ALLOCATED( evc ) ) DEALLOCATE( evc )
  !
END SUBROUTINE deallocate_wfc
     !
!=----------------------------------------------------------------------------=!
   END MODULE wavefunctions
!=----------------------------------------------------------------------------=!

