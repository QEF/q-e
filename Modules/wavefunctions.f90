!
! Copyright (C) 2002-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#if defined(__CUDA)
#define PINMEM ,PINNED 
#else
#define PINMEM
#endif

!=----------------------------------------------------------------------------=!
   MODULE wavefunctions
!=----------------------------------------------------------------------------=!
     !! Wavefunction arrays.
     !
     USE kinds, ONLY :  DP
#if defined (__CUDA)
     USE cudafor
#endif

     IMPLICIT NONE
     SAVE

     !
     COMPLEX(DP), ALLOCATABLE, TARGET :: &
       evc(:,:)
       !! wavefunctions in the PW basis set.  
       !! noncolinear case: first index is a combined PW + spin index
       !
#if defined(__CUDA)
       attributes(PINNED) :: evc
#endif
     !
     COMPLEX(DP) , ALLOCATABLE, TARGET :: psic(:)
     !! additional memory for FFT
     COMPLEX(DP) , ALLOCATABLE, TARGET :: psic_nc(:,:)
     !! additional memory for FFT for the noncolinear case
     !
     !
     ! electronic wave functions, CPV code
     ! distributed over gvector and bands
     !
!dir$ attributes align: 4096 :: c0_bgrp, cm_bgrp, phi
     COMPLEX(DP), ALLOCATABLE PINMEM :: c0_bgrp(:,:)  ! wave functions at time t
     COMPLEX(DP), ALLOCATABLE PINMEM :: cm_bgrp(:,:)  ! wave functions at time t-delta t
     COMPLEX(DP), ALLOCATABLE :: phi(:,:) ! |phi> = s'|c0> = |c0> + sum q_ij |i><j|c0>
     COMPLEX(DP), ALLOCATABLE :: c0_d(:,:)  ! wave functions at time t
     COMPLEX(DP), ALLOCATABLE :: cm_d(:,:)  ! wave functions at time t-delta t
#if defined (__CUDA)
     ATTRIBUTES(DEVICE) :: c0_d, cm_d, phi
#endif
     ! for hybrid functionals in CP with Wannier functions
     COMPLEX(DP), ALLOCATABLE :: cv0(:,:) ! Lingzhu Kong

   CONTAINS

      SUBROUTINE deallocate_wavefunctions
       IF( ALLOCATED( cv0) ) DEALLOCATE( cv0)   ! Lingzhu Kong
       IF( ALLOCATED( c0_bgrp ) ) DEALLOCATE( c0_bgrp )
       IF( ALLOCATED( cm_bgrp ) ) DEALLOCATE( cm_bgrp )
       IF( ALLOCATED( phi ) ) DEALLOCATE( phi )
       IF( ALLOCATED( psic_nc ) ) DEALLOCATE( psic_nc )
       IF( ALLOCATED( psic ) ) DEALLOCATE( psic )
       IF( ALLOCATED( evc ) ) DEALLOCATE( evc )
#if defined (__CUDA)
       IF( ALLOCATED( c0_d ) ) DEALLOCATE( c0_d )
       IF( ALLOCATED( cm_d ) ) DEALLOCATE( cm_d )
#endif
     END SUBROUTINE deallocate_wavefunctions

     SUBROUTINE allocate_cp_wavefunctions( ngw, nbspx, vnbsp, lwfpbe0nscf )
       INTEGER, INTENT(IN) :: ngw, nbspx, vnbsp
       LOGICAL, INTENT(IN) :: lwfpbe0nscf
       INTEGER :: ierr
       ALLOCATE( c0_bgrp( ngw, nbspx ), STAT=ierr )
       IF( ierr /= 0 ) &
         CALL errore( ' allocate_cp_wavefunctions ', ' allocating on CPU ', ABS( ierr ) )
       c0_bgrp = (0_DP,0_DP)
       ALLOCATE( cm_bgrp( ngw, nbspx ), STAT=ierr )
       IF( ierr /= 0 ) &
         CALL errore( ' allocate_cp_wavefunctions ', ' allocating on CPU ', ABS( ierr ) )
       cm_bgrp = (0_DP,0_DP)
       ALLOCATE( phi( ngw, nbspx ), STAT=ierr )
       IF( ierr /= 0 ) &
         CALL errore( ' allocate_cp_wavefunctions ', ' allocating on CPU ', ABS( ierr ) )
       phi = (0_DP,0_DP)
       IF(lwfpbe0nscf) THEN
         ALLOCATE(cv0( ngw, vnbsp ), STAT=ierr )   ! Lingzhu Kong
         IF( ierr /= 0 ) &
           CALL errore( ' allocate_cp_wavefunctions ', ' allocating on CPU ', ABS( ierr ) )
         cv0 = (0_DP,0_DP)
       END IF
#if defined (__CUDA)
       ALLOCATE( c0_d( ngw, nbspx ), STAT=ierr )
       IF( ierr /= 0 ) &
         CALL errore( ' allocate_cp_wavefunctions ', ' allocating on GPU ', ABS( ierr ) )
       ALLOCATE( cm_d( ngw, nbspx ), STAT=ierr )
       IF( ierr /= 0 ) &
         CALL errore( ' allocate_cp_wavefunctions ', ' allocating on GPU ', ABS( ierr ) )
#endif
     END SUBROUTINE

!=----------------------------------------------------------------------------=!
   END MODULE wavefunctions
!=----------------------------------------------------------------------------=!
