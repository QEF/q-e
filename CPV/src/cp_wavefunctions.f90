!
! Copyright (C) 2002-2024 Quantum ESPRESSO Foundation
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
   MODULE cp_wavefunctions
!=----------------------------------------------------------------------------=!
     !
     !! electronic wave functions, CPV code
     !! distributed over G-vectors and bands
     !
     USE kinds, ONLY :  DP
     IMPLICIT NONE
     SAVE
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
     !
   CONTAINS
      !
      SUBROUTINE deallocate_cp_wavefunctions
       IF( ALLOCATED( cv0) ) DEALLOCATE( cv0)   ! Lingzhu Kong
       IF( ALLOCATED( c0_bgrp ) ) DEALLOCATE( c0_bgrp )
       IF( ALLOCATED( cm_bgrp ) ) DEALLOCATE( cm_bgrp )
       IF( ALLOCATED( phi ) ) DEALLOCATE( phi )
#if defined (__CUDA)
       IF( ALLOCATED( c0_d ) ) DEALLOCATE( c0_d )
       IF( ALLOCATED( cm_d ) ) DEALLOCATE( cm_d )
#endif
     END SUBROUTINE deallocate_cp_wavefunctions
     !
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
   END MODULE cp_wavefunctions
!=----------------------------------------------------------------------------=!

