!
! Copyright (C) 2001-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE data_structure( gamma_only )
  !-----------------------------------------------------------------------
  ! this routine sets the data structure for the fft arrays
  ! (both the smooth and the dense grid)
  ! In the parallel case, it distributes columns to processes, too
  !
  USE kinds,      ONLY : DP
  USE mp,         ONLY : mp_max
  USE mp_bands,   ONLY : me_bgrp, nproc_bgrp, root_bgrp, intra_bgrp_comm, &
                         ntask_groups
  USE mp_pools,   ONLY : inter_pool_comm
  USE fft_base,   ONLY : dfftp, dffts, dtgs, fft_base_info
  USE fft_base,   ONLY : smap
  USE fft_types,  ONLY : fft_type_init
  USE cell_base,  ONLY : at, bg, tpiba
  USE klist,      ONLY : xk, nks
  USE gvect,      ONLY : gcutm, gvect_init
  USE gvecs,      ONLY : gcutms, gvecs_init
  USE gvecw,      ONLY : gcutw, gkcut
  USE io_global,  ONLY : stdout, ionode
  USE task_groups,ONLY : task_groups_init
  !
  IMPLICIT NONE
  LOGICAL, INTENT(in) :: gamma_only
  INTEGER :: ik, ngm_, ngs_
#if defined (__MPI) && !defined (__USE_3D_FFT)
  LOGICAL :: lpara = .true.
#else
  LOGICAL :: lpara = .false.
#endif
  !
  ! ... calculate gkcut = max |k+G|^2, in (2pi/a)^2 units
  !
  IF (nks == 0) THEN
     !
     ! if k-points are automatically generated (which happens later)
     ! use max(bg)/2 as an estimate of the largest k-point
     !
     gkcut = 0.5d0 * max ( &
        sqrt (sum(bg (1:3, 1)**2) ), &
        sqrt (sum(bg (1:3, 2)**2) ), &
        sqrt (sum(bg (1:3, 3)**2) ) )
  ELSE
     gkcut = 0.0d0
     DO ik = 1, nks
        gkcut = max (gkcut, sqrt ( sum(xk (1:3, ik)**2) ) )
     ENDDO
  ENDIF
  gkcut = (sqrt (gcutw) + gkcut)**2
  !
  ! ... find maximum value among all the processors
  !
  CALL mp_max (gkcut, inter_pool_comm )
  !
  ! ... set up fft descriptors, including parallel stuff: sticks, planes, etc.
  !
  CALL fft_type_init( dffts, smap, "wave", gamma_only, lpara, intra_bgrp_comm, at, bg, gkcut, gcutms/gkcut, ntask_groups )
  CALL fft_type_init( dfftp, smap, "rho", gamma_only, lpara, intra_bgrp_comm, at, bg,  gcutm, ntask_groups=ntask_groups )
  CALL task_groups_init( dffts, dtgs, ntask_groups )
  CALL fft_base_info( ionode, stdout )
  ngs_ = dffts%ngl( dffts%mype + 1 )
  ngm_ = dfftp%ngl( dfftp%mype + 1 )
  IF( gamma_only ) THEN
     ngs_ = (ngs_ + 1)/2
     ngm_ = (ngm_ + 1)/2
  END IF
  !
  !     on output, ngm_ and ngs_ contain the local number of G-vectors
  !     for the two grids. Initialize local and global number of G-vectors
  !
  call gvect_init ( ngm_ , intra_bgrp_comm )
  call gvecs_init ( ngs_ , intra_bgrp_comm )
  !

END SUBROUTINE data_structure

