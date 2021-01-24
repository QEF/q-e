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
  !! This routine sets the data structure for the FFT arrays (both the
  !! smooth and the dense grid).  
  !! In the parallel case, it distributes columns to processes too.  
  !! BEWARE: to compute \(\text{gkcut}\), \(\text{nks}\) and the list of
  !! k-points or \(\text{nks}=0\) and the primitive lattice vectors, 
  !! \(\text{bg}\) are needed.
  !
  USE kinds,      ONLY : DP
  USE mp,         ONLY : mp_max
  USE mp_bands,   ONLY : nproc_bgrp, intra_bgrp_comm, nyfft, ntask_groups
  USE mp_pools,   ONLY : inter_pool_comm
  USE fft_base,   ONLY : dfftp, dffts, fft_base_info, smap
  USE fft_types,  ONLY : fft_type_init
  USE cell_base,  ONLY : at, bg, tpiba
  USE klist,      ONLY : xk, nks
  USE gvect,      ONLY : gcutm, gvect_init
  USE gvecs,      ONLY : gcutms, gvecs_init, doublegrid
  USE gvecw,      ONLY : gcutw, gkcut
  USE io_global,  ONLY : stdout, ionode
  ! FIXME: find a better way to transmit these three variables, or remove them
  USE realus,     ONLY : real_space
  USE symm_base,  ONLY : fft_fact
  USE command_line_options, ONLY: pencil_decomposition_
  USE command_line_options, ONLY : nmany_
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: gamma_only
  !
  ! ... local variables
  !
  INTEGER :: ik, ngm_, ngs_
  LOGICAL :: lpara
  !
  lpara = ( nproc_bgrp > 1 )
  !
  ! ... calculate gkcut = max |k+G|^2, in (2pi/a)^2 units
  !
  IF (nks == 0) THEN
     !
     ! k-point list not available:
     ! use max(bg)/2 as an estimate of the largest k-point
     !
     gkcut = 0.5d0 * MAX( SQRT( SUM(bg(1:3,1)**2) ), &
                          SQRT( SUM(bg(1:3,2)**2) ), &
                          SQRT( SUM(bg(1:3,3)**2) )  )
  ELSE
     gkcut = 0.0d0
     DO ik = 1, nks
        gkcut = MAX(gkcut, SQRT( SUM(xk(1:3,ik)**2) ) )
     ENDDO
  ENDIF
  gkcut = (SQRT(gcutw) + gkcut)**2
  !
  ! ... find maximum value among all the processors
  !
  CALL mp_max( gkcut, inter_pool_comm )
  !
  ! ... set up fft descriptors, including parallel stuff: sticks, planes, etc.
  !
  ! task group are disabled if real_space calculation of calbec is used
  dffts%has_task_groups = (ntask_groups >1) .AND. .NOT. real_space
  CALL fft_type_init( dffts, smap, "wave", gamma_only, lpara, intra_bgrp_comm, &
       at, bg, gkcut, gcutms/gkcut, fft_fact=fft_fact, nyfft=nyfft, nmany=nmany_, use_pd=pencil_decomposition_  )
  CALL fft_type_init( dfftp, smap, "rho" , gamma_only, lpara, intra_bgrp_comm, &
       at, bg, gcutm , 4.d0, fft_fact=fft_fact, nyfft=nyfft, nmany=nmany_, use_pd=pencil_decomposition_ )
  ! define the clock labels ( this enables the corresponding fft too ! )
  dffts%rho_clock_label = 'ffts' ; dffts%wave_clock_label = 'fftw'
  dfftp%rho_clock_label = 'fft'
  ! this makes so that interpolation is just a copy.
  IF (.NOT.doublegrid) dfftp%grid_id = dffts%grid_id
  !
  CALL fft_base_info( ionode, stdout )
  ngs_ = dffts%ngl( dffts%mype + 1 )
  ngm_ = dfftp%ngl( dfftp%mype + 1 )
  !
  IF( gamma_only ) THEN
     ngs_ = (ngs_ + 1)/2
     ngm_ = (ngm_ + 1)/2
  ENDIF
  !
  ! ... on output, ngm_ and ngs_ contain the local number of G-vectors
  ! for the two grids. Initialize local and global number of G-vectors
  !
  CALL gvect_init( ngm_, intra_bgrp_comm )
  CALL gvecs_init( ngs_, intra_bgrp_comm )
  !
END SUBROUTINE data_structure

