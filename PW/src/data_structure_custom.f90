!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE data_structure_custom(fc, smap_exx, gamma_only)
  !-----------------------------------------------------------------------
  ! this routine sets the data structure for the custom fft array
  ! In the parallel case, it distributes columns to processes, too
  !
  USE kinds,      ONLY : DP
  USE cell_base,  ONLY : at, bg, tpiba, tpiba2
  USE klist,      ONLY : xk, nks
  USE mp,         ONLY : mp_sum, mp_max,mp_barrier
  USE mp_exx,     ONLY : me_egrp, negrp, nproc_egrp, inter_egrp_comm, &
                         intra_egrp_comm, root_egrp, ntask_groups 
  USE mp_bands,   ONLY : intra_bgrp_comm
  USE fft_types,  ONLY : fft_type_init
  USE stick_base, ONLY : sticks_map

  USE fft_custom, ONLY : fft_cus, gvec_init
  USE fft_base,   ONLY : dfftp, smap
  USE gvect,      ONLY : gcutm
  !
  !
  IMPLICIT NONE
  
  TYPE(fft_cus) :: fc
  LOGICAL :: gamma_only
  REAL (DP) :: gkcut
  INTEGER :: ik, ngm_, ngs_, ngw_
  INTEGER :: me, nproc, inter_comm, intra_comm, root
  TYPE (sticks_map) :: smap_exx ! Stick map descriptor

  INTEGER :: kpoint
#if defined (__MPI)
  LOGICAL :: lpara = .true.
#else
  LOGICAL :: lpara = .false.
#endif

  ! sticks coordinates
  
  !
  !  Subroutine body
  !

  !
  ! compute gkcut calling an internal procedure
  !

  me = me_egrp
  nproc = nproc_egrp
  inter_comm = inter_egrp_comm
  intra_comm = intra_egrp_comm
  root = root_egrp

  IF (nks == 0) THEN
     !
     ! if k-points are automatically generated (which happens later)
     ! use max(bg)/2 as an estimate of the largest k-point
     !
     gkcut = 0.5d0 * MAX ( &
            &SQRT (SUM(bg (1:3, 1)**2) ), &
            &SQRT (SUM(bg (1:3, 2)**2) ), &
            &SQRT (SUM(bg (1:3, 3)**2) ) )
  ELSE
     gkcut = 0.0d0
     DO kpoint = 1, nks
        gkcut = MAX (gkcut, SQRT ( SUM(xk (1:3, kpoint)**2) ) )
     ENDDO
  ENDIF
  gkcut = (SQRT (fc%ecutt) / tpiba + gkcut)**2
  !
  ! ... find maximum value among all the processors
  !
  CALL mp_max (gkcut, inter_comm )
  !
  ! ... set up fft descriptors, including parallel stuff: sticks, planes, etc.
  !
  IF( negrp == 1 )THEN
     CALL fft_type_init( fc%dfftt, smap, "rho", gamma_only, lpara, &
                         intra_bgrp_comm, at, bg, fc%gcutmt, fc%gcutmt/gkcut )
  ELSE
     CALL fft_type_init( fc%dfftt, smap_exx, "rho", gamma_only, lpara, &
                         intra_comm, at, bg, fc%gcutmt, fc%gcutmt/gkcut )
  END IF
  ngs_ = fc%dfftt%ngl( fc%dfftt%mype + 1 )
  IF( gamma_only ) THEN
     ngs_ = (ngs_ + 1)/2
  END IF
  !
  !     on output, ngm_ and ngs_ contain the local number of G-vectors
  !     for the two grids. Initialize local and global number of G-vectors
  !
  CALL gvec_init (fc, ngs_ , intra_comm )

  
END SUBROUTINE data_structure_custom
