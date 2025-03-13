!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE data_structure_3drism(dfft, gvec, gamma_only, mp_task)
  !---------------------------------------------------------------------------
  !
  ! ... this routine sets the data structure for the 3D-RISM's fft array
  ! ... In the parallel case, it distributes columns to processes, too
  !
  USE cell_base,            ONLY : at, bg
  USE command_line_options, ONLY : nmany_
  USE fft_base,             ONLY : smap
  USE fft_types,            ONLY : fft_type_descriptor, fft_type_init
  USE gvec_3drism,          ONLY : gvec_type, gvec_init_3drism
  USE kinds,                ONLY : DP
  USE mp_bands,             ONLY : nproc_bgrp, nyfft
  USE mp_rism,              ONLY : mp_rism_task
  !USE symm_base,            ONLY : fft_fact
  !
  IMPLICIT NONE
  !
  TYPE(fft_type_descriptor), INTENT(INOUT) :: dfft
  TYPE(gvec_type),           INTENT(INOUT) :: gvec
  LOGICAL,                   INTENT(IN)    :: gamma_only
  TYPE(mp_rism_task),        INTENT(IN)    :: mp_task  ! must be same as intra_bgrp_comm
  !
  INTEGER :: ngm_
  INTEGER :: intra_comm
  LOGICAL :: lpara
  !
  lpara = (nproc_bgrp > 1)
  !
  intra_comm = mp_task%itask_comm
  !
  ! ... set up fft descriptors, including parallel stuff: sticks, planes, etc.
  !
  !CALL fft_type_init(dfft, smap, "rho", gamma_only, lpara, intra_comm, &
  !     at, bg, gvec%gcutm, 1.0_DP, fft_fact=fft_fact, nyfft=nyfft, nmany=nmany_)
  CALL fft_type_init(dfft, smap, "rho", gamma_only, lpara, intra_comm, &
       at, bg, gvec%gcutm, 1.0_DP, nyfft=nyfft, nmany=nmany_)  ! do not use symmetric mesh
  !
  dfft%rho_clock_label = 'fftr'  ! this is the label of FFT for 3D-RISM
  !
  ngm_ = dfft%ngl(dfft%mype + 1)
  IF (gamma_only) THEN
    ngm_ = (ngm_ + 1) / 2
  END IF
  !
  ! ... initialize local and global number of G-vectors
  !
  CALL gvec_init_3drism(gvec, ngm_, intra_comm)
  !
END SUBROUTINE data_structure_3drism
