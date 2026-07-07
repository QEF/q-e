!
! Copyright (C) 2025 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE direct_diag_k( ngk_g, nbnd, hmat, evc_g, e, &
                          use_para_diag, notconv )
  !----------------------------------------------------------------------------
  !
  ! ... Direct diagonalization of an explicit dense Hermitian matrix.
  ! ... Solves H v = e v with the standard (S = I) eigenvalue problem,
  ! ... dispatching between serial (diagh) and parallel (pdiagh) LAXlib
  ! ... eigensolvers based on use_para_diag, with band-group broadcast
  ! ... when do_distr_diag_inside_bgrp is enabled.
  !
  ! ... Differs from iterative KS solvers (cegterg, crmmdiagg, ...) in that
  ! ... the matrix is supplied explicitly: there is no h_psi_ptr argument
  ! ... and no convergence loop. Operates entirely in the GLOBAL plane-wave
  ! ... basis (matrix is replicated on every process); the caller is
  ! ... responsible for any global <-> local layout conversion.
  !
  USE util_param,    ONLY : DP, stdout
  USE mp,            ONLY : mp_bcast
  USE mp_bands_util, ONLY : intra_bgrp_comm, me_bgrp, root_bgrp, &
                            my_bgrp_id, root_bgrp_id, nbgrp, inter_bgrp_comm
  USE LAXlib,        ONLY : diagh, pdiagh
  !
  IMPLICIT NONE
  !
  include 'laxlib.fh'
  !
  ! ... I/O variables
  !
  INTEGER,     INTENT(IN)    :: ngk_g
  !! global matrix dimension (number of plane waves at this k-point summed over procs)
  INTEGER,     INTENT(IN)    :: nbnd
  !! number of eigenpairs to keep
  COMPLEX(DP), INTENT(INOUT) :: hmat(ngk_g, ngk_g)
  !! Hermitian matrix to diagonalize (replicated; unchanged on exit per
  !! the diagh/pdiagh contract, but declared INOUT to match those interfaces)
  COMPLEX(DP), INTENT(OUT)   :: evc_g(ngk_g, nbnd)
  !! global eigenvectors (lowest nbnd, replicated on every process)
  REAL(DP),    INTENT(OUT)   :: e(nbnd)
  !! eigenvalues (lowest nbnd)
  LOGICAL,     INTENT(IN)    :: use_para_diag
  !! if .TRUE. use distributed pdiagh, else serial diagh
  INTEGER,     INTENT(OUT)   :: notconv
  !! number of unconverged roots (always 0 for direct diagonalization;
  !! kept for interface uniformity with iterative solvers)
  !
  ! ... local variables
  !
  REAL(DP),    ALLOCATABLE :: et_tmp(:)
  !! full-length eigenvalue buffer (LAXlib eigensolvers fill all ngk_g)
  COMPLEX(DP), ALLOCATABLE :: vmat(:,:)
  !! full-length eigenvector buffer:
  !! (ngk_g, nbnd) for serial diagh, (ngk_g, ngk_g) for pdiagh
  INTEGER :: neig
  !! number of eigenvector columns to allocate
  !
  ! ... LAXlib descriptor variables for the parallel path
  !
  INTEGER :: idesc(LAX_DESC_SIZE)
  INTEGER, ALLOCATABLE :: rank_ip(:,:)
  INTEGER, ALLOCATABLE :: idesc_ip(:,:,:)
  INTEGER :: nx
  LOGICAL :: la_proc
  LOGICAL :: do_distr_diag_inside_bgrp
  !
  CALL start_clock('direct_diag')
  !
  IF ( use_para_diag ) THEN
     neig = ngk_g
     WRITE(stdout,'(5X,"Diagonalizing ",I0,"x",I0," matrix (parallel)")') ngk_g, ngk_g
  ELSE
     neig = nbnd
     WRITE(stdout,'(5X,"Diagonalizing ",I0,"x",I0," matrix (serial)")') ngk_g, ngk_g
  END IF
  !
  ALLOCATE( et_tmp(ngk_g) )
  ALLOCATE( vmat(ngk_g, neig) )
  !
  IF ( use_para_diag ) THEN
     !
     ! ... Parallel diagonalization (pdiagh).
     ! ... Pattern follows KS_Solvers/DENSE/rotate_wfc_k.f90 for the
     ! ... do_distr_diag_inside_bgrp band-group broadcast.
     !
     CALL laxlib_getval( do_distr_diag_inside_bgrp = do_distr_diag_inside_bgrp )
     !
     CALL desc_init( ngk_g, nx, la_proc, idesc, rank_ip, idesc_ip )
     !
#if defined(__VERBOSE)
     WRITE(stdout,'(5X,"Parallel diagonalization descriptor:")')
     WRITE(stdout,'(5X,"  Matrix dimension (N)      =",I8)') idesc(LAX_DESC_N)
     WRITE(stdout,'(5X,"  Local block dimension     =",I8)') nx
     WRITE(stdout,'(5X,"  Processor grid            =",I4," x",I4)') &
          idesc(LAX_DESC_NPR), idesc(LAX_DESC_NPC)
     WRITE(stdout,'(5X,"  This processor position   = (",I3,",",I3,")")') &
          idesc(LAX_DESC_MYR), idesc(LAX_DESC_MYC)
     WRITE(stdout,'(5X,"  Active in diag?           =",L2)') la_proc
     WRITE(stdout,'(5X,"  do_distr_diag_inside_bgrp =",L2)') do_distr_diag_inside_bgrp
#endif
     !
     IF ( do_distr_diag_inside_bgrp ) THEN
        IF ( my_bgrp_id == root_bgrp_id ) THEN
           CALL pdiagh( ngk_g, hmat, et_tmp, vmat, idesc )
        END IF
        IF ( nbgrp > 1 ) THEN
           CALL mp_bcast( vmat,   root_bgrp_id, inter_bgrp_comm )
           CALL mp_bcast( et_tmp, root_bgrp_id, inter_bgrp_comm )
        END IF
     ELSE
        CALL pdiagh( ngk_g, hmat, et_tmp, vmat, idesc )
     END IF
     !
     IF ( ALLOCATED(rank_ip)  )  DEALLOCATE( rank_ip )
     IF ( ALLOCATED(idesc_ip) )  DEALLOCATE( idesc_ip )
     !
  ELSE
     !
     ! ... Serial Hermitian eigensolver (S = I assumed for the standard problem).
     !
     CALL diagh( ngk_g, nbnd, hmat, et_tmp, vmat, &
                 me_bgrp, root_bgrp, intra_bgrp_comm )
     !
  END IF
  !
  ! ... Copy out the lowest nbnd eigenpairs.
  !
  e(1:nbnd)        = et_tmp(1:nbnd)
  evc_g(:,1:nbnd)  = vmat(:,1:nbnd)
  notconv          = 0
  !
  DEALLOCATE( et_tmp, vmat )
  !
  CALL stop_clock('direct_diag')
  !
END SUBROUTINE direct_diag_k
