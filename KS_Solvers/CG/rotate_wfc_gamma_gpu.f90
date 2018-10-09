!
! Copyright (C) 2003-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE rotate_wfc_gamma_gpu( h_psi_gpu, s_psi_gpu, overlap, &
                             npwx, npw, nstart, nbnd, psi_d, evc_d, e_d )
  !----------------------------------------------------------------------------
  !
  ! ... Serial version of rotate_wfc for Gamma-only calculations
  ! ... This version assumes real wavefunctions (k=0) with only
  ! ... half plane waves stored: psi_d(-G)=psi_d*(G), except G=0
  !
  USE LAXlib,       ONLY : diaghg
#if defined(__CUDA)
  USE cudafor
  USE cublas
#endif
  USE cg_param,      ONLY : DP
  USE mp_bands_util, ONLY : intra_bgrp_comm, inter_bgrp_comm, root_bgrp_id, &
          nbgrp, my_bgrp_id
  USE mp_bands_util, ONLY : gstart ! index of the first nonzero G
  USE mp,            ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER :: npw, npwx, nstart, nbnd
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix psi_d, as declared in the calling pgm unit
    ! input number of states
    ! output number of states
  LOGICAL :: overlap
    ! if .FALSE. : S|psi_d> not needed
  COMPLEX(DP) :: psi_d(npwx,nstart), evc_d(npwx,nbnd)
    ! input and output eigenvectors (may overlap)
  REAL(DP) :: e_d(nbnd)
    ! eigenvalues
#if defined(__CUDA)
  attributes(DEVICE) :: psi_d, evc_d, e_d
#endif
  !
  ! ... local variables
  !
  INTEGER                  :: npw2, npwx2
  COMPLEX(DP), ALLOCATABLE :: aux_d(:,:)
  REAL(DP),    ALLOCATABLE :: hr_d(:,:), sr_d(:,:), vr_d(:,:)
  REAL(DP),    ALLOCATABLE :: en_d(:)
#if defined(__CUDA)
  attributes(DEVICE) :: aux_d, hr_d, sr_d, vr_d, en_d
#endif
  INTEGER :: n_start, n_end, my_n, i, j
  !
  EXTERNAL  h_psi_gpu,    s_psi_gpu
    ! h_psi(npwx,npw,nvec,psi_d,hpsi)
    !     calculates H|psi_d>
    ! s_psi(npwx,npw,nvec,spsi)
    !     calculates S|psi_d> (if needed)
    !     Vectors psi_d,hpsi,spsi are dimensioned (npwx,npol,nvec)

  npw2  = 2 * npw
  npwx2 = 2 * npwx

  IF ( gstart == -1 ) CALL errore( 'regter', 'gstart variable not initialized', 1 )

  !
  ALLOCATE( aux_d(  npwx, nstart ) )
  ALLOCATE( hr_d( nstart, nstart ) )
  ALLOCATE( sr_d( nstart, nstart ) )
  ALLOCATE( vr_d( nstart, nstart ) )
  ALLOCATE( en_d( nstart ) )
  call start_clock('rotwfcg'); !write(*,*) 'start rotwfcg' ; FLUSH(6)
  !
  ! ... Set up the Hamiltonian and Overlap matrix on the subspace :
  !
  ! ...      H_ij = <psi_i| H |psi_j>     S_ij = <psi_i| S |psi_j>
  !
  ! ... set Im[ psi_d(G=0) ] -  needed for numerical stability
  !
  IF ( gstart == 2 ) THEN
     !$cuf kernel do(1)
     DO i=1,nstart
        psi_d(1,i) = CMPLX( DBLE( psi_d(1,i) ), 0.D0,kind=DP)
     END DO
  END IF
  !
  call start_clock('rotwfcg:hpsi'); !write(*,*) 'start rotwfcg:hpsi' ; FLUSH(6)
  CALL h_psi_gpu( npwx, npw, nstart, psi_d, aux_d )
  call stop_clock('rotwfcg:hpsi'); !write(*,*) 'stop rotwfcg:hpsi' ; FLUSH(6)
  !
  call start_clock('rotwfcg:hc'); !write(*,*) 'start rotwfcg:hc' ; FLUSH(6)
  hr_d=0.D0
  CALL divide(inter_bgrp_comm,nstart,n_start,n_end)
  my_n = n_end - n_start + 1; !write (*,*) nstart,n_start,n_end
  if (n_start .le. n_end) &
  CALL cublasDGEMM( 'T','N', nstart, my_n, npw2, 2.D0, psi_d, npwx2, aux_d(1,n_start), npwx2, 0.D0, hr_d(1,n_start), nstart )
  IF ( gstart == 2 ) call CGcudaDGER( nstart, my_n, -1.D0, psi_d, npwx2, aux_d(1,n_start), npwx2, hr_d(1,n_start), nstart )
  CALL mp_sum( hr_d, inter_bgrp_comm )
  !
  CALL mp_sum( hr_d, intra_bgrp_comm )
  !
  sr_d=0.D0
  IF ( overlap ) THEN
     !
     CALL s_psi_gpu( npwx, npw, nstart, psi_d, aux_d )
     !
     if (n_start .le. n_end) &
     CALL cublasDGEMM( 'T','N', nstart, my_n, npw2, 2.D0, psi_d, npwx2, aux_d(1,n_start), npwx2, 0.D0, sr_d(1,n_start), nstart )
     IF ( gstart == 2 ) CALL CGcudaDGER( nstart, my_n, -1.D0, psi_d, npwx2, aux_d(1,n_start), npwx2, sr_d(1,n_start), nstart )
     !
  ELSE
     !
     if (n_start .le. n_end) &
     CALL cublasDGEMM( 'T','N', nstart, my_n, npw2, 2.D0, psi_d, npwx2, psi_d(1,n_start), npwx2, 0.D0, sr_d(1,n_start), nstart )
     IF ( gstart == 2 ) CALL CGcudaDGER( nstart, my_n, -1.D0, psi_d, npwx2, psi_d(1,n_start), npwx2, sr_d(1,n_start), nstart )
     !
  END IF
  CALL mp_sum( sr_d, inter_bgrp_comm )
  !
  CALL mp_sum( sr_d, intra_bgrp_comm )
  call stop_clock('rotwfcg:hc'); !write(*,*) 'stop rotwfcg:hc' ; FLUSH(6)
  !
  ! ... Diagonalize
  !
  call start_clock('rotwfcg:diag'); !write(*,*) 'start rotwfcg:diag' ; FLUSH(6)
  CALL diaghg( nstart, nbnd, hr_d, sr_d, nstart, en_d, vr_d )
  call stop_clock('rotwfcg:diag'); !write(*,*) 'stop rotwfcg:diag' ; FLUSH(6)
  call start_clock('rotwfcg:evc_d'); !write(*,*) 'start rotwfcg:evc_d' ; FLUSH(6)
  !
  !$cuf kernel do(1)
  DO i=1, nbnd
     e_d(i) = en_d(i)
  END DO
  !
  ! ... update the basis set
  !
  aux_d=(0.D0,0.D0)
  if (n_start .le. n_end) &
  CALL cublasDGEMM( 'N','N', npw2, nbnd, my_n, 1.D0, psi_d(1,n_start), npwx2, vr_d(n_start,1), nstart, 0.D0, aux_d, npwx2 )
  CALL mp_sum( aux_d, inter_bgrp_comm )
  !
  !$cuf kernel do(2)
  DO i=1, nbnd
     DO j=1, npwx
        evc_d(j,i) = aux_d(j,i)
     END DO
  END DO
  call stop_clock('rotwfcg:evc_d'); !write(*,*) 'stop rotwfcg:evc_d' ; FLUSH(6)
  !
  DEALLOCATE( en_d )
  DEALLOCATE( vr_d )
  DEALLOCATE( sr_d )
  DEALLOCATE( hr_d )
  DEALLOCATE( aux_d )
  call stop_clock('rotwfcg'); !write(*,*) 'stop rotwfcg' ; FLUSH(6)
  !call print_clock('rotwfcg')
  !call print_clock('rotwfcg:hpsi')
  !call print_clock('rotwfcg:hc')
  !call print_clock('rotwfcg:diag')
  !call print_clock('rotwfcg:evc_d')
  !
  RETURN
  !
END SUBROUTINE rotate_wfc_gamma_gpu

! In principle this can go away .......
SUBROUTINE CGcudaDGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
#if defined(__CUDA)
    use cudafor
    use cublas
#endif
!     .. Scalar Arguments ..
    DOUBLE PRECISION ::  ALPHA
    INTEGER          ::   INCX, INCY, LDA, M, N
!     .. Array Arguments ..
    DOUBLE PRECISION :: A( LDA, * ), X( * ), Y( * )
#if defined(__CUDA)
    attributes(device) :: A, X, Y
#endif
    CALL DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )

END SUBROUTINE CGcudaDGER
