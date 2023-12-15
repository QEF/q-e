!
! Copyright (C) 2019 National Institute of Advanced Industrial Science and Technology (AIST)
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE rotate_xpsi_gamma_gpu( h_psi_ptr, s_psi_ptr, overlap, &
                              npwx, npw, nstart, nbnd, psi_d, evc_d, hevc_d, sevc_d, e_d )
  !----------------------------------------------------------------------------
  !
  ! ... Serial version of rotate_xpsi for Gamma-only calculations
  ! ... This version assumes real wavefunctions (k=0) with only
  ! ... half plane waves stored: psi_d(-G)=psi_d*(G), except G=0
  !
#if defined(__CUDA)
  USE cudafor
  USE cublas
#else
#define cublasDGEMM dgemm
#endif
  USE util_param,    ONLY : DP
  USE mp_bands_util, ONLY : intra_bgrp_comm, inter_bgrp_comm, root_bgrp_id, nbgrp, my_bgrp_id, &
                            me_bgrp, root_bgrp
  USE mp_bands_util, ONLY : gstart ! index of the first nonzero G 
  USE mp,            ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  include 'laxlib.fh'
  !
  ! ... I/O variables
  !
  INTEGER, INTENT(IN) :: npw, npwx, nstart, nbnd
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix psi, as declared in the calling pgm unit
    ! input number of states
    ! output number of states
    ! first G with nonzero norm
  LOGICAL, INTENT(IN) :: overlap
    ! if .FALSE. : S|psi> not needed
  COMPLEX(DP), INTENT(INOUT) :: psi_d(npwx,nstart)
  COMPLEX(DP), INTENT(OUT)   :: evc_d(npwx,nbnd)
    ! input and output eigenvectors (may overlap)
  COMPLEX(DP), INTENT(OUT) :: hevc_d(npwx,nbnd), sevc_d(npwx,nbnd)
    ! H|psi> and S|psi>
  REAL(DP), INTENT(OUT) :: e_d(nbnd)
    ! eigenvalues
#if defined(__CUDA)
  attributes(DEVICE) :: psi_d, evc_d, hevc_d, sevc_d, e_d
#endif
  !
  ! ... local variables:
  !
  INTEGER                  :: npw2, npwx2
  INTEGER                  :: n_start, n_end, my_n, i
  REAL(DP),    ALLOCATABLE :: hr_d(:,:), sr_d(:,:), vr_d(:,:)
  COMPLEX(DP), ALLOCATABLE :: tpsi_d(:,:), hpsi_d(:,:), spsi_d(:,:)
  REAL(DP),    ALLOCATABLE :: en_d(:)
#if defined(__CUDA)
  attributes(DEVICE) :: hr_d, sr_d, vr_d, tpsi_d, hpsi_d, spsi_d, en_d
#endif
  !
  EXTERNAL :: h_psi_ptr, s_psi_ptr
    ! h_psi_ptr(npwx,npw,nbnd,psi,hpsi)
    !     calculates H|psi>
    ! s_psi_ptr(npwx,npw,nbnd,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx,npol,nbnd)


  npw2  = 2 * npw
  npwx2 = 2 * npwx
  !
  IF ( gstart == -1 ) CALL errore( 'rotxpsig', 'gstart variable not initialized', 1 )
  !
  ALLOCATE( tpsi_d( npwx, nstart ) )
  ALLOCATE( hpsi_d( npwx, nstart ) )
  IF ( overlap ) &
  ALLOCATE( spsi_d( npwx, nstart ) )
  ALLOCATE( hr_d( nstart, nstart ) )    
  ALLOCATE( sr_d( nstart, nstart ) )    
  ALLOCATE( vr_d( nstart, nstart ) )    
  ALLOCATE( en_d( nstart ) )
  !
  CALL start_clock('rotxpsig')
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
        psi_d(1,i) = CMPLX( DBLE( psi_d(1,i) ), 0.D0, kind=DP)
     END DO
  END IF
  !
  CALL start_clock('rotxpsig:hpsi')
  !
  CALL h_psi_ptr( npwx, npw, nstart, psi_d, hpsi_d)
  !
  CALL stop_clock('rotxpsig:hpsi')
  !
  IF ( overlap ) THEN
     !
     CALL start_clock('rotxpsig:spsi')
     !
     CALL s_psi_ptr( npwx, npw, nstart, psi_d, spsi_d )
     !
     CALL stop_clock('rotxpsig:spsi')
     !
  END IF
  !
  CALL start_clock('rotxpsig:hc') 
  hr_d = 0.D0
  CALL divide(inter_bgrp_comm, nstart, n_start, n_end)
  my_n = n_end - n_start + 1
  !
  IF ( n_start .le. n_end ) &
  CALL cublasDGEMM( 'T', 'N', nstart, my_n, npw2, 2.D0, psi_d, &
                    npwx2, hpsi_d(1,n_start), npwx2, 0.D0, hr_d(1,n_start), nstart )
  !
  IF ( gstart == 2 ) &
  CALL CGcudaDGER( nstart, my_n, -1.D0, psi_d, npwx2, hpsi_d(1,n_start), npwx2, hr_d(1,n_start), nstart )
  !
  CALL mp_sum( hr_d , inter_bgrp_comm )
  !
  CALL mp_sum( hr_d , intra_bgrp_comm )
  !
  CALL stop_clock('rotxpsig:hc')
  !     
  CALL start_clock('rotxpsig:sc')
  !
  sr_d = 0.D0
  !
  IF ( overlap ) THEN 
     !
     IF ( n_start .le. n_end ) &
     CALL cublasDGEMM( 'T', 'N', nstart, my_n, npw2, 2.D0, psi_d, &
                        npwx2, spsi_d(1,n_start), npwx2, 0.D0, sr_d(1,n_start), nstart )
     !
     IF ( gstart == 2 ) &
     CALL CGcudaDGER( nstart, my_n, -1.D0, psi_d, npwx2, spsi_d(1,n_start), npwx2, sr_d(1,n_start), nstart )
     !              
  ELSE
     !
     IF ( n_start .le. n_end ) &
     CALL cublasDGEMM( 'T', 'N', nstart, my_n, npw2, 2.D0, psi_d, &
                        npwx2, psi_d(1,n_start), npwx2, 0.D0, sr_d(1,n_start), nstart )
     !
     IF ( gstart == 2 ) &
     CALL CGcudaDGER( nstart, my_n, -1.D0, psi_d, npwx2, psi_d(1,n_start), npwx2, sr_d(1,n_start), nstart )
     !
  END IF
  !
  CALL mp_sum( sr_d , inter_bgrp_comm )
  !
  CALL mp_sum( sr_d , intra_bgrp_comm )
  !
  CALL stop_clock('rotxpsig:sc')
  !
  ! ... Diagonalize
  !
  CALL start_clock('rotxpsig:diag')
  !
  CALL diaghg( nstart, nbnd, hr_d, sr_d, nstart, en_d, vr_d, me_bgrp, root_bgrp, intra_bgrp_comm )
  !
  CALL stop_clock('rotxpsig:diag')
  !
  CALL start_clock('rotxpsig:evc')
  !
  !$cuf kernel do(1) <<<*,*>>>
  DO i=1, nbnd
     e_d(i) = en_d(i)
  END DO
  !
  ! ... update the basis set
  !
  tpsi_d = psi_d
  !
  evc_d  = (0.D0, 0.D0)
  hevc_d = (0.D0, 0.D0)
  IF ( overlap ) &
  sevc_d = (0.D0, 0.D0)
  !
  IF ( n_start .le. n_end ) THEN
     !
     CALL cublasDGEMM( 'N', 'N', npw2, nbnd, my_n, 1.D0, &
                 tpsi_d(1,n_start), npwx2,  vr_d(n_start,1), nstart, 0.D0, evc_d,  npwx2 )
     !
     CALL cublasDGEMM( 'N', 'N', npw2, nbnd, my_n, 1.D0, &
                 hpsi_d(1,n_start), npwx2,  vr_d(n_start,1), nstart, 0.D0, hevc_d, npwx2 )
     !
     IF ( overlap ) &
     CALL cublasDGEMM( 'N', 'N', npw2, nbnd, my_n, 1.D0, &
                 spsi_d(1,n_start), npwx2,  vr_d(n_start,1), nstart, 0.D0, sevc_d, npwx2 )
     !
  END IF
  !
  CALL mp_sum( evc_d,  inter_bgrp_comm )
  CALL mp_sum( hevc_d, inter_bgrp_comm )
  IF ( overlap ) &
  CALL mp_sum( sevc_d, inter_bgrp_comm )
  !
  CALL stop_clock('rotxpsig:evc')
  !
  DEALLOCATE( en_d )
  DEALLOCATE( vr_d )
  DEALLOCATE( sr_d )
  DEALLOCATE( hr_d )
  IF ( overlap ) &
  DEALLOCATE( spsi_d )
  DEALLOCATE( hpsi_d )
  DEALLOCATE( tpsi_d )
  !
  CALL stop_clock('rotxpsig')
  !
  !CALL print_clock('rotxpsig')
  !CALL print_clock('rotxpsig:hpsi')
  !CALL print_clock('rotxpsig:spsi')
  !CALL print_clock('rotxpsig:hc')
  !CALL print_clock('rotxpsig:sc')
  !CALL print_clock('rotxpsig:diag')
  !CALL print_clock('rotxpsig:evc')
  !
  RETURN
  !
END SUBROUTINE rotate_xpsi_gamma_gpu
