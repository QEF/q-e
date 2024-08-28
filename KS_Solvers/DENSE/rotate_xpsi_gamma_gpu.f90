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
                              npwx, npw, nstart, nbnd, psi, evc, hevc, sevc, e_d )
  !----------------------------------------------------------------------------
  !
  ! ... Serial version of rotate_xpsi for Gamma-only calculations
  ! ... This version assumes real wavefunctions (k=0) with only
  ! ... half plane waves stored: psi(-G)=psi*(G), except G=0
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
  COMPLEX(DP), INTENT(INOUT) :: psi(npwx,nstart)
  COMPLEX(DP), INTENT(OUT)   :: evc(npwx,nbnd)
    ! input and output eigenvectors (may overlap)
  COMPLEX(DP), INTENT(OUT) :: hevc(npwx,nbnd), sevc(npwx,nbnd)
    ! H|psi> and S|psi>
  REAL(DP), INTENT(OUT) :: e_d(nbnd)
    ! eigenvalues
#if defined(__CUDA)
  attributes(DEVICE) :: e_d
#endif
  !
  ! ... local variables:
  !
  INTEGER                  :: npw2, npwx2
  INTEGER                  :: n_start, n_end, my_n, i
  REAL(DP),    ALLOCATABLE :: hr_d(:,:), sr_d(:,:), vr_d(:,:)
  COMPLEX(DP), ALLOCATABLE :: tpsi_d(:,:), hpsi(:,:), spsi(:,:)
  REAL(DP),    ALLOCATABLE :: en_d(:)
#if defined(__CUDA)
  attributes(DEVICE) :: hr_d, sr_d, vr_d, tpsi_d, en_d
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
  ALLOCATE( hpsi( npwx, nstart ) )
  !$acc enter data create(hpsi)
  IF ( overlap ) THEN
     ALLOCATE( spsi( npwx, nstart ) )
     !$acc enter data create(spsi)
  ENDIF
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
  ! ... set Im[ psi(G=0) ] -  needed for numerical stability
  !
  IF ( gstart == 2 ) THEN
     !$acc kernels      
     DO i=1,nstart
        psi(1,i) = CMPLX( DBLE( psi(1,i) ), 0.D0, kind=DP)
     END DO
     !$acc end kernels
  END IF
  !
  CALL start_clock('rotxpsig:hpsi')
  !
  CALL h_psi_ptr( npwx, npw, nstart, psi, hpsi)
  !
  CALL stop_clock('rotxpsig:hpsi')
  !
  IF ( overlap ) THEN
     !
     CALL start_clock('rotxpsig:spsi')
     !
     CALL s_psi_ptr( npwx, npw, nstart, psi, spsi )
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
  !$acc host_data use_device(psi, hpsi)
  IF ( n_start .le. n_end ) &
  CALL cublasDGEMM( 'T', 'N', nstart, my_n, npw2, 2.D0, psi, &
                    npwx2, hpsi(1,n_start), npwx2, 0.D0, hr_d(1,n_start), nstart )
  !
  IF ( gstart == 2 ) &
  CALL MYDGER( nstart, my_n, -1.D0, psi, npwx2, hpsi(1,n_start), npwx2, hr_d(1,n_start), nstart )
  !$acc end host_data
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
  !$acc host_data use_device(psi, spsi, hpsi)
  IF ( overlap ) THEN 
     !
     IF ( n_start .le. n_end ) &
     CALL cublasDGEMM( 'T', 'N', nstart, my_n, npw2, 2.D0, psi, &
                        npwx2, spsi(1,n_start), npwx2, 0.D0, sr_d(1,n_start), nstart )
     !
     IF ( gstart == 2 ) &
     CALL MYDGER( nstart, my_n, -1.D0, psi, npwx2, spsi(1,n_start), npwx2, sr_d(1,n_start), nstart )
     !              
  ELSE
     !
     IF ( n_start .le. n_end ) &
     CALL cublasDGEMM( 'T', 'N', nstart, my_n, npw2, 2.D0, psi, &
                        npwx2, psi(1,n_start), npwx2, 0.D0, sr_d(1,n_start), nstart )
     !
     IF ( gstart == 2 ) &
     CALL MYDGER( nstart, my_n, -1.D0, psi, npwx2, psi(1,n_start), npwx2, sr_d(1,n_start), nstart )
     !
  END IF
  !$acc end host_data
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
  !$acc kernels
  tpsi_d = psi
  !$acc end kernels
  !
  !$acc kernels
  evc  = (0.D0, 0.D0)
  hevc = (0.D0, 0.D0)
  IF ( overlap ) &
  sevc = (0.D0, 0.D0)
  !$acc end kernels
  !
  IF ( n_start .le. n_end ) THEN
     !
     !$acc host_data use_device(hevc, sevc, hpsi, spsi, evc)
     CALL cublasDGEMM( 'N', 'N', npw2, nbnd, my_n, 1.D0, &
                 tpsi_d(1,n_start), npwx2,  vr_d(n_start,1), nstart, 0.D0, evc,  npwx2 )
     !
     CALL cublasDGEMM( 'N', 'N', npw2, nbnd, my_n, 1.D0, &
                 hpsi(1,n_start), npwx2,  vr_d(n_start,1), nstart, 0.D0, hevc, npwx2 )
     !
     IF ( overlap ) &
     CALL cublasDGEMM( 'N', 'N', npw2, nbnd, my_n, 1.D0, &
                 spsi(1,n_start), npwx2,  vr_d(n_start,1), nstart, 0.D0, sevc, npwx2 )
     !$acc end host_data
     !
  END IF
  !
  !$acc host_data use_device(hevc, sevc, evc)
  CALL mp_sum( evc,  inter_bgrp_comm )
  CALL mp_sum( hevc, inter_bgrp_comm )
  IF ( overlap ) &
  CALL mp_sum( sevc, inter_bgrp_comm )
  !$acc end host_data
  !
  CALL stop_clock('rotxpsig:evc')
  !
  DEALLOCATE( en_d )
  DEALLOCATE( vr_d )
  DEALLOCATE( sr_d )
  DEALLOCATE( hr_d )
  IF ( overlap ) THEN
     !$acc exit data delete(spsi)
     DEALLOCATE( spsi )
  ENDIF
  !$acc exit data delete(hpsi)
  DEALLOCATE( hpsi )
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
