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
                              npwx, npw, nstart, nbnd, psi, evc, hevc, sevc, e )
  !----------------------------------------------------------------------------
  !
  ! ... Serial version of rotate_xpsi for Gamma-only calculations
  ! ... This version assumes real wavefunctions (k=0) with only
  ! ... half plane waves stored: psi(-G)=psi*(G), except G=0
  !
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
  REAL(DP), INTENT(OUT) :: e(nbnd)
    ! eigenvalues
  !
  ! ... local variables:
  !
  INTEGER                  :: npw2, npwx2
  INTEGER                  :: n_start, n_end, my_n, i
  REAL(DP),    ALLOCATABLE :: hr(:,:), sr(:,:), vr(:,:)
  COMPLEX(DP), ALLOCATABLE :: tpsi(:,:), hpsi(:,:), spsi(:,:)
  REAL(DP),    ALLOCATABLE :: en(:)
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
  ALLOCATE( tpsi( npwx, nstart ) )
  ALLOCATE( hpsi( npwx, nstart ) )
  ALLOCATE( hr( nstart, nstart ) )    
  ALLOCATE( sr( nstart, nstart ) )    
  ALLOCATE( vr( nstart, nstart ) )    
  ALLOCATE( en( nstart ) )
  !$acc enter data create(hpsi, hr, sr, vr, tpsi, en )
  IF ( overlap ) THEN
     ALLOCATE( spsi( npwx, nstart ) )
     !$acc enter data create(spsi)
  ENDIF
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
  !$acc kernels
  hr = 0.D0
  !$acc end kernels
  CALL divide(inter_bgrp_comm, nstart, n_start, n_end)
  my_n = n_end - n_start + 1
  !
  !$acc host_data use_device(psi, hpsi, hr)
  IF ( n_start .le. n_end ) &
  CALL MYDGEMM( 'T', 'N', nstart, my_n, npw2, 2.D0, psi, &
                    npwx2, hpsi(1,n_start), npwx2, 0.D0, hr(1,n_start), nstart )
  !
  IF ( gstart == 2 ) &
  CALL MYDGER( nstart, my_n, -1.D0, psi, npwx2, hpsi(1,n_start), npwx2, hr(1,n_start), nstart )
  !
  CALL mp_sum( hr, inter_bgrp_comm )
  !
  CALL mp_sum( hr, intra_bgrp_comm )
  !$acc end host_data
  !
  CALL stop_clock('rotxpsig:hc')
  !     
  CALL start_clock('rotxpsig:sc')
  !
  !$acc kernels
  sr = 0.D0
  !$acc end kernels
  !
  !$acc host_data use_device(psi, spsi, hpsi, sr)
  IF ( overlap ) THEN 
     !
     IF ( n_start .le. n_end ) &
     CALL MYDGEMM( 'T', 'N', nstart, my_n, npw2, 2.D0, psi, &
                        npwx2, spsi(1,n_start), npwx2, 0.D0, sr(1,n_start), nstart )
     !
     IF ( gstart == 2 ) &
     CALL MYDGER( nstart, my_n, -1.D0, psi, npwx2, spsi(1,n_start), npwx2, sr(1,n_start), nstart )
     !              
  ELSE
     !
     IF ( n_start .le. n_end ) &
     CALL MYDGEMM( 'T', 'N', nstart, my_n, npw2, 2.D0, psi, &
                        npwx2, psi(1,n_start), npwx2, 0.D0, sr(1,n_start), nstart )
     !
     IF ( gstart == 2 ) &
     CALL MYDGER( nstart, my_n, -1.D0, psi, npwx2, psi(1,n_start), npwx2, sr(1,n_start), nstart )
     !
  END IF
  !
  CALL mp_sum( sr, inter_bgrp_comm )
  !
  CALL mp_sum( sr, intra_bgrp_comm )
  !$acc end host_data
  !
  CALL stop_clock('rotxpsig:sc')
  !
  ! ... Diagonalize
  !
  CALL start_clock('rotxpsig:diag')
  !
  !$acc host_data use_device(hr, sr, en, vr)
  CALL diaghg( nstart, nbnd, hr, sr, nstart, en, vr, me_bgrp, root_bgrp, intra_bgrp_comm )
  !$acc end host_data
  !
  CALL stop_clock('rotxpsig:diag')
  !
  CALL start_clock('rotxpsig:evc')
  !
  !$acc kernels
  DO i=1, nbnd
     e(i) = en(i)
  END DO
  !$acc end kernels
  !
  ! ... update the basis set
  !
  !$acc kernels
  tpsi = psi
  evc  = (0.D0, 0.D0)
  hevc = (0.D0, 0.D0)
  IF ( overlap ) &
  sevc = (0.D0, 0.D0)
  !$acc end kernels
  !
  IF ( n_start .le. n_end ) THEN
     !
     !$acc host_data use_device(hevc, sevc, hpsi, spsi, evc, tpsi, vr)
     CALL MYDGEMM( 'N', 'N', npw2, nbnd, my_n, 1.D0, &
                 tpsi(1,n_start), npwx2,  vr(n_start,1), nstart, 0.D0, evc,  npwx2 )
     !
     CALL MYDGEMM( 'N', 'N', npw2, nbnd, my_n, 1.D0, &
                 hpsi(1,n_start), npwx2,  vr(n_start,1), nstart, 0.D0, hevc, npwx2 )
     !
     IF ( overlap ) &
     CALL MYDGEMM( 'N', 'N', npw2, nbnd, my_n, 1.D0, &
                 spsi(1,n_start), npwx2,  vr(n_start,1), nstart, 0.D0, sevc, npwx2 )
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
  IF ( overlap ) THEN
     !$acc exit data delete(spsi)
     DEALLOCATE( spsi )
  ENDIF
  !$acc exit data delete(hpsi, en, vr, sr, hr, hpsi, tpsi )
  DEALLOCATE( en, vr, sr, hr, hpsi, tpsi )
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
