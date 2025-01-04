!
! Copyright (C) 2019 National Institute of Advanced Industrial Science and Technology (AIST)
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE rotate_xpsi_k_gpu( h_psi_ptr, s_psi_ptr, overlap, &
                          npwx, npw, nstart, nbnd, npol, psi, evc, hevc, sevc, e )
  !----------------------------------------------------------------------------
  !
  ! ... Serial version of rotate_xpsi for colinear, k-point calculations
  !
#if defined(__CUDA)
  USE cudafor
  USE cublas
#endif
  USE util_param,    ONLY : DP
  USE mp_bands_util, ONLY : intra_bgrp_comm, inter_bgrp_comm, nbgrp, my_bgrp_id, &
                            me_bgrp, root_bgrp
  USE mp,            ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  include 'laxlib.fh'
  !
  ! ... I/O variables
  !
  INTEGER, INTENT(IN) :: npw, npwx, nstart, nbnd, npol
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix psi, as declared in the calling pgm unit
    ! input number of states
    ! output number of states
    ! number of spin polarizations
  LOGICAL, INTENT(IN) :: overlap
    ! if .FALSE. : S|psi> not needed
  COMPLEX(DP), INTENT(INOUT) :: psi(npwx*npol,nstart)
  COMPLEX(DP), INTENT(OUT)   :: evc(npwx*npol,nbnd) !(intent inout?)
    ! input and output eigenvectors (may overlap)
  COMPLEX(DP), INTENT(OUT)   :: hevc(npwx*npol,nbnd), sevc(npwx*npol,nbnd)
    ! H|psi> and S|psi>
  REAL(DP), INTENT(OUT)      :: e(nbnd)
    ! eigenvalues
  !
  ! ... local variables
  !
  INTEGER                  :: kdim, kdmx
  INTEGER                  :: n_start, n_end, my_n, i, j
  COMPLEX(DP), ALLOCATABLE :: hc_d(:,:), sc_d(:,:), vc_d(:,:)
  COMPLEX(DP), ALLOCATABLE :: tpsi_d(:,:), hpsi(:,:), spsi(:,:)
  REAL(DP),    ALLOCATABLE :: en_d(:)
#if defined(__CUDA)
  attributes(DEVICE) :: hc_d, sc_d, vc_d, tpsi_d, en_d
#endif
  !
  EXTERNAL :: h_psi_ptr, s_psi_ptr
  ! h_psi_ptr(npwx,npw,nvec,psi,hpsi)
  !     calculates H|psi>
  ! s_psi_ptr(npwx,npw,nvec,spsi)
  !     calculates S|psi> (if needed)
  !     Vectors psi,hpsi,spsi are dimensioned (npwx,npol,nvec)

  IF ( npol == 1 ) THEN
     !
     kdim = npw
     kdmx = npwx
     !
  ELSE
     !
     kdim = npwx*npol
     kdmx = npwx*npol
     !
  END IF
  !
  !
  ALLOCATE( tpsi_d( kdmx, nstart ) )
  ALLOCATE( hpsi( kdmx, nstart ) )
  !$acc enter data create(hpsi)
  IF ( overlap ) THEN
    ALLOCATE( spsi(kdmx, nstart ) )
    !$acc enter data create(spsi)
  ENDIF
  ALLOCATE( hc_d( nstart, nstart) )    
  ALLOCATE( sc_d( nstart, nstart) )    
  ALLOCATE( vc_d( nstart, nstart) )    
  ALLOCATE( en_d( nstart ) )
  !
  CALL start_clock('rotxpsik')
  !
  ! ... Set up the Hamiltonian and Overlap matrix on the subspace :
  !
  ! ...      H_ij = <psi_i| H |psi_j>     S_ij = <psi_i| S |psi_j>
  !
  CALL start_clock('rotxpsik:hpsi')
  !
  CALL h_psi_ptr( npwx, npw, nstart, psi, hpsi )
  !
  CALL stop_clock('rotxpsik:hpsi')
  !
  IF ( overlap ) THEN
     !
     CALL start_clock('rotxpsik:spsi')
     !
     CALL s_psi_ptr( npwx, npw, nstart, psi, spsi )
     !
     CALL stop_clock('rotxpsik:spsi')
     !
  END IF
  !
  CALL divide(inter_bgrp_comm, nstart, n_start, n_end)
  my_n = n_end - n_start + 1
  !
  CALL start_clock('rotxpsik:hc')
  !
  hc_d = (0.D0, 0.D0)
  !
  !$acc host_data use_device(psi, hpsi)
  IF ( n_start .le. n_end ) &
        CALL ZGEMM( 'C', 'N', nstart, my_n, kdim, (1.D0, 0.D0), psi, &
                    kdmx, hpsi(1,n_start), kdmx, (0.D0, 0.D0), hc_d(1,n_start), nstart )
  !$acc end host_data
  !
  CALL mp_sum( hc_d, inter_bgrp_comm )
  !
  CALL mp_sum( hc_d, intra_bgrp_comm )
  !
  CALL stop_clock('rotxpsik:hc')
  !
  CALL start_clock('rotxpsik:sc')
  !
  sc_d = (0.D0, 0.D0)
  !
  IF ( overlap ) THEN
     !
     !$acc host_data use_device(psi, spsi)
     IF ( n_start .le. n_end ) &
     CALL ZGEMM( 'C', 'N', nstart, my_n, kdim, (1.D0, 0.D0), psi, &
                 kdmx, spsi(1,n_start), kdmx, (0.D0, 0.D0), sc_d(1,n_start), nstart )
     !$acc end host_data
     !
  ELSE
     !
     !$acc host_data use_device(psi)
     IF ( n_start .le. n_end ) &
     CALL ZGEMM( 'C', 'N', nstart, my_n, kdim, (1.D0, 0.D0), psi, &
                 kdmx, psi(1,n_start), kdmx, (0.D0, 0.D0), sc_d(1,n_start), nstart )
     !$acc end host_data
     !
  END IF
  !
  CALL mp_sum( sc_d, inter_bgrp_comm )
  !
  CALL mp_sum( sc_d, intra_bgrp_comm )
  !
  CALL stop_clock('rotxpsik:sc')
  !
  ! ... Diagonalize
  !
  CALL start_clock('rotxpsik:diag')
  !
  CALL diaghg( nstart, nbnd, hc_d, sc_d, nstart, en_d, vc_d, me_bgrp, root_bgrp, intra_bgrp_comm )
  !
  CALL stop_clock('rotxpsik:diag')
  CALL start_clock('rotxpsik:evc')
  !
  !$acc kernels
  DO i=1,nbnd
        e(i) = en_d(i)
  END DO
  !$acc end kernels
  !
  ! ... update the basis set
  !
  !$acc parallel loop collapse(2)
  DO i=1, nbnd
     DO j=1, kdmx
        tpsi_d(j,i) = psi(j,i)
     END DO
  END DO
  !
  !$acc kernels
  evc  = (0.D0, 0.D0)
  hevc = (0.D0, 0.D0)
  !
  IF ( overlap ) sevc = (0.D0, 0.D0)
  !$acc end kernels
  !
  IF ( n_start .le. n_end ) THEN
     !
     !$acc host_data use_device(evc, hpsi, spsi, hevc, sevc)
     CALL ZGEMM( 'N', 'N', kdim, nbnd, my_n, (1.D0, 0.D0), &
                 tpsi_d(1,n_start), kdmx, vc_d(n_start,1), nstart, (0.D0, 0.D0), evc,  kdmx )
     !
     CALL ZGEMM( 'N', 'N', kdim, nbnd, my_n, (1.D0, 0.D0), &
                 hpsi(1,n_start), kdmx, vc_d(n_start,1), nstart, (0.D0, 0.D0), hevc, kdmx )
     !
     IF ( overlap ) &
     CALL ZGEMM( 'N', 'N', kdim, nbnd, my_n, (1.D0, 0.D0), &
                 spsi(1,n_start), kdmx, vc_d(n_start,1), nstart, (0.D0, 0.D0), sevc, kdmx )
     !$acc end host_data
     !
  END IF
  !
  !$acc host_data use_device(evc, hevc, sevc)
  CALL mp_sum( evc,  inter_bgrp_comm )
  CALL mp_sum( hevc, inter_bgrp_comm )
  IF ( overlap ) &
  CALL mp_sum( sevc, inter_bgrp_comm )
  !$acc end host_data
  !
  CALL stop_clock('rotxpsik:evc')
  !
  DEALLOCATE( en_d )
  DEALLOCATE( vc_d )
  DEALLOCATE( sc_d )
  DEALLOCATE( hc_d )
  IF ( overlap ) THEN
    !$acc exit data delete(spsi)
    DEALLOCATE( spsi )
  ENDIF
  !$acc exit data delete(hpsi)
  DEALLOCATE( hpsi )
  DEALLOCATE( tpsi_d )
  !
  CALL stop_clock('rotxpsik')
  !
  !CALL print_clock('rotxpsik')
  !CALL print_clock('rotxpsik:hpsi')
  !CALL print_clock('rotxpsik:spsi')
  !CALL print_clock('rotxpsik:hc')
  !CALL print_clock('rotxpsik:sc')
  !CALL print_clock('rotxpsik:diag')
  !CALL print_clock('rotxpsik:evc')
  !
  RETURN
  !
END SUBROUTINE rotate_xpsi_k_gpu
!
!
!----------------------------------------------------------------------------
