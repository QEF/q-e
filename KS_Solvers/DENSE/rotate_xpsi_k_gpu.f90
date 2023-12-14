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
                          npwx, npw, nstart, nbnd, npol, psi_d, evc_d, hevc_d, sevc_d, e_d )
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
  COMPLEX(DP), INTENT(INOUT) :: psi_d(npwx*npol,nstart)
  COMPLEX(DP), INTENT(OUT)   :: evc_d(npwx*npol,nbnd) !(intent inout?)
    ! input and output eigenvectors (may overlap)
  COMPLEX(DP), INTENT(OUT)   :: hevc_d(npwx*npol,nbnd), sevc_d(npwx*npol,nbnd)
    ! H|psi> and S|psi>
  REAL(DP), INTENT(OUT)      :: e_d(nbnd)
    ! eigenvalues
#if defined(__CUDA)
  attributes(DEVICE) :: psi_d, evc_d, hevc_d, sevc_d, e_d
#endif
  !
  ! ... local variables
  !
  INTEGER                  :: kdim, kdmx
  INTEGER                  :: n_start, n_end, my_n, i, j
  COMPLEX(DP), ALLOCATABLE :: hc_d(:,:), sc_d(:,:), vc_d(:,:)
  COMPLEX(DP), ALLOCATABLE :: tpsi_d(:,:), hpsi_d(:,:), spsi_d(:,:)
  REAL(DP),    ALLOCATABLE :: en_d(:)
#if defined(__CUDA)
  attributes(DEVICE) :: hc_d, sc_d, vc_d, tpsi_d, hpsi_d, spsi_d, en_d
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
  ALLOCATE( hpsi_d( kdmx, nstart ) )
  IF ( overlap ) &
  ALLOCATE( spsi_d(kdmx, nstart ) )
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
  CALL h_psi_ptr( npwx, npw, nstart, psi_d, hpsi_d )
  !
  CALL stop_clock('rotxpsik:hpsi')
  !
  IF ( overlap ) THEN
     !
     CALL start_clock('rotxpsik:spsi')
     !
     CALL s_psi_ptr( npwx, npw, nstart, psi_d, spsi_d )
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
  IF ( n_start .le. n_end ) &
        CALL ZGEMM( 'C', 'N', nstart, my_n, kdim, (1.D0, 0.D0), psi_d, &
                    kdmx, hpsi_d(1,n_start), kdmx, (0.D0, 0.D0), hc_d(1,n_start), nstart )
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
     IF ( n_start .le. n_end ) &
     CALL ZGEMM( 'C', 'N', nstart, my_n, kdim, (1.D0, 0.D0), psi_d, &
                 kdmx, spsi_d(1,n_start), kdmx, (0.D0, 0.D0), sc_d(1,n_start), nstart )
     !
  ELSE
     !
     IF ( n_start .le. n_end ) &
     CALL ZGEMM( 'C', 'N', nstart, my_n, kdim, (1.D0, 0.D0), psi_d, &
                 kdmx, psi_d(1,n_start), kdmx, (0.D0, 0.D0), sc_d(1,n_start), nstart )
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
  !$cuf kernel do(1) <<<*,*>>>
  DO i=1,nbnd
        e_d(i) = en_d(i)
  END DO
  !
  ! ... update the basis set
  !
  !$cuf kernel do(1) <<<*,*>>>
  DO i=1, nbnd
     DO j=1, kdmx
        tpsi_d(j,i) = psi_d(j,i)
     END DO
  END DO
  !
  evc_d  = (0.D0, 0.D0)
  hevc_d = (0.D0, 0.D0)
  !
  IF ( overlap ) sevc_d = (0.D0, 0.D0)
  !
  IF ( n_start .le. n_end ) THEN
     !
     CALL ZGEMM( 'N', 'N', kdim, nbnd, my_n, (1.D0, 0.D0), &
                 tpsi_d(1,n_start), kdmx, vc_d(n_start,1), nstart, (0.D0, 0.D0), evc_d,  kdmx )
     !
     CALL ZGEMM( 'N', 'N', kdim, nbnd, my_n, (1.D0, 0.D0), &
                 hpsi_d(1,n_start), kdmx, vc_d(n_start,1), nstart, (0.D0, 0.D0), hevc_d, kdmx )
     !
     IF ( overlap ) &
     CALL ZGEMM( 'N', 'N', kdim, nbnd, my_n, (1.D0, 0.D0), &
                 spsi_d(1,n_start), kdmx, vc_d(n_start,1), nstart, (0.D0, 0.D0), sevc_d, kdmx )
     !
  END IF
  !
  CALL mp_sum( evc_d,  inter_bgrp_comm )
  CALL mp_sum( hevc_d, inter_bgrp_comm )
  IF ( overlap ) &
  CALL mp_sum( sevc_d, inter_bgrp_comm )
  !
  CALL stop_clock('rotxpsik:evc')
  !
  DEALLOCATE( en_d )
  DEALLOCATE( vc_d )
  DEALLOCATE( sc_d )
  DEALLOCATE( hc_d )
  IF ( overlap ) DEALLOCATE( spsi_d )
  DEALLOCATE( hpsi_d )
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
