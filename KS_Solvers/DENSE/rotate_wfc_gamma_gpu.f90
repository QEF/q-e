!
! Copyright (C) 2003-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE rotate_wfc_gamma_gpu( h_psi_ptr, s_psi_ptr, overlap, &
                             npwx, npw, nstart, nbnd, psi, evc_d, e_d )
  !----------------------------------------------------------------------------
  !
  ! ... Serial version of rotate_wfc for Gamma-only calculations
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
  USE mp_bands_util, ONLY : intra_bgrp_comm, inter_bgrp_comm, root_bgrp_id, &
          nbgrp, my_bgrp_id, me_bgrp, root_bgrp
  USE mp_bands_util, ONLY : gstart ! index of the first nonzero G 
  USE mp,            ONLY : mp_sum 
  !
  IMPLICIT NONE
  !
  INCLUDE 'laxlib.fh'
  !
  ! ... I/O variables
  !
  INTEGER :: npw, npwx, nstart, nbnd
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix psi, as declared in the calling pgm unit
    ! input number of states
    ! output number of states
  LOGICAL :: overlap
    ! if .FALSE. : S|psi> not needed
  COMPLEX(DP) :: psi(npwx,nstart), evc_d(npwx,nbnd)
    ! input and output eigenvectors (may overlap)
  REAL(DP) :: e_d(nbnd)
    ! eigenvalues
#if defined(__CUDA)
  attributes(DEVICE) :: evc_d, e_d
#endif
  !
  ! ... local variables
  !
  INTEGER                  :: npw2, npwx2
  COMPLEX(DP), ALLOCATABLE :: aux(:,:)
  REAL(DP),    ALLOCATABLE :: hr_d(:,:), sr_d(:,:), vr_d(:,:)
  REAL(DP),    ALLOCATABLE :: en_d(:)
#if defined(__CUDA)
  attributes(DEVICE) :: hr_d, sr_d, vr_d, en_d
#endif
  INTEGER :: n_start, n_end, my_n, i, j
  !
  EXTERNAL  h_psi_ptr,    s_psi_ptr
    ! h_psi_ptr(npwx,npw,nvec,psi,hpsi)
    !     calculates H|psi>
    ! s_psi_ptr(npwx,npw,nvec,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx,npol,nvec)

  npw2  = 2 * npw
  npwx2 = 2 * npwx

  IF ( gstart == -1 ) CALL errore( 'regter', 'gstart variable not initialized', 1 )

  !
  ALLOCATE( aux(  npwx, nstart ) )
  ALLOCATE( hr_d( nstart, nstart ) )
  ALLOCATE( sr_d( nstart, nstart ) )
  ALLOCATE( vr_d( nstart, nstart ) )
  ALLOCATE( en_d( nstart ) )
  !$acc enter data create(aux)
  !
  call start_clock('rotwfcg'); !write(*,*) 'start rotwfcg' ; FLUSH(6)
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
        psi(1,i) = CMPLX( DBLE( psi(1,i) ), 0.D0,kind=DP)
     END DO
     !$acc end kernels
  END IF
  !
  call start_clock('rotwfcg:hpsi'); !write(*,*) 'start rotwfcg:hpsi' ; FLUSH(6)
  CALL h_psi_ptr( npwx, npw, nstart, psi, aux )
  call stop_clock('rotwfcg:hpsi'); !write(*,*) 'stop rotwfcg:hpsi' ; FLUSH(6)
  !
  call start_clock('rotwfcg:hc'); !write(*,*) 'start rotwfcg:hc' ; FLUSH(6)
  hr_d=0.D0
  CALL divide(inter_bgrp_comm,nstart,n_start,n_end)
  my_n = n_end - n_start + 1; !write (*,*) nstart,n_start,n_end
  !$acc host_data use_device(psi, aux)
  if (n_start .le. n_end) &
  CALL cublasDGEMM( 'T','N', nstart, my_n, npw2, 2.D0, psi, &
                    npwx2, aux(1,n_start), npwx2, 0.D0, hr_d(1,n_start), nstart )
  IF ( gstart == 2 ) call MYDGER( nstart, my_n, -1.D0, psi, &
                                      npwx2, aux(1,n_start), npwx2, hr_d(1,n_start), nstart )
  !$acc end host_data
  CALL mp_sum( hr_d, inter_bgrp_comm )
  !
  CALL mp_sum( hr_d, intra_bgrp_comm )
  !
  sr_d=0.D0
  IF ( overlap ) THEN
     !
     CALL s_psi_ptr( npwx, npw, nstart, psi, aux )
     !
     !$acc host_data use_device(psi, aux)
     if (n_start .le. n_end) &
     CALL cublasDGEMM( 'T','N', nstart, my_n, npw2, 2.D0, psi, &
                       npwx2, aux(1,n_start), npwx2, 0.D0, sr_d(1,n_start), nstart )
     IF ( gstart == 2 ) CALL MYDGER( nstart, my_n, -1.D0, psi, &
                                         npwx2, aux(1,n_start), npwx2, sr_d(1,n_start), nstart )
     !$acc end host_data
     !
  ELSE
     !
     !$acc host_data use_device(psi)
     if (n_start .le. n_end) &
     CALL cublasDGEMM( 'T','N', nstart, my_n, npw2, 2.D0, psi, &
                       npwx2, psi(1,n_start), npwx2, 0.D0, sr_d(1,n_start), nstart )
     IF ( gstart == 2 ) CALL MYDGER( nstart, my_n, -1.D0, psi, &
                                         npwx2, psi(1,n_start), npwx2, sr_d(1,n_start), nstart )
     !$acc end host_data
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
  CALL diaghg( nstart, nbnd, hr_d, sr_d, nstart, en_d, vr_d, me_bgrp, root_bgrp, intra_bgrp_comm )
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
  !$acc kernels
  aux=(0.D0,0.D0)
  !$acc end kernels
  !$acc host_data use_device(psi, aux)
  if (n_start .le. n_end) &
  CALL cublasDGEMM( 'N','N', npw2, nbnd, my_n, 1.D0, psi(1,n_start), &
                     npwx2, vr_d(n_start,1), nstart, 0.D0, aux, npwx2 )
  CALL mp_sum( aux, inter_bgrp_comm )
  !$acc end host_data
  !
  !$acc parallel loop collapse(2)
  DO i=1, nbnd
     DO j=1, npwx
        evc_d(j,i) = aux(j,i)
     END DO
  END DO
  !$acc end parallel 
  call stop_clock('rotwfcg:evc_d'); !write(*,*) 'stop rotwfcg:evc_d' ; FLUSH(6)
  !
  !$acc exit data delete(aux)
  DEALLOCATE( en_d )
  DEALLOCATE( vr_d )
  DEALLOCATE( sr_d )
  DEALLOCATE( hr_d )
  DEALLOCATE( aux )
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
