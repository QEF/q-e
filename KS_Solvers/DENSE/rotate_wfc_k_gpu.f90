!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE rotate_wfc_k_gpu( h_psi_ptr, s_psi_ptr, overlap, &
                             npwx, npw, nstart, nbnd, npol, psi, evc_d, e_d )
  !----------------------------------------------------------------------------
  !
  ! ... Serial version of rotate_wfc for colinear, k-point calculations
  !
#if defined(__CUDA)
  USE cudafor
  USE cublas
#endif
  USE util_param,    ONLY : DP
  USE mp_bands_util, ONLY : intra_bgrp_comm, inter_bgrp_comm, root_bgrp_id,&
          nbgrp, my_bgrp_id, me_bgrp, root_bgrp
  USE mp,            ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INCLUDE 'laxlib.fh'
  !
  ! ... I/O variables
  !
  INTEGER, INTENT(IN) :: npw, npwx, nstart, nbnd, npol
  ! dimension of the matrix to be diagonalized
  ! leading dimension of matrix psi, as declared in the calling pgm unit
  ! input number of states
  ! output number of states
  ! number of spin polarizations
  LOGICAL :: overlap
  ! if .FALSE. : S|psi> not needed
  COMPLEX(DP) :: psi(npwx*npol,nstart), evc_d(npwx*npol,nbnd)
  ! input and output eigenvectors (may overlap)
  REAL(DP)    :: e_d(nbnd)
  ! eigenvalues
#if defined(__CUDA)
  attributes(DEVICE) :: evc_d, e_d
#endif
  !
  ! ... local variables
  !
  INTEGER :: kdim, kdmx
  COMPLEX(DP), ALLOCATABLE :: aux(:,:)
  COMPLEX(DP), ALLOCATABLE :: hc_d(:,:), sc_d(:,:), vc_d(:,:)
  REAL(DP),    ALLOCATABLE :: en_d(:)
#if defined(__CUDA)
  attributes(DEVICE) :: hc_d, sc_d, vc_d, en_d
#endif
  INTEGER :: n_start, n_end, my_n, i, j
  !
  EXTERNAL  h_psi_ptr,    s_psi_ptr
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
  ALLOCATE( aux(kdmx, nstart ) )
  ALLOCATE( hc_d( nstart, nstart) )
  ALLOCATE( sc_d( nstart, nstart) )
  ALLOCATE( vc_d( nstart, nstart) )
  ALLOCATE( en_d( nstart ) )
  !$acc enter data create(aux)
  !
  call start_clock('rotwfck'); !write(*,*) 'start rotwfck';FLUSH(6)
  !
  ! ... Set up the Hamiltonian and Overlap matrix on the subspace :
  !
  ! ...      H_ij = <psi_i| H |psi_j>     S_ij = <psi_i| S |psi_j>
  !
  call start_clock('rotwfck:hpsi'); !write(*,*) 'start rotwfck:hpsi';FLUSH(6)
  CALL h_psi_ptr( npwx, npw, nstart, psi, aux )
  call stop_clock('rotwfck:hpsi') ; !write(*,*) 'stop rotwfck:hpsi';FLUSH(6)
  !
  call start_clock('rotwfck:hc'); !write(*,*) 'start rotwfck:hc';FLUSH(6)
  hc_d=(0.D0,0.D0)
  CALL divide(inter_bgrp_comm,nstart,n_start,n_end)
  my_n = n_end - n_start + 1; !write (*,*) nstart,n_start,n_end
  IF (n_start .le. n_end) THEN
  !$acc host_data use_device(psi, aux)
    call ZGEMM( 'C','N', nstart, my_n, kdim, (1.D0,0.D0), psi, &
                 kdmx, aux(1,n_start), kdmx, (0.D0,0.D0), hc_d(1,n_start), nstart )
  !$acc end host_data
  END IF
  CALL mp_sum( hc_d, inter_bgrp_comm )
  !
  CALL mp_sum( hc_d, intra_bgrp_comm )
  !
  sc_d=(0.D0,0.D0)
  IF ( overlap ) THEN
     !
     CALL s_psi_ptr( npwx, npw, nstart, psi, aux )
     IF (n_start .le. n_end) THEN
       !$acc host_data use_device(psi, aux)
       CALL ZGEMM( 'C','N', nstart, my_n, kdim, (1.D0,0.D0), psi, &
                   kdmx, aux(1,n_start), kdmx, (0.D0,0.D0), sc_d(1,n_start), nstart )
       !$acc end host_data
     END IF
     !
  ELSE
     !
     IF (n_start .le. n_end) THEN
       !$acc host_data use_device(psi)
       CALL ZGEMM( 'C','N', nstart, my_n, kdim, (1.D0,0.D0), psi, &
                   kdmx, psi(1,n_start), kdmx, (0.D0,0.D0), sc_d(1,n_start), nstart )
       !$acc end host_data
     END IF
     !
  END IF
  CALL mp_sum( sc_d, inter_bgrp_comm )
  !
  CALL mp_sum( sc_d, intra_bgrp_comm )
  call stop_clock('rotwfck:hc'); !write(*,*) 'stop rotwfck:hc';FLUSH(6)
  !
  ! ... Diagonalize
  !
  call start_clock('rotwfck:diag');  !write(*,*) 'start rotwfck:diag';FLUSH(6)
  !! cdiaghg on GPU. See interface from LAXlib module
  CALL diaghg( nstart, nbnd, hc_d, sc_d, nstart, en_d, vc_d, me_bgrp, root_bgrp, intra_bgrp_comm )
  call stop_clock('rotwfck:diag');  !write(*,*) 'stop rotwfck:diag';FLUSH(6)
  call start_clock('rotwfck:evc'); !write(*,*) 'start rotwfck:evc';FLUSH(6)
  !
  !$cuf kernel do(1) <<<*,*>>>
  DO i=1,nbnd
    e_d(i) = en_d(i)
  END DO
  !
  ! ...  update the basis set
  !
  !$acc kernels
  aux=(0.D0,0.D0)
  !$acc end kernels
  IF (n_start .le. n_end) THEN
    !$acc host_data use_device(psi, aux)
    CALL ZGEMM( 'N','N', kdim, nbnd, my_n, (1.D0,0.D0), psi(1,n_start), &
                kdmx, vc_d(n_start,1), nstart, (0.D0,0.D0), aux, kdmx )
    !$acc end host_data
  END IF
  !$acc host_data use_device(aux)
  CALL mp_sum( aux, inter_bgrp_comm )
  !$acc end host_data
  !
  !
  !$acc parallel loop collapse(2) 
  DO i=1, nbnd
     DO j=1, kdmx
       evc_d(j,i) = aux(j,i)
     END DO
  END DO
  !$acc end parallel
  !
  call stop_clock('rotwfck:evc') ! ; write(*,*) 'stop rotwfck;evc';FLUSH(6)
  !
  !$acc exit data delete(aux)
  DEALLOCATE( en_d )
  DEALLOCATE( vc_d )
  DEALLOCATE( sc_d )
  DEALLOCATE( hc_d )
  DEALLOCATE( aux )
  call stop_clock('rotwfck'); !write(*,*) 'stop rotwfck';FLUSH(6)
  !call print_clock('rotwfck')
  !call print_clock('rotwfck:hpsi')
  !call print_clock('rotwfck:hc')
  !call print_clock('rotwfck:diag')
  !call print_clock('rotwfck:evc')
  !
  RETURN
  !
END SUBROUTINE rotate_wfc_k_gpu
!
