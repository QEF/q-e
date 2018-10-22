!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE rotate_wfc_k_gpu( h_psi_gpu, s_psi_gpu, overlap, &
                             npwx, npw, nstart, nbnd, npol, psi_d, evc_d, e_d )
  !----------------------------------------------------------------------------
  !
  ! ... Serial version of rotate_wfc for colinear, k-point calculations
  !
  USE LAXlib
#if defined(__CUDA)
  USE cudafor
  USE cublas
#endif
  USE cg_param,      ONLY : DP
  USE mp_bands_util, ONLY : intra_bgrp_comm, inter_bgrp_comm, root_bgrp_id,&
       nbgrp, my_bgrp_id
  USE mp,            ONLY : mp_sum
  !
  IMPLICIT NONE
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
  COMPLEX(DP) :: psi_d(npwx*npol,nstart), evc_d(npwx*npol,nbnd)
  ! input and output eigenvectors (may overlap)
  REAL(DP)    :: e_d(nbnd)
  ! eigenvalues
#if defined(__CUDA)
  attributes(DEVICE) :: psi_d, evc_d, e_d
#endif
  !
  ! ... local variables
  !
  INTEGER :: kdim, kdmx
  COMPLEX(DP), ALLOCATABLE :: aux_d(:,:)
  COMPLEX(DP), ALLOCATABLE :: hc_d(:,:), sc_d(:,:), vc_d(:,:)
  REAL(DP),    ALLOCATABLE :: en_d(:)
#if defined(__CUDA)
  attributes(DEVICE) :: aux_d, hc_d, sc_d, vc_d, en_d
#endif
  INTEGER :: n_start, n_end, my_n, i, j
  !
  EXTERNAL  h_psi_gpu,    s_psi_gpu
  ! h_psi(npwx,npw,nvec,psi,hpsi)
  !     calculates H|psi>
  ! s_psi(npwx,npw,nvec,spsi)
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
  ALLOCATE( aux_d(kdmx, nstart ) )
  ALLOCATE( hc_d( nstart, nstart) )
  ALLOCATE( sc_d( nstart, nstart) )
  ALLOCATE( vc_d( nstart, nstart) )
  ALLOCATE( en_d( nstart ) )
  call start_clock('rotwfck'); !write(*,*) 'start rotwfck';FLUSH(6)
  !
  ! ... Set up the Hamiltonian and Overlap matrix on the subspace :
  !
  ! ...      H_ij = <psi_i| H |psi_j>     S_ij = <psi_i| S |psi_j>
  !
  call start_clock('rotwfck:hpsi'); !write(*,*) 'start rotwfck:hpsi';FLUSH(6)
  CALL h_psi_gpu( npwx, npw, nstart, psi_d, aux_d )
  call stop_clock('rotwfck:hpsi') ; !write(*,*) 'stop rotwfck:hpsi';FLUSH(6)
  !
  call start_clock('rotwfck:hc'); !write(*,*) 'start rotwfck:hc';FLUSH(6)
  hc_d=(0.D0,0.D0)
  CALL divide(inter_bgrp_comm,nstart,n_start,n_end)
  my_n = n_end - n_start + 1; !write (*,*) nstart,n_start,n_end
  if (n_start .le. n_end) &
       call ZGEMM( 'C','N', nstart, my_n, kdim, (1.D0,0.D0), psi_d, &
                    kdmx, aux_d(1,n_start), kdmx, (0.D0,0.D0), hc_d(1,n_start), nstart )
  CALL mp_sum( hc_d, inter_bgrp_comm )
  !
  CALL mp_sum( hc_d, intra_bgrp_comm )
  !
  sc_d=(0.D0,0.D0)
  IF ( overlap ) THEN
     !
     CALL s_psi_gpu( npwx, npw, nstart, psi_d, aux_d )
     if (n_start .le. n_end) &
          CALL ZGEMM( 'C','N', nstart, my_n, kdim, (1.D0,0.D0), psi_d, &
                      kdmx, aux_d(1,n_start), kdmx, (0.D0,0.D0), sc_d(1,n_start), nstart )
     !
  ELSE
     !
     if (n_start .le. n_end) &
          CALL ZGEMM( 'C','N', nstart, my_n, kdim, (1.D0,0.D0), psi_d, &
                      kdmx, psi_d(1,n_start), kdmx, (0.D0,0.D0), sc_d(1,n_start), nstart )
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
  CALL diaghg( nstart, nbnd, hc_d, sc_d, nstart, en_d, vc_d )
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
  aux_d=(0.D0,0.D0)
  if (n_start .le. n_end) &
       CALL ZGEMM( 'N','N', kdim, nbnd, my_n, (1.D0,0.D0), psi_d(1,n_start), &
                   kdmx, vc_d(n_start,1), nstart, (0.D0,0.D0), aux_d, kdmx )
  CALL mp_sum( aux_d, inter_bgrp_comm )
  !
  !
  !$cuf kernel do(1) <<<*,*>>>
  DO i=1, nbnd
     DO j=1, kdmx
       evc_d(j,i) = aux_d(j,i)
     END DO
  END DO
  !
  call stop_clock('rotwfck:evc') ! ; write(*,*) 'stop rotwfck;evc';FLUSH(6)
  !
  DEALLOCATE( en_d )
  DEALLOCATE( vc_d )
  DEALLOCATE( sc_d )
  DEALLOCATE( hc_d )
  DEALLOCATE( aux_d )
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
