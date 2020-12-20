!
! Copyright (C) 2003-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
#define ZERO ( 0.D0, 0.D0 )
!----------------------------------------------------------------------------
SUBROUTINE rotate_HSpsi_k_gpu( npwx, npw, nstart, nbnd, npol, psi_d, hpsi_d, overlap, spsi_d, e_d )
  !----------------------------------------------------------------------------
  !
  ! ... Serial version of rotate_wfc for colinear, k-point calculations
  !
#if defined (__CUDA)
  USE cudafor
#endif
  USE util_param,    ONLY : DP
  USE mp_bands_util, ONLY : intra_bgrp_comm, inter_bgrp_comm, root_bgrp_id, nbgrp, my_bgrp_id, &
                            me_bgrp, root_bgrp
  USE mp,            ONLY : mp_sum, mp_barrier, mp_allgather, mp_type_create_column_section, mp_type_free
  USE device_memcpy_m,    ONLY: dev_memcpy, dev_memset
  !
  IMPLICIT NONE
  !
  INCLUDE 'laxlib.fh'
  !
  ! ... I/O variables
  !
  INTEGER, INTENT(IN) :: &
    npw,                 &         ! dimension of the matrices (psi,Hpsi,Spsi) to be rotated
    npwx,                &         ! leading dimension of the wavefunction-related matrices
    nstart,              &         ! input number of states 
    nbnd,                &         ! output number of states
    npol                           ! number of spin polarizations
  COMPLEX(DP), INTENT(INOUT) :: psi_d(npwx*npol,nstart), hpsi_d(npwx*npol,nstart) ! input and output psi, Hpsi,
  COMPLEX(DP), INTENT(INOUT), OPTIONAL :: spsi_d(npwx*npol,nstart)              ! ...   and optionnally Spsi
  LOGICAL, INTENT(IN) :: overlap   ! if .FALSE. : spsi is not needed (and not used)
  REAL(DP), INTENT(OUT) :: e_d(nbnd) ! eigenvalues of the reduced H matrix
  !
  ! ... local variables
  !
  INTEGER :: kdim, kdmx
  INTEGER :: n_start, n_end, my_n, recv_counts(nbgrp), displs(nbgrp), column_type
  INTEGER :: ii, jj  ! indexes for cuf kernel loops 
  !
  ! ... device variables
  !
  COMPLEX(DP), ALLOCATABLE :: aux_d(:,:)
  COMPLEX(DP), ALLOCATABLE :: hh_d(:,:), ss_d(:,:), vv_d(:,:)
  REAL(DP),    ALLOCATABLE :: en_d(:)
#if defined (__CUDA)
  attributes(device) :: aux_d, psi_d, hpsi_d, spsi_d 
  attributes(device) :: hh_d, ss_d, vv_d, en_d, e_d
#endif
  !
  IF ( overlap .AND..NOT.present(spsi_d) ) call errore( 'rotHSw','spsi_d array needed with overlap=.TRUE.',1)
  !
  call start_clock('rotHSw'); !write(*,*) 'start rotHSw' ; FLUSH(6)
  !
  if (npol == 2 .and. npw < npwx ) then ! pack wfcs so that pw's are contiguous 
    call start_clock('rotHSw:move'); !write(*,*) 'start rotHSw:move' ; FLUSH(6)
    ALLOCATE ( aux_d ( npwx, nstart ) )
!$cuf kernel do(2)
    DO ii = 1, npw
      DO jj = 1, nstart 
        aux_d(ii,jj) = psi_d (npwx+ii,jj)
      END DO 
    END DO 
!$cuf kernel do(2)
    DO ii = 1, npw
      DO jj = 1, nstart 
        psi_d (npw+ii,jj) = aux_d(ii,jj)
      END DO 
    END DO 
!$cuf kernel do(2)
    DO ii = 1, npw
      DO jj = 1, nstart 
        aux_d(ii,jj) = hpsi_d(npwx+ii,jj)
      END DO 
    END DO 
!$cuf kernel do(2)
    DO ii = 1, npw
      DO jj = 1, nstart 
        hpsi_d(npw+ii,jj) = aux_d(ii,jj)
      END DO 
    END DO 
!$cuf kernel do(2)
    DO ii = 1, npw
      DO jj = 1, nstart 
        aux_d(ii,jj) = spsi_d(npwx+ii,jj)
      END DO 
    END DO 
!$cuf kernel do(2)
    DO ii = 1, npw
      DO jj = 1, nstart 
        spsi_d(npw+ii,jj) = aux_d(ii,jj)
      END DO 
    END DO 
    DEALLOCATE( aux_d )
    call stop_clock('rotHSw:move'); !write(*,*) 'stop rotHSw:move' ; FLUSH(6)
  end if

  kdim = npol * npw 
  kdmx = npol * npwx
  !
  ALLOCATE( hh_d( nstart, nstart ) )
  ALLOCATE( ss_d( nstart, nstart ) )
  ALLOCATE( vv_d( nstart, nstart ) )
  ALLOCATE( en_d( nstart ) )
  !
  ! ... Set up the Hamiltonian and Overlap matrix on the subspace :
  !
  ! ...      H_ij = <psi_i| H |psi_j>     S_ij = <psi_i| S |psi_j>
  !
  call start_clock('rotHSw:hc'); !write(*,*) 'start rotHSw:hc' ; FLUSH(6)

  CALL mp_type_create_column_section(hh_d(1,1), 0, nstart, nstart, column_type)
  CALL divide_all(inter_bgrp_comm,nstart,n_start,n_end,recv_counts,displs)

  my_n = n_end - n_start + 1; !write (*,*) nstart,n_start,n_end
  if (n_start .le. n_end) &
    CALL gpu_ZGEMM( 'C','N', nstart, my_n, kdim, (1.D0,0.D0), psi_d, kdmx, hpsi_d(1,n_start), kdmx, (0.D0,0.D0), &
                                                                                         hh_d(1,n_start), nstart )
  call start_clock('rotHSw:hc:s1')
  CALL mp_sum( hh_d(:,n_start:n_end), intra_bgrp_comm ) ! this section only needs to be collected inside bgrp
  call stop_clock('rotHSw:hc:s1')
  !call start_clock('rotHSw:hc:b1'); CALL mp_barrier( inter_bgrp_comm ); call stop_clock('rotHSw:hc:b1')
  call start_clock('rotHSw:hc:s2')
  CALL mp_allgather(hh_d, column_type, recv_counts, displs, inter_bgrp_comm)
  call stop_clock('rotHSw:hc:s2')
  !
  IF ( overlap ) THEN
     !
     if (n_start .le. n_end) &
       CALL gpu_ZGEMM( 'C','N', nstart, my_n, kdim, (1.D0,0.D0), psi_d, kdmx, spsi_d(1,n_start), kdmx, &
                                                                          (0.D0,0.D0), ss_d(1,n_start), nstart )
     !
  ELSE
     !
     if (n_start .le. n_end) &
       CALL gpu_ZGEMM( 'C','N', nstart, my_n, kdim, (1.D0,0.D0), psi_d, kdmx, psi_d(1,n_start), kdmx, &
                                                                           (0.D0,0.D0), ss_d(1,n_start), nstart )
     !
  END IF
  call start_clock('rotHSw:hc:s3')
  CALL mp_sum( ss_d(:,n_start:n_end), intra_bgrp_comm ) ! this section only needs to be collected inside bgrp
  call stop_clock('rotHSw:hc:s3')
  !call start_clock('rotHSw:hc:b2'); CALL mp_barrier( inter_bgrp_comm ); call stop_clock('rotHSw:hc:b2')
  call start_clock('rotHSw:hc:s4')
  CALL mp_allgather(ss_d, column_type, recv_counts, displs, inter_bgrp_comm)
  call stop_clock('rotHSw:hc:s4')

  CALL mp_type_free( column_type )

  call stop_clock('rotHSw:hc'); !write(*,*) 'stop rotHSw:hc' ; FLUSH(6)
  !
  ! ... Diagonalize
  !
  call start_clock('rotHSw:diag'); !write(*,*) 'start rotHSw:diag' ; FLUSH(6)
  CALL diaghg( nstart, nbnd, hh_d, ss_d, nstart, en_d, vv_d, me_bgrp, root_bgrp, intra_bgrp_comm )
  CALL dev_memcpy(e_d, en_d, [1,nbnd])
  call stop_clock('rotHSw:diag'); !write(*,*) 'stop rotHSw:diag' ; FLUSH(6)
  !
  ! ... update the basis set
  !
  call start_clock('rotHSw:evc'); !write(*,*) 'start rotHSw:evc' ; FLUSH(6)

  CALL mp_type_create_column_section(psi_d(1,1), 0, npwx, npwx, column_type)
  CALL divide_all(inter_bgrp_comm,nbnd,n_start,n_end,recv_counts,displs)

  ALLOCATE( aux_d ( kdmx, nbnd ) )

  my_n = n_end - n_start + 1; !write (*,*) nstart,n_start,n_end
  if (n_start .le. n_end) &
    CALL gpu_ZGEMM( 'N','N', kdim, my_n, nstart, (1.D0,0.D0), psi_d, kdmx, vv_d(1,n_start), nstart, &
                                                                      (0.D0,0.D0), aux_d(1,n_start), kdmx )

  CALL dev_memcpy(psi_d, aux_d, [1, kdmx], 1, [n_start,n_end])
  !call start_clock('rotHSw:ev:b3'); CALL mp_barrier( inter_bgrp_comm ); call stop_clock('rotHSw:ev:b3')
  call start_clock('rotHSw:ev:s5')
  CALL mp_allgather(psi_d(:,1:nbnd), column_type, recv_counts, displs, inter_bgrp_comm)
  call stop_clock('rotHSw:ev:s5')

  if (n_start .le. n_end) &
    CALL gpu_ZGEMM( 'N','N', kdim, my_n, nstart, (1.D0,0.D0), hpsi_d, kdmx, vv_d(1,n_start), nstart, &
                                                                       (0.D0,0.D0), aux_d(1,n_start), kdmx )
  CALL dev_memcpy(hpsi_d, aux_d, [1, kdmx], 1, [n_start,n_end])  !call start_clock('rotHSw:ev:b4'); CALL mp_barrier( inter_bgrp_comm ); call stop_clock('rotHSw:ev:b4')
  call start_clock('rotHSw:ev:s6')
  CALL mp_allgather(hpsi_d(:,1:nbnd), column_type, recv_counts, displs, inter_bgrp_comm)
  call stop_clock('rotHSw:ev:s6')

  IF (overlap) THEN
     if (n_start .le. n_end) &
       CALL gpu_ZGEMM( 'N','N', kdim, my_n, nstart, (1.D0,0.D0), spsi_d, kdmx, vv_d(1,n_start), &
                                                          nstart, (0.D0,0.D0), aux_d(1,n_start), kdmx )
     CALL dev_memcpy(spsi_d, aux_d, [1, kdmx], 1, [n_start,n_end])     !call start_clock('rotHSw:ev:b5'); CALL mp_barrier( inter_bgrp_comm ); call stop_clock('rotHSw:ev:b5')
     call start_clock('rotHSw:ev:s7')
     CALL mp_allgather(spsi_d(:,1:nbnd), column_type, recv_counts, displs, inter_bgrp_comm)
     call stop_clock('rotHSw:ev:s7')

  ELSE IF (present(spsi_d)) THEN

     CALL dev_memcpy(spsi_d, psi_d, [1, kdmx], 1, [n_start,n_end])
  END IF

  DEALLOCATE( aux_d )
  CALL mp_type_free( column_type )
  call stop_clock('rotHSw:evc'); !write(*,*) 'stop rotHSw:evc' ; FLUSH(6)
  !
  DEALLOCATE( vv_d )
  DEALLOCATE( ss_d )
  DEALLOCATE( hh_d )
  DEALLOCATE( en_d )
  !
  call stop_clock('rotHSw'); !write(*,*) 'stop rotHSw' ; FLUSH(6)
  !call print_clock('rotHSw')
  !call print_clock('rotHSw:hc')
  !call print_clock('rotHSw:diag')
  !call print_clock('rotHSw:evc')
  !
  if (npol== 2 .and. npw < npwx ) then   ! unpack wfcs to its original order
    call start_clock('rotHSw:move'); !write(*,*) 'start rotHSw:move' ; FLUSH(6)
    ALLOCATE ( aux_d ( npwx, nbnd ) )
!$cuf kernel do(2)
    DO ii = 1, npw
      DO jj = 1, nbnd 
        aux_d(ii,jj) = psi_d (npw+ii,jj)
      END DO 
    END DO 
!$cuf kernel do(2)
    DO ii = 1, npw
      DO jj = 1, nbnd 
        psi_d(npwx+ii,jj) = aux_d(ii,jj)
      END DO 
    END DO 
!$cuf kernel do(2)
    DO ii = 1, npw
      DO jj = 1, nbnd 
        aux_d(ii,jj) = hpsi_d(npw+ii,jj)
      END DO 
    END DO 
!$cuf kernel do(2)
    DO ii = 1, npw
      DO jj = 1, nbnd 
        hpsi_d(npwx+ii,jj) = aux_d(ii,jj)
      END DO 
    END DO 
!$cuf kernel do(2)
    DO ii = 1, npw
      DO jj = 1, nbnd 
        aux_d(ii,jj) = spsi_d(npw+ii,jj)
      END DO 
    END DO 
!$cuf kernel do(2)
    DO ii = 1, npw
      DO jj = 1, nbnd 
        spsi_d(npwx+ii,jj) = aux_d(ii,jj)
      END DO 
    END DO 
!$cuf kernel do(2)
    DO ii = npw+1, npwx
      DO jj = 1, nbnd 
        psi_d(ii,jj) = ZERO
      END DO 
    END DO 
!$cuf kernel do(2)
    DO ii = npw+1, npwx
      DO jj = 1, nbnd 
        hpsi_d(ii,jj) = ZERO 
      END DO 
    END DO 
!$cuf kernel do(2)
    DO ii = npw+1, npwx
      DO jj = 1, nbnd 
        spsi_d(ii,jj) = ZERO
      END DO 
    END DO 
    DEALLOCATE( aux_d )
    call stop_clock('rotHSw:move'); !write(*,*) 'stop rotHSw:move' ; FLUSH(6)
  end if
  !
  RETURN
  !
END SUBROUTINE rotate_HSpsi_k_gpu
