!
! Copyright (C) 2019 National Institute of Advanced Industrial Science and Technology (AIST)
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE rotate_xpsi_gamma( h_psi_ptr, s_psi_ptr, overlap, &
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
  INTEGER                  :: n_start, n_end, my_n
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
  CALL start_clock('rotxpsig')
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
  ! ... Set up the Hamiltonian and Overlap matrix on the subspace :
  !
  ! ...      H_ij = <psi_i| H |psi_j>     S_ij = <psi_i| S |psi_j>
  !
  ! ... set Im[ psi(G=0) ] -  needed for numerical stability
  !
  IF ( gstart == 2 ) THEN
     !$acc kernels      
     psi(1,1:nstart) = CMPLX( DBLE( psi(1,1:nstart) ), 0.D0, kind=DP)
     !$acc end kernels
  END IF
  !
  CALL start_clock('rotxpsig:hpsi')
  !
  CALL h_psi_ptr( npwx, npw, nstart, psi, hpsi )
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
  CALL divide(inter_bgrp_comm, nstart, n_start, n_end)
  my_n = n_end - n_start + 1
  !
  CALL start_clock('rotxpsig:hc')
  !
  !$acc kernels
  hr = 0.D0
  !$acc end kernels
  !
  !$acc host_data use_device(psi, hpsi, hr)
  IF ( n_start <= n_end ) &
  CALL MYDGEMM( 'T', 'N', nstart, my_n, npw2, 2.D0, &
              psi, npwx2, hpsi(1,n_start), npwx2, 0.D0, hr(1,n_start), nstart )
  !
  IF ( gstart == 2 ) &
  CALL MYDGER( nstart, my_n, -1.D0, psi, npwx2, hpsi(1,n_start), npwx2, hr(1,n_start), nstart )
  !
  CALL mp_sum( hr , inter_bgrp_comm )
  !
  CALL mp_sum( hr , intra_bgrp_comm )
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
     IF ( n_start <= n_end ) &
     CALL MYDGEMM( 'T', 'N', nstart, my_n, npw2, 2.D0, &
                 psi, npwx2, spsi(1,n_start), npwx2, 0.D0, sr(1,n_start), nstart )
     !
     IF ( gstart == 2 ) &
     CALL MYDGER( nstart, my_n, -1.D0, psi, npwx2, spsi(1,n_start), npwx2, sr(1,n_start), nstart )
     !              
  ELSE
     !
     IF ( n_start <= n_end ) &
     CALL MYDGEMM( 'T', 'N', nstart, my_n, npw2, 2.D0, &
                 psi, npwx2, psi(1,n_start), npwx2, 0.D0, sr(1,n_start), nstart )
     !
     IF ( gstart == 2 ) &
     CALL MYDGER( nstart, my_n, -1.D0, psi, npwx2, psi(1,n_start), npwx2, sr(1,n_start), nstart )
     !
  END IF
  !
  CALL mp_sum( sr , inter_bgrp_comm )
  !
  CALL mp_sum( sr , intra_bgrp_comm )
  !$acc end host_data
  !
  CALL stop_clock('rotxpsig:sc')
  !
  ! ... Diagonalize
  !
  CALL start_clock('rotxpsig:diag')
  !
  !$acc host_data use_device(hr, sr, vr, en)
  CALL diaghg( nstart, nbnd, hr, sr, nstart, en, vr, me_bgrp, root_bgrp, intra_bgrp_comm )
  !$acc end host_data
  !
  !$acc kernels
  e(:) = en(1:nbnd)
  !$acc end kernels
  !
  CALL stop_clock('rotxpsig:diag')
  !
  ! ... update the basis set
  !
  CALL start_clock('rotxpsig:evc')
  !
  !$acc kernels
  tpsi = psi
  evc  = (0.D0, 0.D0)
  hevc = (0.D0, 0.D0)
  IF ( overlap ) &
  sevc = (0.D0, 0.D0)
  !$acc end kernels
  !
  IF ( n_start <= n_end ) THEN
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
  !$acc exit data delete(hpsi, vr, sr, hr, tpsi, en )
  DEALLOCATE( vr, sr, hr, hpsi, tpsi, en )
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
END SUBROUTINE rotate_xpsi_gamma
!
!
!----------------------------------------------------------------------------
SUBROUTINE protate_xpsi_gamma( h_psi_ptr, s_psi_ptr, overlap, &
                               npwx, npw, nstart, nbnd, psi, evc, hevc, sevc, e )
  !----------------------------------------------------------------------------
  !
  ! ... Parallel version of rotate_xpsi for Gamma-only calculations
  ! ... Subroutine with distributed matrices, written by Carlo Cavazzoni
  ! ... This version assumes real wavefunctions (k=0) with only
  ! ... half plane waves stored: psi(-G)=psi*(G), except G=0
  !
  USE util_param,       ONLY : DP
  USE mp_bands_util,    ONLY : intra_bgrp_comm, inter_bgrp_comm, &
                               nbgrp, root_bgrp_id, my_bgrp_id
  USE mp_bands_util,    ONLY : gstart ! index of the first nonzero G 
  USE mp,               ONLY : mp_bcast, mp_root_sum, mp_sum, mp_barrier
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
  ! ... local variables
  !
  INTEGER                  :: npw2, npwx2
  REAL(DP),    ALLOCATABLE :: hr(:,:), sr(:,:), vr(:,:)
  COMPLEX(DP), ALLOCATABLE :: tpsi(:,:), hpsi(:,:), spsi(:,:)
  REAL(DP),    ALLOCATABLE :: en(:)
  !
  INTEGER :: idesc(LAX_DESC_SIZE)
    ! matrix distribution descriptors
  INTEGER :: nx
    ! maximum local block dimension
  LOGICAL :: la_proc
    ! flag to distinguish procs involved in linear algebra
  LOGICAL :: do_distr_diag_inside_bgrp
  INTEGER :: ortho_parent_comm
  INTEGER, ALLOCATABLE :: idesc_ip( :, :, : )
  INTEGER, ALLOCATABLE :: rank_ip( :, : )
  !
  EXTERNAL :: h_psi_ptr, s_psi_ptr
    ! h_psi_ptr(npwx,npw,nvec,psi,hpsi)
    !     calculates H|psi>
    ! s_psi_ptr(npwx,npw,nvec,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx,npol,nvec)

  CALL start_clock('protxpsig')
  !
  CALL laxlib_getval( do_distr_diag_inside_bgrp = do_distr_diag_inside_bgrp, &
       ortho_parent_comm = ortho_parent_comm )
  CALL desc_init( nstart, nx, la_proc, idesc, rank_ip, idesc_ip )
  !
  npw2  = 2 * npw
  npwx2 = 2 * npwx
  !
  IF ( gstart == -1 ) CALL errore( 'protxpsig', 'gstart variable not initialized', 1 )
  !
  ALLOCATE( tpsi( npwx, nstart ) )
  ALLOCATE( hpsi( npwx, nstart ) )
  IF ( overlap ) &
  ALLOCATE( spsi( npwx, nstart ) )
  ALLOCATE( hr( nx, nx ) )    
  ALLOCATE( sr( nx, nx ) )    
  ALLOCATE( vr( nx, nx ) )    
  ALLOCATE( en( nstart ) )
  !
  ! ... Set up the Hamiltonian and Overlap matrix on the subspace :
  !
  ! ...      H_ij = <psi_i| H |psi_j>     S_ij = <psi_i| S |psi_j>
  !
  ! ... set Im[ psi(G=0) ] -  needed for numerical stability
  !
  IF ( gstart == 2 ) &
     psi(1,1:nstart) = CMPLX( DBLE( psi(1,1:nstart) ), 0.D0, kind=DP)
  !
  CALL start_clock('protxpsig:hpsi')
  !
  CALL h_psi_ptr( npwx, npw, nstart, psi, hpsi )
  !
  CALL stop_clock('protxpsig:hpsi')
  !
  IF ( overlap ) THEN
     !
     CALL start_clock('protxpsig:spsi')
     !
     CALL s_psi_ptr( npwx, npw, nstart, psi, spsi )
     !
     CALL stop_clock('protxpsig:spsi')
     !
  END IF
  !
  CALL start_clock('protxpsig:hc')
  !
  CALL compute_distmat( hr, psi, hpsi )
  !
  CALL stop_clock('protxpsig:hc')
  !
  CALL start_clock('protxpsig:sc')
  !
  IF ( overlap ) THEN
     !
     CALL compute_distmat( sr, psi, spsi )
     !
  ELSE
     !
     CALL compute_distmat( sr, psi, psi )
     !
  END IF
  !
  CALL stop_clock('protxpsig:sc')
  !
  ! ... Diagonalize
  !
  CALL start_clock('protxpsig:diag')
  !
  IF ( do_distr_diag_inside_bgrp ) THEN ! NB on output of pdiaghg en and vr are the same across ortho_parent_comm
     ! only the first bgrp performs the diagonalization
     IF( my_bgrp_id == root_bgrp_id ) CALL pdiaghg( nstart, hr, sr, nx, en, vr, idesc )
     IF( nbgrp > 1 ) THEN ! results must be brodcast to the other band groups
       CALL mp_bcast( vr, root_bgrp_id, inter_bgrp_comm )
       CALL mp_bcast( en, root_bgrp_id, inter_bgrp_comm )
     ENDIF
  ELSE
     CALL pdiaghg( nstart, hr, sr, nx, en, vr, idesc )
  END IF
  !
  e(:) = en(1:nbnd)
  !
  CALL stop_clock('protxpsig:diag')
  !
  ! ... update the basis set
  !
  CALL start_clock('protxpsig:evc')
  !
  tpsi = psi
  !
  CALL refresh_evc()
  !
  CALL stop_clock('protxpsig:evc')
  !
  DEALLOCATE( idesc_ip )
  DEALLOCATE( rank_ip )
  DEALLOCATE( en )
  DEALLOCATE( vr )
  DEALLOCATE( sr )
  DEALLOCATE( hr )
  IF ( overlap ) &
  DEALLOCATE( spsi )
  DEALLOCATE( hpsi )
  DEALLOCATE( tpsi )
  !
  CALL stop_clock('protxpsig')
  !
  !CALL print_clock('protxpsig')
  !CALL print_clock('protxpsig:hpsi')
  !CALL print_clock('protxpsig:spsi')
  !CALL print_clock('protxpsig:hc')
  !CALL print_clock('protxpsig:sc')
  !CALL print_clock('protxpsig:diag')
  !CALL print_clock('protxpsig:evc')
  !
  RETURN
  !
CONTAINS
  !
  SUBROUTINE compute_distmat( dm, v, w )
     !
     !  This subroutine compute <vi|wj> and store the
     !  result in distributed matrix dm 
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     REAL(DP), INTENT(OUT) :: dm( :, : )
     COMPLEX(DP) :: v(:,:), w(:,:)
     REAL(DP), ALLOCATABLE :: work( :, : )
     !
     ALLOCATE( work( nx, nx ) )
     !
     work = 0.0d0
     !
     DO ipc = 1, idesc(LAX_DESC_NPC) !  loop on column procs 
        !
        nc = idesc_ip( LAX_DESC_NC, 1, ipc )
        ic = idesc_ip( LAX_DESC_IC, 1, ipc )
        !
        DO ipr = 1, ipc ! use symmetry for the loop on row procs
           !
           nr = idesc_ip( LAX_DESC_NR, ipr, ipc )
           ir = idesc_ip( LAX_DESC_IR, ipr, ipc )
           !
           !  rank of the processor for which this block (ipr,ipc) is destinated
           !
           root = rank_ip( ipr, ipc )

           ! use blas subs. on the matrix block

           CALL DGEMM( 'T', 'N', nr, nc, npw2, 2.D0, v(1,ir), npwx2, w(1,ic), npwx2, 0.D0, work, nx )

           IF ( gstart == 2 ) &
              CALL DGER( nr, nc, -1.D0, v(1,ir), npwx2, w(1,ic), npwx2, work, nx )

           ! accumulate result on dm of root proc.

           CALL mp_root_sum( work, dm, root, ortho_parent_comm )

        END DO
        !
     END DO

     if (ortho_parent_comm.ne.intra_bgrp_comm .and. nbgrp > 1) dm = dm/nbgrp
     !
     CALL laxlib_dsqmsym( nstart, dm, nx, idesc )
     !
     DEALLOCATE( work )
     !
     RETURN
  END SUBROUTINE compute_distmat
  !
  !
  SUBROUTINE refresh_evc( )
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     REAL(DP), ALLOCATABLE :: vtmp( :, : )
     REAL(DP) :: beta

     ALLOCATE( vtmp( nx, nx ) )
     !
     DO ipc = 1, idesc(LAX_DESC_NPC) !  loop on column procs 
        !
        nc = idesc_ip( LAX_DESC_NC, 1, ipc )
        ic = idesc_ip( LAX_DESC_IC, 1, ipc )
        !
        IF( ic <= nbnd ) THEN
           !
           nc = min( nc, nbnd - ic + 1 )
           !
           beta = 0.0d0

           DO ipr = 1, idesc(LAX_DESC_NPR)
              !
              nr = idesc_ip( LAX_DESC_NR, ipr, ipc )
              ir = idesc_ip( LAX_DESC_IR, ipr, ipc )
              !
              root = rank_ip( ipr, ipc )

              IF( ipr-1 == idesc(LAX_DESC_MYR) .AND. ipc-1 == idesc(LAX_DESC_MYC) .AND. la_proc ) THEN
                 !
                 !  this proc sends his block
                 ! 
                 CALL mp_bcast( vr(:,1:nc), root, ortho_parent_comm )
                 !
                 CALL DGEMM( 'N', 'N', npw2, nc, nr, 1.D0,  &
                             tpsi(1,ir), npwx2, vr, nx, beta, evc(1,ic),  npwx2 )
                 !
                 CALL DGEMM( 'N', 'N', npw2, nc, nr, 1.D0,  &
                             hpsi(1,ir), npwx2, vr, nx, beta, hevc(1,ic), npwx2 )
                 !
                 IF ( overlap ) &
                 CALL DGEMM( 'N', 'N', npw2, nc, nr, 1.D0,  &
                             spsi(1,ir), npwx2, vr, nx, beta, sevc(1,ic), npwx2 )
                 !
              ELSE
                 !
                 !  all other procs receive
                 ! 
                 CALL mp_bcast( vtmp(:,1:nc), root, ortho_parent_comm )
                 !
                 CALL DGEMM( 'N', 'N', npw2, nc, nr, 1.D0,  &
                             tpsi(1,ir), npwx2, vtmp, nx, beta, evc(1,ic),  npwx2 )
                 !
                 CALL DGEMM( 'N', 'N', npw2, nc, nr, 1.D0,  &
                             hpsi(1,ir), npwx2, vtmp, nx, beta, hevc(1,ic), npwx2 )
                 !
                 IF ( overlap ) &
                 CALL DGEMM( 'N', 'N', npw2, nc, nr, 1.D0,  &
                             spsi(1,ir), npwx2, vtmp, nx, beta, sevc(1,ic), npwx2 )
                 !
              END IF
              ! 

              beta = 1.0d0

           END DO
           !
        END IF
        !
     END DO
     !
     DEALLOCATE( vtmp )

     RETURN
  END SUBROUTINE refresh_evc
  !
END SUBROUTINE protate_xpsi_gamma
