!
! Copyright (C) 2019 National Institute of Advanced Industrial Science and Technology (AIST)
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE rotate_xpsi_k( h_psi_ptr, s_psi_ptr, overlap, &
                          npwx, npw, nstart, nbnd, npol, psi, evc, hevc, sevc, e )
  !----------------------------------------------------------------------------
  !
  ! ... Serial version of rotate_xpsi for colinear, k-point calculations
  !
  USE util_param,     ONLY : DP
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
  COMPLEX(DP), INTENT(OUT)   :: evc(npwx*npol,nbnd)
    ! input and output eigenvectors (may overlap)
  COMPLEX(DP), INTENT(OUT) :: hevc(npwx*npol,nbnd), sevc(npwx*npol,nbnd)
    ! H|psi> and S|psi>
  REAL(DP), INTENT(OUT) :: e(nbnd)
    ! eigenvalues
  !
  ! ... local variables
  !
  INTEGER                  :: kdim, kdmx
  INTEGER                  :: n_start, n_end, my_n
  COMPLEX(DP), ALLOCATABLE :: hc(:,:), sc(:,:), vc(:,:)
  COMPLEX(DP), ALLOCATABLE :: tpsi(:,:), hpsi(:,:), spsi(:,:)
  REAL(DP),    ALLOCATABLE :: en(:)
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
  CALL start_clock('rotxpsik')
  !
  ALLOCATE( tpsi( kdmx, nstart ) )
  ALLOCATE( hpsi( kdmx, nstart ) )
  ALLOCATE( hc( nstart, nstart) )    
  ALLOCATE( sc( nstart, nstart) )    
  ALLOCATE( vc( nstart, nstart) )    
  ALLOCATE( en( nstart ) )
  !$acc enter data create(hpsi, hc, sc, vc, tpsi, en )
  IF ( overlap ) THEN
    ALLOCATE( spsi(kdmx, nstart ) )
    !$acc enter data create(spsi)
  ENDIF
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
  !$acc kernels
  hc = (0.D0, 0.D0)
  !$acc end kernels
  !
  !$acc host_data use_device(psi, hpsi, hc)
  IF ( n_start <= n_end ) &
  CALL MYZGEMM( 'C', 'N', nstart, my_n, kdim, (1.D0, 0.D0), &
              psi, kdmx, hpsi(1,n_start), kdmx, (0.D0, 0.D0), hc(1,n_start), nstart )
  !
  CALL mp_sum( hc, inter_bgrp_comm )
  !
  CALL mp_sum( hc, intra_bgrp_comm )
  !$acc end host_data
  !
  CALL stop_clock('rotxpsik:hc')
  !
  CALL start_clock('rotxpsik:sc')
  !
  !$acc kernels
  sc = (0.D0, 0.D0)
  !$acc end kernels
  !
  IF ( overlap ) THEN
     !
     !$acc host_data use_device(psi, spsi, sc)
     IF ( n_start <= n_end ) &
     CALL MYZGEMM( 'C', 'N', nstart, my_n, kdim, (1.D0, 0.D0), &
                 psi, kdmx, spsi(1,n_start), kdmx, (0.D0, 0.D0), sc(1,n_start), nstart )
     !$acc end host_data
     !
  ELSE
     !
     !$acc host_data use_device(psi, sc)
     IF ( n_start <= n_end ) &
     CALL MYZGEMM( 'C', 'N', nstart, my_n, kdim, (1.D0, 0.D0), &
                 psi, kdmx, psi(1,n_start), kdmx, (0.D0, 0.D0), sc(1,n_start), nstart )
     !$acc end host_data
     !
  END IF
  !
  !$acc host_data use_device(sc)
  CALL mp_sum( sc, inter_bgrp_comm )
  !
  CALL mp_sum( sc, intra_bgrp_comm )
  !$acc end host_data
  !
  CALL stop_clock('rotxpsik:sc')
  !
  ! ... Diagonalize
  !
  CALL start_clock('rotxpsik:diag')
  !
  !$acc host_data use_device(hc, sc, en, vc)
  CALL diaghg( nstart, nbnd, hc, sc, nstart, en, vc, me_bgrp, root_bgrp, intra_bgrp_comm )
  !$acc end host_data
  !
  !$acc kernels
  e(:) = en(1:nbnd)
  !$acc end kernels
  !
  CALL stop_clock('rotxpsik:diag')
  !
  ! ... update the basis set
  !
  CALL start_clock('rotxpsik:evc')
  !
  !$acc kernels
  tpsi = psi
  !
  evc  = (0.D0, 0.D0)
  hevc = (0.D0, 0.D0)
  IF ( overlap ) &
  sevc = (0.D0, 0.D0)
  !$acc end kernels
  !
  IF ( n_start <= n_end ) THEN
     !
     !$acc host_data use_device(evc, hpsi, spsi, hevc, sevc, tpsi, vc)
     CALL MYZGEMM( 'N', 'N', kdim, nbnd, my_n, (1.D0, 0.D0), &
                 tpsi(1,n_start), kdmx, vc(n_start,1), nstart, (0.D0, 0.D0), evc,  kdmx )
     !
     CALL MYZGEMM( 'N', 'N', kdim, nbnd, my_n, (1.D0, 0.D0), &
                 hpsi(1,n_start), kdmx, vc(n_start,1), nstart, (0.D0, 0.D0), hevc, kdmx )
     !
     IF ( overlap ) &
     CALL MYZGEMM( 'N', 'N', kdim, nbnd, my_n, (1.D0, 0.D0), &
                 spsi(1,n_start), kdmx, vc(n_start,1), nstart, (0.D0, 0.D0), sevc, kdmx )
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
  IF ( overlap ) THEN
    !$acc exit data delete(spsi)
    DEALLOCATE( spsi )
  ENDIF
  !$acc exit data delete(hpsi, en, vc, sc, hc, tpsi )
  DEALLOCATE( en, vc, sc, hc, hpsi, tpsi )
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
END SUBROUTINE rotate_xpsi_k
!
!
!----------------------------------------------------------------------------
SUBROUTINE protate_xpsi_k( h_psi_ptr, s_psi_ptr, overlap, &
                           npwx, npw, nstart, nbnd, npol, psi, evc, hevc, sevc, e )
  !----------------------------------------------------------------------------
  !
  ! ... Parallel version of rotate_xpsi for colinear, k-point calculations
  ! ... Subroutine with distributed matrices, written by Carlo Cavazzoni
  !
  USE util_param,        ONLY : DP
  USE mp_bands_util,    ONLY : intra_bgrp_comm, inter_bgrp_comm, &
                               nbgrp, root_bgrp_id, my_bgrp_id
  USE mp,               ONLY : mp_bcast, mp_root_sum, mp_sum, mp_barrier
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
  COMPLEX(DP), INTENT(OUT)   :: evc(npwx*npol,nbnd)
    ! input and output eigenvectors (may overlap)
  COMPLEX(DP), INTENT(OUT) :: hevc(npwx*npol,nbnd), sevc(npwx*npol,nbnd)
    ! H|psi> and S|psi>
  REAL(DP), INTENT(OUT) :: e(nbnd)
    ! eigenvalues
  !
  ! ... local variables
  !
  INTEGER :: kdim, kdmx
  COMPLEX(DP), ALLOCATABLE :: hc(:,:), sc(:,:), vc(:,:)
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

  CALL start_clock('protxpsik')
  !
  CALL laxlib_getval( do_distr_diag_inside_bgrp = do_distr_diag_inside_bgrp, &
       ortho_parent_comm = ortho_parent_comm )
  CALL desc_init( nstart, nx, la_proc, idesc, rank_ip, idesc_ip )
  !
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
  ALLOCATE( tpsi( kdmx, nstart ) )
  ALLOCATE( hpsi( kdmx, nstart ) )
  IF ( overlap ) &
  ALLOCATE( spsi( kdmx, nstart ) )
  ALLOCATE( hc( nx, nx) )    
  ALLOCATE( sc( nx, nx) )    
  ALLOCATE( vc( nx, nx) )    
  ALLOCATE( en( nstart ) )
  !
  ! ... Set up the Hamiltonian and Overlap matrix on the subspace :
  !
  ! ...      H_ij = <psi_i| H |psi_j>     S_ij = <psi_i| S |psi_j>
  !
  CALL start_clock('protxpsik:hpsi')
  !
  CALL h_psi_ptr( npwx, npw, nstart, psi, hpsi )
  !
  CALL stop_clock('protxpsik:hpsi')
  !
  IF ( overlap ) THEN
     !
     CALL start_clock('protxpsik:spsi')
     !
     CALL s_psi_ptr( npwx, npw, nstart, psi, spsi )
     !
     CALL stop_clock('protxpsik:spsi')
     !
  END IF
  !
  CALL start_clock('protxpsik:hc')
  !
  CALL compute_distmat( hc, psi, hpsi )
  !
  CALL stop_clock('protxpsik:hc')
  !
  CALL start_clock('protxpsik:sc')
  !            
  IF ( overlap ) THEN
     !
     CALL compute_distmat( sc, psi, spsi )
     !
  ELSE
     !
     CALL compute_distmat( sc, psi, psi )
     ! 
  END IF
  !
  CALL stop_clock('protxpsik:sc')
  !
  ! ... Diagonalize
  !
  CALL start_clock('protxpsik:diag')
  !
  IF ( do_distr_diag_inside_bgrp ) THEN ! NB on output of pdiaghg en and vc are the same across ortho_parent_comm
     ! only the first bgrp performs the diagonalization
     IF( my_bgrp_id == root_bgrp_id ) CALL pdiaghg( nstart, hc, sc, nx, en, vc, idesc )
     IF( nbgrp > 1 ) THEN ! results must be brodcast to the other band groups
       CALL mp_bcast( vc, root_bgrp_id, inter_bgrp_comm )
       CALL mp_bcast( en, root_bgrp_id, inter_bgrp_comm )
     ENDIF
  ELSE
     CALL pdiaghg( nstart, hc, sc, nx, en, vc, idesc )
  END IF
  !
  e(:) = en(1:nbnd)
  !
  CALL stop_clock('protxpsik:diag')
  !
  ! ... update the basis set
  !
  CALL start_clock('protxpsik:evc')
  !
  tpsi = psi
  !
  CALL refresh_evc()
  !
  CALL stop_clock('protxpsik:evc')
  !
  DEALLOCATE( en )
  DEALLOCATE( vc )
  DEALLOCATE( sc )
  DEALLOCATE( hc )
  IF ( overlap ) &
  DEALLOCATE( spsi )
  DEALLOCATE( hpsi )
  DEALLOCATE( tpsi )
  !
  DEALLOCATE( idesc_ip )
  DEALLOCATE( rank_ip )
  !
  CALL stop_clock('protxpsik')
  !
  !CALL print_clock('protxpsik')
  !CALL print_clock('protxpsik:hpsi')
  !CALL print_clock('protxpsik:spsi')
  !CALL print_clock('protxpsik:hc')
  !CALL print_clock('protxpsik:sc')
  !CALL print_clock('protxpsik:diag')
  !CALL print_clock('protxpsik:evc')
  !
  RETURN
  !
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
     COMPLEX(DP), INTENT(OUT) :: dm( :, : )
     COMPLEX(DP) :: v(:,:), w(:,:)
     COMPLEX(DP), ALLOCATABLE :: work( :, : )
     !
     ALLOCATE( work( nx, nx ) )
     !
     work = ( 0.0_DP, 0.0_DP )
     !
     DO ipc = 1, idesc(LAX_DESC_NPC) !  loop on column procs 
        !
        nc = idesc_ip( LAX_DESC_NC, 1, ipc )
        ic = idesc_ip( LAX_DESC_IC, 1, ipc )
        !
        DO ipr = 1, ipc ! desc%npr ! ipc ! use symmetry for the loop on row procs
           !
           nr = idesc_ip( LAX_DESC_NR, ipr, ipc )
           ir = idesc_ip( LAX_DESC_IR, ipr, ipc )
           !
           !  rank of the processor for which this block (ipr,ipc) is destinated
           !
           root = rank_ip( ipr, ipc )

           ! use blas subs. on the matrix block

           CALL ZGEMM( 'C', 'N', nr, nc, kdim, ( 1.D0, 0.D0 ),  v(1,ir), kdmx, w(1,ic), kdmx, ( 0.D0, 0.D0 ), work, nx )

           ! accumulate result on dm of root proc.
           CALL mp_root_sum( work, dm, root, ortho_parent_comm )

        END DO
        !
     END DO
     if (ortho_parent_comm.ne.intra_bgrp_comm .and. nbgrp > 1) dm = dm/nbgrp
     !
     CALL laxlib_zsqmher( nstart, dm, nx, idesc )
     !
     DEALLOCATE( work )
     !
     RETURN
  END SUBROUTINE compute_distmat


  SUBROUTINE refresh_evc( )
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     COMPLEX(DP), ALLOCATABLE :: vtmp( :, : )
     COMPLEX(DP) :: beta

     ALLOCATE( vtmp( nx, nx ) )
     !
     DO ipc = 1, idesc(LAX_DESC_NPC)
        !
        nc = idesc_ip( LAX_DESC_NC, 1, ipc )
        ic = idesc_ip( LAX_DESC_IC, 1, ipc )
        !
        IF( ic <= nbnd ) THEN
           !
           nc = min( nc, nbnd - ic + 1 )
           !
           beta = ( 0.D0, 0.D0 )

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
                 CALL mp_bcast( vc(:,1:nc), root, ortho_parent_comm )
                 !
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ( 1.D0, 0.D0 ),  &
                             tpsi(1,ir), kdmx, vc, nx, beta, evc(1,ic),  kdmx )
                 !
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ( 1.D0, 0.D0 ),  &
                             hpsi(1,ir), kdmx, vc, nx, beta, hevc(1,ic), kdmx )
                 !
                 IF ( overlap ) &
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ( 1.D0, 0.D0 ),  &
                             spsi(1,ir), kdmx, vc, nx, beta, sevc(1,ic), kdmx )
                 !
              ELSE
                 !
                 !  all other procs receive
                 ! 
                 CALL mp_bcast( vtmp(:,1:nc), root, ortho_parent_comm )
                 !
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ( 1.D0, 0.D0 ),  &
                             tpsi(1,ir), kdmx, vtmp, nx, beta, evc(1,ic),  kdmx )
                 !
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ( 1.D0, 0.D0 ),  &
                             hpsi(1,ir), kdmx, vtmp, nx, beta, hevc(1,ic), kdmx )
                 !
                 IF ( overlap ) &
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ( 1.D0, 0.D0 ),  &
                             spsi(1,ir), kdmx, vtmp, nx, beta, sevc(1,ic), kdmx )
                 !
              END IF
              ! 

              beta = ( 1.D0, 0.D0 )

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
END SUBROUTINE protate_xpsi_k
