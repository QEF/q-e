!
! Copyright (C) 2003-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE rotate_wfc_gamma( h_psi, s_psi, overlap, &
                             npwx, npw, nstart, nbnd, psi, evc, e )
  !----------------------------------------------------------------------------
  !
  ! ... Serial version of rotate_wfc for Gamma-only calculations
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
  COMPLEX(DP) :: psi(npwx,nstart), evc(npwx,nbnd)
    ! input and output eigenvectors (may overlap)
  REAL(DP) :: e(nbnd)
    ! eigenvalues
  !
  ! ... local variables
  !
  INTEGER                  :: npw2, npwx2
  COMPLEX(DP), ALLOCATABLE :: aux(:,:)
  REAL(DP),    ALLOCATABLE :: hr(:,:), sr(:,:), vr(:,:)
  REAL(DP),    ALLOCATABLE :: en(:)
  INTEGER :: n_start, n_end, my_n
  !
  EXTERNAL  h_psi,    s_psi
    ! h_psi(npwx,npw,nvec,psi,hpsi)
    !     calculates H|psi>
    ! s_psi(npwx,npw,nvec,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx,npol,nvec)

  npw2  = 2 * npw
  npwx2 = 2 * npwx

  IF ( gstart == -1 ) CALL errore( 'regter', 'gstart variable not initialized', 1 )

  !
  ALLOCATE( aux(  npwx, nstart ) )    
  ALLOCATE( hr( nstart, nstart ) )    
  ALLOCATE( sr( nstart, nstart ) )    
  ALLOCATE( vr( nstart, nstart ) )    
  ALLOCATE( en( nstart ) )
  call start_clock('rotwfcg'); !write(*,*) 'start rotwfcg' ; FLUSH(6)
  !
  ! ... Set up the Hamiltonian and Overlap matrix on the subspace :
  !
  ! ...      H_ij = <psi_i| H |psi_j>     S_ij = <psi_i| S |psi_j>
  !
  ! ... set Im[ psi(G=0) ] -  needed for numerical stability
  !
  IF ( gstart == 2 ) &
     psi(1,1:nstart) = CMPLX( DBLE( psi(1,1:nstart) ), 0.D0,kind=DP)
  !
  call start_clock('rotwfcg:hpsi'); !write(*,*) 'start rotwfcg:hpsi' ; FLUSH(6)
  CALL h_psi( npwx, npw, nstart, psi, aux )
  call stop_clock('rotwfcg:hpsi'); !write(*,*) 'stop rotwfcg:hpsi' ; FLUSH(6)
  !
  call start_clock('rotwfcg:hc'); !write(*,*) 'start rotwfcg:hc' ; FLUSH(6)
  hr=0.D0
  CALL divide(inter_bgrp_comm,nstart,n_start,n_end)
  my_n = n_end - n_start + 1; !write (*,*) nstart,n_start,n_end
  if (n_start .le. n_end) &
  CALL DGEMM( 'T','N', nstart, my_n, npw2, 2.D0, psi, npwx2, aux(1,n_start), npwx2, 0.D0, hr(1,n_start), nstart )
  IF ( gstart == 2 ) call DGER( nstart, my_n, -1.D0, psi, npwx2, aux(1,n_start), npwx2, hr(1,n_start), nstart )
  CALL mp_sum( hr, inter_bgrp_comm )
  !     
  CALL mp_sum( hr, intra_bgrp_comm )
  !     
  sr=0.D0
  IF ( overlap ) THEN 
     ! 
     CALL s_psi( npwx, npw, nstart, psi, aux )
     !
     if (n_start .le. n_end) &
     CALL DGEMM( 'T','N', nstart, my_n, npw2, 2.D0, psi, npwx2, aux(1,n_start), npwx2, 0.D0, sr(1,n_start), nstart )
     IF ( gstart == 2 ) CALL DGER( nstart, my_n, -1.D0, psi, npwx2, aux(1,n_start), npwx2, sr(1,n_start), nstart )
     !              
  ELSE
     !
     if (n_start .le. n_end) &
     CALL DGEMM( 'T','N', nstart, my_n, npw2, 2.D0, psi, npwx2, psi(1,n_start), npwx2, 0.D0, sr(1,n_start), nstart )
     IF ( gstart == 2 ) CALL DGER( nstart, my_n, -1.D0, psi, npwx2, psi(1,n_start), npwx2, sr(1,n_start), nstart )
     !
  END IF
  CALL mp_sum( sr, inter_bgrp_comm )
  !
  CALL mp_sum( sr, intra_bgrp_comm )
  call stop_clock('rotwfcg:hc'); !write(*,*) 'stop rotwfcg:hc' ; FLUSH(6)
  !
  ! ... Diagonalize
  !
  call start_clock('rotwfcg:diag'); !write(*,*) 'start rotwfcg:diag' ; FLUSH(6)
  CALL diaghg( nstart, nbnd, hr, sr, nstart, en, vr, me_bgrp, root_bgrp, intra_bgrp_comm )
  call stop_clock('rotwfcg:diag'); !write(*,*) 'stop rotwfcg:diag' ; FLUSH(6)
  call start_clock('rotwfcg:evc'); !write(*,*) 'start rotwfcg:evc' ; FLUSH(6)
  !
  e(:) = en(1:nbnd)
  !
  ! ... update the basis set
  !
  aux=(0.D0,0.D0)
  if (n_start .le. n_end) &
  CALL DGEMM( 'N','N', npw2, nbnd, my_n, 1.D0, psi(1,n_start), npwx2, vr(n_start,1), nstart, 0.D0, aux, npwx2 )
  CALL mp_sum( aux, inter_bgrp_comm )
  !   
  evc(:,:) = aux(:,1:nbnd)
  call stop_clock('rotwfcg:evc'); !write(*,*) 'stop rotwfcg:evc' ; FLUSH(6)
  !
  DEALLOCATE( en )
  DEALLOCATE( vr )
  DEALLOCATE( sr )
  DEALLOCATE( hr )
  DEALLOCATE( aux )
  call stop_clock('rotwfcg'); !write(*,*) 'stop rotwfcg' ; FLUSH(6)
  !call print_clock('rotwfcg')
  !call print_clock('rotwfcg:hpsi')
  !call print_clock('rotwfcg:hc')
  !call print_clock('rotwfcg:diag')
  !call print_clock('rotwfcg:evc')
  !
  RETURN
  !
END SUBROUTINE rotate_wfc_gamma
!
!
!----------------------------------------------------------------------------
SUBROUTINE protate_wfc_gamma( h_psi, s_psi, overlap, &
                              npwx, npw, nstart, nbnd, psi, evc, e )
  !----------------------------------------------------------------------------
  !
  ! ... Parallel version of rotate_wfc for Gamma-only calculations
  ! ... Subroutine with distributed matrices, written by Carlo Cavazzoni
  ! ... This version assumes real wavefunctions (k=0) with only
  ! ... half plane waves stored: psi(-G)=psi*(G), except G=0
  !
  USE util_param,       ONLY : DP
  USE mp_bands_util,    ONLY : intra_bgrp_comm, inter_bgrp_comm, root_bgrp_id, nbgrp, my_bgrp_id
  USE mp_bands_util,    ONLY : gstart ! index of the first nonzero G 
  USE mp,               ONLY : mp_bcast, mp_root_sum, mp_sum, mp_barrier
  !
  IMPLICIT NONE
  !
  include 'laxlib.fh'
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
  COMPLEX(DP) :: psi(npwx,nstart), evc(npwx,nbnd)
    ! input and output eigenvectors (may overlap)
  REAL(DP) :: e(nbnd)
    ! eigenvalues
  !
  ! ... local variables:
  !
  INTEGER             :: npw2, npwx2
  COMPLEX(DP), ALLOCATABLE :: aux(:,:)
  REAL(DP),    ALLOCATABLE :: hr(:,:), sr(:,:), vr(:,:)
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
  EXTERNAL  h_psi,    s_psi
    ! h_psi(npwx,npw,nvec,psi,hpsi)
    !     calculates H|psi>
    ! s_psi(npwx,npw,nvec,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx,npol,nvec)

  call start_clock('protwfcg'); !write(*,*) 'start protwfcg' ; FLUSH(6)
  !
  CALL laxlib_getval( do_distr_diag_inside_bgrp = do_distr_diag_inside_bgrp, &
       ortho_parent_comm = ortho_parent_comm )
  CALL desc_init( nstart, nx, la_proc, idesc, rank_ip, idesc_ip )
  !
  npw2  = 2 * npw
  npwx2 = 2 * npwx

  IF ( gstart == -1 ) CALL errore( 'regter', 'gstart variable not initialized', 1 )
  !
  ALLOCATE( aux(  npwx, nstart ) )    
  ALLOCATE( hr( nx, nx ) )    
  ALLOCATE( sr( nx, nx ) )    
  ALLOCATE( vr( nx, nx ) )    
  ALLOCATE( en( nstart ) )

  aux=(0.0_DP,0.0_DP) 
  !
  ! ... Set up the Hamiltonian and Overlap matrix on the subspace :
  !
  ! ...      H_ij = <psi_i| H |psi_j>     S_ij = <psi_i| S |psi_j>
  !
  ! ... set Im[ psi(G=0) ] -  needed for numerical stability
  IF ( gstart == 2 ) &
     psi(1,1:nstart) = CMPLX( DBLE( psi(1,1:nstart) ), 0.D0, kind=DP)
  !
  call start_clock('protwfcg:hpsi'); !write(*,*) 'start protwfcg:hpsi' ; FLUSH(6)
  CALL h_psi( npwx, npw, nstart, psi, aux )
  call stop_clock('protwfcg:hpsi'); !write(*,*) 'stop protwfcg:hpsi' ; FLUSH(6)
  !
  call start_clock('protwfcg:hc'); !write(*,*) 'start protwfcg:hc' ; FLUSH(6)
  CALL compute_distmat( hr, psi, aux )
  !
  IF ( overlap ) THEN
     !
     CALL s_psi( npwx, npw, nstart, psi, aux )
     CALL compute_distmat( sr, psi, aux )
     !              
  ELSE
     !
     CALL compute_distmat( sr, psi, psi )
     !
  END IF
  call stop_clock('protwfcg:hc'); !write(*,*) 'stop protwfcg:hc' ; FLUSH(6)
  !
  ! ... Diagonalize
  !
  call start_clock('protwfcg:diag'); !write(*,*) 'start protwfcg:diag' ; FLUSH(6)
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
  call stop_clock('protwfcg:diag'); !write(*,*) 'stop protwfcg:diag' ; FLUSH(6)
  !
  e(:) = en(1:nbnd)
  !
  ! ... update the basis set
  !
  call start_clock('protwfcg:evc'); !write(*,*) 'start protwfcg:evc' ; FLUSH(6)
  CALL refresh_evc( )
  !   
  evc(:,:) = aux(:,1:nbnd)
  call stop_clock('protwfcg:evc'); !write(*,*) 'stop protwfcg:evc' ; FLUSH(6)
  !
  DEALLOCATE( en )
  DEALLOCATE( vr )
  DEALLOCATE( sr )
  DEALLOCATE( hr )
  DEALLOCATE( aux )
  !
  DEALLOCATE( idesc_ip )
  DEALLOCATE( rank_ip )
  call stop_clock('protwfcg'); !write(*,*) 'stop protwfcg' ; FLUSH(6)
  !call print_clock('protwfcg')
  !call print_clock('protwfcg:hpsi')
  !call print_clock('protwfcg:hc')
  !call print_clock('protwfcg:diag')
  !call print_clock('protwfcg:evc')

  !
  RETURN
  !
CONTAINS
  !
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
        nc = idesc_ip(LAX_DESC_NC, 1, ipc )
        ic = idesc_ip(LAX_DESC_IC, 1, ipc ) 
        !
        DO ipr = 1, ipc ! use symmetry for the loop on row procs
           !
           nr = idesc_ip(LAX_DESC_NR, ipr, ipc ) 
           ir = idesc_ip(LAX_DESC_IR, ipr, ipc )
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
        nc = idesc_ip(LAX_DESC_NC, 1, ipc )
        ic = idesc_ip(LAX_DESC_IC, 1, ipc ) 
        !
        IF( ic <= nbnd ) THEN
           !
           nc = min( nc, nbnd - ic + 1 )
           !
           beta = 0.0d0

           DO ipr = 1, idesc(LAX_DESC_NPR)
              !
              nr = idesc_ip(LAX_DESC_NR, ipr, ipc ) 
              ir = idesc_ip(LAX_DESC_IR, ipr, ipc )
              !
              root = rank_ip( ipr, ipc )

              IF( ipr-1 == idesc(LAX_DESC_MYR) .AND. ipc-1 == idesc(LAX_DESC_MYC) .AND. la_proc ) THEN
                 !
                 !  this proc sends his block
                 ! 
                 CALL mp_bcast( vr(:,1:nc), root, ortho_parent_comm )
                 CALL DGEMM( 'N', 'N', npw2, nc, nr, 1.D0, psi(1,ir), npwx2, vr, nx, beta, aux(1,ic), npwx2 )
              ELSE
                 !
                 !  all other procs receive
                 ! 
                 CALL mp_bcast( vtmp(:,1:nc), root, ortho_parent_comm )
                 CALL DGEMM( 'N', 'N', npw2, nc, nr, 1.D0, psi(1,ir), npwx2, vtmp, nx, beta, aux(1,ic), npwx2 )
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
END SUBROUTINE protate_wfc_gamma
