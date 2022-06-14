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
SUBROUTINE rotate_HSpsi_gamma( npwx, npw, nstart, nbnd, psi, hpsi, overlap, spsi, e )
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
  USE mp,            ONLY : mp_sum, mp_barrier, mp_allgather, mp_type_create_column_section, mp_type_free
  USE device_memcpy_m,    ONLY: dev_memcpy
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
    nbnd                           ! output number of states
  COMPLEX(DP), INTENT(INOUT) :: psi(npwx,nstart), hpsi(npwx,nstart) ! input and output psi, Hpsi,
  COMPLEX(DP), INTENT(INOUT), OPTIONAL :: spsi(npwx,nstart)         ! ...   and optionnally Spsi
  LOGICAL, INTENT(IN) :: overlap   ! if .FALSE. : spsi is not needed (and not used)
  REAL(DP), INTENT(OUT) :: e(nbnd) ! eigenvalues of the reduced H matrix
  !
  ! ... local variables
  !
  INTEGER                  :: kdim, kdmx
  COMPLEX(DP), ALLOCATABLE :: aux(:,:)
  !$acc declare device_resident(aux)
  REAL(DP),    ALLOCATABLE :: hh(:,:), ss(:,:), vv(:,:)
  !$acc declare device_resident(hh, ss, vv)
  REAL(DP),    ALLOCATABLE :: en(:)
  !$acc declare device_resident(en)
  INTEGER :: n_start, n_end, my_n, recv_counts(nbgrp), displs(nbgrp), column_type
  !
  IF ( gstart == -1 ) CALL errore( 'rotHSw', 'gstart variable not initialized', 1 )
  IF ( overlap .AND..NOT.present(spsi) ) call errore( 'rotHSw','spsi array needed with overlap=.TRUE.',1)
  !
  call start_clock('rotHSw'); !write(*,*) 'start rotHSw' ; FLUSH(6)
  !
  ! ... set Im[ psi(G=0) ] etc -  needed for numerical stability
  !
  IF ( gstart == 2 ) then
     !$acc kernels
     psi (1,1:nstart) = CMPLX( DBLE( psi (1,1:nstart) ), 0.D0,kind=DP)
     hpsi(1,1:nstart) = CMPLX( DBLE( hpsi(1,1:nstart) ), 0.D0,kind=DP)
     if (overlap) spsi(1,1:nstart) = CMPLX( DBLE( spsi(1,1:nstart) ), 0.D0,kind=DP)
     !$acc end kernels
  END IF
  
  kdim = 2 * npw
  kdmx = 2 * npwx
  !
  ALLOCATE( hh( nstart, nstart ) )
  ALLOCATE( ss( nstart, nstart ) )
  ALLOCATE( vv( nstart, nstart ) )
  ALLOCATE( en( nstart ) )
  !
  !$acc host_data use_device(psi, hpsi, spsi, hh, ss, vv, en)
  !
  ! ... Set up the Hamiltonian and Overlap matrix on the subspace :
  !
  ! ...      H_ij = <psi_i| H |psi_j>     S_ij = <psi_i| S |psi_j>
  !
  call start_clock('rotHSw:hc'); !write(*,*) 'start rotHSw:hc' ; FLUSH(6)

  CALL mp_type_create_column_section(hh(1,1), 0, nstart, nstart, column_type)
  CALL divide_all(inter_bgrp_comm,nstart,n_start,n_end,recv_counts,displs)

  my_n = n_end - n_start + 1; !write (*,*) nstart,n_start,n_end
  if (n_start .le. n_end) &
  CALL MYDGEMM( 'T','N', nstart,my_n,kdim, 2.D0, psi,kdmx, hpsi(1,n_start),kdmx, 0.D0, hh(1,n_start),nstart )
  IF ( gstart == 2 ) call MYDGER( nstart, my_n, -1.D0, psi,kdmx, hpsi(1,n_start),kdmx, hh(1,n_start),nstart )
  call start_clock('rotHSw:hc:s1')
  CALL mp_sum( hh(:,n_start:n_end), intra_bgrp_comm ) ! this section only needs to be collected inside bgrp
  call stop_clock('rotHSw:hc:s1')
!  call start_clock('rotHSw:hc:b1'); CALL mp_barrier( inter_bgrp_comm ); call stop_clock('rotHSw:hc:b1')
  call start_clock('rotHSw:hc:s2')
  CALL mp_allgather(hh, column_type, recv_counts, displs, inter_bgrp_comm)
  call stop_clock('rotHSw:hc:s2')
  !
  IF ( overlap ) THEN
     !
     if (n_start .le. n_end) &
     CALL MYDGEMM('T','N', nstart,my_n,kdim, 2.D0, psi,kdmx, spsi(1,n_start),kdmx, 0.D0, ss(1,n_start),nstart)
     IF ( gstart == 2 ) CALL MYDGER(nstart, my_n, -1.D0, psi,kdmx, spsi(1,n_start),kdmx, ss(1,n_start),nstart)
     !
  ELSE
     !
     if (n_start .le. n_end) &
     CALL MYDGEMM('T','N', nstart,my_n,kdim, 2.D0, psi,kdmx, psi(1,n_start),kdmx, 0.D0, ss(1,n_start),nstart)
     IF ( gstart == 2 ) CALL MYDGER(nstart, my_n, -1.D0, psi,kdmx, psi(1,n_start),kdmx, ss(1,n_start),nstart)
     !
  END IF
  call start_clock('rotHSw:hc:s3')
  CALL mp_sum( ss(:,n_start:n_end), intra_bgrp_comm ) ! this section only needs to be collected inside bgrp
  call stop_clock('rotHSw:hc:s3')
!  call start_clock('rotHSw:hc:b2'); CALL mp_barrier( inter_bgrp_comm ); call stop_clock('rotHSw:hc:b2')
  call start_clock('rotHSw:hc:s4')
  CALL mp_allgather(ss, column_type, recv_counts, displs, inter_bgrp_comm)
  call stop_clock('rotHSw:hc:s4')

  CALL mp_type_free( column_type )

  call stop_clock('rotHSw:hc'); !write(*,*) 'stop rotHSw:hc' ; FLUSH(6)
  !
  ! ... Diagonalize
  !
  call start_clock('rotHSw:diag'); !write(*,*) 'start rotHSw:diag' ; FLUSH(6)
  CALL diaghg( nstart, nbnd, hh, ss, nstart, en, vv, me_bgrp, root_bgrp, intra_bgrp_comm )
  !$acc data deviceptr(e)
  CALL dev_memcpy(e, en, [1,nbnd])
  !$acc end data
  call stop_clock('rotHSw:diag'); !write(*,*) 'stop rotHSw:diag' ; FLUSH(6)
  !
  ! ... update the basis set
  !
  call start_clock('rotHSw:evc'); !write(*,*) 'start rotHSw:evc' ; FLUSH(6)

  CALL mp_type_create_column_section(psi(1,1), 0, npwx, npwx, column_type)
  CALL divide_all(inter_bgrp_comm,nbnd,n_start,n_end,recv_counts,displs)
  !$acc end host_data

  ALLOCATE( aux ( npwx, nbnd ) )
  !$acc host_data use_device(psi, hpsi, spsi, hh, ss, vv, en, aux)

  my_n = n_end - n_start + 1; !write (*,*) nstart,n_start,n_end
  if (n_start .le. n_end) &
  CALL MYDGEMM( 'N','N', kdim,my_n,nstart, 1.D0, psi,kdmx,vv(1,n_start),nstart, 0.D0, aux(1,n_start),kdmx )
  CALL dev_memcpy(psi, aux, [1, npwx], 1, [n_start,n_end])
!  call start_clock('rotHSw:ev:b3'); CALL mp_barrier( inter_bgrp_comm ); call stop_clock('rotHSw:ev:b3')
  call start_clock('rotHSw:ev:s5')
  CALL mp_allgather(psi(:,1:nbnd), column_type, recv_counts, displs, inter_bgrp_comm)
  call stop_clock('rotHSw:ev:s5')

  if (n_start .le. n_end) &
  CALL MYDGEMM( 'N','N', kdim,my_n,nstart, 1.D0,hpsi,kdmx,vv(1,n_start),nstart, 0.D0, aux(1,n_start),kdmx )
  CALL dev_memcpy(hpsi, aux, [1, npwx], 1, [n_start,n_end])  
!  call start_clock('rotHSw:ev:b4'); CALL mp_barrier( inter_bgrp_comm ); call stop_clock('rotHSw:ev:b4')
  call start_clock('rotHSw:ev:s6')
  CALL mp_allgather(hpsi(:,1:nbnd), column_type, recv_counts, displs, inter_bgrp_comm)
  call stop_clock('rotHSw:ev:s6')

  IF (overlap) THEN

     if (n_start .le. n_end) &
     CALL MYDGEMM( 'N','N', kdim,my_n,nstart, 1.D0,spsi,kdmx,vv(1,n_start),nstart, 0.D0, aux(1,n_start),kdmx )
     CALL dev_memcpy(spsi, aux, [1, npwx], 1, [n_start,n_end])     
!     call start_clock('rotHSw:ev:b5'); CALL mp_barrier( inter_bgrp_comm ); call stop_clock('rotHSw:ev:b5')
     call start_clock('rotHSw:ev:s7')
     CALL mp_allgather(spsi(:,1:nbnd), column_type, recv_counts, displs, inter_bgrp_comm)
     call stop_clock('rotHSw:ev:s7')

  ELSE IF (present(spsi)) THEN

     CALL dev_memcpy(spsi, psi, [1, npwx], 1, [1,nbnd])

  END IF
  !
  !$acc end host_data
  !
  DEALLOCATE( aux )
  CALL mp_type_free( column_type )

  call stop_clock('rotHSw:evc'); !write(*,*) 'stop rotHSw:evc' ; FLUSH(6)
  !
  DEALLOCATE( vv )
  DEALLOCATE( ss )
  DEALLOCATE( hh )
  DEALLOCATE( en )
  !
  call stop_clock('rotHSw'); !write(*,*) 'stop rotHSw' ; FLUSH(6)
  !call print_clock('rotHSw')
  !call print_clock('rotHSw:hc')
  !call print_clock('rotHSw:diag')
  !call print_clock('rotHSw:evc')
  !
  RETURN
  !
END SUBROUTINE rotate_HSpsi_gamma
!
!
!----------------------------------------------------------------------------
SUBROUTINE protate_HSpsi_gamma( npwx, npw, nstart, nbnd, psi, hpsi, overlap, spsi, e )
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
  USE mp,               ONLY : mp_bcast, mp_root_sum, mp_sum, mp_barrier, mp_rank
  USE mp,               ONLY : mp_allgather, mp_type_create_column_section, mp_type_free
  !
  IMPLICIT NONE
  !
  include 'laxlib.fh'
  !
  ! ... I/O variables
  !
  INTEGER, INTENT(IN) :: &
    npw,                 &         ! dimension of the matrices (psi,Hpsi,Spsi) to be rotated
    npwx,                &         ! leading dimension of the wavefunction-related matrices
    nstart,              &         ! input number of states
    nbnd                           ! output number of states
  COMPLEX(DP), INTENT(INOUT) :: psi(npwx,nstart), hpsi(npwx,nstart) ! input and output psi, Hpsi,
  COMPLEX(DP), INTENT(INOUT), OPTIONAL :: spsi(npwx,nstart)         ! ...   and optionnally Spsi
  LOGICAL, INTENT(IN) :: overlap   ! if .FALSE. : spsi is not needed (and not used)
  REAL(DP), INTENT(OUT) :: e(nbnd) ! eigenvalues of the reduced H matrix
  !
  ! ... local variables
  !
  INTEGER                  :: kdim, kdmx
  COMPLEX(DP), ALLOCATABLE :: aux(:,:)
  REAL(DP),    ALLOCATABLE :: hh(:,:), ss(:,:), vv(:,:), v(:,:)
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
  INTEGER :: npw_s, npw_e, my_kdim, ib_s, ib_e, my_nbnd, recv_counts(nbgrp), displs(nbgrp), column_type
  !
  IF ( gstart == -1 ) CALL errore( 'regter', 'gstart variable not initialized', 1 )
  IF ( overlap .AND..NOT.present(spsi) ) call errore( 'rotHSw','spsi array needed with overlap=.TRUE.',1)
  !
  call start_clock('rotHSw'); !write(*,*) 'start rotHSw' ; FLUSH(6)
  !
  CALL laxlib_getval( do_distr_diag_inside_bgrp = do_distr_diag_inside_bgrp, &
       ortho_parent_comm = ortho_parent_comm )
  CALL desc_init( nstart, nx, la_proc, idesc, rank_ip, idesc_ip )
  !
  ! ... set Im[ psi(G=0) ] -  needed for numerical stability
  !
  IF ( gstart == 2 )  psi(1,1:nstart) = CMPLX( DBLE(  psi(1,1:nstart) ), 0.D0, kind=DP)
  IF ( gstart == 2 ) hpsi(1,1:nstart) = CMPLX( DBLE( hpsi(1,1:nstart) ), 0.D0, kind=DP)
  IF ( gstart == 2 ) spsi(1,1:nstart) = CMPLX( DBLE( spsi(1,1:nstart) ), 0.D0, kind=DP)
  !
  kdim = 2 * npw
  kdmx = 2 * npwx
  !
  npw_s = 1; npw_e = npw ;
  IF ( ortho_parent_comm .ne. intra_bgrp_comm ) CALL divide_all( inter_bgrp_comm, npw, npw_s, npw_e, recv_counts, displs )
  my_kdim = 2 * ( npw_e - npw_s + 1 )
  !WRITE(6,*) 'ortho_parent_comm, intra_bgrp_comm, inter_bgrp_comm ',ortho_parent_comm, intra_bgrp_comm, inter_bgrp_comm
  !WRITE(6,*) 'kdim, npw_s, npw_e, my_kdim ', kdim, npw_s, npw_e, my_kdim
  !
  ALLOCATE( vv( nx, nx ) ) ! needed later

  ALLOCATE( hh( nx, nx ), ss( nx, nx ), en( nstart ) )

  !
  ! ... Set up the Hamiltonian and Overlap matrix on the subspace :
  !
  ! ...      H_ij = <psi_i| H |psi_j>     S_ij = <psi_i| S |psi_j>
  !
  call start_clock('rotHSw:hc'); !write(*,*) 'start rotHSw:hc' ; FLUSH(6)
!  call start_clock('rotHSw:hc:b0'); CALL mp_barrier( ortho_parent_comm ); call stop_clock('rotHSw:hc:b0')

  CALL compute_distmat( hh, psi, hpsi )
  !
  IF ( overlap ) THEN
     !
     CALL compute_distmat( ss, psi, spsi )
     !
  ELSE
     !
     CALL compute_distmat( ss, psi, psi )
     !
  END IF
  call stop_clock('rotHSw:hc'); !write(*,*) 'stop rotHSw:hc' ; FLUSH(6)
  !
  ! ... Diagonalize
  !
  call start_clock('rotHSw:diag'); !write(*,*) 'start rotHSw:diag' ; FLUSH(6)
  IF ( do_distr_diag_inside_bgrp ) THEN ! NB on output of pdiaghg en and vv are the same across ortho_parent_comm
     ! only the first bgrp performs the diagonalization
     IF( my_bgrp_id == root_bgrp_id ) CALL pdiaghg( nstart, hh, ss, nx, en, vv, idesc )
     IF( nbgrp > 1 ) THEN ! results must be brodcast to the other band groups
       CALL mp_bcast( vv, root_bgrp_id, inter_bgrp_comm )
       CALL mp_bcast( en, root_bgrp_id, inter_bgrp_comm )
     ENDIF
  ELSE
     CALL pdiaghg( nstart, hh, ss, nx, en, vv, idesc )
  END IF
  call stop_clock('rotHSw:diag'); !write(*,*) 'stop rotHSw:diag' ; FLUSH(6)
  !
  e(:) = en(1:nbnd)
  !
  DEALLOCATE( en , ss , hh )
  !
  ! ... update the basis set
  !
  call start_clock('rotHSw:evc'); !write(*,*) 'start rotHSw:evc' ; FLUSH(6)
!  call start_clock('rotHSw:ev:b0'); CALL mp_barrier( ortho_parent_comm ); call stop_clock('rotHSw:ev:b0')
  !
  ib_s = 1 ; ib_e = nbnd ;
  IF ( ortho_parent_comm .ne. intra_bgrp_comm ) CALL divide_all( inter_bgrp_comm, nbnd, ib_s, ib_e, recv_counts, displs )
  my_nbnd = ib_e - ib_s + 1
  !WRITE(6,*) 'ortho_parent_comm, intra_bgrp_comm, inter_bgrp_comm ',ortho_parent_comm, intra_bgrp_comm, inter_bgrp_comm
  !WRITE(6,*) 'nbnd, ib_s, ib_e, my_nbnd ', nbnd, ib_s, ib_e, my_nbnd
  !

  ALLOCATE ( v( nstart, my_nbnd ) )

  CALL collect_v (  v )

  CALL mp_type_create_column_section(psi(1,1), 0, npwx, npwx, column_type)

  ALLOCATE( aux( npwx, my_nbnd ) )

  CALL refresh_evc( psi , v )
  CALL refresh_evc( hpsi, v )
  IF ( overlap ) THEN
     CALL refresh_evc( spsi, v )
  ELSE IF (present(spsi)) THEN
     spsi(:,1:nbnd) = psi(:,1:nbnd)
  END IF

  DEALLOCATE( aux )

  DEALLOCATE( v, vv )
  !
  CALL mp_type_free( column_type )
  !
  call stop_clock('rotHSw:evc'); !write(*,*) 'stop rotHSw:evc' ; FLUSH(6)
  !
  DEALLOCATE( idesc_ip )
  DEALLOCATE( rank_ip )
  call stop_clock('rotHSw'); !write(*,*) 'stop rotHSw' ; FLUSH(6)

  !call print_clock('rotHSw')
  !call print_clock('rotHSw:hc')
  !call print_clock('rotHSw:diag')
  !call print_clock('rotHSw:evc')

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
     INTEGER :: nr, nc, ir, ic, root, ix
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
           work = 0.0d0
           !
           nr = idesc_ip(LAX_DESC_NR, ipr, ipc )
           ir = idesc_ip(LAX_DESC_IR, ipr, ipc )
           !
           !  rank of the processor for which this block (ipr,ipc) is destinated
           !
           root = rank_ip( ipr, ipc )

           ! use blas subs. on the matrix block

           call start_clock('rotHSw:hc:comp')
           CALL DGEMM( 'T','N', nr,nc,my_kdim, 2.D0, v(npw_s,ir),kdmx,w(npw_s,ic),kdmx, 0.D0, work, nx )
           IF ( gstart == 2 .AND. npw_e == npw ) CALL DGER( nr,nc, -1.D0, v(1,ir),kdmx,w(1,ic),kdmx, work, nx )
           call stop_clock('rotHSw:hc:comp')

!           call start_clock('rotHSw:hc:b1'); CALL mp_barrier( ortho_parent_comm ); call stop_clock('rotHSw:hc:b1')
           ! accumulate result on dm of root proc.
           call start_clock('rotHSw:hc:rs')
           CALL mp_root_sum( work, dm, root, ortho_parent_comm )
           call stop_clock('rotHSw:hc:rs')

!           call start_clock('rotHSw:hc:b2'); CALL mp_barrier( ortho_parent_comm ); call stop_clock('rotHSw:hc:b2')

           !write (6,*) 'nx: ', nx
           !write (6,*) 'ipc, ipr, root', ipc, ipr, root
           !write (6,*) 'ic, nc,  ir, nr', ic, nc, ir, nr 
           !do ix=1,nx
           !   write (6,'(16f12.8)')  work(1:nx,ix)
           !enddo

        END DO
        !
     END DO
     !
     call start_clock('rotHSw:hc:sy')
     CALL laxlib_dsqmsym( nstart, dm, nx, idesc )
     call stop_clock('rotHSw:hc:sy')
!     call start_clock('rotHSw:hc:b3'); CALL mp_barrier( ortho_parent_comm ); call stop_clock('rotHSw:hc:b3')
     !
     DEALLOCATE( work )
     !
     RETURN
  END SUBROUTINE compute_distmat
  !
  SUBROUTINE collect_v ( v )
     !
     REAL(DP) :: v( nstart, * )
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, fc, root
     INTEGER :: i_s, ii_s, i_e, ii_e
     REAL(DP), ALLOCATABLE :: vtmp( :, : )

     ALLOCATE( vtmp( nx, nx ) )
     !
!     call start_clock('rotHSw:hc:b3'); CALL mp_barrier( ortho_parent_comm ); call stop_clock('rotHSw:hc:b3')
     call start_clock('rotHSw:ev:bc')
     DO ipc = 1, idesc(LAX_DESC_NPC) !  loop on column procs
        !
        nc = idesc_ip(LAX_DESC_NC, 1, ipc )
        ic = idesc_ip(LAX_DESC_IC, 1, ipc )
        !
        IF( ic <= nbnd ) THEN
           !
           nc = min( nc, nbnd - ic + 1 )
           fc = ic + nc - 1
           !
           i_s =  max ( ib_s, ic ) - ib_s + 1; ii_s = max ( ib_s, ic ) - ic + 1
           i_e =  min ( ib_e, fc ) - ib_s + 1; ii_e = min ( ib_e, fc ) - ic + 1

           IF (i_e-i_s /= ii_e - ii_s) THEN
              WRITE (*,*) ' SOMETHING WRONG '
              write (*,*) ib_s,ib_e, my_nbnd
              write (*,*) ic, fc, nc
              write (*,*) i_s, i_e
              write (*,*) ii_s, ii_e
              FLUSH(6)
           END IF

           DO ipr = 1, idesc(LAX_DESC_NPR)
              !
              nr = idesc_ip(LAX_DESC_NR, ipr, ipc )
              ir = idesc_ip(LAX_DESC_IR, ipr, ipc )
              !
              root = rank_ip( ipr, ipc )                             ! this proc has the needed block 
              IF( ipr-1 == idesc(LAX_DESC_MYR) .AND. ipc-1 == idesc(LAX_DESC_MYC) .AND. la_proc ) vtmp(:,1:nc) = vv(:,1:nc)
              CALL mp_bcast( vtmp(:,1:nc), root, ortho_parent_comm ) ! the other ones will receive it
              !
              IF (ic > ib_e .OR. fc < ib_s ) CYCLE
              !
              v( ir:ir+nr-1, i_s:i_e ) = vtmp( 1:nr, ii_s:ii_e )

           END DO
           !
        END IF
        !
     END DO
     DEALLOCATE( vtmp )
     call stop_clock('rotHSw:ev:bc')

  END SUBROUTINE collect_v
  !
  SUBROUTINE refresh_evc( psi, v )
     !
     COMPLEX(DP) :: psi(npwx,* )
     REAL(DP) :: v ( nstart, * )

     !
!     call start_clock('rotHSw:ev:b4'); CALL mp_barrier( ortho_parent_comm ); call stop_clock('rotHSw:ev:b4')
     call start_clock('rotHSw:ev:comp')
     CALL DGEMM( 'N', 'N', kdim, my_nbnd, nstart, 1.D0, psi, kdmx, v, nstart, 0.D0, aux, kdmx )
     call stop_clock('rotHSw:ev:comp')
     !
     psi(:,ib_s:ib_e) = aux(:,1:my_nbnd)
!     call start_clock('rotHSw:ev:b5'); CALL mp_barrier( ortho_parent_comm ); call stop_clock('rotHSw:ev:b5')

     call start_clock('rotHSw:ev:sum')
     IF (ortho_parent_comm.ne.intra_bgrp_comm .AND. nbgrp > 1) THEN
        CALL mp_allgather(psi(:,1:nbnd), column_type, recv_counts, displs, inter_bgrp_comm)
     END IF
     call stop_clock('rotHSw:ev:sum')

!     call start_clock('rotHSw:ev:b6'); CALL mp_barrier( ortho_parent_comm ); call stop_clock('rotHSw:ev:b6')

     RETURN
  END SUBROUTINE refresh_evc
  !
END SUBROUTINE protate_HSpsi_gamma
