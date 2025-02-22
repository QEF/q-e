!
! Copyright (C) 2001-2024 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! NOTE (Ivan Carnimeo, May, 05th, 2022): 
!   cegterg and regterg have been ported to GPU with OpenACC, 
!   the previous CUF versions (cegterg_gpu and regterg_gpu) have been removed, 
!   and now cegterg and regterg are used for both CPU and GPU execution.
!   If you want to see the previous code checkout to commit: df3080b231c5daf52295c23501fbcaa9bfc4bfcc (on Thu Apr 21 06:18:02 2022 +0000)
!
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
#if !defined(__CUDA)
#define cublasDgemm dgemm
! workaround for some old compilers that don't like CUDA fortran code
SUBROUTINE pregterg_gpu( )
end SUBROUTINE pregterg_gpu
#else
!----------------------------------------------------------------------------
!
!  Wrapper for subroutine with distributed matrixes (written by Carlo Cavazzoni)
!
!----------------------------------------------------------------------------
SUBROUTINE pregterg_gpu(h_psi_ptr, s_psi_ptr, uspp, g_psi_ptr, &  
                    npw, npwx, nvec, nvecx, evc_d, ethr, &
                    e_d, btype, notcnv, lrot, dav_iter, nhpsi )
  !----------------------------------------------------------------------------
  !
  ! ... iterative solution of the eigenvalue problem:
  !
  ! ... ( H - e S ) * evc = 0
  !
  ! ... where H is an hermitean operator, e is a real scalar,
  ! ... S is an uspp matrix, evc is a complex vector
  ! ... (real wavefunctions with only half plane waves stored)
  !
  USE util_param,        ONLY : DP, stdout
  USE mp_bands_util,     ONLY : intra_bgrp_comm, inter_bgrp_comm, root_bgrp_id, nbgrp, my_bgrp_id
  USE mp_bands_util,     ONLY : gstart
  USE mp,                ONLY : mp_bcast, mp_root_sum, mp_sum
  !
  IMPLICIT NONE
  !
  include 'laxlib.fh'
  !
  INTEGER, INTENT(IN) :: npw, npwx, nvec, nvecx
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix evc, as declared in the calling pgm unit
    ! integer number of searched low-lying roots
    ! maximum dimension of the reduced basis set
    !    (the basis set is refreshed when its dimension would exceed nvecx)
  COMPLEX(DP), INTENT(INOUT) :: evc_d(npwx,nvec)
    !  evc   contains the  refined estimates of the eigenvectors
  REAL(DP), INTENT(IN) :: ethr
    ! energy threshold for convergence: root improvement is stopped,
    ! when two consecutive estimates of the root differ by less than ethr.
  LOGICAL, INTENT(IN) :: uspp
    ! if .FALSE. : S|psi> not needed
  INTEGER, INTENT(IN) :: btype(nvec)
    ! band type ( 1 = occupied, 0 = empty )
  LOGICAL, INTENT(IN) :: lrot
    ! .TRUE. if the wfc have already been rotated
  REAL(DP), INTENT(OUT) :: e_d(nvec)
#if defined(__CUDA)
  attributes(DEVICE) :: e_d
#endif
    ! contains the estimated roots.
  INTEGER, INTENT(OUT) :: dav_iter, notcnv
    ! integer  number of iterations performed
    ! number of unconverged roots
  INTEGER, INTENT(OUT) :: nhpsi
    ! number of individual Hpsi made
  !
  ! ... LOCAL variables
  !
  !COMPLEX(DP), ALLOCATABLE :: evc(:,:)
  REAL(DP), ALLOCATABLE :: e(:)
  
  INTEGER, PARAMETER :: maxter = 20
    ! maximum number of iterations
  !
  INTEGER :: kter, nbase, np, n, m, nb1, i, j, k
    ! counter on iterations
    ! dimension of the reduced basis
    ! counter on the reduced basis vectors
    ! do-loop counters
  INTEGER :: ierr
  REAL(DP), ALLOCATABLE :: ew(:)
  REAL(DP), ALLOCATABLE :: hl(:,:), sl(:,:), vl(:,:)
    ! Hamiltonian on the reduced basis
    ! S matrix on the reduced basis
    ! eigenvectors of the Hamiltonian
    ! eigenvalues of the reduced hamiltonian
  COMPLEX(DP), ALLOCATABLE :: psi(:,:), hpsi(:,:), spsi(:,:)
    ! work space, contains psi
    ! the product of H and psi
    ! the product of S and psi
  LOGICAL, ALLOCATABLE :: conv(:)
    ! true if the root is converged
  REAL(DP) :: empty_ethr 
    ! threshold for empty bands
  INTEGER :: npw2, npwx2
  INTEGER :: idesc(LAX_DESC_SIZE), idesc_old(LAX_DESC_SIZE)
  INTEGER, ALLOCATABLE :: irc_ip( : )
  INTEGER, ALLOCATABLE :: nrc_ip( : )
  INTEGER, ALLOCATABLE :: rank_ip( :, : )
    ! matrix distribution descriptors
  INTEGER :: nx
    ! maximum local block dimension
  LOGICAL :: la_proc
    ! flag to distinguish procs involved in linear algebra
  INTEGER, ALLOCATABLE :: notcnv_ip( : )
  INTEGER, ALLOCATABLE :: ic_notcnv( : )
  !
  INTEGER :: np_ortho(2), ortho_parent_comm
  LOGICAL :: do_distr_diag_inside_bgrp
  !
  REAL(DP), EXTERNAL :: ddot
  !
  EXTERNAL  h_psi_ptr, s_psi_ptr, g_psi_ptr
    ! h_psi_ptr(npwx,npw,nvec,psi,hpsi)
    !     calculates H|psi> 
    ! s_psi_ptr(npwx,npw,nvec,psi,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx,nvec)
    ! g_psi_ptr(npwx,npw,notcnv,psi,e)
    !    calculates (diag(h)-e)^-1 * psi, diagonal approx. to (h-e)^-1*psi
    !    the first nvec columns contain the trial eigenvectors
  INTEGER :: info
    ! infor from buffer 
  !
  !
  CALL start_clock( 'regterg' )
  ! 
  CALL laxlib_getval( np_ortho = np_ortho, ortho_parent_comm = ortho_parent_comm, &
    do_distr_diag_inside_bgrp = do_distr_diag_inside_bgrp )
  ! 
  IF ( nvec > nvecx / 2 ) CALL errore( 'pregter', 'nvecx is too small', 1 )
  !
  IF ( gstart == -1 ) CALL errore( 'pregter', 'gstart variable not initialized', 1 )
  !
  ! ... threshold for empty bands
  !
  empty_ethr = MAX( ( ethr * 5.D0 ), 1.D-5 )
  !
  ALLOCATE(  e( nvec ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pregterg ',' cannot allocate e (host) ', ABS(ierr) )
  !
  ALLOCATE( psi(  npwx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'pregterg ',' cannot allocate psi ', ABS(ierr) )
  !
  ALLOCATE( hpsi( npwx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'pregterg ',' cannot allocate hpsi ', ABS(ierr) )
  !
  IF ( uspp ) THEN
     ALLOCATE( spsi( npwx, nvecx ), STAT=ierr )
     IF( ierr /= 0 ) &
          CALL errore( 'pregterg ',' cannot allocate spsi ', ABS(ierr) )
     !$acc enter data create(spsi)
  END IF
  !
  ! ... Initialize the matrix descriptor
  !
  ALLOCATE( ic_notcnv( np_ortho(2) ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'pregterg ',' cannot allocate ic_notcnv ', ABS(ierr) )
  !
  ALLOCATE( notcnv_ip( np_ortho(2) ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'pregterg ',' cannot allocate notcnv_ip ', ABS(ierr) )
  !
  ALLOCATE( irc_ip( np_ortho(1) ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'pregterg ',' cannot allocate irc_ip ', ABS(ierr) )
  !
  ALLOCATE( nrc_ip( np_ortho(1) ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'pregterg ',' cannot allocate nrc_ip ', ABS(ierr) )
  !
  ALLOCATE( rank_ip( np_ortho(1), np_ortho(2) ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'pregterg ',' cannot allocate rank_ip ', ABS(ierr) )
  !
  CALL desc_init( nvec, nx, la_proc, idesc, rank_ip, irc_ip, nrc_ip  )
  !
  IF( la_proc ) THEN
     !
     ! only procs involved in the diagonalization need to allocate local 
     ! matrix block.
     !
     ALLOCATE( vl( nx , nx ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( 'pregterg ',' cannot allocate vl ', ABS(ierr) )
     !
     ALLOCATE( sl( nx , nx ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( 'pregterg ',' cannot allocate sl ', ABS(ierr) )
     !
     ALLOCATE( hl( nx , nx ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( 'pregterg ',' cannot allocate hl ', ABS(ierr) )
     !
  ELSE
     !
     ALLOCATE( vl( 1 , 1 ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( 'pregterg ',' cannot allocate vl ', ABS(ierr) )
     !
     ALLOCATE( sl( 1 , 1 ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( 'pregterg ',' cannot allocate sl ', ABS(ierr) )
     !
     ALLOCATE( hl( 1 , 1 ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( 'pregterg ',' cannot allocate hl ', ABS(ierr) )
     !
  END IF
  !
  ALLOCATE( ew( nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'pregterg ',' cannot allocate ew ', ABS(ierr) )
  !$acc enter data create(ew)
  !
  ALLOCATE( conv( nvec ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'pregterg ',' cannot allocate conv ', ABS(ierr) )
  !
  npw2  = 2*npw
  npwx2  = 2*npwx
  notcnv = nvec
  nbase  = nvec
  conv   = .FALSE.
  !
  !$acc enter data create(psi, hpsi)
  !$acc update host( evc_d)
  !$acc kernels
  psi(:,1:nvec) = evc_d(:,1:nvec)
  !$acc end kernels
  ! ... set Im[ psi(G=0) ] -  needed for numerical stability
  !$acc kernels
  IF ( gstart == 2 ) psi(1,1:nvec) = CMPLX( DBLE( psi(1,1:nvec) ), 0.D0 ,kind=DP)
  !$acc end kernels
  !$acc update host( psi)
  !
  ! ... hpsi contains h times the basis vectors
  !
  CALL h_psi_ptr( npwx, npw, nvec, psi, hpsi )  ; nhpsi = nvec
  !$acc update host(hpsi)
  !
  IF ( uspp ) THEN
     CALL s_psi_ptr( npwx, npw, nvec, psi, spsi )
     !$acc update host(spsi)
  END IF
  !
  ! ... hl contains the projection of the hamiltonian onto the reduced
  ! ... space, vl contains the eigenvectors of hl. Remember hl, vl and sl
  ! ... are all distributed across processors, global replicated matrixes
  ! ... here are never allocated
  !
  CALL start_clock( 'regterg:init' )

  CALL compute_distmat( hl, psi, hpsi )
  !
  IF ( uspp ) THEN
     !
     CALL compute_distmat( sl, psi, spsi )
     !
  ELSE
     !
     CALL compute_distmat( sl, psi, psi )
     !
  END IF
  CALL stop_clock( 'regterg:init' )
  !
  IF ( lrot ) THEN
     !
     CALL set_e_from_h()
     e_d = e
     !
     CALL set_to_identity( vl, idesc )
     !
  ELSE
     !
     ! ... diagonalize the reduced hamiltonian
     !     Calling block parallel algorithm
     !
     CALL start_clock( 'regterg:diag' )
     IF ( do_distr_diag_inside_bgrp ) THEN ! NB on output of pdiaghg ew and vl are the same across ortho_parent_comm
        ! only the first bgrp performs the diagonalization
        IF( my_bgrp_id == root_bgrp_id ) CALL pdiaghg( nbase, hl, sl, nx, ew, vl, idesc )
        IF( nbgrp > 1 ) THEN ! results must be brodcast to the other band groups
           CALL mp_bcast( vl, root_bgrp_id, inter_bgrp_comm )
           CALL mp_bcast( ew, root_bgrp_id, inter_bgrp_comm )
        ENDIF
     ELSE
        CALL pdiaghg( nbase, hl, sl, nx, ew, vl, idesc )
     END IF
     CALL stop_clock( 'regterg:diag' )
     !
     e(1:nvec) = ew(1:nvec)
     e_d(1:nvec) = ew(1:nvec)
     !
  END IF
  !
  ! ... iterate
  !
  iterate: DO kter = 1, maxter
     !
     dav_iter = kter
     !
     CALL start_clock( 'regterg:update' )
     !
     CALL reorder_v()
     !
     nb1 = nbase + 1
     !
     ! ... expand the basis set with new basis vectors ( H - e*S )|psi> ...
     !
     CALL hpsi_dot_v()
     !
     CALL stop_clock( 'regterg:update' )
     !
     ! ... approximate inverse iteration
     !
     !$acc update device(ew,psi)
     CALL g_psi_ptr( npwx, npw, notcnv, 1, psi(1,nb1), ew(nb1) )
     !$acc update host(psi)
     !
     ! ... "normalize" correction vectors psi(:,nb1:nbase+notcnv) in 
     ! ... order to improve numerical stability of subspace diagonalization 
     ! ... (rdiaghg) ew is used as work array :
     !
     ! ...         ew = <psi_i|psi_i>,  i = nbase + 1, nbase + notcnv
     !
     DO n = 1, notcnv
        !
        ew(n) = 2.D0 * ddot( npw2, psi(1,nbase+n), 1, psi(1,nbase+n), 1 )
        !
        IF ( gstart == 2 ) ew(n) = ew(n) - psi(1,nbase+n) * psi(1,nbase+n)
        !
     END DO
     !
     CALL mp_sum( ew( 1:notcnv ), intra_bgrp_comm )
     !
     DO n = 1, notcnv
        !
        psi(:,nbase+n) = psi(:,nbase+n) / SQRT( ew(n) )
        ! ... set Im[ psi(G=0) ] -  needed for numerical stability
        IF ( gstart == 2 ) psi(1,nbase+n) = CMPLX( DBLE(psi(1,nbase+n)), 0.D0 ,kind=DP)
        !
     END DO
     !
     ! ... here compute the hpsi and spsi of the new functions
     !
     !$acc update device(psi)
     CALL h_psi_ptr( npwx, npw, notcnv, psi(1,nb1), hpsi(1,nb1) ) ; nhpsi = nhpsi + notcnv
     !$acc update host(hpsi)
     !
     IF ( uspp ) THEN
        CALL s_psi_ptr( npwx, npw, notcnv, psi(1,nb1), spsi(1,nb1) )
        !$acc update host(spsi)
     END IF
     !
     ! ... update the reduced hamiltonian
     !
     CALL start_clock( 'regterg:overlap' )
     !
     ! we need to save the old descriptor in order to redistribute matrices 
     !
     idesc_old = idesc
     !
     ! ... RE-Initialize the matrix descriptor
     !
     CALL desc_init( nbase+notcnv, nx, la_proc, idesc, rank_ip, irc_ip, nrc_ip  )
     !
     IF( la_proc ) THEN

        !  redistribute hl and sl (see dsqmred), since the dimension of the subspace has changed
        !
        vl = hl
        DEALLOCATE( hl )
        ALLOCATE( hl( nx , nx ), STAT=ierr )
        IF( ierr /= 0 ) &
           CALL errore( 'pregterg ',' cannot allocate hl ', ABS(ierr) )

        CALL laxlib_dsqmred( nbase, vl, idesc_old(LAX_DESC_NRCX), idesc_old, nbase+notcnv, hl, nx, idesc )

        vl = sl
        DEALLOCATE( sl )
        ALLOCATE( sl( nx , nx ), STAT=ierr )
        IF( ierr /= 0 ) &
           CALL errore( 'pregterg ',' cannot allocate sl ', ABS(ierr) )

        CALL laxlib_dsqmred( nbase, vl, idesc_old(LAX_DESC_NRCX), idesc_old, nbase+notcnv, sl, nx, idesc )

        DEALLOCATE( vl )
        ALLOCATE( vl( nx , nx ), STAT=ierr )
        IF( ierr /= 0 ) &
           CALL errore( 'pregterg ',' cannot allocate vl ', ABS(ierr) )

     END IF
     !
     !
     CALL update_distmat( hl, psi, hpsi )
     !
     IF ( uspp ) THEN
        !
        CALL update_distmat( sl, psi, spsi )
        !
     ELSE
        !
        CALL update_distmat( sl, psi, psi )
        !
     END IF
     !
     CALL stop_clock( 'regterg:overlap' )
     !
     nbase = nbase + notcnv
     !
     ! ... diagonalize the reduced hamiltonian
     !     Call block parallel algorithm
     !
     CALL start_clock( 'regterg:diag' )
     IF ( do_distr_diag_inside_bgrp ) THEN ! NB on output of pdiaghg ew and vl are the same across ortho_parent_comm
        ! only the first bgrp performs the diagonalization
        IF( my_bgrp_id == root_bgrp_id ) CALL pdiaghg( nbase, hl, sl, nx, ew, vl, idesc )
        IF( nbgrp > 1 ) THEN ! results must be brodcast to the other bnd groups
           CALL mp_bcast( vl, root_bgrp_id, inter_bgrp_comm )
           CALL mp_bcast( ew, root_bgrp_id, inter_bgrp_comm )
        ENDIF
     ELSE
        CALL pdiaghg( nbase, hl, sl, nx, ew, vl, idesc )
     END IF
     CALL stop_clock( 'regterg:diag' )
     !
     ! ... test for convergence
     !
     WHERE( btype(1:nvec) == 1 )
        !
        conv(1:nvec) = ( ( ABS( ew(1:nvec) - e(1:nvec) ) < ethr ) )
        !
     ELSEWHERE
        !
        conv(1:nvec) = ( ( ABS( ew(1:nvec) - e(1:nvec) ) < empty_ethr ) )
        !
     END WHERE
     ! ... next line useful for band parallelization of exact exchange
     IF ( nbgrp > 1 ) CALL mp_bcast(conv,root_bgrp_id,inter_bgrp_comm)
     !
     notcnv = COUNT( .NOT. conv(:) )
     !
     e(1:nvec) = ew(1:nvec)
     e_d(1:nvec) = e(1:nvec)
     !
     ! ... if overall convergence has been achieved, or the dimension of
     ! ... the reduced basis set is becoming too large, or in any case if
     ! ... we are at the last iteration refresh the basis set. i.e. replace
     ! ... the first nvec elements with the current estimate of the
     ! ... eigenvectors;  set the basis dimension to nvec.
     !
     IF ( notcnv == 0 .OR. nbase+notcnv > nvecx .OR. dav_iter == maxter ) THEN
        !
        CALL start_clock( 'regterg:last' )
        !
        CALL refresh_evc()
        !$acc update device(evc_d)
        !
        IF ( notcnv == 0 ) THEN
           !
           ! ... all roots converged: return
           !
           CALL stop_clock( 'regterg:last' )
           !
           EXIT iterate
           !
        ELSE IF ( dav_iter == maxter ) THEN
           !
           ! ... last iteration, some roots not converged: return
           !
           !!!WRITE( stdout, '(5X,"WARNING: ",I5, &
           !!!     &   " eigenvalues not converged in regterg")' ) notcnv
           !
           CALL stop_clock( 'regterg:last' )
           !
           EXIT iterate
           !
        END IF
        !
        ! ... refresh psi, H*psi and S*psi
        !
        psi(:,1:nvec) = evc_d(:,1:nvec)
        !$acc update device(psi)
        !
        IF ( uspp ) THEN
           CALL refresh_spsi()
           !$acc update device(spsi)
        END IF
        !
        CALL refresh_hpsi()
        !$acc update device(hpsi)
        !
        ! ... refresh the reduced hamiltonian
        !
        nbase = nvec
        !
        CALL desc_init( nvec, nx, la_proc, idesc, rank_ip, irc_ip, nrc_ip  )
        !
        IF( la_proc ) THEN
           !
           ! note that nx has been changed by desc_init
           ! we need to re-alloc with the new size.
           !
           DEALLOCATE( vl, hl, sl )
           ALLOCATE( vl( nx, nx ), STAT=ierr )
           IF( ierr /= 0 ) &
              CALL errore( 'pregterg ',' cannot allocate vl ', ABS(ierr) )
           ALLOCATE( hl( nx, nx ), STAT=ierr )
           IF( ierr /= 0 ) &
              CALL errore( 'pregterg ',' cannot allocate hl ', ABS(ierr) )
           ALLOCATE( sl( nx, nx ), STAT=ierr )
           IF( ierr /= 0 ) &
              CALL errore( 'pregterg ',' cannot allocate sl ', ABS(ierr) )
           !
        END IF
        !
        CALL set_h_from_e( )
        !
        CALL set_to_identity( vl, idesc )
        CALL set_to_identity( sl, idesc )
        !
        CALL stop_clock( 'regterg:last' )
        !
     END IF
     !
  END DO iterate
  !
  DEALLOCATE( vl, hl, sl )
  !
  DEALLOCATE( rank_ip )
  DEALLOCATE( ic_notcnv )
  DEALLOCATE( irc_ip )
  DEALLOCATE( nrc_ip )
  DEALLOCATE( notcnv_ip )
  DEALLOCATE( conv )
  DEALLOCATE( e )
  !$acc exit data delete(ew)
  DEALLOCATE( ew )
  !
  IF ( uspp ) THEN
     !$acc exit data delete(spsi)
     DEALLOCATE( spsi )
  END IF
  !$acc exit data delete (psi, hpsi)
  DEALLOCATE( hpsi )
  DEALLOCATE( psi )  
  !
  CALL stop_clock( 'regterg' )
  !call print_clock( 'regterg' )
  !call print_clock( 'regterg:init' )
  !call print_clock( 'regterg:diag' )
  !call print_clock( 'regterg:update' )
  !call print_clock( 'regterg:overlap' )
  !call print_clock( 'regterg:last' )

  !
  RETURN
  !
  !
CONTAINS
  !
  !
  SUBROUTINE set_to_identity( distmat, idesc )
     INTEGER, INTENT(IN)  :: idesc(LAX_DESC_SIZE)
     REAL(DP), INTENT(OUT) :: distmat(:,:)
     INTEGER :: i
     distmat = 0_DP
     IF( idesc(LAX_DESC_MYC) == idesc(LAX_DESC_MYR) .AND. idesc(LAX_DESC_ACTIVE_NODE) > 0 ) THEN
        DO i = 1, idesc(LAX_DESC_NC)
           distmat( i, i ) = 1_DP
        END DO
     END IF 
     RETURN
  END SUBROUTINE set_to_identity
  !
  !
  SUBROUTINE reorder_v()
     !
     INTEGER :: ipc, ipr
     INTEGER :: nc, ic
     INTEGER :: nl, npl
     !
     np = 0
     !
     notcnv_ip = 0
     !
     n = 0
     !
     DO ipc = 1, idesc(LAX_DESC_NPC)
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        npl = 0
        !
        IF( ic <= nvec ) THEN
           !
           DO nl = 1, min( nvec - ic + 1, nc )
              !
              n  = n  + 1
              !
              IF ( .NOT. conv(n) ) THEN
                 !
                 ! ... this root not yet converged ... 
                 !
                 np  = np  + 1
                 npl = npl + 1
                 IF( npl == 1 ) ic_notcnv( ipc ) = np
                 !
                 ! ... reorder eigenvectors so that coefficients for unconverged
                 ! ... roots come first. This allows to use quick matrix-matrix 
                 ! ... multiplications to set a new basis vector (see below)
                 !
                 notcnv_ip( ipc ) = notcnv_ip( ipc ) + 1
                 !
                 IF ( npl /= nl ) THEN
                    IF( la_proc .AND. idesc(LAX_DESC_MYC) == ipc-1 ) THEN
                       vl( :, npl) = vl( :, nl )
                    END IF
                 END IF
                 !
                 ! ... for use in g_psi_ptr
                 !
                 ew(nbase+np) = e(n)
                 !   
              END IF
              !
           END DO
           !
        END IF
        !
     END DO
     !
  END SUBROUTINE reorder_v
  !
  !
  SUBROUTINE hpsi_dot_v()
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, notcl, root, np
     REAL(DP), ALLOCATABLE :: vtmp( :, : )
     COMPLEX(DP), ALLOCATABLE :: ptmp( :, : )
     REAL(DP) :: beta

     ALLOCATE( vtmp( nx, nx ) )
     ALLOCATE( ptmp( npwx, nx ) )

     DO ipc = 1, idesc(LAX_DESC_NPC)
        !
        IF( notcnv_ip( ipc ) > 0 ) THEN

           notcl = notcnv_ip( ipc )
           ic    = ic_notcnv( ipc ) 

           beta = 0.0d0

           DO ipr = 1, idesc(LAX_DESC_NPR)
              !
              nr = nrc_ip( ipr )
              ir = irc_ip( ipr )
              !
              root = rank_ip( ipr, ipc )

              IF( ipr-1 == idesc(LAX_DESC_MYR) .AND. ipc-1 == idesc(LAX_DESC_MYC) .AND. la_proc ) THEN
                 vtmp(:,1:notcl) = vl(:,1:notcl)
              END IF

              CALL mp_bcast( vtmp(:,1:notcl), root, ortho_parent_comm )
              ! 
              IF ( uspp ) THEN
                 !
                 CALL DGEMM( 'N', 'N', npw2, notcl, nr, 1.D0, &
                    spsi( 1, ir ), npwx2, vtmp, nx, beta, psi(1,nb1+ic-1), npwx2 )
                 !
              ELSE
                 !
                 CALL DGEMM( 'N', 'N', npw2, notcl, nr, 1.D0, &
                    psi( 1, ir ), npwx2, vtmp, nx, beta, psi(1,nb1+ic-1), npwx2 )
                 !
              END IF
              !
              CALL DGEMM( 'N', 'N', npw2, notcl, nr, 1.D0, &
                      hpsi( 1, ir ), npwx2, vtmp, nx, beta, ptmp, npwx2 )

              beta = 1.0d0

           END DO

           DO np = 1, notcl
              !
              psi(1:npw,nbase+np+ic-1) = ptmp(1:npw,np) - ew(nbase+np+ic-1) * psi(1:npw,nbase+np+ic-1)
              !
           END DO
           !
        END IF
        !
     END DO

     
     DEALLOCATE( vtmp )
     DEALLOCATE( ptmp )

     RETURN
  END SUBROUTINE hpsi_dot_v
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
     DO ipc = 1, idesc(LAX_DESC_NPC)
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        IF( ic <= nvec ) THEN
           !
           nc = min( nc, nvec - ic + 1 )
           !
           beta = 0.0d0

           DO ipr = 1, idesc(LAX_DESC_NPR)
              !
              nr = nrc_ip( ipr )
              ir = irc_ip( ipr )
              !
              root = rank_ip( ipr, ipc )

              IF( ipr-1 == idesc(LAX_DESC_MYR) .AND. ipc-1 == idesc(LAX_DESC_MYC) .AND. la_proc ) THEN
                 !
                 !  this proc sends his block
                 ! 
                 CALL mp_bcast( vl(:,1:nc), root, ortho_parent_comm )
                 CALL DGEMM( 'N', 'N', npw2, nc, nr, 1.D0, &
                          psi(1,ir), npwx2, vl, nx, beta, evc_d(1,ic), npwx2 )
              ELSE
                 !
                 !  all other procs receive
                 ! 
                 CALL mp_bcast( vtmp(:,1:nc), root, ortho_parent_comm )
                 CALL DGEMM( 'N', 'N', npw2, nc, nr, 1.D0, &
                          psi(1,ir), npwx2, vtmp, nx, beta, evc_d(1,ic), npwx2 )
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
  !
  SUBROUTINE refresh_spsi( )
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     REAL(DP), ALLOCATABLE :: vtmp( :, : )
     REAL(DP) :: beta

     ALLOCATE( vtmp( nx, nx ) )
     !
     DO ipc = 1, idesc(LAX_DESC_NPC)
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        IF( ic <= nvec ) THEN
           !
           nc = min( nc, nvec - ic + 1 )
           !
           beta = 0_DP
           !
           DO ipr = 1, idesc(LAX_DESC_NPR)
              !
              nr = nrc_ip( ipr )
              ir = irc_ip( ipr )
              !
              root = rank_ip( ipr, ipc )

              IF( ipr-1 == idesc(LAX_DESC_MYR) .AND. ipc-1 == idesc(LAX_DESC_MYC) .AND. la_proc ) THEN
                 !
                 !  this proc sends his block
                 ! 
                 CALL mp_bcast( vl(:,1:nc), root, ortho_parent_comm )
                 CALL DGEMM( 'N', 'N', npw2, nc, nr, 1.D0, &
                          spsi(1,ir), npwx2, vl, nx, beta, psi(1,nvec+ic), npwx2 )
              ELSE
                 !
                 !  all other procs receive
                 ! 
                 CALL mp_bcast( vtmp(:,1:nc), root, ortho_parent_comm )
                 CALL DGEMM( 'N', 'N', npw2, nc, nr, 1.D0, &
                          spsi(1,ir), npwx2, vtmp, nx, beta, psi(1,nvec+ic), npwx2 )
              END IF
              ! 
              beta = 1_DP

           END DO
           !
        END IF
        !
     END DO
     !
     spsi(:,1:nvec) = psi(:,nvec+1:nvec+nvec)
     !
     DEALLOCATE( vtmp )

     RETURN
  END SUBROUTINE refresh_spsi
  !
  !
  !
  SUBROUTINE refresh_hpsi( )
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     REAL(DP), ALLOCATABLE :: vtmp( :, : )
     REAL(DP) :: beta

     ALLOCATE( vtmp( nx, nx ) )
     !
     DO ipc = 1, idesc(LAX_DESC_NPC)
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        IF( ic <= nvec ) THEN
           !
           nc = min( nc, nvec - ic + 1 )
           !
           beta = 0.0d0
           !
           DO ipr = 1, idesc(LAX_DESC_NPR)
              !
              nr = nrc_ip( ipr )
              ir = irc_ip( ipr )
              !
              root = rank_ip( ipr, ipc )

              IF( ipr-1 == idesc(LAX_DESC_MYR) .AND. ipc-1 == idesc(LAX_DESC_MYC) .AND. la_proc ) THEN
                 !
                 !  this proc sends his block
                 ! 
                 CALL mp_bcast( vl(:,1:nc), root, ortho_parent_comm )
                 CALL DGEMM( 'N', 'N', npw2, nc, nr, 1.D0, &
                          hpsi(1,ir), npwx2, vl, nx, beta, psi(1,nvec+ic), npwx2 )
              ELSE
                 !
                 !  all other procs receive
                 ! 
                 CALL mp_bcast( vtmp(:,1:nc), root, ortho_parent_comm )
                 CALL DGEMM( 'N', 'N', npw2, nc, nr, 1.D0, &
                          hpsi(1,ir), npwx2, vtmp, nx, beta, psi(1,nvec+ic), npwx2 )
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

     hpsi(:,1:nvec) = psi(:,nvec+1:nvec+nvec)

     RETURN
  END SUBROUTINE refresh_hpsi
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
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        DO ipr = 1, ipc ! use symmetry for the loop on row procs
           !
           nr = nrc_ip( ipr )
           ir = irc_ip( ipr )
           !
           !  rank of the processor for which this block (ipr,ipc) is destinated
           !
           root = rank_ip( ipr, ipc )

           ! use blas subs. on the matrix block

           CALL DGEMM( 'T', 'N', nr, nc, npw2, 2.D0 , &
                       v(1,ir), npwx2, w(1,ic), npwx2, 0.D0, work, nx )

           IF ( gstart == 2 ) &
              CALL DGER( nr, nc, -1.D0, v(1,ir), npwx2, w(1,ic), npwx2, work, nx )

           ! accumulate result on dm of root proc.

           CALL mp_root_sum( work, dm, root, ortho_parent_comm )

        END DO
        !
     END DO
     IF (ortho_parent_comm.ne.intra_bgrp_comm .and. nbgrp > 1) dm = dm/nbgrp
     !
     CALL laxlib_dsqmsym( nbase, dm, nx, idesc )
     !
     DEALLOCATE( work )
     !
     RETURN
  END SUBROUTINE compute_distmat
  !
  !
  SUBROUTINE update_distmat( dm, v, w )
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root, icc, ii
     REAL(DP)    :: dm( :, : )
     COMPLEX(DP) :: v(:,:), w(:,:)
     REAL(DP), ALLOCATABLE :: vtmp( :, : )

     ALLOCATE( vtmp( nx, nx ) )
     !
     vtmp = 0.0d0
     !
     DO ipc = 1, idesc(LAX_DESC_NPC)
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        IF( ic+nc-1 >= nb1 ) THEN

           nc = MIN( nc, ic+nc-1 - nb1 + 1 )
           IF( ic >= nb1 ) THEN
              ii = ic
              icc = 1
           ELSE
              ii = nb1
              icc = nb1-ic+1
           END IF

           DO ipr = 1, ipc ! desc%npr use symmetry
              !
              nr = nrc_ip( ipr )
              ir = irc_ip( ipr )
              !
              root = rank_ip( ipr, ipc )

              CALL DGEMM( 'T', 'N', nr, nc, npw2, 2.D0, v( 1, ir ), &
                          npwx2, w(1,ii), npwx2, 0.D0, vtmp, nx )
              !
              IF ( gstart == 2 ) &
                 CALL DGER( nr, nc, -1.D0, v( 1, ir ), npwx2, w(1,ii), npwx2, vtmp, nx )
              IF (ortho_parent_comm.ne.intra_bgrp_comm .and. nbgrp > 1) vtmp = vtmp/nbgrp

              IF(  (idesc(LAX_DESC_ACTIVE_NODE) > 0) .AND. &
                   (ipr-1 == idesc(LAX_DESC_MYR)) .AND. (ipc-1 == idesc(LAX_DESC_MYC)) ) THEN
                 CALL mp_root_sum( vtmp(:,1:nc), dm(:,icc:icc+nc-1), root, ortho_parent_comm )
              ELSE
                 CALL mp_root_sum( vtmp(:,1:nc), dm, root, ortho_parent_comm )
              END IF


           END DO
           !
        END IF
        !
     END DO
     !
     CALL laxlib_dsqmsym( nbase+notcnv, dm, nx, idesc )
     !
     DEALLOCATE( vtmp )
     RETURN
  END SUBROUTINE update_distmat
  !
  !
  !
  SUBROUTINE set_e_from_h()
     INTEGER :: nc, ic, i
     e(1:nbase) = 0.0d0
     IF( idesc(LAX_DESC_MYC) == idesc(LAX_DESC_MYR) .AND. la_proc ) THEN
        nc = idesc(LAX_DESC_NC)
        ic = idesc(LAX_DESC_IC)
        DO i = 1, nc
           e( i + ic - 1 ) = hl( i, i )
        END DO
     END IF
     CALL mp_sum( e(1:nbase), ortho_parent_comm )
     RETURN
  END SUBROUTINE set_e_from_h
  !
  SUBROUTINE set_h_from_e()
     INTEGER :: nc, ic, i
     IF( la_proc ) THEN
        hl = 0.0d0
        IF( idesc(LAX_DESC_MYC) == idesc(LAX_DESC_MYR) ) THEN
           nc = idesc(LAX_DESC_NC)
           ic = idesc(LAX_DESC_IC)
           DO i = 1, nc
              hl(i,i) = e( i + ic - 1 )
           END DO
        END IF
     END IF
     RETURN
  END SUBROUTINE set_h_from_e
  !
END SUBROUTINE pregterg_gpu

#endif
