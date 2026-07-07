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
!
!
!----------------------------------------------------------------------------
SUBROUTINE regterg(  h_psi_ptr, s_psi_ptr, uspp, g_psi_ptr, &
                    npw, npwx, nvec, nvecx, evc, ethr, &
                    e, btype, notcnv, lrot, dav_iter, nhpsi )
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
#if defined(__CUDA)
  use cublas
#endif
  USE util_param,    ONLY : DP, stdout
  USE mp_bands_util, ONLY : intra_bgrp_comm, inter_bgrp_comm, root_bgrp_id, &
          nbgrp, my_bgrp_id, me_bgrp, root_bgrp
  USE mp_bands_util, ONLY : gstart
  USE mp,            ONLY : mp_sum, mp_bcast
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
  COMPLEX(DP), INTENT(INOUT) :: evc(npwx,nvec)
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
  REAL(DP), INTENT(OUT) :: e(nvec)
    ! contains the estimated roots.
  INTEGER, INTENT(OUT) :: dav_iter, notcnv
    ! integer  number of iterations performed
    ! number of unconverged roots
  INTEGER, INTENT(OUT) :: nhpsi
    ! number of individual Hpsi made
  !
  ! ... LOCAL variables
  !
  INTEGER, PARAMETER :: maxter = 20
    ! maximum number of iterations
  !
  INTEGER :: kter, nbase, np, npw2, npwx2, n, m, nb1, nbn
    ! counter on iterations
    ! dimension of the reduced basis
    ! counter on the reduced basis vectors
    ! do-loop counters
    ! counter on the bands
  INTEGER :: n_start, n_end, my_n
  INTEGER :: ierr
  REAL(DP), ALLOCATABLE :: hr(:,:), sr(:,:), vr(:,:), ew(:)
  !$acc declare device_resident(hr, sr, vr)
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
  INTEGER :: i,j,k
  !
  REAL(DP), EXTERNAL :: MYDDOT_VECTOR_GPU 
  !$acc routine(MYDDOT_VECTOR_GPU) vector
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
  !
  CALL start_clock( 'regterg' ) !; write(6,*) 'enter regterg' ; FLUSH(6)
  ! 
  IF ( nvec > nvecx / 2 ) CALL errore( 'regter', 'nvecx is too small', 1 )
  !
  IF ( gstart == -1 ) CALL errore( 'regter', 'gstart variable not initialized', 1 )
  !
  ! ... threshold for empty bands
  !
  empty_ethr = MAX( ( ethr * 5.D0 ), 1.D-5 )
  !
  npw2  = 2*npw
  npwx2  = 2*npwx
  !
  ALLOCATE( psi(  npwx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'regterg ',' cannot allocate psi ', ABS(ierr) )
  ALLOCATE( hpsi( npwx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'regterg ',' cannot allocate hpsi ', ABS(ierr) )
  !$acc enter data create(psi, hpsi)
  !
  IF ( uspp ) THEN
     ALLOCATE( spsi( npwx, nvecx ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' regterg ',' cannot allocate spsi ', ABS(ierr) )
     !$acc enter data create(spsi)
  END IF
  !
  ALLOCATE( sr( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'regterg ',' cannot allocate sr ', ABS(ierr) )
  ALLOCATE( hr( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'regterg ',' cannot allocate hr ', ABS(ierr) )
  ALLOCATE( vr( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'regterg ',' cannot allocate vr ', ABS(ierr) )
  ALLOCATE( ew( nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'regterg ',' cannot allocate ew ', ABS(ierr) )
  !$acc enter data create(ew)
  ALLOCATE( conv( nvec ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'regterg ',' cannot allocate conv ', ABS(ierr) )
  !
  notcnv = nvec
  nbase  = nvec
  conv   = .FALSE.
  !  
  !$acc kernels
  hpsi = ZERO
  psi  = ZERO
  IF ( uspp ) spsi = ZERO
  !$acc end kernels
  !
  !$acc parallel vector_length(64) 
  !$acc loop gang independent
  DO k=1,nvec
     psi(1,k) = evc(1,k)
  ! ... set Im[ psi(G=0) ] -  needed for numerical stability
     IF (gstart == 2) psi(1,k) = CMPLX( DBLE( psi(1,k) ), 0.D0 ,kind=DP)
     !$acc loop vector 
     DO i=2,npwx
         psi(i,k) = evc(i,k)
     END DO
  END DO
  !$acc end parallel
  !
  ! ... hpsi contains h times the basis vectors
  !
  CALL h_psi_ptr( npwx, npw, nvec, psi, hpsi )  ; nhpsi = nvec
  !
  ! ... spsi contains s times the basis vectors
  !
  IF ( uspp ) CALL s_psi_ptr( npwx, npw, nvec, psi, spsi )
  !
  ! ... hr contains the projection of the hamiltonian onto the reduced
  ! ... space vr contains the eigenvectors of hr
  !
  CALL start_clock( 'regterg:init' )           
  !$acc kernels
  hr(:,:) = 0.D0
  sr(:,:) = 0.D0
  vr(:,:) = 0.D0
  !$acc end kernels
  !
  !$acc host_data use_device(psi, hpsi, spsi, hr, sr)
  CALL divide(inter_bgrp_comm,nbase,n_start,n_end)
  my_n = n_end - n_start + 1; !write (*,*) nbase,n_start,n_end
  if (n_start .le. n_end) &
  CALL DGEMM( 'T','N', nbase, my_n, npw2, 2.D0 , psi, npwx2, hpsi(1,n_start), npwx2, 0.D0, hr(1,n_start), nvecx )
  IF ( gstart == 2 ) CALL MYDGER( nbase, my_n, -1.D0, psi, npwx2, hpsi(1,n_start), npwx2, hr(1,n_start), nvecx )
  CALL mp_sum( hr( :, 1:nbase ), inter_bgrp_comm )
  !
  CALL mp_sum( hr( :, 1:nbase ), intra_bgrp_comm )
  !
  IF ( uspp ) THEN
     !
     if (n_start .le. n_end) &
     CALL DGEMM( 'T','N', nbase, my_n, npw2, 2.D0, psi, npwx2, spsi(1,n_start), npwx2, 0.D0, sr(1,n_start), nvecx )
     IF ( gstart == 2 ) CALL MYDGER( nbase, my_n, -1.D0, psi, npwx2, spsi(1,n_start), npwx2, sr(1,n_start), nvecx )
     !
  ELSE
     !
     if (n_start .le. n_end) &
     CALL DGEMM( 'T','N', nbase, my_n, npw2, 2.D0, psi, npwx2, psi(1,n_start), npwx2, 0.D0, sr(1,n_start), nvecx )
     IF ( gstart == 2 ) CALL MYDGER( nbase, my_n, -1.D0, psi, npwx2, psi(1,n_start), npwx2, sr(1,n_start), nvecx )
     !
  END IF
  CALL mp_sum( sr( :, 1:nbase ), inter_bgrp_comm )
  !
  CALL mp_sum( sr( :, 1:nbase ), intra_bgrp_comm )
  !$acc end host_data
  !
  CALL stop_clock( 'regterg:init' )
  !
  IF ( lrot ) THEN
     !
     !$acc parallel loop 
     DO n = 1, nbase
        !
        e(n) = hr(n,n)
        vr(n,n) = 1.D0
        !
     END DO
     !
  ELSE
     !
     ! ... diagonalize the reduced hamiltonian
     !
     !$acc host_data use_device(hr, sr, vr, ew)
     CALL start_clock( 'regterg:diag' )
     IF( my_bgrp_id == root_bgrp_id ) THEN
        CALL diaghg( nbase, nvec, hr, sr, nvecx, ew, vr, me_bgrp, root_bgrp, intra_bgrp_comm )
     END IF
     IF( nbgrp > 1 ) THEN
        CALL mp_bcast( vr, root_bgrp_id, inter_bgrp_comm )
        CALL mp_bcast( ew, root_bgrp_id, inter_bgrp_comm )
     ENDIF
     CALL stop_clock( 'regterg:diag' ) 
     !$acc end host_data
     !
     !$acc parallel loop 
     DO i = 1, nvec
        e(i) = ew(i)
     END DO
     !
  END IF
  !
  ! ... iterate
  !
  iterate: DO kter = 1, maxter
     !
     dav_iter = kter ; !write(*,*) kter, notcnv, conv
     !
     CALL start_clock( 'regterg:update' )
     !
     np = 0
     ! 
     DO n = 1, nvec
        !
        IF ( .NOT. conv(n) ) THEN
           !
           ! ... this root not yet converged ... 
           !
           np = np + 1
           !
           ! ... reorder eigenvectors so that coefficients for unconverged
           ! ... roots come first. This allows to use quick matrix-matrix 
           ! ... multiplications to set a new basis vector (see below)
           !
           IF ( np /= n ) THEN 
             !$acc parallel loop 
             DO i = 1, nvecx
               vr(i,np) = vr(i,n)
             END DO 
           END IF
           !
           ! ... for use in g_psi_ptr
           !
           !$acc kernels 
           ew(nbase+np) = e(n)
           !$acc end kernels
           !
        END IF
        !
     END DO
     !
     nb1 = nbase + 1
     !
     ! ... expand the basis set with new basis vectors ( H - e*S )|psi> ...
     !
     CALL divide(inter_bgrp_comm,nbase,n_start,n_end)
     my_n = n_end - n_start + 1; !write (*,*) nbase,n_start,n_end
     !$acc parallel loop collapse(2)
     DO i=1, notcnv
        DO k=1,npwx
           psi(k,nbase+i)=ZERO
        END DO
     END DO
     !$acc host_data use_device(psi, spsi, vr)
     IF ( uspp ) THEN
        !
        if (n_start .le. n_end) &
        CALL DGEMM( 'N','N', npw2, notcnv, my_n, 1.D0, spsi(1,n_start), npwx2, vr(n_start,1), nvecx, 0.D0, psi(1,nb1), npwx2 )
        !     
     ELSE
        !
        if (n_start .le. n_end) &
        CALL DGEMM( 'N','N', npw2, notcnv, my_n, 1.D0, psi(1,n_start), npwx2, vr(n_start,1), nvecx, 0.D0, psi(1,nb1), npwx2 )
        !
     END IF
     !$acc end host_data
! NB: must not call mp_sum over inter_bgrp_comm here because it is done later to the full correction
     !
     !$acc parallel loop collapse(2) 
     DO np=1,notcnv
        DO k=1,npwx
          psi(k,nbase+np) = - ew(nbase+np) * psi(k,nbase+np)
        END DO
     END DO
     !
     !$acc host_data use_device(psi, hpsi, vr)
     if (n_start .le. n_end) &
     CALL DGEMM( 'N','N', npw2, notcnv, my_n, 1.D0, hpsi(1,n_start), npwx2, vr(n_start,1), nvecx, 1.D0, psi(1,nb1), npwx2 )
     CALL mp_sum( psi(:,nb1:nbase+notcnv), inter_bgrp_comm )
     !$acc end host_data
     !
     CALL stop_clock( 'regterg:update' )
     !
     ! ... approximate inverse iteration
     !
     CALL g_psi_ptr( npwx, npw, notcnv, 1, psi(1,nb1), ew(nb1) )
     !
     ! ... "normalize" correction vectors psi(:,nb1:nbase+notcnv) in 
     ! ... order to improve numerical stability of subspace diagonalization 
     ! ... (rdiaghg) ew is used as work array :
     !
     ! ...         ew = <psi_i|psi_i>,  i = nbase + 1, nbase + notcnv
     !
     !$acc parallel vector_length(96) 
     !$acc loop gang private(nbn) 
     DO n = 1, notcnv
        !
        nbn = nbase + n
        !
        ew(n) = 2.D0 * MYDDOT_VECTOR_GPU( npw2, psi(1,nbn), psi(1,nbn) )
        IF (gstart == 2) ew(n) = ew(n) - DBLE(psi(1,nbn) * psi(1,nbn)) ! psi(1,nbn) * psi(1,nbn)
        !
     END DO
     !$acc end parallel
     !
     !$acc host_data use_device(ew)
     CALL mp_sum( ew( 1:notcnv ), intra_bgrp_comm )
     !$acc end host_data
     !
     !$acc parallel vector_length(96) 
     !$acc loop gang  
     DO i = 1,notcnv
        psi(1,nbase+i) = psi(1,nbase+i)/SQRT( ew(i) )
        ! ... set Im[ psi(G=0) ] -  needed for numerical stability
        IF (gstart == 2) psi(1,nbase+i) = CMPLX( DBLE(psi(1,nbase+i)), 0.D0 ,kind=DP)
        !$acc loop vector 
        DO k=2,npwx
           psi(k,nbase+i) = psi(k,nbase+i)/SQRT( ew(i) )
        END DO
     END DO
     !$acc end parallel
     !
     ! ... here compute the hpsi and spsi of the new functions
     !
     CALL h_psi_ptr( npwx, npw, notcnv, psi(1,nb1), hpsi(1,nb1) ) ; nhpsi = nhpsi + notcnv
     !
     IF ( uspp ) CALL s_psi_ptr( npwx, npw, notcnv, psi(1,nb1), spsi(1,nb1) )
     !
     ! ... update the reduced hamiltonian
     !
     CALL start_clock( 'regterg:overlap' )
     !
     !$acc parallel loop collapse(2)
     DO i=0,notcnv-1
        DO j=1, nvecx
          hr( j, nb1+i )=0.d0
        END DO
     END DO
     !
     !$acc host_data use_device(psi, hpsi, hr)
     CALL divide(inter_bgrp_comm,nbase+notcnv,n_start,n_end)
     my_n = n_end - n_start + 1; !write (*,*) nbase+notcnv,n_start,n_end
     CALL DGEMM( 'T','N', my_n, notcnv, npw2, 2.D0, psi(1,n_start), npwx2, hpsi(1,nb1), npwx2, 0.D0, hr(n_start,nb1), nvecx )
     IF ( gstart == 2 ) CALL MYDGER( my_n, notcnv, -1.D0, psi(1,n_start), npwx2, hpsi(1,nb1), npwx2, hr(n_start,nb1), nvecx )
     CALL mp_sum( hr( :, nb1:nb1+notcnv-1 ), inter_bgrp_comm )
     !
     CALL mp_sum( hr( :, nb1:nb1+notcnv-1 ), intra_bgrp_comm )
     !$acc end host_data
     !
     !$acc parallel loop collapse(2)
     DO i=0,notcnv-1
        DO j=1, nvecx
          sr( j, nb1+i )=0.d0
        END DO
     END DO
     !
     !$acc host_data use_device(psi, spsi, sr)
     CALL divide(inter_bgrp_comm,nbase+notcnv,n_start,n_end)
     my_n = n_end - n_start + 1; !write (*,*) nbase+notcnv,n_start,n_end
     IF ( uspp ) THEN
        !
        CALL DGEMM( 'T','N', my_n, notcnv, npw2, 2.D0, psi(1,n_start), npwx2, spsi(1,nb1), npwx2, 0.D0, sr(n_start,nb1), nvecx )
        IF ( gstart == 2 ) CALL MYDGER( my_n, notcnv, -1.D0, psi(1,n_start), npwx2, spsi(1,nb1), npwx2, sr(n_start,nb1), nvecx )
        !
     ELSE
        !
        CALL DGEMM( 'T','N', my_n, notcnv, npw2, 2.D0, psi(1,n_start), npwx2, psi(1,nb1), npwx2, 0.D0, sr(n_start,nb1) , nvecx )
        IF ( gstart == 2 ) CALL MYDGER( my_n, notcnv, -1.D0, psi(1,n_start), npwx2, psi(1,nb1), npwx2, sr(n_start,nb1), nvecx )
        !
     END IF
     CALL mp_sum( sr( :, nb1:nb1+notcnv-1 ), inter_bgrp_comm )
     !
     CALL mp_sum( sr( :, nb1:nb1+notcnv-1 ), intra_bgrp_comm  )
     !$acc end host_data
     !
     CALL stop_clock( 'regterg:overlap' )
     !
     nbase = nbase + notcnv
     !
     !$acc parallel 
     !$acc loop gang
     DO n = 1, nbase
        !
        !$acc loop vector 
        DO m = n + 1, nbase
           !
           hr(m,n) = hr(n,m)
           sr(m,n) = sr(n,m)
           !
        END DO
        !
     END DO
     !$acc end parallel 
     !
     ! ... diagonalize the reduced hamiltonian
     !
     CALL start_clock( 'regterg:diag' )
     !$acc host_data use_device(hr, sr, ew, vr)
     IF( my_bgrp_id == root_bgrp_id ) THEN
        CALL diaghg( nbase, nvec, hr, sr, nvecx, ew, vr, me_bgrp, root_bgrp, intra_bgrp_comm )
     END IF
     IF( nbgrp > 1 ) THEN
        CALL mp_bcast( vr, root_bgrp_id, inter_bgrp_comm )
        CALL mp_bcast( ew, root_bgrp_id, inter_bgrp_comm )
     ENDIF
     !$acc end host_data
     CALL stop_clock( 'regterg:diag' )
     !
     ! ... test for convergence
     !
     !$acc parallel loop copy(conv(1:nvec)) copyin(btype(1:nvec))
     DO i = 1, nvec
       IF(btype(i) == 1) THEN
         conv(i) = ( ( ABS( ew(i) - e(i) ) < ethr ) )
       ELSE
         conv(i) = ( ( ABS( ew(i) - e(i) ) < empty_ethr ) )
       END IF 
     END DO 
     !
     ! ... next line useful for band parallelization of exact exchange
     IF ( nbgrp > 1 ) CALL mp_bcast(conv,root_bgrp_id,inter_bgrp_comm)
     !
     notcnv = COUNT( .NOT. conv(:) )
     !
     !$acc parallel loop 
     DO i=1,nvec
       e(i) = ew(i)
     END DO
     !
     ! ... if overall convergence has been achieved, or the dimension of
     ! ... the reduced basis set is becoming too large, or in any case if
     ! ... we are at the last iteration refresh the basis set. i.e. replace
     ! ... the first nvec elements with the current estimate of the
     ! ... eigenvectors;  set the basis dimension to nvec.
     !
     IF ( notcnv == 0 .OR. &
          nbase+notcnv > nvecx .OR. dav_iter == maxter ) THEN
        !
        CALL start_clock( 'regterg:last' )
        !
        !$acc parallel loop collapse(2)  
        DO k=1,nvec
           DO i=1,npwx
              evc(i,k) = ZERO
           END DO
        END DO
        !
        CALL divide(inter_bgrp_comm,nbase,n_start,n_end)
        my_n = n_end - n_start + 1; !write (*,*) nbase,n_start,n_end
        !$acc host_data use_device(evc, psi, vr)
        CALL DGEMM( 'N','N', npw2, nvec, my_n, 1.D0, psi(1,n_start), npwx2, vr(n_start,1), nvecx, 0.D0, evc, npwx2 )
        CALL mp_sum( evc, inter_bgrp_comm )
        !$acc end host_data
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
        !$acc parallel loop collapse(2) 
        DO i=1,nvec
           DO k=1,npwx
              psi(k,i) = evc(k,i)
           END DO
        END DO
        !
        IF ( uspp ) THEN
           !
           !$acc parallel loop collapse(2)
           DO i = 1, npwx
             DO j = nvec+1, nvec+nvec
               psi(i,j) = ZERO
             END DO 
           END DO 
           !
           !$acc host_data use_device(psi, spsi, vr)
           CALL DGEMM( 'N','N', npw2, nvec, my_n, 1.D0, spsi(1,n_start), npwx2, vr(n_start,1), nvecx, 0.D0, psi(1,nvec+1), npwx2 )
           CALL mp_sum( psi(:,nvec+1:nvec+nvec), inter_bgrp_comm )
           !$acc end host_data
           !
           !$acc parallel loop collapse(2) 
           DO i=1,nvec
              DO k=1,npwx
                 spsi(k,i) = psi(k,i+nvec)
              END DO
           END DO
           !
        END IF
        !
        !$acc kernels
        psi(:,nvec+1:nvec+nvec) = ZERO
        !$acc end kernels
        !$acc host_data use_device(psi, hpsi, vr)
        CALL DGEMM( 'N','N', npw2, nvec, my_n, 1.D0, hpsi(1,n_start), npwx2, vr(n_start,1), nvecx, 0.D0, psi(1,nvec+1), npwx2 )
        CALL mp_sum( psi(:,nvec+1:nvec+nvec), inter_bgrp_comm )
        !$acc end host_data
        !
        !$acc parallel loop collapse(2) 
        DO i=1,nvec
           DO k=1, npwx
              hpsi(k,i) = psi(k,i+nvec)
           END DO
        END DO
        !
        ! ... refresh the reduced hamiltonian
        !
        nbase = nvec
        !
        !$acc parallel loop collapse(2) 
        DO i = 1, nvecx
          DO j = 1, nbase
            hr(i,j) = 0.D0
            sr(i,j) = 0.D0
            vr(i,j) = 0.D0
          END DO 
        END DO 
        !
        !$acc parallel loop 
        DO j = 1, nbase
          hr(j,j) = e(j)
          sr(j,j) = 1.D0
          vr(j,j) = 1.D0
        END DO 
        !
        CALL stop_clock( 'regterg:last' )
        !
     END IF
     !
  END DO iterate
  !
  DEALLOCATE( conv )
  DEALLOCATE( vr )
  DEALLOCATE( hr )
  DEALLOCATE( sr )
  !$acc exit data delete(ew)
  DEALLOCATE( ew )
  !
  IF ( uspp ) THEN
     !$acc exit data delete(spsi)
     DEALLOCATE( spsi )
  END IF
  !
  !$acc exit data delete(psi, hpsi)
  DEALLOCATE( hpsi )
  DEALLOCATE( psi )  
  !
  !$acc update host(e)
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
END SUBROUTINE regterg
