!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
#if defined(__CUDA)
!----------------------------------------------------------------------------
SUBROUTINE cegterg_gpu( h_psi_gpu, s_psi_gpu, uspp, g_psi_gpu, &
                    npw, npwx, nvec, nvecx, npol, evc_d, ethr, &
                    e_d, btype, notcnv, lrot, dav_iter )
  !----------------------------------------------------------------------------
  ! PB : 17/8/18, restored original algorithm since aligned version
  !      in commit 354a86b is much slower, probably due to MPI and/or
  !      load unbalance. More careful analysis needed.
  !
  ! ... iterative solution of the eigenvalue problem:
  !
  ! ... ( H - e S ) * evc = 0
  !
  ! ... where H is an hermitean operator, e is a real scalar,
  ! ... S is an overlap matrix, evc is a complex vector
  !
  use cudafor
  use cublas
  USE LAXlib,        ONLY : diaghg
  USE david_param,   ONLY : DP
  USE mp_bands_util, ONLY : intra_bgrp_comm, inter_bgrp_comm, root_bgrp_id,&
          nbgrp, my_bgrp_id
  USE mp,            ONLY : mp_sum, mp_bcast
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: npw, npwx, nvec, nvecx, npol
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix evc, as declared in the calling pgm unit
    ! integer number of searched low-lying roots
    ! maximum dimension of the reduced basis set :
    !    (the basis set is refreshed when its dimension would exceed nvecx)
    ! umber of spin polarizations
  COMPLEX(DP), DEVICE, INTENT(INOUT) :: evc_d(npwx,npol,nvec)
    !  evc contains the  refined estimates of the eigenvectors  
  REAL(DP), INTENT(IN) :: ethr
    ! energy threshold for convergence :
    !   root improvement is stopped, when two consecutive estimates of the root
    !   differ by less than ethr.
  LOGICAL, INTENT(IN) :: uspp
    ! if .FALSE. : do not calculate S|psi>
  INTEGER, INTENT(IN) :: btype(nvec)
    ! band type ( 1 = occupied, 0 = empty )
  LOGICAL, INTENT(IN) :: lrot
    ! .TRUE. if the wfc have already been rotated
  REAL(DP), DEVICE, INTENT(OUT) :: e_d(nvec)
    ! contains the estimated roots.
  INTEGER, INTENT(OUT) :: dav_iter, notcnv
    ! integer number of iterations performed
    ! number of unconverged roots
  !
  ! ... LOCAL variables
  !@njs: hc, sc, vc, psi, hpsi, spsi
  !
  INTEGER, PARAMETER :: maxter = 20
    ! maximum number of iterations
  !
  INTEGER :: kter, nbase, np, kdim, kdmx, n, m, nb1, nbn
    ! counter on iterations
    ! dimension of the reduced basis
    ! counter on the reduced basis vectors
    ! adapted npw and npwx
    ! do-loop counters
  INTEGER :: n_start, n_end, my_n
  INTEGER :: ierr
  COMPLEX(DP), DEVICE, ALLOCATABLE :: hc_d(:,:), sc_d(:,:), vc_d(:,:)
    ! Hamiltonian on the reduced basis
    ! S matrix on the reduced basis
    ! the eigenvectors of the Hamiltonian
  REAL(DP), DEVICE, ALLOCATABLE :: ew_d(:)
    ! eigenvalues of the reduced hamiltonian
  COMPLEX(DP), DEVICE, ALLOCATABLE :: psi_d(:,:,:), hpsi_d(:,:,:), spsi_d(:,:,:)
    ! work space, contains psi
    ! the product of H and psi
    ! the product of S and psi
  LOGICAL, ALLOCATABLE  :: conv(:)
    ! true if the root is converged
  REAL(DP) :: empty_ethr 
    ! threshold for empty bands
  REAL(DP), ALLOCATABLE :: ew_host(:)
  REAL(DP), ALLOCATABLE :: e_host(:)
    ! auxiliary variables for performing dot product
  INTEGER :: i,j,k
  !
  !
  EXTERNAL  h_psi_gpu,    s_psi_gpu,    g_psi_gpu
    ! h_psi(npwx,npw,nvec,psi,hpsi)
    !     calculates H|psi>
    ! s_psi(npwx,npw,nvec,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx,npol,nvec)
    ! g_psi(npwx,npw,notcnv,psi,e)
    !    calculates (diag(h)-e)^-1 * psi, diagonal approx. to (h-e)^-1*psi
    !    the first nvec columns contain the trial eigenvectors
  !
  CALL start_clock( 'cegterg' ); !write(*,*) 'start cegterg' ; FLUSH(6)
  !
  IF ( nvec > nvecx / 2 ) CALL errore( 'cegterg', 'nvecx is too small', 1 )
  !
  ! ... threshold for empty bands
  !
  empty_ethr = MAX( ( ethr * 5.D0 ), 1.D-5 )
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
  ALLOCATE(  psi_d( npwx, npol, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate psi ', ABS(ierr) )
  ALLOCATE( hpsi_d( npwx, npol, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate hpsi ', ABS(ierr) )
  !
  IF ( uspp ) THEN
     ALLOCATE( spsi_d( npwx, npol, nvecx ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' cegterg ',' cannot allocate spsi ', ABS(ierr) )
  END IF
  !
  ALLOCATE( sc_d( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate sc_d ', ABS(ierr) )
  ALLOCATE( hc_d( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate hc_d ', ABS(ierr) )
  ALLOCATE( vc_d( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate vc_d ', ABS(ierr) )
  ALLOCATE( ew_d( nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate ew_d ', ABS(ierr) )
  ALLOCATE( ew_host( nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate ew_host ', ABS(ierr) )
  ALLOCATE( e_host( nvec ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate e_host ', ABS(ierr) )
  ALLOCATE( conv( nvec ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate conv ', ABS(ierr) )
  !
  notcnv = nvec
  nbase  = nvec
  conv   = .FALSE.
  !
  IF ( uspp ) spsi_d = ZERO
  !
  hpsi_d = ZERO
  psi_d  = ZERO
  !$cuf kernel do(3) <<<*,*>>>
  DO k=1,nvec
     DO j=1,npol
        DO i=1,npwx
           psi_d(i,j,k) = evc_d(i,j,k)
      	END DO
     END DO
  END DO
  !
  ! ... hpsi contains h times the basis vectors
  !
  CALL h_psi_gpu( npwx, npw, nvec, psi_d, hpsi_d )
  !
  ! ... spsi contains s times the basis vectors
  !
  IF ( uspp ) CALL s_psi_gpu( npwx, npw, nvec, psi_d, spsi_d )
  !
  ! ... hc contains the projection of the hamiltonian onto the reduced 
  ! ... space vc contains the eigenvectors of hc
  !
  CALL start_clock( 'cegterg:init' )
  hc_d(:,:) = ZERO
  sc_d(:,:) = ZERO
  vc_d(:,:) = ZERO
  !
  CALL divide(inter_bgrp_comm,nbase,n_start,n_end)
  my_n = n_end - n_start + 1; !write (*,*) nbase,n_start,n_end
  if (n_start .le. n_end) &
  CALL ZGEMM( 'C','N', nbase, my_n, kdim, ONE, psi_d, kdmx, hpsi_d(1,1,n_start), kdmx, ZERO, hc_d(1,n_start), nvecx )
  CALL mp_sum( hc_d( :, 1:nbase ), inter_bgrp_comm )
  !
  CALL mp_sum( hc_d( :, 1:nbase ), intra_bgrp_comm )
  !
  IF ( uspp ) THEN
     !
     if (n_start .le. n_end) &
     CALL ZGEMM( 'C','N', nbase, my_n, kdim, ONE, psi_d, kdmx, spsi_d(1,1,n_start), kdmx, &
                 ZERO, sc_d(1,n_start), nvecx )
     !     
  ELSE
     !
     if (n_start .le. n_end) &
     CALL ZGEMM( 'C','N', nbase, my_n, kdim, ONE, psi_d, kdmx, psi_d(1,1,n_start), kdmx, &
                 ZERO, sc_d(1,n_start), nvecx )
     !
  END IF
  CALL mp_sum( sc_d( :, 1:nbase ), inter_bgrp_comm )
  !
  CALL mp_sum( sc_d( :, 1:nbase ), intra_bgrp_comm )

  CALL stop_clock( 'cegterg:init' )
  !
  IF ( lrot ) THEN
     !
!$cuf kernel do(1) <<<*,*>>>
     DO n = 1, nbase
        !
        e_d(n) = REAL( hc_d(n,n) )
        !
        vc_d(n,n) = ONE
        !
     END DO
     !
  ELSE
     !
     ! ... diagonalize the reduced hamiltonian
     !
     CALL start_clock( 'cegterg:diag' )
     IF( my_bgrp_id == root_bgrp_id ) THEN
        CALL diaghg( nbase, nvec, hc_d, sc_d, nvecx, ew_d, vc_d )
     END IF
     IF( nbgrp > 1 ) THEN
        CALL mp_bcast( vc_d, root_bgrp_id, inter_bgrp_comm )
        CALL mp_bcast( ew_d, root_bgrp_id, inter_bgrp_comm )
     ENDIF
     CALL stop_clock( 'cegterg:diag' )
     !
!$cuf kernel do(1) <<<*,*>>>
     DO i = 1, nvec
        e_d(i) = ew_d(i)
     END DO
     !
  END IF
  !
  ! ... iterate
  !
  iterate: DO kter = 1, maxter
     !
     dav_iter = kter
     !
     CALL start_clock( 'cegterg:update' )
     !
     !  ======== FROM HERE =====
     !np = 0
     !
     !DO n = 1, nvec
     !   !
     !   IF ( .NOT. conv(n) ) THEN
     !      !
     !      ! ... this root not yet converged ... 
     !      !
     !      np = np + 1
     !      !
     !      ! ... reorder eigenvectors so that coefficients for unconverged
     !      ! ... roots come first. This allows to use quick matrix-matrix 
     !      ! ... multiplications to set a new basis vector (see below)
     !      !
     !      IF ( np /= n ) vc_d(:,np) = vc_d(:,n)
     !      !
     !      ! ... for use in g_psi
     !      !
     !      ew_d(nbase+np) = e_d(n)
     !      !
     !   END IF
     !   !
     !END DO
     ! ========= TO HERE, REPLACED BY =======

     CALL reorder_evals_cevecs(nbase, nvec, nvecx, conv, e_d, ew_d, vc_d)
     !
     nb1 = nbase + 1
     !
     ! ... expand the basis set with new basis vectors ( H - e*S )|psi> ...
     !
     CALL divide(inter_bgrp_comm,nbase,n_start,n_end)
     my_n = n_end - n_start + 1; !write (*,*) nbase,n_start,n_end
!$cuf kernel do(3) <<<*,*>>>
     DO i=1, notcnv
        DO j=1,npol
           DO k=npw,npwx         ! pn;y cleanup what needs to be cleaned up
              psi_d(k,j,nbase+i)=ZERO
           END DO
        END DO
     END DO
     IF ( uspp ) THEN
        !
        if (n_start .le. n_end) &
        CALL ZGEMM( 'N','N', kdim, notcnv, my_n, ONE, spsi_d(1,1,n_start), kdmx, vc_d(n_start,1), nvecx, &
                    ZERO, psi_d(1,1,nb1), kdmx )
        !     
     ELSE
        !
        if (n_start .le. n_end) &
        CALL ZGEMM( 'N','N', kdim, notcnv, my_n, ONE, psi_d(1,1,n_start), kdmx, vc_d(n_start,1), nvecx, &
                    ZERO, psi_d(1,1,nb1), kdmx )
        !
     END IF
! NB: must not call mp_sum over inter_bgrp_comm here because it is done later to the full correction
     !

!$cuf kernel do(3) <<<*,*>>>
     DO np=1,notcnv
        DO j=1,npol
           DO k=1,npwx
             psi_d(k,j,nbase+np) = - ew_d(nbase+np)*psi_d(k,j,nbase+np)
           END DO
        END DO
     END DO
     !
     if (n_start .le. n_end) &
     CALL ZGEMM( 'N','N', kdim, notcnv, my_n, ONE, hpsi_d(1,1,n_start), kdmx, vc_d(n_start,1), nvecx, &
                 ONE, psi_d(1,1,nb1), kdmx )
     CALL mp_sum( psi_d(:,:,nb1:nbase+notcnv), inter_bgrp_comm )
     !
     CALL stop_clock( 'cegterg:update' )
     !
     ! ... approximate inverse iteration
     !
     CALL g_psi_gpu( npwx, npw, notcnv, npol, psi_d(1,1,nb1), ew_d(nb1) )
     !
     ! ... "normalize" correction vectors psi(:,nb1:nbase+notcnv) in
     ! ... order to improve numerical stability of subspace diagonalization
     ! ... (cdiaghg) ew is used as work array :
     !
     ! ...         ew = <psi_i|psi_i>,  i = nbase + 1, nbase + notcnv
     !
     !!! == OPTIMIZE HERE ==
     DO n = 1, notcnv
        !
        nbn = nbase + n
        !
        IF ( npol == 1 ) THEN
           !
           ew_host(n) = DBLE(cublasZdotc( npw, psi_d(1,1,nbn), 1, psi_d(1,1,nbn), 1 )) ! DBLE(aux1)
           !
        ELSE
           !
           ew_host(n) = DBLE(cublasZdotc( npw, psi_d(1,1,nbn), 1, psi_d(1,1,nbn), 1 ) + &
                        cublasZdotc( npw, psi_d(1,2,nbn), 1, psi_d(1,2,nbn), 1 ))
           !
        END IF
        !
     END DO
     !
     CALL mp_sum( ew_host( 1:notcnv ), intra_bgrp_comm )
     ew_d(1:notcnv) = ew_host(1:notcnv)
     !
!$cuf kernel do(3) <<<*,*>>>
     DO i = 1,notcnv
        DO j=1,npol
           DO k=1,npwx
             psi_d(k,j,nbase+i) = psi_d(k,j,nbase+i)/SQRT( ew_d(i) )
           END DO
        END DO
     END DO
     !
     ! ... here compute the hpsi and spsi of the new functions
     !
     CALL h_psi_gpu( npwx, npw, notcnv, psi_d(:,:,nb1), hpsi_d(:,:,nb1) )
     !
     IF ( uspp ) CALL s_psi_gpu( npwx, npw, notcnv, psi_d(1,1,nb1), spsi_d(1,1,nb1) )
     !
     ! ... update the reduced hamiltonian
     !
     CALL start_clock( 'cegterg:overlap' )
     !
!$cuf kernel do(2) <<<*,*>>>
     DO i=0,notcnv-1
        DO j=1, nvecx
          hc_d( j, nb1+i )=ZERO
          sc_d( j, nb1+i )=ZERO
        END DO
     END DO
     CALL divide(inter_bgrp_comm,nbase+notcnv,n_start,n_end)
     my_n = n_end - n_start + 1; !write (*,*) nbase+notcnv,n_start,n_end
     CALL ZGEMM( 'C','N', my_n, notcnv, kdim, ONE, psi_d(1,1,n_start), kdmx, hpsi_d(1,1,nb1), kdmx, &
                 ZERO, hc_d(n_start,nb1), nvecx )
     CALL mp_sum( hc_d( :, nb1:nb1+notcnv-1 ), inter_bgrp_comm )
     CALL mp_sum( hc_d( :, nb1:nb1+notcnv-1 ), intra_bgrp_comm )
     !
     IF ( uspp ) THEN
        !
        CALL ZGEMM( 'C','N', my_n, notcnv, kdim, ONE, psi_d(1,1,n_start), kdmx, spsi_d(1,1,nb1), kdmx, &
                    ZERO, sc_d(n_start,nb1), nvecx )
        !     
     ELSE
        !
        CALL ZGEMM( 'C','N', my_n, notcnv, kdim, ONE, psi_d(1,1,n_start), kdmx, psi_d(1,1,nb1), kdmx, &
                     ZERO, sc_d(n_start,nb1), nvecx )
        !
     END IF
     CALL mp_sum( sc_d( :, nb1:nb1+notcnv-1 ), inter_bgrp_comm )
     CALL mp_sum( sc_d( :, nb1:nb1+notcnv-1 ), intra_bgrp_comm )
     !
     CALL stop_clock( 'cegterg:overlap' )
     !
     nbase = nbase + notcnv
     !
!$cuf kernel do(1) <<<*,*>>>
     DO n = 1, nbase
        !
        ! ... the diagonal of hc and sc must be strictly real 
        !
        hc_d(n,n) = CMPLX( REAL( hc_d(n,n) ), 0.D0 ,kind=DP)
        sc_d(n,n) = CMPLX( REAL( sc_d(n,n) ), 0.D0 ,kind=DP)
        !
        DO m = n + 1, nbase
           !
           hc_d(m,n) = CONJG( hc_d(n,m) )
           sc_d(m,n) = CONJG( sc_d(n,m) )
           !
        END DO
        !
     END DO
     !
     ! ... diagonalize the reduced hamiltonian
     !
     CALL start_clock( 'cegterg:diag' )
     IF( my_bgrp_id == root_bgrp_id ) THEN
        CALL diaghg( nbase, nvec, hc_d, sc_d, nvecx, ew_d, vc_d )
     END IF
     IF( nbgrp > 1 ) THEN
        CALL mp_bcast( vc_d, root_bgrp_id, inter_bgrp_comm )
        CALL mp_bcast( ew_d, root_bgrp_id, inter_bgrp_comm )
     ENDIF
     CALL stop_clock( 'cegterg:diag' )
     !
     ! ... test for convergence (on the CPU)
     !
     ew_host(1:nvec) = ew_d(1:nvec)
     e_host(1:nvec) = e_d(1:nvec)
     WHERE( btype(1:nvec) == 1 )
        !
        conv(1:nvec) = ( ( ABS( ew_host(1:nvec) - e_host(1:nvec) ) < ethr ) )
        !
     ELSEWHERE
        !
        conv(1:nvec) = ( ( ABS( ew_host(1:nvec) - e_host(1:nvec) ) < empty_ethr ) )
        !
     END WHERE
     ! ... next line useful for band parallelization of exact exchange
     IF ( nbgrp > 1 ) CALL mp_bcast(conv,root_bgrp_id,inter_bgrp_comm)
     !
     notcnv = COUNT( .NOT. conv(:) )
     !
!$cuf kernel do(1) <<<*,*>>>
     DO i=1,nvec
       e_d(i) = ew_d(i)
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
        CALL start_clock( 'cegterg:last' )
        !
        ! Only reset what needs to be reset
        !$cuf kernel do(3) <<<*,*>>>
        DO k=1,nvec
           DO j=1,npol
              DO i=npw,npwx
                 evc_d(i,j,k) = ZERO
              END DO
           END DO
        END DO
        CALL divide(inter_bgrp_comm,nbase,n_start,n_end)
        my_n = n_end - n_start + 1; !write (*,*) nbase,n_start,n_end
        CALL ZGEMM( 'N','N', kdim, nvec, my_n, ONE, psi_d(1,1,n_start), kdmx, vc_d(n_start,1), nvecx, &
                    ZERO, evc_d, kdmx )
        CALL mp_sum( evc_d, inter_bgrp_comm )
        !
        IF ( notcnv == 0 ) THEN
           !
           ! ... all roots converged: return
           !
           CALL stop_clock( 'cegterg:last' )
           !
           EXIT iterate
           !
        ELSE IF ( dav_iter == maxter ) THEN
           !
           ! ... last iteration, some roots not converged: return
           !
           !!!WRITE( stdout, '(5X,"WARNING: ",I5, &
           !!!     &   " eigenvalues not converged")' ) notcnv
           !
           CALL stop_clock( 'cegterg:last' )
           !
           EXIT iterate
           !
        END IF
        !
        ! ... refresh psi, H*psi and S*psi
        !
        !$cuf kernel do(3) <<<*,*>>>
        DO i=1,nvec
           DO j=1,npol
              DO k=1,npwx
                 psi_d(k,j,i) = evc_d(k,j,i)
              END DO
           END DO
        END DO
        !
        IF ( uspp ) THEN
           !
           !psi_d(:,:,nvec+1:nvec+nvec) = ZERO (only clean what needs to be cleaned)
           !$cuf kernel do(3) <<<*,*>>>
           DO i=1, nvec
              DO j=1, npol
                 DO k=npw, npwx
                    psi_d(k,j,nvec+i) = ZERO
                 END DO
              END DO
           END DO
           CALL ZGEMM( 'N','N', kdim, nvec, my_n, ONE, spsi_d(1,1,n_start), kdmx, vc_d(n_start,1), nvecx, &
                       ZERO, psi_d(1,1,nvec+1), kdmx)
           CALL mp_sum( psi_d(:,:,nvec+1:nvec+nvec), inter_bgrp_comm )
           !
           !spsi_d(:,:,1:nvec) = psi_d(:,:,nvec+1:nvec+nvec)
           !$cuf kernel do(3) <<<*,*>>>
           DO i=1,nvec
              DO j=1, npol
                 DO k=1, npwx
                    spsi_d(k,j,i) = psi_d(k,j,i+nvec)
                 END DO
              END DO
           END DO
           !
        END IF
        !
        !psi_d(:,:,nvec+1:nvec+nvec) = ZERO (only clean what needs to be cleaned)
        !$cuf kernel do(3) <<<*,*>>>
        DO i=1, nvec
           DO j=1, npol
              DO k=npw, npwx
                 psi_d(k,j,nvec+i) = ZERO
              END DO
           END DO
        END DO
        CALL ZGEMM( 'N','N', kdim, nvec, my_n, ONE, hpsi_d(1,1,n_start), kdmx, vc_d(n_start,1), nvecx, &
                    ZERO, psi_d(1,1,nvec+1), kdmx )
        CALL mp_sum( psi_d(:,:,nvec+1:nvec+nvec), inter_bgrp_comm )
        !
        !hpsi_d(:,:,1:nvec) = psi_d(:,:,nvec+1:nvec+nvec)
        !$cuf kernel do(3) <<<*,*>>>
        DO i=1,nvec
           DO j=1, npol
              DO k=1, npwx
                 hpsi_d(k,j,i) = psi_d(k,j,i+nvec)
              END DO
           END DO
        END DO
        !
        ! ... refresh the reduced hamiltonian 
        !
        nbase = nvec
        !
        hc_d(:,1:nbase) = ZERO
        sc_d(:,1:nbase) = ZERO
        vc_d(:,1:nbase) = ZERO
        !
        !$cuf kernel do(1) <<<*,*>>>
        DO n = 1, nbase
           !
!           hc(n,n) = REAL( e(n) )
           hc_d(n,n) = CMPLX( e_d(n), 0.0_DP ,kind=DP)
           !
           sc_d(n,n) = ONE
           vc_d(n,n) = ONE
           !
        END DO
        !
        CALL stop_clock( 'cegterg:last' )
        !
     END IF
     !
  END DO iterate
  !
  DEALLOCATE( conv )
  DEALLOCATE( e_host, ew_host, ew_d )
  DEALLOCATE( vc_d )
  DEALLOCATE( hc_d )
  DEALLOCATE( sc_d )
  !
  IF ( uspp ) DEALLOCATE( spsi_d )
  !
  DEALLOCATE( hpsi_d )
  DEALLOCATE( psi_d )
  !
  CALL stop_clock( 'cegterg' ); !write(*,*) 'stop cegterg' ; FLUSH(6)
  !call print_clock( 'cegterg' )
  !call print_clock( 'cegterg:init' )
  !call print_clock( 'cegterg:diag' )
  !call print_clock( 'cegterg:update' )
  !call print_clock( 'cegterg:overlap' )
  !call print_clock( 'cegterg:last' )
  !
  RETURN
  !
END SUBROUTINE cegterg_gpu

SUBROUTINE reorder_evals_cevecs(nbase, nvec, nvecx, conv, e_d, ew_d, v_d)
   USE david_param,   ONLY : DP
   USE gbuffers,  ONLY : buffer => dev_buf
   implicit none
   INTEGER, INTENT(IN) :: nbase, nvec, nvecx
   LOGICAL, INTENT(IN) :: conv(nvec)
   REAL(DP), DEVICE :: e_d(nvecx), ew_d(nvecx)
   COMPLEX(DP), DEVICE :: v_d(nvecx,nvecx)
   !
   INTEGER :: j, k, n, np, info
   INTEGER, ALLOCATABLE :: conv_idx(:)
   INTEGER, DEVICE, POINTER :: conv_idx_d(:)
   COMPLEX(DP), DEVICE, POINTER :: vtmp_d(:,:)
   !
   np = 0
   ALLOCATE(conv_idx(nvec))
   DO n = 1, nvec
      conv_idx(n) = -1
      IF ( .NOT. conv(n) ) THEN
         np = np + 1
         conv_idx(n) = np
      END IF
   END DO

   CALL buffer%lock_buffer(conv_idx_d, nvec, info)
   CALL buffer%lock_buffer(vtmp_d, (/nvecx, nvecx/), info)

   conv_idx_d(1:nvec) = conv_idx(1:nvec)

!$cuf kernel do(2) <<<*,*>>>
   DO j=1,nvec
      DO k=1,nvecx
         vtmp_d(k,j) = v_d(k,j)
      END DO
   END DO

!$cuf kernel do(2) <<<*,*>>>
   DO j=1,nvec
      DO k=1,nvecx
         IF(conv_idx_d(j) /= -1) THEN
           v_d(k,conv_idx_d(j)) = vtmp_d(k,j)
           IF(k==1) ew_d(nbase+conv_idx_d(j)) = e_d(j)
         END IF
      END DO
   END DO
   !
   CALL buffer%release_buffer(conv_idx_d, info)
   CALL buffer%release_buffer(vtmp_d, info)
   !
   DEALLOCATE(conv_idx)
END SUBROUTINE reorder_evals_cevecs

!
!  Wrapper for subroutine with distributed matrixes (written by Carlo Cavazzoni)
!
!----------------------------------------------------------------------------
SUBROUTINE pcegterg_gpu(h_psi_gpu, s_psi_gpu, uspp, g_psi_gpu, &  
                    npw, npwx, nvec, nvecx, npol, evc_d, ethr, &
                    e_d, btype, notcnv, lrot, dav_iter )
  !----------------------------------------------------------------------------
  !
  ! ... iterative solution of the eigenvalue problem:
  !
  ! ... ( H - e S ) * evc = 0
  !
  ! ... where H is an hermitean operator, e is a real scalar,
  ! ... S is an uspp matrix, evc is a complex vector
  !
  USE david_param,      ONLY : DP, stdout
  USE mp_bands_util,    ONLY : intra_bgrp_comm, inter_bgrp_comm, root_bgrp_id, nbgrp, my_bgrp_id
  USE mp_diag,          ONLY : ortho_comm, np_ortho, me_ortho, ortho_comm_id, leg_ortho, &
                                 ortho_parent_comm, ortho_cntx, do_distr_diag_inside_bgrp
  USE descriptors,      ONLY : la_descriptor, descla_init , descla_local_dims
  USE parallel_toolkit, ONLY : zsqmred, zsqmher, zsqmdst
  USE mp,               ONLY : mp_bcast, mp_root_sum, mp_sum, mp_barrier
  USE gbuffers,         ONLY : buffer => dev_buf
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: npw, npwx, nvec, nvecx, npol
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix evc, as declared in the calling pgm unit
    ! integer number of searched low-lying roots
    ! maximum dimension of the reduced basis set
    !    (the basis set is refreshed when its dimension would exceed nvecx)
    ! number of spin polarizations
  INTEGER, PARAMETER :: blocksize = 256
  INTEGER :: numblock
    ! chunking parameters
  COMPLEX(DP), DEVICE, INTENT(INOUT) :: evc_d(npwx,npol,nvec)
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
  REAL(DP), DEVICE, INTENT(OUT) :: e_d(nvec)
    ! contains the estimated roots.
  INTEGER, INTENT(OUT) :: dav_iter, notcnv
    ! integer  number of iterations performed
    ! number of unconverged roots
  !
  ! ... LOCAL variables
  !
  COMPLEX(DP), ALLOCATABLE :: evc(:,:,:)
  REAL(DP), ALLOCATABLE :: e(:)
  
  INTEGER, PARAMETER :: maxter = 20
    ! maximum number of iterations
  !
  INTEGER :: kter, nbase, np, kdim, kdmx, n, m, ipol, nb1, nbn
    ! counter on iterations
    ! dimension of the reduced basis
    ! counter on the reduced basis vectors
    ! do-loop counters
  INTEGER :: i, j, k, ierr
  REAL(DP), ALLOCATABLE :: ew(:)
  REAL(DP), DEVICE, POINTER :: ew_d(:)
  COMPLEX(DP), ALLOCATABLE :: hl(:,:), sl(:,:), vl(:,:)
    ! Hamiltonian on the reduced basis
    ! S matrix on the reduced basis
    ! eigenvectors of the Hamiltonian
    ! eigenvalues of the reduced hamiltonian
  COMPLEX(DP), ALLOCATABLE :: psi(:,:,:), hpsi(:,:,:), spsi(:,:,:)
  COMPLEX(DP), DEVICE, POINTER :: psi_d(:,:,:), hpsi_d(:,:,:), spsi_d(:,:,:)
    ! work space, contains psi
    ! the product of H and psi
    ! the product of S and psi
  LOGICAL, ALLOCATABLE :: conv(:)
    ! true if the root is converged
  REAL(DP) :: empty_ethr 
    ! threshold for empty bands
  TYPE(la_descriptor) :: desc, desc_old
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
  REAL(DP), EXTERNAL :: ddot
  !
  EXTERNAL  h_psi_gpu, s_psi_gpu, g_psi_gpu
    ! h_psi(npwx,npw,nvec,psi,hpsi)
    !     calculates H|psi> 
    ! s_psi(npwx,npw,nvec,psi,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx,nvec)
    ! g_psi(npwx,npw,notcnv,psi,e)
    !    calculates (diag(h)-e)^-1 * psi, diagonal approx. to (h-e)^-1*psi
    !    the first nvec columns contain the trial eigenvectors
  !
  !
  CALL start_clock( 'cegterg' )
  !
  IF ( nvec > nvecx / 2 ) CALL errore( 'pcegterg', 'nvecx is too small', 1 )
  !
  ! ... threshold for empty bands
  !
  empty_ethr = MAX( ( ethr * 5.D0 ), 1.D-5 )
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
  ! compute the number of chuncks
  numblock  = (npw+blocksize-1)/blocksize

  ALLOCATE(  evc( npwx, npol, nvec ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate evc (host) ', ABS(ierr) )
  !
  ALLOCATE(  e( nvec ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate e (host) ', ABS(ierr) )
  !
  ALLOCATE(  psi( npwx, npol, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate psi ', ABS(ierr) )
  !
  ALLOCATE( hpsi( npwx, npol, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate hpsi ', ABS(ierr) )
  !
  IF ( uspp ) THEN
     ALLOCATE( spsi( npwx, npol, nvecx ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' pcegterg ',' cannot allocate spsi ', ABS(ierr) )
  END IF
  !
  ! ... Initialize the matrix descriptor
  !
  ALLOCATE( ic_notcnv( np_ortho(2) ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate ic_notcnv ', ABS(ierr) )
  !
  ALLOCATE( notcnv_ip( np_ortho(2) ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate notcnv_ip ', ABS(ierr) )
  !
  ALLOCATE( irc_ip( np_ortho(1) ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate irc_ip ', ABS(ierr) )
  !
  ALLOCATE( nrc_ip( np_ortho(1) ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate nrc_ip ', ABS(ierr) )
  !
  ALLOCATE( rank_ip( np_ortho(1), np_ortho(2) ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate rank_ip ', ABS(ierr) )
  !
  CALL desc_init( nvec, desc, irc_ip, nrc_ip )
  !
  IF( la_proc ) THEN
     !
     ! only procs involved in the diagonalization need to allocate local 
     ! matrix block.
     !
     ALLOCATE( vl( nx , nx ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' pcegterg ',' cannot allocate vl ', ABS(ierr) )
     !
     ALLOCATE( sl( nx , nx ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' pcegterg ',' cannot allocate sl ', ABS(ierr) )
     !
     ALLOCATE( hl( nx , nx ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' pcegterg ',' cannot allocate hl ', ABS(ierr) )
     !
  ELSE
     !
     ALLOCATE( vl( 1 , 1 ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' pcegterg ',' cannot allocate vl ', ABS(ierr) )
     !
     ALLOCATE( sl( 1 , 1 ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' pcegterg ',' cannot allocate sl ', ABS(ierr) )
     !
     ALLOCATE( hl( 1 , 1 ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' pcegterg ',' cannot allocate hl ', ABS(ierr) )
     !
  END IF
  !
  ALLOCATE( ew( nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate ew ', ABS(ierr) )
  !
  ALLOCATE( conv( nvec ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate conv ', ABS(ierr) )
  !
  notcnv = nvec
  nbase  = nvec
  conv   = .FALSE.
  !
  IF ( uspp ) spsi = ZERO
  !
  hpsi = ZERO
  psi  = ZERO
  IF (.not. buffer%is_initialized) CALL buffer%init(4, ierr)
  CALL buffer%lock_buffer(psi_d, (/npwx, npol, nvecx/), ierr)
  CALL buffer%lock_buffer(hpsi_d, (/npwx, npol, nvecx/), ierr)
  CALL buffer%lock_buffer(spsi_d, (/npwx, npol, nvecx/), ierr)
  CALL buffer%lock_buffer(ew_d, nvecx, ierr)
  
  
  evc(:,:,1:nvec) = evc_d(:,:,1:nvec)
  psi(:,:,1:nvec) = evc(:,:,1:nvec)
  !$cuf kernel do(3) <<<*,*>>>
  do k=1,nvec
    do j=1,npol
      do i=1,npwx
        psi_d(i,j,k) = evc_d(i,j,k)
      end do
    end do
  end do
  !
  ! ... hpsi contains h times the basis vectors
  !
  CALL h_psi_gpu( npwx, npw, nvec, psi_d, hpsi_d )
  hpsi(1:npwx, 1:npol, 1:nvec) = hpsi_d(1:npwx, 1:npol, 1:nvec)
  !
  IF ( uspp ) CALL s_psi_gpu( npwx, npw, nvec, psi_d, spsi_d )
  IF ( uspp ) spsi(1:npwx, 1:npol, 1:nvec) = spsi_d(1:npwx, 1:npol, 1:nvec)
  !
  ! ... hl contains the projection of the hamiltonian onto the reduced
  ! ... space, vl contains the eigenvectors of hl. Remember hl, vl and sl
  ! ... are all distributed across processors, global replicated matrixes
  ! ... here are never allocated
  !
  CALL start_clock( 'cegterg:init' )

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
  CALL stop_clock( 'cegterg:init' )
  !
  IF ( lrot ) THEN
     !
     CALL set_e_from_h()
     e_d = e
     !
     CALL set_to_identity( vl, desc )
     !
  ELSE
     !
     ! ... diagonalize the reduced hamiltonian
     !     Calling block parallel algorithm
     !
     CALL start_clock( 'cegterg:diag' )
     IF ( do_distr_diag_inside_bgrp ) THEN ! NB on output of pcdiaghg ew and vl are the same across ortho_parent_comm
        ! only the first bgrp performs the diagonalization
        IF( my_bgrp_id == root_bgrp_id ) CALL pcdiaghg( nbase, hl, sl, nx, ew, vl, desc )
        IF( nbgrp > 1 ) THEN ! results must be brodcast to the other band groups
           CALL mp_bcast( vl, root_bgrp_id, inter_bgrp_comm )
           CALL mp_bcast( ew, root_bgrp_id, inter_bgrp_comm )
        ENDIF
     ELSE
        CALL pcdiaghg( nbase, hl, sl, nx, ew, vl, desc )
     END IF
     CALL stop_clock( 'cegterg:diag' )
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
     CALL start_clock( 'cegterg:update' )
     !
     CALL reorder_v()
     !
     nb1 = nbase + 1
     !
     ! ... expand the basis set with new basis vectors ( H - e*S )|psi> ...
     !
     CALL hpsi_dot_v()
     !
     CALL stop_clock( 'cegterg:update' )
     !
     ! ... approximate inverse iteration
     !
     ew_d = ew
     psi_d(1:npwx, 1:npol, nb1:nb1+notcnv) = psi(1:npwx, 1:npol, nb1:nb1+notcnv)
     CALL g_psi_gpu( npwx, npw, notcnv, npol, psi_d(1,1,nb1), ew_d(nb1) )
     psi(1:npwx, 1:npol, nb1:nb1+notcnv) = psi_d(1:npwx, 1:npol, nb1:nb1+notcnv)
     !
     ! ... "normalize" correction vectors psi(:,nb1:nbase+notcnv) in 
     ! ... order to improve numerical stability of subspace diagonalization 
     ! ... (cdiaghg) ew is used as work array :
     !
     ! ...         ew = <psi_i|psi_i>,  i = nbase + 1, nbase + notcnv
     !
     DO n = 1, notcnv
        !
        nbn = nbase + n
        !
        IF ( npol == 1 ) THEN
           !
           ew(n) = ddot( 2*npw, psi(1,1,nbn), 1, psi(1,1,nbn), 1 )
           !
        ELSE
           !
           ew(n) = ddot( 2*npw, psi(1,1,nbn), 1, psi(1,1,nbn), 1 ) + &
                   ddot( 2*npw, psi(1,2,nbn), 1, psi(1,2,nbn), 1 )
           !
        END IF
        !
     END DO
     !
     CALL mp_sum( ew( 1:notcnv ), intra_bgrp_comm )
     !
     !$omp parallel do collapse(3)
     DO n = 1, notcnv
        DO ipol = 1, npol
           DO m = 1, numblock
              psi((m-1)*blocksize+1:MIN(npw, m*blocksize),ipol,nbase+n) = &
                 psi((m-1)*blocksize+1:MIN(npw, m*blocksize),ipol,nbase+n) / SQRT( ew(n) )
           END DO
        END DO
     END DO
     !$omp end parallel do
     !
     ! ... here compute the hpsi and spsi of the new functions
     !
     psi_d(1:npwx, 1:npol, nb1:nb1+notcnv) = psi(1:npwx, 1:npol, nb1:nb1+notcnv)
     CALL h_psi_gpu( npwx, npw, notcnv, psi_d(1,1,nb1), hpsi_d(1,1,nb1) )
     hpsi(1:npwx, 1:npol, nb1:nb1+notcnv) = hpsi_d(1:npwx, 1:npol, nb1:nb1+notcnv)
     !
     IF ( uspp ) CALL s_psi_gpu( npwx, npw, notcnv, psi_d(1,1,nb1), spsi_d(1,1,nb1) )
     IF ( uspp ) spsi = spsi_d
     !
     ! ... update the reduced hamiltonian
     !
     CALL start_clock( 'cegterg:overlap' )
     !
     ! we need to save the old descriptor in order to redistribute matrices 
     !
     desc_old = desc
     !
     ! ... RE-Initialize the matrix descriptor
     !
     CALL desc_init( nbase+notcnv, desc, irc_ip, nrc_ip )
     !
     IF( la_proc ) THEN

        !  redistribute hl and sl (see dsqmred), since the dimension of the subspace has changed
        !
        vl = hl
        DEALLOCATE( hl )
        ALLOCATE( hl( nx , nx ), STAT=ierr )
        IF( ierr /= 0 ) &
           CALL errore( ' pcegterg ',' cannot allocate hl ', ABS(ierr) )

        CALL zsqmred( nbase, vl, desc_old%nrcx, desc_old, nbase+notcnv, hl, nx, desc )

        vl = sl
        DEALLOCATE( sl )
        ALLOCATE( sl( nx , nx ), STAT=ierr )
        IF( ierr /= 0 ) &
           CALL errore( ' pcegterg ',' cannot allocate sl ', ABS(ierr) )

        CALL zsqmred( nbase, vl, desc_old%nrcx, desc_old, nbase+notcnv, sl, nx, desc )

        DEALLOCATE( vl )
        ALLOCATE( vl( nx , nx ), STAT=ierr )
        IF( ierr /= 0 ) &
           CALL errore( ' pcegterg ',' cannot allocate vl ', ABS(ierr) )

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
     CALL stop_clock( 'cegterg:overlap' )
     !
     nbase = nbase + notcnv
     !
     ! ... diagonalize the reduced hamiltonian
     !     Call block parallel algorithm
     !
     CALL start_clock( 'cegterg:diag' )
     IF ( do_distr_diag_inside_bgrp ) THEN ! NB on output of pcdiaghg ew and vl are the same across ortho_parent_comm
        ! only the first bgrp performs the diagonalization
        IF( my_bgrp_id == root_bgrp_id ) CALL pcdiaghg( nbase, hl, sl, nx, ew, vl, desc )
        IF( nbgrp > 1 ) THEN ! results must be brodcast to the other band groups
           CALL mp_bcast( vl, root_bgrp_id, inter_bgrp_comm )
           CALL mp_bcast( ew, root_bgrp_id, inter_bgrp_comm )
        ENDIF
     ELSE
        CALL pcdiaghg( nbase, hl, sl, nx, ew, vl, desc )
     END IF
     CALL stop_clock( 'cegterg:diag' )
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
        CALL start_clock( 'cegterg:last' )
        !
        CALL refresh_evc()
        evc_d = evc       
        !
        IF ( notcnv == 0 ) THEN
           !
           ! ... all roots converged: return
           !
           CALL stop_clock( 'cegterg:last' )
           !
           EXIT iterate
           !
        ELSE IF ( dav_iter == maxter ) THEN
           !
           ! ... last iteration, some roots not converged: return
           !
           !!!WRITE( stdout, '(5X,"WARNING: ",I5, &
           !!!     &   " eigenvalues not converged")' ) notcnv
           !
           CALL stop_clock( 'cegterg:last' )
           !
           EXIT iterate
           !
        END IF
        !
        ! ... refresh psi, H*psi and S*psi
        !
        CALL threaded_memcpy(psi, evc, nvec*npol*npwx*2)
        !
        IF ( uspp ) THEN
           !
           CALL refresh_spsi()
           ! 
        END IF
        !
        CALL refresh_hpsi()
        !
        ! ... refresh the reduced hamiltonian
        !
        nbase = nvec
        !
        CALL desc_init( nvec, desc, irc_ip, nrc_ip )
        !
        IF( la_proc ) THEN
           !
           ! note that nx has been changed by desc_init
           ! we need to re-alloc with the new size.
           !
           DEALLOCATE( vl, hl, sl )
           ALLOCATE( vl( nx, nx ), STAT=ierr )
           IF( ierr /= 0 ) &
              CALL errore( ' pcegterg ',' cannot allocate vl ', ABS(ierr) )
           ALLOCATE( hl( nx, nx ), STAT=ierr )
           IF( ierr /= 0 ) &
              CALL errore( ' pcegterg ',' cannot allocate hl ', ABS(ierr) )
           ALLOCATE( sl( nx, nx ), STAT=ierr )
           IF( ierr /= 0 ) &
              CALL errore( ' pcegterg ',' cannot allocate sl ', ABS(ierr) )
           !
        END IF
        !
        CALL set_h_from_e( )
        !
        CALL set_to_identity( vl, desc )
        CALL set_to_identity( sl, desc )
        !
        CALL stop_clock( 'cegterg:last' )
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
  DEALLOCATE( ew )
  DEALLOCATE( evc )
  DEALLOCATE( e )
  
  CALL buffer%release_buffer(psi_d, ierr)
  CALL buffer%release_buffer(hpsi_d, ierr)
  CALL buffer%release_buffer(spsi_d, ierr)
  CALL buffer%release_buffer(ew_d, ierr)
  !
  IF ( uspp ) DEALLOCATE( spsi )
  !
  DEALLOCATE( hpsi )
  DEALLOCATE( psi )  
  !
  CALL stop_clock( 'cegterg' )
  !call print_clock( 'cegterg' )
  !call print_clock( 'cegterg:init' )
  !call print_clock( 'cegterg:diag' )
  !call print_clock( 'cegterg:update' )
  !call print_clock( 'cegterg:overlap' )
  !call print_clock( 'cegterg:last' )
  !
  RETURN
  !
  !
CONTAINS
  !
  !
  SUBROUTINE desc_init( nsiz, desc, irc_ip, nrc_ip )
     !
     INTEGER, INTENT(IN)  :: nsiz
     TYPE(la_descriptor), INTENT(OUT) :: desc
     INTEGER, INTENT(OUT) :: irc_ip(:) 
     INTEGER, INTENT(OUT) :: nrc_ip(:) 
     INTEGER :: i, j, rank
     !
     CALL descla_init( desc, nsiz, nsiz, np_ortho, me_ortho, ortho_comm, ortho_cntx, ortho_comm_id )
     !
     nx = desc%nrcx
     !
     DO j = 0, desc%npc - 1
        CALL descla_local_dims( irc_ip( j + 1 ), nrc_ip( j + 1 ), desc%n, desc%nx, np_ortho(1), j )
        DO i = 0, desc%npr - 1
           CALL GRID2D_RANK( 'R', desc%npr, desc%npc, i, j, rank )
           rank_ip( i+1, j+1 ) = rank * leg_ortho
        END DO
     END DO
     !
     la_proc = .FALSE.
     IF( desc%active_node > 0 ) la_proc = .TRUE.
     !
     RETURN
  END SUBROUTINE desc_init
  !
  !
  SUBROUTINE set_to_identity( distmat, desc )
     TYPE(la_descriptor), INTENT(IN)  :: desc
     COMPLEX(DP), INTENT(OUT) :: distmat(:,:)
     INTEGER :: i
     distmat = ( 0_DP , 0_DP )
     IF( desc%myc == desc%myr .AND. desc%active_node > 0 ) THEN
        DO i = 1, desc%nc
           distmat( i, i ) = ( 1_DP , 0_DP )
        END DO
     END IF 
     RETURN
  END SUBROUTINE set_to_identity
  !
  !
  SUBROUTINE reorder_v()
     !
     INTEGER :: ipc
     INTEGER :: nc, ic
     INTEGER :: nl, npl
     !
     np = 0
     !
     notcnv_ip = 0
     !
     n = 0
     !
     DO ipc = 1, desc%npc
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
                    IF( la_proc .AND. desc%myc == ipc-1 ) THEN
                       vl( :, npl) = vl( :, nl )
                    END IF
                 END IF
                 !
                 ! ... for use in g_psi
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
     INTEGER :: nr, ir, ic, notcl, root, np, ipol, ib
     COMPLEX(DP), ALLOCATABLE :: vtmp( :, : )
     COMPLEX(DP), ALLOCATABLE :: ptmp( :, :, : )
     COMPLEX(DP) :: beta

     ALLOCATE( vtmp( nx, nx ) )
     ALLOCATE( ptmp( npwx, npol, nx ) )

     DO ipc = 1, desc%npc
        !
        IF( notcnv_ip( ipc ) > 0 ) THEN

           notcl = notcnv_ip( ipc )
           ic    = ic_notcnv( ipc )

           beta = ZERO

           DO ipr = 1, desc%npr
              !
              nr = nrc_ip( ipr )
              ir = irc_ip( ipr )
              !
              root = rank_ip( ipr, ipc )

              IF( ipr-1 == desc%myr .AND. ipc-1 == desc%myc .AND. la_proc ) THEN
                 vtmp(:,1:notcl) = vl(:,1:notcl)
              END IF

              CALL mp_bcast( vtmp(:,1:notcl), root, ortho_parent_comm )
              !
              IF ( uspp ) THEN
                 !
                 CALL ZGEMM( 'N', 'N', kdim, notcl, nr, ONE, &
                    spsi( 1, 1, ir ), kdmx, vtmp, nx, beta, psi(1,1,nb1+ic-1), kdmx )
                 !
              ELSE
                 !
                 CALL ZGEMM( 'N', 'N', kdim, notcl, nr, ONE, &
                    psi( 1, 1, ir ), kdmx, vtmp, nx, beta, psi(1,1,nb1+ic-1), kdmx )
                 !
              END IF
              !
              CALL ZGEMM( 'N', 'N', kdim, notcl, nr, ONE, &
                      hpsi( 1, 1, ir ), kdmx, vtmp, nx, beta, ptmp, kdmx )

              beta = ONE

           END DO

           !$omp parallel do collapse(3)
           DO np = 1, notcl
              DO ipol = 1, npol
                 DO ib = 1, numblock
                    !
                    psi((ib-1)*blocksize+1:MIN(npw, ib*blocksize),ipol,nbase+np+ic-1) = &
                       ptmp((ib-1)*blocksize+1:MIN(npw, ib*blocksize),ipol,np) &
                       - ew(nbase+np+ic-1) * psi((ib-1)*blocksize+1:MIN(npw, ib*blocksize),ipol,nbase+np+ic-1)
                    !
                 END DO
              END DO
           END DO
           !$omp end parallel do
           !
           ! clean up garbage if there is any
           IF (npw < npwx) psi(npw+1:npwx,:,nbase+ic:nbase+notcl+ic-1) = ZERO
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
     COMPLEX(DP), ALLOCATABLE :: vtmp( :, : )
     COMPLEX(DP) :: beta

     ALLOCATE( vtmp( nx, nx ) )
     !
     DO ipc = 1, desc%npc
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        IF( ic <= nvec ) THEN
           !
           nc = min( nc, nvec - ic + 1 )
           !
           beta = ZERO

           DO ipr = 1, desc%npr
              !
              nr = nrc_ip( ipr )
              ir = irc_ip( ipr )
              !
              root = rank_ip( ipr, ipc )

              IF( ipr-1 == desc%myr .AND. ipc-1 == desc%myc .AND. la_proc ) THEN
                 !
                 !  this proc sends his block
                 ! 
                 CALL mp_bcast( vl(:,1:nc), root, ortho_parent_comm )
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ONE, &
                          psi(1,1,ir), kdmx, vl, nx, beta, evc(1,1,ic), kdmx )
              ELSE
                 !
                 !  all other procs receive
                 ! 
                 CALL mp_bcast( vtmp(:,1:nc), root, ortho_parent_comm )
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ONE, &
                          psi(1,1,ir), kdmx, vtmp, nx, beta, evc(1,1,ic), kdmx )
              END IF
              ! 

              beta = ONE

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
     COMPLEX(DP), ALLOCATABLE :: vtmp( :, : )
     COMPLEX(DP) :: beta

     ALLOCATE( vtmp( nx, nx ) )
     !
     DO ipc = 1, desc%npc
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        IF( ic <= nvec ) THEN
           !
           nc = min( nc, nvec - ic + 1 )
           !
           beta = ZERO
           !
           DO ipr = 1, desc%npr
              !
              nr = nrc_ip( ipr )
              ir = irc_ip( ipr )
              !
              root = rank_ip( ipr, ipc )

              IF( ipr-1 == desc%myr .AND. ipc-1 == desc%myc .AND. la_proc ) THEN
                 !
                 !  this proc sends his block
                 ! 
                 CALL mp_bcast( vl(:,1:nc), root, ortho_parent_comm )
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ONE, &
                          spsi(1,1,ir), kdmx, vl, nx, beta, psi(1,1,nvec+ic), kdmx )
              ELSE
                 !
                 !  all other procs receive
                 ! 
                 CALL mp_bcast( vtmp(:,1:nc), root, ortho_parent_comm )
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ONE, &
                          spsi(1,1,ir), kdmx, vtmp, nx, beta, psi(1,1,nvec+ic), kdmx )
              END IF
              ! 
              beta = ONE

           END DO
           !
        END IF
        !
     END DO
     !
     CALL threaded_memcpy(spsi, psi(1,1,nvec+1), nvec*npol*npwx*2)
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
     COMPLEX(DP), ALLOCATABLE :: vtmp( :, : )
     COMPLEX(DP) :: beta

     ALLOCATE( vtmp( nx, nx ) )
     !
     DO ipc = 1, desc%npc
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        IF( ic <= nvec ) THEN
           !
           nc = min( nc, nvec - ic + 1 )
           !
           beta = ZERO
           !
           DO ipr = 1, desc%npr
              !
              nr = nrc_ip( ipr )
              ir = irc_ip( ipr )
              !
              root = rank_ip( ipr, ipc )

              IF( ipr-1 == desc%myr .AND. ipc-1 == desc%myc .AND. la_proc ) THEN
                 !
                 !  this proc sends his block
                 ! 
                 CALL mp_bcast( vl(:,1:nc), root, ortho_parent_comm )
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ONE, &
                          hpsi(1,1,ir), kdmx, vl, nx, beta, psi(1,1,nvec+ic), kdmx )
              ELSE
                 !
                 !  all other procs receive
                 ! 
                 CALL mp_bcast( vtmp(:,1:nc), root, ortho_parent_comm )
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ONE, &
                          hpsi(1,1,ir), kdmx, vtmp, nx, beta, psi(1,1,nvec+ic), kdmx )
              END IF
              ! 
              beta = ONE

           END DO
           !
        END IF
        !
     END DO
     !
     DEALLOCATE( vtmp )
     !
     CALL threaded_memcpy(hpsi, psi(1,1,nvec+1), nvec*npol*npwx*2)
     !
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
     COMPLEX(DP), INTENT(OUT) :: dm( :, : )
     COMPLEX(DP) :: v(:,:,:), w(:,:,:)
     COMPLEX(DP), ALLOCATABLE :: work( :, : )
     !
     ALLOCATE( work( nx, nx ) )
     !
     work = ZERO
     !
     !  Only upper triangle is computed, then the matrix is hermitianized
     !
     DO ipc = 1, desc%npc !  loop on column procs 
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        DO ipr = 1, ipc ! desc%npr ! ipc ! use symmetry for the loop on row procs
           !
           nr = nrc_ip( ipr )
           ir = irc_ip( ipr )
           !
           !  rank of the processor for which this block (ipr,ipc) is destinated
           !
           root = rank_ip( ipr, ipc )

           ! use blas subs. on the matrix block

           CALL ZGEMM( 'C', 'N', nr, nc, kdim, ONE , &
                       v(1,1,ir), kdmx, w(1,1,ic), kdmx, ZERO, work, nx )

           ! accumulate result on dm of root proc.

           CALL mp_root_sum( work, dm, root, ortho_parent_comm )

        END DO
        !
     END DO
     if (ortho_parent_comm.ne.intra_bgrp_comm .and. nbgrp > 1) dm = dm/nbgrp
     !
     !  The matrix is hermitianized using upper triangle
     !
     CALL zsqmher( nbase, dm, nx, desc )
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
     COMPLEX(DP) :: dm( :, : )
     COMPLEX(DP) :: v(:,:,:), w(:,:,:)
     COMPLEX(DP), ALLOCATABLE :: vtmp( :, : )

     ALLOCATE( vtmp( nx, nx ) )
     !
     vtmp = ZERO
     !
     DO ipc = 1, desc%npc
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        IF( ic+nc-1 >= nb1 ) THEN
           !
           nc = MIN( nc, ic+nc-1 - nb1 + 1 )
           IF( ic >= nb1 ) THEN
              ii = ic
              icc = 1
           ELSE
              ii = nb1
              icc = nb1-ic+1
           END IF
           !
           ! icc to nc is the local index of the unconverged bands
           ! ii is the global index of the first unconverged bands
           !
           DO ipr = 1, ipc ! desc%npr use symmetry
              !
              nr = nrc_ip( ipr )
              ir = irc_ip( ipr )
              !
              root = rank_ip( ipr, ipc )

              CALL ZGEMM( 'C', 'N', nr, nc, kdim, ONE, v( 1, 1, ir ), &
                          kdmx, w(1,1,ii), kdmx, ZERO, vtmp, nx )
              IF (ortho_parent_comm.ne.intra_bgrp_comm .and. nbgrp > 1) vtmp = vtmp/nbgrp
              !
              IF(  (desc%active_node > 0) .AND. (ipr-1 == desc%myr) .AND. (ipc-1 == desc%myc) ) THEN
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
     CALL zsqmher( nbase+notcnv, dm, nx, desc )
     !
     DEALLOCATE( vtmp )
     RETURN
  END SUBROUTINE update_distmat
  !
  !
  SUBROUTINE set_e_from_h()
     INTEGER :: nc, ic, i
     e(1:nbase) = 0_DP
     IF( desc%myc == desc%myr .AND. la_proc ) THEN
        nc = desc%nc
        ic = desc%ic
        DO i = 1, nc
           e( i + ic - 1 ) = REAL( hl( i, i ) )
        END DO
     END IF
     CALL mp_sum( e(1:nbase), ortho_parent_comm )
     RETURN
  END SUBROUTINE set_e_from_h
  !
  SUBROUTINE set_h_from_e()
     INTEGER :: nc, ic, i
     IF( la_proc ) THEN
        hl = ZERO
        IF( desc%myc == desc%myr ) THEN
           nc = desc%nc
           ic = desc%ic
           DO i = 1, nc
              hl(i,i) = CMPLX( e( i + ic - 1 ), 0_DP ,kind=DP)
           END DO
        END IF
     END IF
     RETURN
  END SUBROUTINE set_h_from_e
  !
END SUBROUTINE pcegterg_gpu

#endif
