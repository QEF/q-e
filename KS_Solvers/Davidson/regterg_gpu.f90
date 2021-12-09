!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
#if !defined(__CUDA)
#define cublasDgemm dgemm
! workaround for some old compilers that don't like CUDA fortran code
SUBROUTINE regterg_gpu( )
end SUBROUTINE regterg_gpu
SUBROUTINE pregterg_gpu( )
end SUBROUTINE pregterg_gpu
#else
!----------------------------------------------------------------------------
SUBROUTINE regterg_gpu( h_psi_gpu, s_psi_gpu, uspp, g_psi_gpu, &
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
#if defined(__CUDA)
  use cudafor
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
  COMPLEX(DP), INTENT(INOUT) :: evc_d(npwx,nvec)
#if defined(__CUDA)
  attributes(DEVICE) :: evc_d
#endif
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
  REAL(DP), ALLOCATABLE :: hr_d(:,:), sr_d(:,:), vr_d(:,:), ew_d(:)
#if defined(__CUDA)
  attributes(DEVICE) :: hr_d, sr_d, vr_d, ew_d
#endif
    ! Hamiltonian on the reduced basis
    ! S matrix on the reduced basis
    ! eigenvectors of the Hamiltonian
    ! eigenvalues of the reduced hamiltonian
  COMPLEX(DP), ALLOCATABLE :: psi_d(:,:), hpsi_d(:,:), spsi_d(:,:)
#if defined(__CUDA)
  attributes(DEVICE) :: psi_d, hpsi_d, spsi_d
#endif
    ! work space, contains psi
    ! the product of H and psi
    ! the product of S and psi
  LOGICAL, ALLOCATABLE :: conv(:)
    ! true if the root is converged
  REAL(DP) :: empty_ethr 
    ! threshold for empty bands
  REAL(DP), ALLOCATABLE :: ew_host(:)
  REAL(DP), ALLOCATABLE :: e_host(:)
    ! auxiliary variables for performing dot product
  INTEGER :: i,j,k
  REAL(DP):: aux
  !
  !
  EXTERNAL  h_psi_gpu, s_psi_gpu, g_psi_gpu
    ! h_psi_gpu(npwx,npw,nvec,psi,hpsi)
    !     calculates H|psi> 
    ! s_psi_gpu(npwx,npw,nvec,psi,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx,nvec)
    ! g_psi_gpu(npwx,npw,notcnv,psi,e)
    !    calculates (diag(h)-e)^-1 * psi, diagonal approx. to (h-e)^-1*psi
    !    the first nvec columns contain the trial eigenvectors
  !
  CALL start_clock( 'regterg' )
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
  ALLOCATE(  psi_d( npwx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'regterg ',' cannot allocate psi ', ABS(ierr) )
  ALLOCATE( hpsi_d( npwx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'regterg ',' cannot allocate hpsi ', ABS(ierr) )
  !
  IF ( uspp ) THEN
     ALLOCATE( spsi_d( npwx, nvecx ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' regterg ',' cannot allocate spsi ', ABS(ierr) )
  END IF
  !
  ALLOCATE( sr_d( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' regterg ',' cannot allocate sr_d ', ABS(ierr) )
  ALLOCATE( hr_d( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' regterg ',' cannot allocate hr_d ', ABS(ierr) )
  ALLOCATE( vr_d( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' regterg ',' cannot allocate vr_d ', ABS(ierr) )
  ALLOCATE( ew_d( nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' regterg ',' cannot allocate ew_d ', ABS(ierr) )
  ALLOCATE( e_host( nvec ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' regterg ',' cannot allocate e_host ', ABS(ierr) )
  ALLOCATE( ew_host( nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' regterg ',' cannot allocate ew_host ', ABS(ierr) )
  ALLOCATE( conv( nvec ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'regterg ',' cannot allocate conv ', ABS(ierr) )
  !
  notcnv = nvec
  nbase  = nvec
  conv   = .FALSE.
  !
  IF ( uspp ) spsi_d = ZERO
  !
  hpsi_d = ZERO
  psi_d  = ZERO
  !$cuf kernel do(2) <<<*,*>>>
  DO k=1,nvec
     DO i=1,npwx
         psi_d(i,k) = evc_d(i,k)
         IF ( (gstart == 2) .and. (i == 1) ) psi_d(1,k) = CMPLX( DBLE( psi_d(1,k) ), 0.D0 ,kind=DP)
     END DO
  END DO
  !
  ! ... hpsi contains h times the basis vectors
  !
  CALL h_psi_gpu( npwx, npw, nvec, psi_d, hpsi_d )  ; nhpsi = nvec
  !
  ! ... spsi contains s times the basis vectors
  !
  IF ( uspp ) CALL s_psi_gpu( npwx, npw, nvec, psi_d, spsi_d )
  !
  ! ... hr contains the projection of the hamiltonian onto the reduced 
  ! ... space vr contains the eigenvectors of hr
  !
  CALL start_clock( 'regterg:init' )
  hr_d(:,:) = 0.D0
  sr_d(:,:) = 0.D0
  vr_d(:,:) = 0.D0
  !
  CALL divide(inter_bgrp_comm,nbase,n_start,n_end)
  my_n = n_end - n_start + 1; !write (*,*) nbase,n_start,n_end
  if (n_start .le. n_end) &
  CALL cublasDgemm( 'T','N', nbase, my_n, npw2, 2.D0 , psi_d, npwx2, hpsi_d(1,n_start), npwx2, 0.D0, hr_d(1,n_start), nvecx )
  IF ( gstart == 2 ) CALL KScudaDGER(nbase, my_n, -1.D0, psi_d, npwx2, hpsi_d(1,n_start), npwx2, hr_d(1,n_start), nvecx )
  CALL mp_sum( hr_d( :, 1:nbase ), inter_bgrp_comm )
  !
  CALL mp_sum( hr_d( :, 1:nbase ), intra_bgrp_comm )
  !
  IF ( uspp ) THEN
     !
     if (n_start .le. n_end) &
     CALL cublasDgemm( 'T','N', nbase, my_n, npw2, 2.D0, psi_d, npwx2, spsi_d(1,n_start), npwx2, 0.D0, sr_d(1,n_start), nvecx )
     IF ( gstart == 2 ) CALL KScudaDGER(nbase, my_n, -1.D0, psi_d, npwx2, spsi_d(1,n_start), npwx2, sr_d(1,n_start), nvecx )
     !
  ELSE
     !
     if (n_start .le. n_end) &
     CALL cublasDgemm( 'T','N', nbase, my_n, npw2, 2.D0, psi_d, npwx2, psi_d(1,n_start), npwx2, 0.D0, sr_d(1,n_start), nvecx )
     IF ( gstart == 2 ) CALL KScudaDGER(nbase, my_n, -1.D0, psi_d, npwx2, psi_d(1,n_start), npwx2, sr_d(1,n_start), nvecx )
     !
  END IF
  CALL mp_sum( sr_d( :, 1:nbase ), inter_bgrp_comm )
  !
  CALL mp_sum( sr_d( :, 1:nbase ), intra_bgrp_comm )

  CALL stop_clock( 'regterg:init' )
  !
  IF ( lrot ) THEN
     !
!$cuf kernel do(1) <<<*,*>>>
     DO n = 1, nbase
        !
        e_d(n) = hr_d(n,n)
        vr_d(n,n) = 1.D0
        !
     END DO
     !
  ELSE
     !
     ! ... diagonalize the reduced hamiltonian
     !
     CALL start_clock( 'regterg:diag' )
     IF( my_bgrp_id == root_bgrp_id ) THEN
        CALL diaghg( nbase, nvec, hr_d, sr_d, nvecx, ew_d, vr_d, me_bgrp, root_bgrp, intra_bgrp_comm )
     END IF
     IF( nbgrp > 1 ) THEN
        CALL mp_bcast( vr_d, root_bgrp_id, inter_bgrp_comm )
        CALL mp_bcast( ew_d, root_bgrp_id, inter_bgrp_comm )
     ENDIF
     CALL stop_clock( 'regterg:diag' ) 
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
     CALL start_clock( 'regterg:update' )
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

     CALL reorder_evals_revecs(nbase, nvec, nvecx, conv, e_d, ew_d, vr_d)
     !
     nb1 = nbase + 1
     !
     ! ... expand the basis set with new basis vectors ( H - e*S )|psi> ...
     !
     CALL divide(inter_bgrp_comm,nbase,n_start,n_end)
     my_n = n_end - n_start + 1; !write (*,*) nbase,n_start,n_end
!$cuf kernel do(2) <<<*,*>>>
     DO i=1, notcnv
        DO k=1,npwx
           psi_d(k,nbase+i)=ZERO
        END DO
     END DO
     IF ( uspp ) THEN
        !
        if (n_start .le. n_end) &
        CALL cublasDgemm( 'N','N', npw2, notcnv, my_n, 1.D0, spsi_d(1,n_start), npwx2, vr_d(n_start,1), nvecx, 0.D0, psi_d(1,nb1), npwx2 )
        !     
     ELSE
        !
        if (n_start .le. n_end) &
        CALL cublasDgemm( 'N','N', npw2, notcnv, my_n, 1.D0, psi_d(1,n_start), npwx2, vr_d(n_start,1), nvecx, 0.D0, psi_d(1,nb1), npwx2 )
        !
     END IF
! NB: must not call mp_sum over inter_bgrp_comm here because it is done later to the full correction
     !
!$cuf kernel do(2) <<<*,*>>>
     DO np=1,notcnv
        DO k=1,npwx
          psi_d(k,nbase+np) = - ew_d(nbase+np) * psi_d(k,nbase+np)
        END DO
     END DO
     !
     if (n_start .le. n_end) &
     CALL cublasDgemm( 'N','N', npw2, notcnv, my_n, 1.D0, hpsi_d(1,n_start), npwx2, vr_d(n_start,1), nvecx, 1.D0, psi_d(1,nb1), npwx2 )
     CALL mp_sum( psi_d(:,nb1:nbase+notcnv), inter_bgrp_comm )
     !
     CALL stop_clock( 'regterg:update' )
     !
     ! ... approximate inverse iteration
     !
     CALL g_psi_gpu( npwx, npw, notcnv, 1, psi_d(1,nb1), ew_d(nb1) )
     !
     ! ... "normalize" correction vectors psi(:,nb1:nbase+notcnv) in 
     ! ... order to improve numerical stability of subspace diagonalization 
     ! ... (rdiaghg) ew is used as work array :
     !
     ! ...         ew = <psi_i|psi_i>,  i = nbase + 1, nbase + notcnv
     !
     ! == TO BE OPTIMIZED == !!!!
     !$cuf kernel do(1) <<<*,*>>>
     DO n = 1, notcnv
        !
        nbn = nbase + n
        aux = 0.0_DP
        DO i = 1, npw
          aux = aux + 2.D0 * DCONJG(psi_d(i,nbn)) * psi_d(i, nbn)
        END DO
        !
        IF (gstart == 2) aux = aux - DBLE(psi_d(1,nbn) * psi_d(1,nbn)) ! psi_d(1,nbn) * psi_d(1,nbn)
        ew_d(n) = aux
     END DO
     ! == OPTIMIZE ABOVE ==
     CALL mp_sum( ew_d( 1:notcnv ), intra_bgrp_comm )
     !
!$cuf kernel do(2) <<<*,*>>>
     DO i = 1,notcnv
        DO k=1,npwx
           psi_d(k,nbase+i) = psi_d(k,nbase+i)/SQRT( ew_d(i) )
           IF ( (gstart == 2) .and. (k==1)) psi_d(1,nbase+i) = CMPLX( DBLE(psi_d(1,nbase+i)), 0.D0 ,kind=DP)
        END DO
     END DO
     !
     ! ... here compute the hpsi and spsi of the new functions
     !
     CALL h_psi_gpu( npwx, npw, notcnv, psi_d(1,nb1), hpsi_d(1,nb1) ) ; nhpsi = nhpsi + notcnv
     !
     IF ( uspp ) CALL s_psi_gpu( npwx, npw, notcnv, psi_d(1,nb1), spsi_d(1,nb1) )
     !
     ! ... update the reduced hamiltonian
     !
     CALL start_clock( 'regterg:overlap' )
     !
!$cuf kernel do(2) <<<*,*>>>
     DO i=0,notcnv-1
        DO j=1, nvecx
          hr_d( j, nb1+i )=0.d0
        END DO
     END DO
     CALL divide(inter_bgrp_comm,nbase+notcnv,n_start,n_end)
     my_n = n_end - n_start + 1; !write (*,*) nbase+notcnv,n_start,n_end
     CALL cublasDgemm( 'T','N', my_n, notcnv, npw2, 2.D0, psi_d(1,n_start), npwx2, &
                       hpsi_d(1,nb1), npwx2, 0.D0, hr_d(n_start,nb1), nvecx )
     IF ( gstart == 2 ) CALL KScudaDGER( my_n, notcnv, -1.D0, psi_d(1,n_start), npwx2, &
                                         hpsi_d(1,nb1), npwx2, hr_d(n_start,nb1), nvecx )
     CALL mp_sum( hr_d( :, nb1:nb1+notcnv-1 ), inter_bgrp_comm )
     !
     CALL mp_sum( hr_d( :, nb1:nb1+notcnv-1 ), intra_bgrp_comm )
     !
!$cuf kernel do(2) <<<*,*>>>
     DO i=0,notcnv-1
        DO j=1, nvecx
          sr_d( j, nb1+i )=0.d0
        END DO
     END DO
     CALL divide(inter_bgrp_comm,nbase+notcnv,n_start,n_end)
     my_n = n_end - n_start + 1; !write (*,*) nbase+notcnv,n_start,n_end
     IF ( uspp ) THEN
        !
        CALL cublasDgemm( 'T','N', my_n, notcnv, npw2, 2.D0, psi_d(1,n_start), npwx2, &
                           spsi_d(1,nb1), npwx2, 0.D0, sr_d(n_start,nb1), nvecx )
        IF ( gstart == 2 ) CALL KScudaDGER( my_n, notcnv, -1.D0, psi_d(1,n_start), npwx2, &
                                            spsi_d(1,nb1), npwx2, sr_d(n_start,nb1), nvecx )
        !
     ELSE
        !
        CALL cublasDgemm( 'T','N', my_n, notcnv, npw2, 2.D0, psi_d(1,n_start), npwx2, &
                          psi_d(1,nb1), npwx2, 0.D0, sr_d(n_start,nb1), nvecx )
        IF ( gstart == 2 ) CALL KScudaDGER( my_n, notcnv, -1.D0, psi_d(1,n_start), npwx2, &
                                            psi_d(1,nb1), npwx2, sr_d(n_start,nb1), nvecx )

        !
     END IF
     CALL mp_sum( sr_d( :, nb1:nb1+notcnv-1 ), inter_bgrp_comm )
     !
     CALL mp_sum( sr_d( :, nb1:nb1+notcnv-1 ), intra_bgrp_comm )
     !
     CALL stop_clock( 'regterg:overlap' )
     !
     nbase = nbase + notcnv
     !
!$cuf kernel do(1) <<<*,*>>>
     DO n = 1, nbase
        !
        DO m = n + 1, nbase
           !
           hr_d(m,n) = hr_d(n,m)
           sr_d(m,n) = sr_d(n,m)
           !
        END DO
        !
     END DO
     !
     ! ... diagonalize the reduced hamiltonian
     !
     CALL start_clock( 'regterg:diag' )
     IF( my_bgrp_id == root_bgrp_id ) THEN
        CALL diaghg( nbase, nvec, hr_d, sr_d, nvecx, ew_d, vr_d, me_bgrp, root_bgrp, intra_bgrp_comm )
     END IF
     IF( nbgrp > 1 ) THEN
        CALL mp_bcast( vr_d, root_bgrp_id, inter_bgrp_comm )
        CALL mp_bcast( ew_d, root_bgrp_id, inter_bgrp_comm )
     ENDIF
     CALL stop_clock( 'regterg:diag' )
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
        CALL start_clock( 'regterg:last' )
        !
        !$cuf kernel do(2) <<<*,*>>>
        DO k=1,nvec
           DO i=1,npwx
              evc_d(i,k) = ZERO
           END DO
        END DO
        CALL divide(inter_bgrp_comm,nbase,n_start,n_end)
        my_n = n_end - n_start + 1; !write (*,*) nbase,n_start,n_end
        CALL cublasDgemm( 'N','N', npw2, nvec, my_n, 1.D0, psi_d(1,n_start), npwx2, &
                          vr_d(n_start,1), nvecx, 0.D0, evc_d, npwx2 )
        CALL mp_sum( evc_d, inter_bgrp_comm )
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
        !$cuf kernel do(2) <<<*,*>>>
        DO i=1,nvec
           DO k=1,npwx
              psi_d(k,i) = evc_d(k,i)
           END DO
        END DO
        !
        IF ( uspp ) THEN
           !
           psi_d(:,nvec+1:nvec+nvec) = ZERO
           CALL cublasDgemm( 'N','N', npw2, nvec, my_n, 1.D0, spsi_d(1,n_start), npwx2, &
                             vr_d(n_start,1), nvecx, 0.D0, psi_d(1,nvec+1), npwx2 )
           CALL mp_sum( psi_d(:,nvec+1:nvec+nvec), inter_bgrp_comm )
           !
           !$cuf kernel do(2) <<<*,*>>>
           DO i=1,nvec
              DO k=lbound(psi_d,1),ubound(psi_d,1)
                 spsi_d(k,i) = psi_d(k,i+nvec)
              END DO
           END DO
           !
        END IF
        !
        psi_d(:,nvec+1:nvec+nvec) = ZERO
        CALL cublasDgemm( 'N','N', npw2, nvec, my_n, 1.D0, hpsi_d(1,n_start), npwx2, &
                          vr_d(n_start,1), nvecx, 0.D0, psi_d(1,nvec+1), npwx2 )
        CALL mp_sum( psi_d(:,nvec+1:nvec+nvec), inter_bgrp_comm )
        !
        !hpsi_d(:,1:nvec) = psi_d(:,nvec+1:nvec+nvec)
        !$cuf kernel do(2) <<<*,*>>>
        DO i=1,nvec
           DO k=lbound(psi_d,1),ubound(psi_d,1)
              hpsi_d(k,i) = psi_d(k,i+nvec)
           END DO
        END DO
        !
        ! ... refresh the reduced hamiltonian
        !
        nbase = nvec
        !
        hr_d(:,1:nbase) = 0.D0
        sr_d(:,1:nbase) = 0.D0
        vr_d(:,1:nbase) = 0.D0
        !
        !$cuf kernel do(1) <<<*,*>>>
        DO n = 1, nbase
           !
           hr_d(n,n) = e_d(n)
           sr_d(n,n) = 1.D0
           vr_d(n,n) = 1.D0
           !
        END DO
        !
        CALL stop_clock( 'regterg:last' )
        !
     END IF
     !
  END DO iterate
  !
  DEALLOCATE( conv )
  DEALLOCATE( e_host, ew_host, ew_d )
  DEALLOCATE( vr_d )
  DEALLOCATE( hr_d )
  DEALLOCATE( sr_d )
  !
  IF ( uspp ) DEALLOCATE( spsi_d )
  !
  DEALLOCATE( hpsi_d )
  DEALLOCATE( psi_d )
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
END SUBROUTINE regterg_gpu

SUBROUTINE reorder_evals_revecs(nbase, nvec, nvecx, conv, e_d, ew_d, v_d)
  USE util_param,   ONLY : DP
  USE device_fbuff_m,  ONLY : buffer => dev_buf
  implicit none
  INTEGER, INTENT(IN) :: nbase, nvec, nvecx
  LOGICAL, INTENT(IN) :: conv(nvec)
  REAL(DP) :: e_d(nvecx), ew_d(nvecx)
  REAL(DP) :: v_d(nvecx,nvecx)
#if defined(__CUDA)
  attributes(DEVICE) :: e_d, ew_d, v_d
#endif
  !
  INTEGER :: j, k, n, np, info
  INTEGER, ALLOCATABLE :: conv_idx(:)
  INTEGER, POINTER :: conv_idx_d(:)
  REAL(DP), POINTER :: vtmp_d(:,:)
#if defined(__CUDA)
  attributes(DEVICE) :: conv_idx_d, vtmp_d
#endif
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
  IF( info /= 0 ) &
     CALL errore( ' reorder_evals_revecs ',' cannot allocate vtmp_d ', ABS(info) )

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
END SUBROUTINE reorder_evals_revecs


!
!  Wrapper for subroutine with distributed matrixes (written by Carlo Cavazzoni)
!
!----------------------------------------------------------------------------
SUBROUTINE pregterg_gpu(h_psi_gpu, s_psi_gpu, uspp, g_psi_gpu, &  
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
  USE device_fbuff_m,  ONLY : buffer => dev_buf
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
#if defined(__CUDA)
  attributes(DEVICE) :: evc_d
#endif
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
  COMPLEX(DP), ALLOCATABLE :: evc(:,:)
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
  REAL(DP), POINTER :: ew_d(:)
#if defined(__CUDA)
  attributes(DEVICE) :: ew_d
#endif
  REAL(DP), ALLOCATABLE :: hl(:,:), sl(:,:), vl(:,:)
    ! Hamiltonian on the reduced basis
    ! S matrix on the reduced basis
    ! eigenvectors of the Hamiltonian
    ! eigenvalues of the reduced hamiltonian
  COMPLEX(DP), ALLOCATABLE :: psi(:,:), hpsi(:,:), spsi(:,:)
  COMPLEX(DP), POINTER :: psi_d(:,:), hpsi_d(:,:), spsi_d(:,:)
#if defined(__CUDA)
  attributes(DEVICE) :: psi_d, hpsi_d, spsi_d
#endif
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
  EXTERNAL  h_psi_gpu, s_psi_gpu, g_psi_gpu
    ! h_psi(npwx,npw,nvec,psi,hpsi)
    !     calculates H|psi> 
    ! s_psi(npwx,npw,nvec,psi,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx,nvec)
    ! g_psi(npwx,npw,notcnv,psi,e)
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
  ALLOCATE(  evc( npwx, nvec ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pregterg ',' cannot allocate evc (host) ', ABS(ierr) )
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
  IF ( uspp ) spsi = ZERO
  !
  hpsi = ZERO
  psi  = ZERO
  CALL buffer%lock_buffer(psi_d, (/npwx, nvecx/), info)  
  CALL buffer%lock_buffer(hpsi_d, (/npwx, nvecx/), info)  
  CALL buffer%lock_buffer(spsi_d, (/npwx, nvecx/), info)  
  CALL buffer%lock_buffer(ew_d, nvecx, info)

  
  evc(:,1:nvec) = evc_d(:,1:nvec)
  psi(:,1:nvec) = evc(:,1:nvec)
  ! ... set Im[ psi(G=0) ] -  needed for numerical stability
  IF ( gstart == 2 ) psi(1,1:nvec) = CMPLX( DBLE( psi(1,1:nvec) ), 0.D0 ,kind=DP)

  !$cuf kernel do(2) <<<*,*>>>
  do k=1,nvec
     do i=1,npwx
        psi_d(i,k) = evc_d(i,k)
        ! ... set Im[ psi(G=0) ] -  needed for numerical stability
        IF ( (gstart == 2) .and. (i == 1) ) psi_d(1,k) = CMPLX( DBLE( psi_d(1,k) ), 0.D0 ,kind=DP)
     end do
  end do
  !
  ! ... hpsi contains h times the basis vectors
  !
  CALL h_psi_gpu( npwx, npw, nvec, psi_d, hpsi_d )  ; nhpsi = nvec
  hpsi = hpsi_d
  !
  IF ( uspp ) CALL s_psi_gpu( npwx, npw, nvec, psi_d, spsi_d )
  IF ( uspp ) spsi = spsi_d
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
     ew_d = ew
     psi_d = psi
     CALL g_psi_gpu( npwx, npw, notcnv, 1, psi_d(1,nb1), ew_d(nb1) )
     psi = psi_d
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
     psi_d = psi
     CALL h_psi_gpu( npwx, npw, notcnv, psi_d(1,nb1), hpsi_d(1,nb1) ) ; nhpsi = nhpsi + notcnv
     hpsi = hpsi_d
     !
     IF ( uspp ) CALL s_psi_gpu( npwx, npw, notcnv, psi_d(1,nb1), spsi_d(1,nb1) )
     IF ( uspp ) spsi = spsi_d
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
        evc_d = evc
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
        psi(:,1:nvec) = evc(:,1:nvec)
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
  DEALLOCATE( ew )
  DEALLOCATE( evc )
  DEALLOCATE( e )
  
  CALL buffer%release_buffer(psi_d, info)  
  CALL buffer%release_buffer(hpsi_d, info)  
  CALL buffer%release_buffer(spsi_d, info)  
  CALL buffer%release_buffer(ew_d, info)  
  !
  IF ( uspp ) DEALLOCATE( spsi )
  !
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
                          psi(1,ir), npwx2, vl, nx, beta, evc(1,ic), npwx2 )
              ELSE
                 !
                 !  all other procs receive
                 ! 
                 CALL mp_bcast( vtmp(:,1:nc), root, ortho_parent_comm )
                 CALL DGEMM( 'N', 'N', npw2, nc, nr, 1.D0, &
                          psi(1,ir), npwx2, vtmp, nx, beta, evc(1,ic), npwx2 )
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


! In principle this can go away .......
SUBROUTINE KScudaDGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
#if defined(__CUDA)
    use cudafor
    use cublas
#endif
!     .. Scalar Arguments ..
    DOUBLE PRECISION ::  ALPHA
    INTEGER          ::   INCX, INCY, LDA, M, N
!     .. Array Arguments ..
    DOUBLE PRECISION :: A( LDA, * ), X( * ), Y( * )
#if defined(__CUDA)
    attributes(device) :: A, X, Y
#endif
    CALL DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )

END SUBROUTINE KScudaDGER
#endif
