!
! Copyright (C) 2003-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE regterg( ndim, ndmx, nvec, nvecx, evc, ethr, &
                    overlap, gstart, e, btype, notcnv, iter )
  !----------------------------------------------------------------------------
  !
  ! ... iterative solution of the eigenvalue problem:
  !
  ! ... ( H - e S ) * evc = 0
  !
  ! ... where H is an hermitean operator, e is a real scalar,
  ! ... S is an overlap matrix, evc is a complex vector
  ! ... (real wavefunctions with only half plane waves stored)
  !
  USE io_global,  ONLY : stdout
  USE kinds,      ONLY : DP
  !
  IMPLICIT NONE
  !
  ! ... on INPUT
  !
  INTEGER :: ndim, ndmx, nvec, nvecx, gstart
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix evc, as declared in the calling pgm unit
    ! integer number of searched low-lying roots
    ! maximum dimension of the reduced basis set
    !    (the basis set is refreshed when its dimension would exceed nvecx)
  COMPLEX (KIND=DP) :: evc(ndmx,nvec)
    !  evc   contains the  refined estimates of the eigenvectors
  REAL (KIND=DP) :: ethr
    ! energy threshold for convergence: root improvement is stopped,
    ! when two consecutive estimates of the root differ by less than ethr.
  LOGICAL :: overlap
    ! if .FALSE. : S|psi> not needed
  INTEGER, INTENT(IN) :: btype(nvec)
    ! band type ( 1 = occupied, 0 = empty )
  !
  ! ... on OUTPUT
  !
  REAL (KIND=DP) :: e(nvec)
    ! contains the estimated roots.
  INTEGER :: iter, notcnv
    ! integer  number of iterations performed
    ! number of unconverged roots
  !
  ! ... LOCAL variables
  !
  INTEGER, PARAMETER :: maxter = 20
    ! maximum number of iterations
  !
  INTEGER :: kter, nbase, np, n, m
    ! counter on iterations
    ! dimension of the reduced basis
    ! counter on the reduced basis vectors
    ! do-loop counters
  REAL (KIND=DP), ALLOCATABLE :: hr(:,:),  sr(:,:), vr(:,:), ew(:)
    ! Hamiltonian on the reduced basis
    ! S matrix on the reduced basis
    ! eigenvectors of the Hamiltonian
    ! eigenvalues of the reduced hamiltonian
  REAL (KIND=DP), EXTERNAL :: DDOT
  COMPLEX (KIND=DP), ALLOCATABLE :: psi(:,:), hpsi(:,:), spsi(:,:)
    ! work space, contains psi
    ! the product of H and psi
    ! the product of S and psi
  LOGICAL, ALLOCATABLE  :: conv(:)
    ! true if the root is converged
  REAL (KIND=DP) :: empty_ethr 
    ! threshold for empty bands
  INTEGER :: ndim2, ndmx2
  !
  ! ... Called routines:
  !
  EXTERNAL  h_psi, s_psi, g_psi
    ! h_psi(ndmx,ndim,nvec,psi,hpsi)
    !     calculates H|psi> 
    ! s_psi(ndmx,ndim,nvec,psi,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (ndmx,nvec)
    ! g_psi(ndmx,ndim,notcnv,psi,e)
    !    calculates (diag(h)-e)^-1 * psi, diagonal approx. to (h-e)^-1*psi
    !    the first nvec columns contain the trial eigenvectors
  !
  !
  CALL start_clock( 'cegterg' )
  !
  ! ... ALLOCATE the work arrays
  !
  ALLOCATE( psi( ndmx, nvecx ) )    
  ALLOCATE( hpsi( ndmx, nvecx ) )    
  IF ( overlap ) ALLOCATE( spsi( ndmx, nvecx ) )    
  ALLOCATE( sr( nvecx, nvecx ) )    
  ALLOCATE( hr( nvecx, nvecx ) )    
  ALLOCATE( vr( nvecx, nvecx ) )    
  ALLOCATE( ew( nvecx ) )    
  ALLOCATE( conv( nvec ) )    
  !
  IF ( nvec > nvecx / 2 ) CALL errore( 'regter', 'nvecx is too small', 1 )
  !
  !
  ! ... threshold for empty bands
  !
  empty_ethr = MAX( ( ethr * 5.D0 ), 1.D-5 )
  !
  ! ... prepare the hamiltonian for the first iteration
  !
  ndim2  = 2 * ndim
  ndmx2  = 2 * ndmx
  notcnv = nvec
  nbase  = nvec
  conv   = .FALSE.
  !
  IF ( overlap ) spsi = ZERO
  psi  = ZERO
  hpsi = ZERO
  psi(:,1:nvec) = evc(:,1:nvec)
  !
  ! ... hpsi contains h times the basis vectors
  !
  CALL h_psi( ndmx, ndim, nvec, psi, hpsi )
  !
  IF ( overlap ) CALL s_psi( ndmx, ndim, nvec, psi, spsi )
  !
  ! ... hr contains the projection of the hamiltonian onto the reduced space
  ! ... vr contains the eigenvectors of hr
  !
  hr(:,:) = 0.D0
  vr(:,:) = 0.D0
  !
  CALL DGEMM( 'T', 'N', nbase, nbase, ndim2, 2.D0 , &
              psi, ndmx2, hpsi, ndmx2, 0.D0, hr, nvecx )
  !
  IF ( gstart == 2 ) &
     CALL DGER( nbase, nbase, -1.D0, psi, ndmx2, hpsi, ndmx2, hr, nvecx )
  !
  CALL reduce( nbase * nvecx, hr )
  !
  sr(:,:) = 0.D0
  !
  IF ( overlap ) THEN
     !
     CALL DGEMM( 'T', 'N', nbase, nbase, ndim2, 2.D0, &
                 psi, ndmx2, spsi, ndmx2, 0.D0, sr, nvecx )
     !
     IF ( gstart == 2 ) &
        CALL DGER( nbase, nbase, -1.D0, psi, ndmx2, spsi, ndmx2, sr, nvecx )
     !
  ELSE
     !
     CALL DGEMM( 'T', 'N', nbase, nbase, ndim2, 2.D0, &
                 psi, ndmx2, psi, ndmx2, 0.D0, sr, nvecx )
     !
     IF ( gstart == 2 ) &
        CALL DGER( nbase, nbase, -1.D0, psi, ndmx2, psi, ndmx2, sr, nvecx )
     !
  END IF
  !
  CALL reduce( nbase * nvecx, sr )
  !
  FORALL( n = 1 : nbase )
     !
     e(n) = hr(n,n)
     !
     vr(n,n) = 1.D0
     !
  END FORALL
  !
  ! ... iterate
  !
  iterate: DO kter = 1, maxter
     !
     iter = kter
     !
     CALL start_clock( 'update' )
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
           IF ( np /= n ) vr(:,np) = vr(:,n)
           !
           ! ... for use in g_psi
           !
           ew(nbase+np) = e(n)
           !   
        END IF
        !
     END DO
     !
     ! ... expand the basis set with new basis vectors ( H - e*S )|psi> ...
     !
     IF ( overlap ) THEN
        !
        CALL DGEMM( 'N', 'N', ndim2, notcnv, nbase, 1.D0, spsi, &
                    ndmx2, vr, nvecx, 0.D0, psi(1,nbase+1), ndmx2 )
        !
     ELSE
        !
        CALL DGEMM( 'N', 'N', ndim2, notcnv, nbase, 1.D0, psi, &
                    ndmx2, vr, nvecx, 0.D0, psi(1,nbase+1), ndmx2 )
        !
     END IF
     !
! workaround for g95 bug
!     FORALL( np = 1: notcnv ) &
     DO np = 1, notcnv
        psi(:,nbase+np) = - ew(nbase+np) * psi(:,nbase+np)
     END DO
     !
     CALL DGEMM( 'N', 'N', ndim2, notcnv, nbase, 1.D0, hpsi, &
                 ndmx2, vr, nvecx, 1.D0, psi(1,nbase+1), ndmx2 )
     !
     CALL stop_clock( 'update' )
     !
     ! ... approximate inverse iteration
     !
     CALL g_psi( ndmx, ndim, notcnv, psi(1,nbase+1), ew(nbase+1) )
     !
     ! ... "normalize" correction vectors psi(:,nbase+1:nbase+notcnv) in order
     ! ... to improve numerical stability of subspace diagonalization rdiaghg
     ! ... ew is used as work array : ew = <psi_i|psi_i>, i=nbase+1,nbase+notcnv
     !
     DO n = 1, notcnv
        !
        ew(n) = 2.D0 * DDOT( ndim2, psi(1,nbase+n), 1, psi(1,nbase+n), 1 )
        !
        IF ( gstart == 2 ) ew(n) = ew(n) - psi(1,nbase+n) * psi(1,nbase+n)
        !
     END DO
     !
     CALL reduce( notcnv, ew )
     !
! workaround for g95 bug
!     FORALL( n = 1 : notcnv )
     DO n = 1, notcnv
        !
        psi(:,nbase+n) = psi(:,nbase+n) / SQRT( ew(n) )
        !
!     END FORALL
     END DO
     !
     ! ... here compute the hpsi and spsi of the new functions
     !
     CALL h_psi( ndmx, ndim, notcnv, psi(1,nbase+1), hpsi(1,nbase+1) )
     !
     IF ( overlap ) &
        CALL s_psi( ndmx, ndim, notcnv, psi(1,nbase+1), spsi(1,nbase+1) )
     !
     ! ... update the reduced hamiltonian
     !
     CALL start_clock( 'overlap' )
     !
     CALL DGEMM( 'T', 'N', nbase+notcnv, notcnv, ndim2, 2.D0, psi, &
                 ndmx2, hpsi(1,nbase+1), ndmx2, 0.D0, hr(1,nbase+1), nvecx )
     !
     IF ( gstart == 2 ) &
        CALL DGER( nbase+notcnv, notcnv, -1.D0, psi, ndmx2, &
                   hpsi(1,nbase+1), ndmx2, hr(1,nbase+1), nvecx )
     !
     CALL reduce( nvecx * notcnv, hr(1,nbase+1) )
     !
     IF ( overlap ) THEN
        !
        CALL DGEMM( 'T', 'N', nbase+notcnv, notcnv, ndim2, 2.D0, psi, ndmx2, &
                    spsi(1,nbase+1), ndmx2, 0.D0, sr(1,nbase+1), nvecx )
        !
        IF ( gstart == 2 ) &
           CALL DGER( nbase+notcnv, notcnv, -1.D0, psi, ndmx2, &
                      spsi(1,nbase+1), ndmx2, sr(1,nbase+1), nvecx )
        !
     ELSE
        !
        CALL DGEMM( 'T', 'N', nbase+notcnv, notcnv, ndim2, 2.D0, psi, ndmx2, &
                    psi(1,nbase+1), ndmx2, 0.D0, sr(1,nbase+1) , nvecx )
        !
        IF ( gstart == 2 ) &
           CALL DGER( nbase+notcnv, notcnv, -1.D0, psi, ndmx2, &
                      psi(1,nbase+1), ndmx2, sr(1,nbase+1), nvecx )
        !
     END IF
     !
     CALL reduce( nvecx * notcnv, sr(1,nbase+1) )
     !
     CALL stop_clock( 'overlap' )
     !
     nbase = nbase + notcnv
     !
     FORALL( n = 1 : nbase )
        !
        FORALL( m = n + 1 : nbase )
           !
           hr(m,n) = hr(n,m)
           sr(m,n) = sr(n,m)
           !
        END FORALL
        !
     END FORALL
     !
     ! ... diagonalize the reduced hamiltonian
     !
     CALL rdiaghg( nbase, nvec, hr, sr, nvecx, ew, vr )
     !
     ! ... test for convergence
     !
     WHERE( btype(:) == 1 )
        !
        conv(:) = ( ( ABS( ew(:) - e(:) ) < ethr ) )
        !
     ELSEWHERE
        !
        conv(:) = ( ( ABS( ew(:) - e(:) ) < empty_ethr ) )
        !
     END WHERE
     !
     notcnv = COUNT( .NOT. conv(:) )
     !
     e(1:nvec) = ew(1:nvec)
     !
     ! ... if overall convergence has been achieved, OR
     ! ... the dimension of the reduced basis set is becoming too large, OR
     ! ... in any case if we are at the last iteration
     ! ... refresh the basis set. i.e. replace the first nvec elements
     ! ... with the current estimate of the eigenvectors;
     ! ... set the basis dimension to nvec.
     !
     IF ( notcnv == 0 .OR. nbase+notcnv > nvecx .OR. iter == maxter ) THEN
        !
        CALL start_clock( 'last' )
        !
        CALL DGEMM( 'N', 'N', ndim2, nvec, nbase, 1.D0, &
                    psi, ndmx2, vr, nvecx, 0.D0, evc, ndmx2 )
        !
        IF ( notcnv == 0 ) THEN
           !
           ! ... all roots converged: return
           !
           CALL stop_clock( 'last' )
           !
           EXIT iterate
           !
        ELSE IF ( iter == maxter ) THEN
           !
           ! ... last iteration, some roots not converged: return
           !
           WRITE( UNIT = stdout, &
                  FMT = '("   WARNING: ",i5," eigenvalues not converged")' ) &
                notcnv
           !
           CALL stop_clock( 'last' )
           !
           EXIT iterate
           !
        END IF
        !
        ! ... refresh psi, H*psi and S*psi
        !
        psi(:,1:nvec) = evc(:,1:nvec)
        !
        IF ( overlap ) THEN
           !
           CALL DGEMM( 'N', 'N', ndim2, nvec, nbase, 1.D0, spsi, &
                       ndmx2, vr, nvecx, 0.D0, psi(1,nvec+1), ndmx2 )
           !
           spsi(:,1:nvec) = psi(:,nvec+1:2*nvec)
           !
        END IF
        !
        CALL DGEMM( 'N', 'N', ndim2, nvec, nbase, 1.D0, hpsi, &
                    ndmx2, vr, nvecx, 0.D0, psi(1,nvec+1), ndmx2 )
        !
        hpsi(:,1:nvec) = psi(:,nvec+1:2*nvec)
        !
        ! ... refresh the reduced hamiltonian
        !
        nbase = nvec
        !
        hr(:,1:nbase) = 0.D0
        sr(:,1:nbase) = 0.D0
        vr(:,1:nbase) = 0.D0
        !
        FORALL( n = 1 : nbase )
           !
           hr(n,n) = e(n)
           sr(n,n) = 1.D0
           vr(n,n) = 1.D0
           !
        END FORALL
        !
        CALL stop_clock( 'last' )
        !
     END IF
     !
  END DO iterate
  !
  DEALLOCATE( conv )
  DEALLOCATE( ew )
  DEALLOCATE( vr )
  DEALLOCATE( hr )
  DEALLOCATE( sr )
  IF ( overlap ) DEALLOCATE( spsi )
  DEALLOCATE( hpsi )
  DEALLOCATE( psi )
  !
  CALL stop_clock( 'cegterg' )
  !
  RETURN
  !
END SUBROUTINE regterg
