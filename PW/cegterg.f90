!
! Copyright (C) 2001-2006 Quantum-ESPRESSO group
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
SUBROUTINE cegterg( ndim, ndmx, nvec, nvecx, evc, ethr, &
                    uspp, e, btype, notcnv, lrot, dav_iter )
  !----------------------------------------------------------------------------
  !
  ! ... iterative solution of the eigenvalue problem:
  !
  ! ... ( H - e S ) * evc = 0
  !
  ! ... where H is an hermitean operator, e is a real scalar,
  ! ... S is an uspp matrix, evc is a complex vector
  !
  USE kinds,            ONLY : DP
  USE io_global,        ONLY : stdout
  USE noncollin_module, ONLY : noncolin, npol
  USE bp,               ONLY : lelfield
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ndim, ndmx, nvec, nvecx
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix evc, as declared in the calling pgm unit
    ! integer number of searched low-lying roots
    ! maximum dimension of the reduced basis set :
    !    (the basis set is refreshed when its dimension would exceed nvecx)
    ! k-point considered
  COMPLEX(DP), INTENT(INOUT) :: evc(ndmx,npol,nvec)
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
  REAL(DP), INTENT(OUT) :: e(nvec)
    ! contains the estimated roots.
  INTEGER, INTENT(OUT) :: dav_iter, notcnv
    ! integer number of iterations performed
    ! number of unconverged roots
  !
  ! ... LOCAL variables
  !
  INTEGER, PARAMETER :: maxter = 20
    ! maximum number of iterations
  !
  INTEGER :: kter, nbase, np, kdim, kdmx, n, m, nb1, nbn
    ! counter on iterations
    ! dimension of the reduced basis
    ! counter on the reduced basis vectors
    ! adapted ndim and ndmx
    ! do-loop counters
  COMPLEX(DP), ALLOCATABLE :: hc(:,:), sc(:,:), vc(:,:)
    ! Hamiltonian on the reduced basis
    ! S matrix on the reduced basis
    ! the eigenvectors of the Hamiltonian
  COMPLEX(DP), ALLOCATABLE :: psi(:,:,:), hpsi(:,:,:), spsi(:,:,:)
    ! work space, contains psi
    ! the product of H and psi
    ! the product of S and psi
  REAL(DP), ALLOCATABLE :: ew(:)
    ! eigenvalues of the reduced hamiltonian
  LOGICAL, ALLOCATABLE  :: conv(:)
    ! true if the root is converged
  REAL(DP) :: empty_ethr 
    ! threshold for empty bands
  !
  REAL(DP), EXTERNAL :: DDOT
  !
  EXTERNAL  h_psi,    s_psi,    g_psi
  EXTERNAL  h_psi_nc, s_psi_nc, g_psi_nc
    ! h_psi(ndmx,ndim,nvec,psi,hpsi)
    !     calculates H|psi>
    ! s_psi(ndmx,ndim,nvec,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (ndmx,npol,nvec)
    ! g_psi(ndmx,ndim,notcnv,psi,e)
    !    calculates (diag(h)-e)^-1 * psi, diagonal approx. to (h-e)^-1*psi
    !    the first nvec columns contain the trial eigenvectors
  !
  !
  CALL start_clock( 'cegterg' )
  !
  IF ( nvec > nvecx / 2 ) CALL errore( 'cegter', 'nvecx is too small', 1 )
  !
  ! ... threshold for empty bands
  !
  empty_ethr = MAX( ( ethr * 5.D0 ), 1.D-5 )
  !
  IF ( npol == 1 ) THEN
     !
     kdim = ndim
     kdmx = ndmx
     !
  ELSE
     !
     kdim = ndmx*npol
     kdmx = ndmx*npol
     !
  END IF
  !
  ALLOCATE(  psi( ndmx, npol, nvecx ) )
  ALLOCATE( hpsi( ndmx, npol, nvecx ) )
  !
  IF ( uspp ) ALLOCATE( spsi( ndmx, npol, nvecx ) )
  !
  ALLOCATE( sc( nvecx, nvecx ) )
  ALLOCATE( hc( nvecx, nvecx ) )
  ALLOCATE( vc( nvecx, nvecx ) )
  ALLOCATE( ew( nvecx ) )
  ALLOCATE( conv( nvec ) )
  !
  notcnv = nvec
  nbase  = nvec
  conv   = .FALSE.
  !
  IF ( uspp ) spsi = ZERO
  !
  hpsi = ZERO
  psi  = ZERO
  psi(:,:,1:nvec) = evc(:,:,1:nvec)
  !
  ! ... hpsi contains h times the basis vectors
  !
  IF ( noncolin ) THEN
     !
     CALL h_psi_nc( ndmx, ndim, nvec, psi, hpsi )
     !
     IF ( lelfield ) CALL h_epsi_her_apply( ndmx, ndim, nvec, psi, hpsi )
     IF ( uspp ) CALL s_psi_nc( ndmx, ndim, nvec, psi, spsi )
     !
  ELSE
     !
     CALL h_psi( ndmx, ndim, nvec, psi, hpsi )
     !
     IF ( lelfield ) CALL h_epsi_her_apply( ndmx, ndim, nvec, psi, hpsi )
     IF ( uspp ) CALL s_psi( ndmx, ndim, nvec, psi, spsi )
     !
  END IF
  !
  ! ... hc contains the projection of the hamiltonian onto the reduced 
  ! ... space vc contains the eigenvectors of hc
  !
  hc(:,:) = ZERO
  sc(:,:) = ZERO
  vc(:,:) = ZERO
  !
  CALL ZGEMM( 'C', 'N', nbase, nbase, kdim, ONE, &
              psi, kdmx, hpsi, kdmx, ZERO, hc, nvecx )
  !
  CALL reduce( 2*nbase*nvecx, hc )
  !
  IF ( uspp ) THEN
     !
     CALL ZGEMM( 'C', 'N', nbase, nbase, kdim, ONE, &
                 psi, kdmx, spsi, kdmx, ZERO, sc, nvecx )
     !     
  ELSE
     !
     CALL ZGEMM( 'C', 'N', nbase, nbase, kdim, ONE, &
                 psi, kdmx, psi, kdmx, ZERO, sc, nvecx )
     !
  END IF
  !
  CALL reduce( 2*nbase*nvecx, sc )
  !
  IF ( lrot ) THEN
     !
     DO n = 1, nbase
        !
        e(n) = REAL( hc(n,n) )
        !
        vc(n,n) = ONE
        !
     END DO
     !
  ELSE
     !
     ! ... diagonalize the reduced hamiltonian
     !
     CALL cdiaghg( nbase, nvec, hc, sc, nvecx, ew, vc )
     !
     e(1:nvec) = ew(1:nvec)
     !
  END IF
  !
  ! ... iterate
  !
  iterate: DO kter = 1, maxter
     !
     dav_iter = kter
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
           IF ( np /= n ) vc(:,np) = vc(:,n)
           !
           ! ... for use in g_psi
           !
           ew(nbase+np) = e(n)
           !
        END IF
        !
     END DO
     !
     nb1 = nbase + 1
     !
     ! ... expand the basis set with new basis vectors ( H - e*S )|psi> ...
     !
     IF ( uspp ) THEN
        !
        CALL ZGEMM( 'N', 'N', kdim, notcnv, nbase, ONE, spsi, &
                    kdmx, vc, nvecx, ZERO, psi(1,1,nb1), kdmx )
        !     
     ELSE
        !
        CALL ZGEMM( 'N', 'N', kdim, notcnv, nbase, ONE, psi, &
                    kdmx, vc, nvecx, ZERO, psi(1,1,nb1), kdmx )
        !
     END IF
     !
     DO np = 1, notcnv
        !
        psi(:,:,nbase+np) = - ew(nbase+np)*psi(:,:,nbase+np)
        !
     END DO
     !
     CALL ZGEMM( 'N', 'N', kdim, notcnv, nbase, ONE, hpsi, &
                 kdmx, vc, nvecx, ONE, psi(1,1,nb1), kdmx )
     !
     CALL stop_clock( 'update' )
     !
     ! ... approximate inverse iteration
     !
     IF ( noncolin ) THEN
        !
        CALL g_psi_nc( ndmx, ndim, notcnv, npol, psi(1,1,nb1), ew(nb1) )
        !
     ELSE
        !
        CALL g_psi( ndmx, ndim, notcnv, psi(1,1,nb1), ew(nb1) )
        !
     END IF
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
           ew(n) = DDOT( 2*ndim, psi(1,1,nbn), 1, psi(1,1,nbn), 1 )
           !
        ELSE
           !
           ew(n) = DDOT( 2*ndim, psi(1,1,nbn), 1, psi(1,1,nbn), 1 ) + &
                   DDOT( 2*ndim, psi(1,2,nbn), 1, psi(1,2,nbn), 1 )
           !
        END IF
        !
     END DO
     !
     CALL reduce( notcnv, ew )
     !
     DO n = 1, notcnv
        !
        psi(:,:,nbase+n) = psi(:,:,nbase+n) / SQRT( ew(n) )
        !
     END DO
     !
     ! ... here compute the hpsi and spsi of the new functions
     !
     IF ( noncolin ) THEN
        !
        CALL h_psi_nc( ndmx, ndim, notcnv, psi(1,1,nb1), hpsi(1,1,nb1) )
        !
        IF ( lelfield ) &
           CALL h_epsi_her_apply( ndmx, ndim, notcnv, psi(1,1,nb1), hpsi(1,1,nb1) )
        !
        IF ( uspp ) &
           CALL s_psi_nc( ndmx, ndim, notcnv, psi(1,1,nb1), spsi(1,1,nb1) )
        !
     ELSE
        !
        CALL h_psi( ndmx, ndim, notcnv, psi(1,1,nb1), hpsi(1,1,nb1) )
        !
        IF ( lelfield ) &
           CALL h_epsi_her_apply( ndmx, ndim, notcnv, psi(1,1,nb1), hpsi(1,1,nb1) )
        !
        IF ( uspp ) &
           CALL s_psi( ndmx, ndim, notcnv, psi(1,1,nb1), spsi(1,1,nb1) )
        !
     END IF
     !
     ! ... update the reduced hamiltonian
     !
     CALL start_clock( 'uspp' )
     !
     CALL ZGEMM( 'C', 'N', nbase+notcnv, notcnv, kdim, ONE, psi, &
                 kdmx, hpsi(1,1,nb1), kdmx, ZERO, hc(1,nb1), nvecx )
     !
     CALL reduce( 2*nvecx*notcnv, hc(1,nb1) )
     !
     IF ( uspp ) THEN
        !
        CALL ZGEMM( 'C', 'N', nbase+notcnv, notcnv, kdim, ONE, psi, &
                    kdmx, spsi(1,1,nb1), kdmx, ZERO, sc(1,nb1), nvecx )
        !     
     ELSE
        !
        CALL ZGEMM( 'C', 'N', nbase+notcnv, notcnv, kdim, ONE, psi, &
                    kdmx, psi(1,1,nb1), kdmx, ZERO, sc(1,nb1), nvecx )
        !
     END IF
     !
     CALL reduce( 2*nvecx*notcnv, sc(1,nb1) )
     !
     CALL stop_clock( 'uspp' )
     !
     nbase = nbase + notcnv
     !
     DO n = 1, nbase
        !
        ! ... the diagonal of hc and sc must be strictly real 
        !
        hc(n,n) = CMPLX( REAL( hc(n,n) ), 0.D0 )
        sc(n,n) = CMPLX( REAL( sc(n,n) ), 0.D0 )
        !
        DO m = n + 1, nbase
           !
           hc(m,n) = CONJG( hc(n,m) )
           sc(m,n) = CONJG( sc(n,m) )
           !
        END DO
        !
     END DO
     !
     ! ... diagonalize the reduced hamiltonian
     !
     CALL cdiaghg( nbase, nvec, hc, sc, nvecx, ew, vc )
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
     !
     notcnv = COUNT( .NOT. conv(:) )
     !
     e(1:nvec) = ew(1:nvec)
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
        CALL start_clock( 'last' )
        !
        CALL ZGEMM( 'N', 'N', kdim, nvec, nbase, ONE, &
                    psi, kdmx, vc, nvecx, ZERO, evc, kdmx )
        !
        IF ( notcnv == 0 ) THEN
           !
           ! ... all roots converged: return
           !
           CALL stop_clock( 'last' )
           !
           EXIT iterate
           !
        ELSE IF ( dav_iter == maxter ) THEN
           !
           ! ... last iteration, some roots not converged: return
           !
           WRITE( stdout, '(5X,"WARNING: ",I5, &
                &   " eigenvalues not converged")' ) notcnv
           !
           CALL stop_clock( 'last' )
           !
           EXIT iterate
           !
        END IF
        !
        ! ... refresh psi, H*psi and S*psi
        !
        psi(:,:,1:nvec) = evc(:,:,1:nvec)
        !
        IF ( uspp ) THEN
           !
           CALL ZGEMM( 'N', 'N', kdim, nvec, nbase, ONE, spsi, &
                       kdmx, vc, nvecx, ZERO, psi(1,1,nvec+1), kdmx )
           !
           spsi(:,:,1:nvec) = psi(:,:,nvec+1:nvec+nvec)
           !
        END IF
        !
        CALL ZGEMM( 'N', 'N', kdim, nvec, nbase, ONE, hpsi, &
                    kdmx, vc, nvecx, ZERO, psi(1,1,nvec+1), kdmx )
        !
        hpsi(:,:,1:nvec) = psi(:,:,nvec+1:nvec+nvec)
        !
        ! ... refresh the reduced hamiltonian 
        !
        nbase = nvec
        !
        hc(:,1:nbase) = ZERO
        sc(:,1:nbase) = ZERO
        vc(:,1:nbase) = ZERO
        !
        DO n = 1, nbase
           !
           hc(n,n) = REAL( e(n) )
           !
           sc(n,n) = ONE
           vc(n,n) = ONE
           !
        END DO
        !
        CALL stop_clock( 'last' )
        !
     END IF
     !
  END DO iterate
  !
  DEALLOCATE( conv )
  DEALLOCATE( ew )
  DEALLOCATE( vc )
  DEALLOCATE( hc )
  DEALLOCATE( sc )
  !
  IF ( uspp ) DEALLOCATE( spsi )
  !
  DEALLOCATE( hpsi )
  DEALLOCATE( psi )
  !
  CALL stop_clock( 'cegterg' )
  !
  RETURN
  !
END SUBROUTINE cegterg
