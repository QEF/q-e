!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
!
!----------------------------------------------------------------------------
SUBROUTINE cegterg_vdw( ndim, ndmx, nvec, nvecx, evc, &
                    ethr, overlap, e, btype, notcnv, iter )
  !----------------------------------------------------------------------------
  !
  ! ... iterative solution of the eigenvalue problem:
  !
  ! ... ( H - e S ) * evc = 0
  !
  ! ... where H is an hermitean operator, e is a real scalar,
  ! ... S is an overlap matrix, evc is a complex vector
  !
  USE io_global,  ONLY : stdout
  USE kinds,      ONLY : DP
  USE mp_global,  ONLY : intra_pool_comm
  USE mp,         ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  ! ... on INPUT
  !
  INTEGER :: ndim, ndmx, nvec, nvecx
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix evc, as declared in the calling pgm unit
    ! integer number of searched low-lying roots
    ! maximum dimension of the reduced basis set
    !    (the basis set is refreshed when its dimension would exceed nvecx)
  COMPLEX (kind=DP) :: evc(ndmx,nvec)
    !  evc   contains the  refined estimates of the eigenvectors
  REAL (kind=DP) :: ethr
    ! energy threshold for convergence
    !   root improvement is stopped, when two consecutive estimates of the root
    !   differ by less than ethr.
  LOGICAL :: overlap
    ! if .FALSE. : do not calculate S|psi>
  INTEGER, INTENT(in) :: btype(nvec)
    ! band type ( 1 = occupied, 0 = empty )
  !
  ! ... on OUTPUT
  !
  REAL(kind=DP) :: e(nvec)
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
  COMPLEX (kind=DP), ALLOCATABLE :: hc(:,:),  sc(:,:), vc(:,:)
    ! Hamiltonian on the reduced basis
    ! S matrix on the reduced basis
    ! the eigenvectors of the Hamiltonian
  COMPLEX(kind=DP), ALLOCATABLE :: psi(:,:), hpsi(:,:), spsi(:,:)
    ! work space, contains psi
    ! the product of H and psi
    ! the product of S and psi
  COMPLEX(kind=DP), EXTERNAL :: zdotc
    ! auxiliary complex variable
  REAL(kind=DP), ALLOCATABLE :: ew(:)
    ! eigenvalues of the reduced hamiltonian
  LOGICAL, ALLOCATABLE  :: conv(:)
    ! true if the root is converged
  REAL (kind=DP) :: empty_ethr
    ! threshold for empty bands
  !
  ! ... Called routines:
  !
  EXTERNAL  h_psi, s_psi, g_psi
    ! h_psi(ndmx,ndim,nvec,psi,hpsi)
    !     calculates H|psi>
    ! s_psi(ndmx,ndim,nvec,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (ndmx,nvec)
    ! g_psi(ndmx,ndim,notcnv,npol,psi,e)
    !    calculates (diag(h)-e)^-1 * psi, diagonal approx. to (h-e)^-1*psi
    !    the first nvec columns contain the trial eigenvectors
  !
  !
  CALL start_clock( 'cegterg' )
  !
  ! ... allocate the work arrays
  !
  ALLOCATE(  psi( ndmx, nvecx ) )
  ALLOCATE( hpsi( ndmx, nvecx ) )
  IF ( overlap ) ALLOCATE( spsi( ndmx, nvecx ) )
  ALLOCATE( sc( nvecx, nvecx ) )
  ALLOCATE( hc( nvecx, nvecx ) )
  ALLOCATE( vc( nvecx, nvecx ) )
  ALLOCATE( ew( nvecx ) )
  ALLOCATE( conv( nvec ) )
  !
  IF ( nvec > nvecx / 2 ) CALL errore( 'cegter', 'nvecx is too small', 1 )
  !
  ! ... threshold for empty bands
  !
  empty_ethr = max( ( ethr * 5.D0 ), 1.D-5 )
  !
  ! ... prepare the hamiltonian for the first iteration
  !
  notcnv = nvec
  nbase  = nvec
  conv   = .false.
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
  ! ... hc contains the projection of the hamiltonian onto the reduced space
  ! ... vc contains the eigenvectors of hc
  !
  hc(:,:) = ZERO
  sc(:,:) = ZERO
  vc(:,:) = ZERO
  !
  CALL zgemm( 'C', 'N', nbase, nbase, ndim, ONE, &
              psi, ndmx, hpsi, ndmx, ZERO, hc, nvecx )
  !
  CALL mp_sum( hc, intra_pool_comm )
  !
  IF ( overlap ) THEN
     !
     CALL zgemm( 'C', 'N', nbase, nbase, ndim, ONE, &
                 psi, ndmx, spsi, ndmx, ZERO, sc, nvecx )
     !
  ELSE
     !
     CALL zgemm( 'C', 'N', nbase, nbase, ndim, ONE, &
                 psi, ndmx,  psi, ndmx, ZERO, sc, nvecx )
     !
  ENDIF
  !
  CALL mp_sum( sc, intra_pool_comm )
  !
  FORALL( n = 1 : nbase )
     !
     e(n) = hc(n,n)
     !
     vc(n,n) = ONE
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
        IF ( .not. conv(n) ) THEN
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
        ENDIF
        !
     ENDDO
     !
     ! ... expand the basis set with new basis vectors ( H - e*S )|psi> ...
     !
     IF ( overlap ) THEN
        !
        CALL zgemm( 'N', 'N', ndim, notcnv, nbase, ONE, spsi, &
                    ndmx, vc, nvecx, ZERO, psi(1,nbase+1), ndmx )
        !
     ELSE
        !
        CALL zgemm( 'N', 'N', ndim, notcnv, nbase, ONE, psi, &
                    ndmx, vc, nvecx, ZERO, psi(1,nbase+1), ndmx )
        !
     ENDIF
     !
     FORALL( np = 1: notcnv ) &
        psi(:,nbase+np) = - ew(nbase+np) * psi(:,nbase+np)
     !
     CALL zgemm( 'N', 'N', ndim, notcnv, nbase, ONE, hpsi, &
                 ndmx, vc, nvecx, ONE, psi(1,nbase+1), ndmx )
     !
     CALL stop_clock( 'update' )
     !
     ! ... approximate inverse iteration
     !
     CALL g_psi( ndmx, ndim, notcnv, 1, psi(1,nbase+1), ew(nbase+1) )
     !
     ! ... "normalize" correction vectors psi(:,nbase+1:nbase+notcnv) in order
     ! ... to improve numerical stability of subspace diagonalization cdiaghg
     ! ... ew is used as work array : ew = <psi_i|psi_i>, i=nbase+1,nbase+notcnv
     !
     DO n = 1, notcnv
        !
        ew(n) = zdotc( ndim, psi(1,nbase+n), 1, psi (1,nbase+n), 1 )
        !
     ENDDO
     !
     CALL mp_sum( ew, intra_pool_comm )
     !
     FORALL( n = 1 : notcnv )
        !
        psi(:,nbase+n) = psi(:,nbase+n) / sqrt( ew(n) )
        !
     END FORALL
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
     CALL zgemm( 'C', 'N', nbase+notcnv, notcnv, ndim, ONE, psi, ndmx, &
                 hpsi(1,nbase+1), ndmx, ZERO, hc(1,nbase+1), nvecx )
     !
     CALL mp_sum( hc(:,nbase+1:nbase+notcnv), intra_pool_comm )
     !
     IF ( overlap ) THEN
        !
        CALL zgemm( 'C', 'N', nbase+notcnv, notcnv, ndim, ONE, psi, ndmx, &
                    spsi(1,nbase+1), ndmx, ZERO, sc(1,nbase+1), nvecx )
        !
     ELSE
        !
        CALL zgemm( 'C', 'N', nbase+notcnv, notcnv, ndim, ONE, psi, ndmx, &
                    psi(1,nbase+1), ndmx, ZERO, sc(1,nbase+1), nvecx )
        !
     ENDIF
     !
     CALL mp_sum( sc(:,nbase+1:nbase+notcnv), intra_pool_comm )
     !
     CALL stop_clock( 'overlap' )
     !
     nbase = nbase + notcnv
     !
     ! ifort 8.1 has some weird problem with forall
     !
     ! FORALL( n = 1 : nbase )
     DO n = 1, nbase
        !
        ! ... the diagonal of hc and sc must be strictly real
        !
        hc(n,n) = cmplx( REAL( hc(n,n) ), 0.D0 ,kind=DP)
        sc(n,n) = cmplx( REAL( sc(n,n) ), 0.D0 ,kind=DP)
        !
        ! FORALL(  m = n + 1 : nbase )
        DO m = n + 1, nbase
           !
           hc(m,n) = conjg( hc(n,m) )
           sc(m,n) = conjg( sc(n,m) )
           !
        ENDDO !  FORALL
        !
     ENDDO ! FORALL
     !
     ! ... diagonalize the reduced hamiltonian
     !
     CALL cdiaghg( nbase, nvec, hc, sc, nvecx, ew, vc )
     !
     ! ... test for convergence
     !
     WHERE( btype(:) == 1 )
        !
        conv(:) = ( ( abs( ew(1:nvec) - e(:) ) < ethr ) )
        !
     ELSEWHERE
        !
        conv(:) = ( ( abs( ew(1:nvec) - e(:) ) < empty_ethr ) )
        !
     END WHERE
     !
     notcnv = count( .not. conv(:) )
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
     IF ( notcnv == 0 .or. nbase+notcnv > nvecx .or. iter == maxter ) THEN
        !
        CALL start_clock( 'last' )
        !
        CALL zgemm( 'N', 'N', ndim, nvec, nbase, ONE, psi, &
                    ndmx, vc, nvecx, ZERO, evc, ndmx )
        !
        IF ( notcnv == 0 ) THEN
           !
           ! ... all roots converged: return
           !
           CALL stop_clock( 'last' )
           !
           exit iterate
           !
        ELSEIF ( iter == maxter ) THEN
           !
           ! ... last iteration, some roots not converged: return
           !
           WRITE( UNIT = stdout, &
                  FMT = '("   WARNING: ",i5," eigenvalues not converged")' ) &
                notcnv
           !
           CALL stop_clock( 'last' )
           !
           exit iterate
           !
        ENDIF
        !
        ! ... refresh psi, H*psi and S*psi
        !
        psi(:,1:nvec) = evc(:,1:nvec)
        !
        IF ( overlap ) THEN
           !
           CALL zgemm( 'N', 'N', ndim, nvec, nbase, ONE, spsi, &
                       ndmx, vc, nvecx, ZERO, psi(1,nvec+1), ndmx )
           !
           spsi(:,1:nvec) = psi(:,nvec+1:2*nvec)
           !
        ENDIF
        !
        CALL zgemm( 'N', 'N', ndim, nvec, nbase, ONE, hpsi, &
                    ndmx, vc, nvecx, ZERO, psi(1,nvec+1), ndmx )
        !
        hpsi(:,1:nvec) = psi(:,nvec+1:2*nvec)
        !
        ! ... refresh the reduced hamiltonian
        !
        nbase = nvec
        !
        hc(:,1:nbase) = ZERO
        sc(:,1:nbase) = ZERO
        vc(:,1:nbase) = ZERO
        !
        FORALL( n = 1 : nbase )
           !
           hc(n,n) = e(n)
           sc(n,n) = ONE
           vc(n,n) = ONE
           !
        END FORALL
        !
        CALL stop_clock( 'last' )
        !
     ENDIF
      !
  ENDDO iterate
  !
  DEALLOCATE( conv )
  DEALLOCATE( ew )
  DEALLOCATE( vc )
  DEALLOCATE( hc )
  DEALLOCATE( sc )
  IF ( overlap ) DEALLOCATE( spsi )
  DEALLOCATE( hpsi )
  DEALLOCATE( psi )
  !
  CALL stop_clock( 'cegterg' )
  !
  RETURN
  !
END SUBROUTINE cegterg_vdw
