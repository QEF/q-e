!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "machine.h"
!
!----------------------------------------------------------------------------
SUBROUTINE cdiisg( ndim, ndmx, nvec, nvecx, evc, e, ethr, &
                   btype, notcnv, diis_iter, iter )
  !----------------------------------------------------------------------------
  !
  ! ... iterative solution of the eigenvalue problem:
  !
  !     ( H - e S ) * evc = 0
  !
  ! ... where H is an hermitean operator, e is a real scalar,
  ! ... S is an overlap matrix, evc is a complex vector.
  ! ... The band-by-band RMM-DIIS method is used. 
  !
  USE io_global,  ONLY : stdout
  USE kinds, ONLY : DP
  USE g_psi_mod,  ONLY : h_diag, s_diag
  USE us,         ONLY : okvan
  !
  IMPLICIT NONE
  !
  ! ... on INPUT
  !
  INTEGER :: ndim, ndmx, nvec, nvecx, btype(nvec), iter
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix evc, as declared in the calling pgm unit
    ! integer number of searched low-lying roots
    ! maximum dimension of the reduced basis set
    !    (the basis set is refreshed when its dimension would exceed nvecx)
    ! band type (0=occupied, 1=empty)
    ! scf iteration
  REAL(KIND=DP) :: ethr
    ! energy threshold for convergence
    !   root improvement is stopped, when two consecutive estimates of the root
    !   differ by less than ethr.
  !
  ! ... on OUTPUT
  !
  COMPLEX(KIND=DP) :: evc(ndmx,nvec)
    !  evc   contains the  refined estimates of the eigenvectors
  REAL(KIND=DP) :: e(nvec), diis_iter
    ! contains the estimated roots.
    ! average number of iterations performed per band
  INTEGER :: notcnv
    ! number of unconverged roots
  !
  ! ... LOCAL variables
  !
  INTEGER, PARAMETER :: maxter = 20
    ! maximum number of iterations
  INTEGER :: kter, nbase, ib, n, m , np
    ! counter on iterations
    ! lower extreme for iteration loop
    ! dimension of the reduced basis
    ! counter on the reduced basis vectors
    ! do-loop counters
  COMPLEX(KIND=DP), ALLOCATABLE :: rc(:,:), hc(:,:), sc(:,:), vc(:), vcn(:,:)
    ! <res_i|res_j> matrix
    ! H matrix on the reduced basis
    ! S matrix on the reduced basis
    ! the eigenvectors of the Hamiltonian
    ! workspace
  COMPLEX(KIND=DP), ALLOCATABLE :: psi(:,:), hpsi(:,:), spsi(:,:), res(:,:)
    ! work space, contains psi
    ! the product of H and psi
    ! the product of S and psi
    ! residual vector
  COMPLEX(KIND=DP) :: hevc(ndmx,nvec), sevc(ndmx,nvec)
    ! the product of H and the best estimate of the eigenvectors evc
    ! the product of S and the best estimate of the eigenvectors evc
  REAL(KIND=DP), ALLOCATABLE :: ew(:)
    ! eigenvalues of the reduced hamiltonian
  REAL(KIND=DP) :: ec, snorm, snorm0, lam, ew0, denm, x
    ! dummy variable
    ! squared norm of current residual
    ! squared norm of initial residual calculated with the old evc
    ! initial eigenvalue from previous iteration
    ! variables for teter preconditioning (not used)
  LOGICAL :: verb
    ! controlling verbosity of printout
  COMPLEX(KIND=DP), EXTERNAL :: ZDOTC
  !
  !
  CALL start_clock( 'diis' )
  !
  verb = .FALSE.
  !
  ! ... allocate the work arrays
  !
  ALLOCATE(  psi( ndmx, nvecx ) )
  ALLOCATE( hpsi( ndmx, nvecx ) )
  ALLOCATE( spsi( ndmx, nvecx ) )
  ALLOCATE(  res( ndmx, nvecx ) )
  ALLOCATE(   rc( nvecx, nvecx ) )
  ALLOCATE(   hc( nvecx, nvecx ) )
  ALLOCATE(   sc( nvecx, nvecx ) )
  ALLOCATE(   vc( nvecx) )
  ALLOCATE(  vcn( nvecx, nvecx ) )
  ALLOCATE(   ew( nvecx ) )
  !
  ! ... rotate
  !
  CALL rotate_wfc( ndmx, ndim, nvec, nvec, evc, okvan, evc, e )
  !
  notcnv = nvec
  lam    = 1.D0
  !
  ! ... Loop over bands
  !
  bands_loop: DO ib = 1, nvec
     !
     ! ... prepare the hamiltonian for the first iteration
     !
     nbase = 1     
     !
     ew(:)     = 0.D0
     vc(:)     = ( 0.D0, 0.D0 )
     vcn(:, :) = ( 0.D0, 0.D0 )
     psi(:, :) = ( 0.D0, 0.D0 )
     res(:, :) = ( 0.D0, 0.D0 )
     hc(:, :)  = ( 0.D0, 0.D0 )
     sc(:, :)  = ( 0.D0, 0.D0 )
     rc(:, :)  = ( 0.D0, 0.D0 )
     !
     ew0    = 1.D10
     snorm0 = 1.D-20
     !
     ! ... |psi_1> is the best approximated eigenvector
     !
     psi(:,nbase) = evc(:,ib)
     vc(nbase)    = ( 1.D0, 0.D0 )
     !
     ! ... calculate hpsi=H|psi_1> and spsi=S|psi_1> 
     !
     CALL h_1psi( ndmx, ndim, psi(1,nbase), hpsi(1,nbase), spsi(1,nbase) )
     !  
     ! ... calculate the first element of the reduced hamiltonian 
     ! ... and overlap matrices
     ! ... hc(1,1)=<psi_1|H|psi_1>    sc(1,1)=<psi_1|S|psi_1>
     !
     hc(1,1) = ZDOTC( ndim, psi(1,1), 1, hpsi(1,1), 1 )
     sc(1,1) = ZDOTC( ndim, psi(1,1), 1, spsi(1,1), 1 )
#ifdef __PARA
     CALL reduce( 2, hc(1,1) )
     CALL reduce( 2, sc(1,1) )
#endif
     !
     ! ... calculate the residual vector |R>=H|psi> - e S|psi>
     !
     CALL ZGEMM( 'N', 'N', ndim, 1, nbase, ( 1.D0, 0.D0 ), spsi, &
                 ndmx, vc, nvecx, ( 0.D0, 0.D0 ), res(1,nbase), ndmx )
     CALL ZGEMM( 'N', 'N', ndim, 1, nbase, ( 1.D0, 0.D0 ), hpsi, &
                 ndmx, vc, nvecx, DCMPLX( -e(ib), 0.D0 ), res(1,nbase), ndmx )
     !
     ! ... calculate the first element of the <res_i|res_j> matrix
     !
     rc(1,1) = ZDOTC( ndim, res(1,1), 1, res(1,1), 1 )
#ifdef __PARA
     CALL reduce( 2, rc(1,1) )
#endif
     !
     ! ... iterate
     !
     iterate: DO kter = 1, maxter
        !
        ! ... preconditionate the  residual vector |P_n>= K*|R_n>
        ! ... and add it to the basis => |psi_n+1>
        !
        CALL ZGEMM( 'N', 'N', ndim, 1, nbase, ( 1.D0, 0.D0 ), res, &
                    ndmx, vc, nvecx, ( 0.D0, 0.D0 ), psi(1,nbase+1), ndmx )
        !
        CALL g_psi( ndmx, ndim, 1, psi(1,nbase+1), e(ib) )
        ! 
        ! ... add a new vector to the basis: kresse method
        ! ... |psi_n+1> = |psi_n> + \lambda * K * |res_n>
        !
!        CALL ZGEMM( 'N', 'N', ndim, 1, nbase, ( 1.D0, 0.D0 ), psi, &
!                    ndmx, vc , nvecx, DCMPLX( lam, 0.D0 ) , &
!                    psi(1,nbase+1), ndmx )
        !
        nbase = nbase + 1
        !
        ! ... normalize new basis vector
        !
        ec = REAL( ZDOTC( ndim, psi(1,nbase), 1, psi(1,nbase), 1 ) )
#ifdef __PARA
        CALL reduce( 1, ec )
#endif
        !
        CALL ZSCAL( ndim, DCMPLX( 1/SQRT( ec ), 0.D0 ), psi(1,nbase), 1 )
        !
        ! ... new eigenvector, normalize eigenvectors
        !
        vc(nbase) = ( 1.D0, 0.D0 )
        ec        = REAL( ZDOTC( nbase, vc , 1, vc , 1 ) )
        !
        CALL ZSCAL( nvecx, DCMPLX( 1/SQRT( ec ), 0.D0 ), vc, 1 )
        !
        ! ... calculate hpsi=H|psi> and spsi=S|psi> 
        !
        CALL h_1psi( ndmx, ndim, psi(1,nbase), hpsi(1,nbase), spsi(1,nbase) )
        !
        ! ... orthogonalize
        !
        CALL cgramg1( ndmx, nvecx, ndim, 1, nbase, psi, spsi, hpsi )
        !
        ! ... calculate the new elements of the reduced hamiltonian
        ! ... and overlap matrices
        ! ... hc(i,j) =<psi_i|H|psi_j>  and sc(i,j)=<psi_i|S|psi_j>
        !
        CALL ZGEMM( 'C', 'N', nbase, 1, ndim, ( 1.D0, 0.D0 ), psi, &
                    ndmx, hpsi(1,nbase), ndmx, ( 0.D0, 0.D0 ), hc(1,nbase), &
                    nvecx )
        CALL ZGEMM( 'C', 'N', nbase, 1, ndim, ( 1.D0, 0.D0 ), psi, &
                    ndmx, spsi(1,nbase), ndmx, ( 0.D0, 0.D0 ), sc(1,nbase), &
                    nvecx )
#ifdef __PARA
        CALL reduce( 2 * nvecx, hc(1,nbase) )
        CALL reduce( 2 * nvecx, sc(1,nbase) )
#endif
        !
        ! ... calculate the residual vector |R>=H|psi> - e S|psi>
        !
        CALL ZGEMM( 'N', 'N', ndim, 1, nbase, ( 1.D0 , 0.D0 ), spsi, &
                    ndmx, vc, nvecx, ( 0.D0, 0.D0 ), res(1,nbase), ndmx )
        CALL ZGEMM( 'N', 'N', ndim, 1, nbase, ( 1.D0, 0.D0 ), hpsi, &
                    ndmx, vc, nvecx, DCMPLX( -e(ib), 0.D0 ), res(1,nbase), &
                    ndmx )
        !
        ! ... calculate the new elements of the <res_i|res_j> matrix
        !
        CALL ZGEMM( 'C', 'N', nbase, 1, ndim, ( 1.D0, 0.D0 ), res, ndmx, &
                    res(1,nbase), ndmx, ( 0.D0, 0.D0 ), rc(1,nbase), nvecx )
        !
        ew(nbase) = rc(nbase,nbase)
        !
#ifdef __PARA
        CALL reduce( 2 * nvecx, rc(1,nbase) )
#endif
        !
        ! ... rc, hc, and sc are hermitian
        !
        DO n = 1, nbase
           !
           ! ... the diagonal of rc and sc must be strictly real 
           !
           rc(n,n) = DCMPLX( REAL( rc(n,n) ), 0.D0 )
           hc(n,n) = DCMPLX( REAL( hc(n,n) ), 0.D0 )
           sc(n,n) = DCMPLX( REAL( sc(n,n) ), 0.D0 )
           !
           DO m = ( n + 1 ), nbase
              rc(m,n) = CONJG( rc(n,m) )
              hc(m,n) = CONJG( hc(n,m) )
              sc(m,n) = CONJG( sc(n,m) )
           END DO
        END DO
        IF ( verb ) THEN
           WRITE( stdout, * ) 'overlap' 
           WRITE( stdout, * ) ( (m,n,sc(n,m), n = 1, nbase), m = 1, nbase )
           WRITE( stdout, * )  
           WRITE( stdout, * ) 'rc' 
           WRITE( stdout, * ) ( (m,n,rc(n,m), n = 1, nbase ), m = 1, nbase )
           WRITE( stdout, * )  
           WRITE( stdout, * ) 'eigval'
        END IF
        !
        ! ... diagonalize the reduced hamiltonian
        !
        CALL cdiaghg( nbase, 1, rc, sc, nvecx, ew, vcn )
        CALL ZCOPY( nvecx, vcn(1,1), 1, vc, 1 )
        !
        IF ( verb ) THEN
           DO n = 1, nbase
              WRITE( stdout, * ) n,ew(n)
           END DO
           WRITE( stdout, * )  
           DO n = 1, nbase
           END DO
           WRITE( stdout, * ) 'eigvec' 
           DO n = 1, nbase
              WRITE( stdout, * ) n, vc(n)
           END DO
        END IF
        !
        ! ... squared norm of current residual
        !
        snorm = ew(1)
        !
        ! ... calculate new eigenvalues
        !
        vcn(:,:)= (0.d0,0.d0)
        CALL ZGEMM( 'N', 'N', nbase, 1, nbase, ( 1.D0 , 0.D0 ), hc, &
                    nvecx, vc, nvecx, ( 0.D0, 0.D0 ), vcn(1,1), nvecx )
        ec = REAL( ZDOTC( nvecx, vc, 1, vcn(1,1), 1 ) )
        CALL ZGEMM( 'N', 'N', nbase, 1, nbase, ( 1.D0 , 0.D0 ), sc, &
                    nvecx, vc, nvecx, ( 0.D0, 0.D0 ), vcn(1,1), nvecx )
        ec = ec / REAL( ZDOTC( nvecx, vc, 1, vcn(1,1), 1 ) )
        !
        IF ( verb ) WRITE( stdout, * ) 'NORM RES=',snorm,'DELTA EIG=',ec-e(ib)
        !
        ! ... Convergence?
        ! ... Non occupied levels are converged with a lower precision than 
        ! ... occupied ones.
        !
        IF ( btype(ib) == 0 ) THEN 
           IF ( ( snorm < snorm0*0.3D0 .OR. ABS( ec - e(ib) ) < ethr ) .AND. &
                ec <= ew0 ) EXIT iterate
        ELSE
           if ( ( snorm < snorm0*0.3D0 .OR. &
                  ABS( ec - e(ib) ) < MAX( ethr*50.D0, 1.D-4 ) ) .AND. &
                ec <= ew0 ) EXIT iterate
        END IF
        !
        IF ( kter == 1 ) THEN
           snorm0 = snorm
           ew0    = ec
        END IF
        e(ib) = ec
        !
        ! ... Size of reduced basis is exceeded: refresh
        !
        IF ( kter == 1 ) THEN
           snorm0 = snorm
           ew0    = ec
        END IF
        e(ib) = ec
        !
        ! ... Size of reduced basis is exceeded: refresh
        !
        IF ( nbase >= nvecx ) THEN
           !
           CALL ZGEMM( 'N', 'N', ndim, 1, nbase, ( 1.D0, 0.D0 ), psi, &
                       ndmx, vc, nvecx, ( 0.D0, 0.D0 ), evc(1,ib), ndmx )
           !
           IF ( verb ) WRITE( stdout, * ) 'rotate band ',ib
           !
           nbase = 1     
           !
           ew(:)     = 0.D0
           vc(:)     = ( 0.D0, 0.D0 )
           vcn(:, :) = ( 0.D0, 0.D0 )
           psi(:, :) = ( 0.D0, 0.D0 )
           res(:, :) = ( 0.D0, 0.D0 )
           hc(:, :)  = ( 0.D0, 0.D0 )
           sc(:, :)  = ( 0.D0, 0.D0 )
           rc(:, :)  = ( 0.D0, 0.D0 )
           !
           ! ... |psi_1> is the best approximated eigenvector
           !
           psi(:,nbase) = evc(:,ib)
           vc(nbase)    = ( 1.D0, 0.D0 )
           !
           ! ... calculate hpsi=H|psi_1> and spsi=S|psi_1> 
           !
           CALL h_1psi( ndmx, ndim, psi(1,nbase), hpsi(1,nbase), spsi(1,nbase) )
           !  
           ! ... calculate the first element of the reduced hamiltonian 
           ! ... and overlap matrices
           ! ... hc(1,1)=<psi_1|H|psi_1>    sc(1,1)=<psi_1|S|psi_1>
           !
           hc(1,1) = ZDOTC( ndim, psi(1,1), 1, hpsi(1,1), 1 )
           sc(1,1) = ZDOTC( ndim, psi(1,1), 1, spsi(1,1), 1 )
#ifdef __PARA
           CALL reduce( 2, hc(1,1) )
           CALL reduce( 2, sc(1,1) )
#endif
           !
           ! ... calculate the residual vector |R>=H|psi> - e S|psi>
           !
           CALL ZGEMM( 'N', 'N', ndim, 1, nbase, ( 1.D0, 0.D0 ), spsi, &
                       ndmx, vc, nvecx, ( 0.D0, 0.D0 ), res(1,nbase), ndmx )
           CALL ZGEMM( 'N', 'N', ndim, 1, nbase, ( 1.D0, 0.D0 ), hpsi, &
                        ndmx, vc, nvecx, DCMPLX( -e(ib), 0.D0 ), res(1,nbase), &
                        ndmx )
           !
           ! ... calculate the first element of the <res_i|res_j> matrix
           !
           rc(1,1) = ZDOTC( ndim, res(1,1), 1, res(1,1), 1 )
           !
#ifdef __PARA
           CALL reduce( 2, rc(1,1) )
#endif           
           !
        END IF
        !
     END DO iterate
     !
     e(ib) = ec
     diis_iter = diis_iter + kter
     !
     IF ( kter > maxter ) THEN
        WRITE( stdout, '("   WARNING: eigenvalue ",i5," not converged")') ib
     ELSE
        notcnv = notcnv - 1
        IF ( verb ) then
           WRITE( stdout, * ) 'BAND ',ib, ' CONVERGED'
           WRITE( stdout, * )
        END IF
     END IF
     !
     ! ... calculate best approximated wavefunction and corresponding 
     ! ... hpsi and spsi
     !
     CALL ZGEMM( 'N', 'N', ndim, 1, nbase, ( 1.D0, 0.D0 ), psi, &
                 ndmx, vc, nvecx, ( 0.D0, 0.D0 ), evc(1,ib), ndmx )
     CALL ZGEMM( 'N', 'N', ndim, 1, nbase, ( 1.D0, 0.D0 ), hpsi, &
                 ndmx, vc, nvecx, ( 0.D0, 0.D0 ), hevc(1,ib), ndmx )
     CALL ZGEMM( 'N', 'N', ndim, 1, nbase, ( 1.D0, 0.D0 ), spsi, &
                 ndmx, vc, nvecx, ( 0.D0, 0.D0 ), sevc(1,ib), ndmx )

  END DO bands_loop
  !
  diis_iter = diis_iter / nvec
  !
  ! ... deallocation of the work arrays
  !
  DEALLOCATE( psi )
  DEALLOCATE( hpsi )
  DEALLOCATE( spsi )
  DEALLOCATE( res )
  DEALLOCATE( rc )
  DEALLOCATE( hc )
  DEALLOCATE( sc )
  DEALLOCATE( vc )
  DEALLOCATE( vcn )
  DEALLOCATE( ew )
  !
  ! ... orthonormalize bands
  !
!  IF ( MOD( iter, 6 ) == 0 ) &
  CALL cgramg1( ndmx, nvec, ndim, 1, nvec, evc, sevc, hevc )
  !
  CALL stop_clock( 'diis' )
  !
  RETURN
  !
END SUBROUTINE cdiisg
