!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine cdiisg (ndim, ndmx, nvec, nvecx, evc, e, ethr, &
      btype, notcnv, diis_iter, iter)
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  !
  !     iterative solution of the eigenvalue problem:
  !
  !     ( H - e S ) * evc = 0
  !
  !     where H is an hermitean operator, e is a real scalar,
  !     S is an overlap matrix, evc is a complex vector.
  !     The band-by-band RMM-DIIS method is used. 
#include "machine.h"
  use parameters, only : DP
  use g_psi_mod
  use pwcom, only : nelec, lgauss, ltetra, okvan

  implicit none
  ! on INPUT
  integer :: ndim, ndmx, nvec, nvecx, btype(nvec), iter
  ! dimension of the matrix to be diagonalized
  ! leading dimension of matrix evc, as declared in the calling pgm unit
  ! integer number of searched low-lying roots
  ! maximum dimension of the reduced basis set
  !    (the basis set is refreshed when its dimension would exceed nvecx)
  ! band type (0=occupied, 1=empty)
  ! scf iteration
  real(kind=DP) :: ethr
  ! energy threshold for convergence
  !   root improvement is stopped, when two consecutive estimates of the root
  !   differ by less than ethr.
  ! on OUTPUT
  complex(kind=DP) :: evc (ndmx, nvec)
  !  evc   contains the  refined estimates of the eigenvectors
  real(kind=DP) :: e (nvec), diis_iter
  ! contains the estimated roots.
  ! average number of iterations performed per band
  integer :: notcnv
  ! number of unconverged roots

  ! LOCAL variables
  !
  integer, parameter :: maxter=20
  ! maximum number of iterations
  !
  integer :: kter, minter, nbase, ib, n, m , np
  ! counter on iterations
  ! lower extreme for iteration loop
  ! dimension of the reduced basis
  ! counter on the reduced basis vectors
  ! do-loop counters
  complex(kind=DP), allocatable :: rc (:,:),  hc (:,:), sc (:,:), &
       vc (:), vcn(:,:)
  ! <res_i|res_j> matrix
  ! H matrix on the reduced basis
  ! S matrix on the reduced basis
  ! the eigenvectors of the Hamiltonian
  ! workspace
  complex(kind=DP), allocatable :: psi(:,:), hpsi (:,:), spsi (:,:), res(:,:)
  ! work space, contains psi
  ! the product of H and psi
  ! the product of S and psi
  ! residual vector
  complex(kind=DP) :: hevc(ndmx, nvec), sevc(ndmx, nvec)
  ! the product of H and the best estimate of the eigenvectors evc
  ! the product of S and the best estimate of the eigenvectors evc
  real(kind=DP), allocatable :: ew (:)
  ! eigenvalues of the reduced hamiltonian
  real(kind=DP) :: ec, snorm, snorm0, lam, ew0, denm, x
  ! dummy variable
  ! squared norm of current residual
  ! squared norm of initial residual calculated with the old evc
  ! initial eigenvalue from previous iteration
  ! variables for teter preconditioning (not used)
  logical :: verb
  ! controlling verbosity of printout
  complex(kind=DP), external ::  ZDOTC
  external h_1psi, cdiagh

  call start_clock ('diis')

  verb = .false.

  !
  ! allocate the work arrays
  !
  allocate( psi (ndmx, nvecx))
  allocate(hpsi (ndmx, nvecx))
  allocate(spsi (ndmx, nvecx))
  allocate( res (ndmx, nvecx))
  allocate( rc (nvecx, nvecx))
  allocate( hc (nvecx, nvecx))
  allocate( sc (nvecx, nvecx))
  allocate( vc (nvecx))
  allocate( vcn(nvecx, nvecx))
  allocate( ew (nvecx))

  !
  ! rotate
  !
  call rotate_wfc (ndmx, ndim, nvec, nvec, evc, okvan, evc, e)

  notcnv = nvec
  lam = 1.d0
  minter = 1

  ! Loop over bands
  do ib = 1, nvec

10   continue
     !
     !     prepare the hamiltonian for the first iteration
     !
     nbase = 1     
     ew  (:) = 0.d0
     vc  (:) = (0.d0, 0.d0)
     vcn (:, :) = (0.d0, 0.d0)
     psi (:, :) = (0.d0, 0.d0)
     res (:, :) = (0.d0, 0.d0)
     hc  (:, :) = (0.d0, 0.d0)
     sc  (:, :) = (0.d0, 0.d0)
     rc  (:, :) = (0.d0, 0.d0)
     if (minter.eq.1) then
        ew0 = 1.d10
        snorm0 = 1.d-20
     endif
     ! |psi_1> is the best approximated eigenvector
     psi(:,nbase) = evc(:, ib)
     vc(nbase) = (1.d0, 0.d0)
     !
     !     calculate hpsi=H|psi_1> and spsi=S|psi_1> 
     !
     call h_1psi (ndmx, ndim, psi(1,nbase), hpsi(1,nbase), spsi(1,nbase))
     !  
     !     calculate the first element of the reduced hamiltonian 
     !     and overlap matrices
     !     hc(1,1)=<psi_1|H|psi_1>    sc(1,1)=<psi_1|S|psi_1>
     !
     hc (1, 1) = ZDOTC (ndim, psi (1, 1), 1, hpsi (1, 1), 1)
     sc (1, 1) = ZDOTC (ndim, psi (1, 1), 1, spsi (1, 1), 1)
     !
     !   calculate the residual vector |R>=H|psi> - e S|psi>
     !
     call ZGEMM ('n', 'n', ndim, 1, nbase, (1.d0, 0.d0), spsi, &
          ndmx, vc , nvecx, (0.d0, 0.d0), res (1, nbase), ndmx)
     call ZGEMM ('n', 'n', ndim, 1, nbase, (1.d0, 0d0), hpsi, &
          ndmx, vc , nvecx, DCMPLX(- e(ib), 0.d0), res (1, nbase), ndmx)
     !
     !   calculate the first element of the <res_i|res_j> matrix
     rc (1, 1) = ZDOTC (ndim, res (1, 1), 1, res (1, 1), 1)
     !
     !  iterate
     !
     do kter = minter, maxter
     !
     !   preconditionate the  residual vector |P_n>= K*|R_n>
     !   and add it to the basis => |psi_n+1>
     !
        call ZGEMM ('n', 'n', ndim, 1, nbase, (1.d0, 0.d0) , res, &
             ndmx, vc , nvecx, (0.d0, 0.d0) , psi (1, nbase+1), ndmx)

        call g_psi( ndmx, ndim, 1, psi (1, nbase+1), e(ib) )

     ! 
     ! add a new vector to the basis: kresse method
     ! |psi_n+1> = |psi_n> + \lambda * K * |res_n>
!        call ZGEMM ('n', 'n', ndim, 1, nbase, (1.d0, 0.d0), psi, &
!             ndmx, vc , nvecx, DCMPLX (lam, 0.d0) , &
!             psi (1, nbase+1), ndmx)

        nbase = nbase + 1
     !
     ! normalize new basis vector
        ec = ZDOTC (ndim, psi (1, nbase), 1, psi (1, nbase), 1)
        call ZSCAL (ndim, DCMPLX (1/dsqrt(ec), 0.d0), psi (1, nbase), 1)
     ! new eigenvector, normalize eigenvectors
        vc(nbase) = (1.d0, 0.d0)
        ec = DREAL(ZDOTC (nbase, vc , 1, vc , 1))
        call ZSCAL (nvecx, DCMPLX (1/dsqrt(ec), 0.d0), vc, 1)
     !
     !     calculate hpsi=H|psi> and spsi=S|psi> 
     !
        call h_1psi (ndmx, ndim, psi(1,nbase), hpsi(1,nbase), spsi(1,nbase))
     !
     !     orthogonalize
     !
        call cgramg1 (ndmx, nvecx, ndim, 1, nbase, psi, spsi, hpsi)
     !
     !     calculate the new elements of the reduced hamiltonian
     !     and overlap matrices
     !     hc(i,j) =<psi_i|H|psi_j>  and sc(i,j)=<psi_i|S|psi_j>
     !
        call ZGEMM ('c', 'n', nbase, 1, ndim, (1.d0, 0.d0) , psi, &
             ndmx, hpsi (1, nbase) , ndmx, (0.d0, 0.d0) , hc (1, nbase), &
             nvecx)
        call ZGEMM ('c', 'n', nbase, 1, ndim, (1.d0, 0.d0) , psi, &
             ndmx, spsi (1, nbase) , ndmx, (0.d0, 0.d0) , sc (1, nbase), &
             nvecx)
     !
     !   calculate the residual vector |R>=H|psi> - e S|psi>
     !
        call ZGEMM ('n', 'n', ndim, 1, nbase, (1.d0 , 0.d0), spsi, &
             ndmx, vc , nvecx, (0.d0, 0.d0), res (1, nbase), ndmx)
        call ZGEMM ('n', 'n', ndim, 1, nbase, (1.d0, 0d0), hpsi, &
             ndmx, vc , nvecx, DCMPLX (-e(ib), 0.d0), res (1, nbase), ndmx)
     !
     !   calculate the new elements of the <res_i|res_j> matrix
        call ZGEMM ('c', 'n', nbase, 1, ndim, (1.d0, 0.d0) , &
             res , ndmx, res (1, nbase) , ndmx, (0.d0, 0.d0) , &
             rc (1, nbase) , nvecx)
        ew(nbase) = rc (nbase, nbase)
     !
     !     rc, hc, and sc are hermitian
     !
     do n = 1, nbase
        !  the diagonal of rc and sc must be strictly real 
        rc (n, n) = DCMPLX (DREAL (rc (n, n) ), 0.d0)
        hc (n, n) = DCMPLX (DREAL (hc (n, n) ), 0.d0)
        sc (n, n) = DCMPLX (DREAL (sc (n, n) ), 0.d0)
        do m = n + 1, nbase
           rc (m, n) = CONJG (rc (n, m) )
           hc (m, n) = CONJG (hc (n, m) )
           sc (m, n) = CONJG (sc (n, m) )
        enddo
     enddo
     if (verb) then
        write(6,*) 'overlap' 
        write(6,*) ((m,n,sc(n,m), n=1,nbase), m=1,nbase)
        write(6,*)  
        write(6,*) 'rc' 
        write(6,*) ((m,n,rc(n,m), n=1,nbase), m=1,nbase)
        write(6,*)  
        write(6,*) 'eigval'
     endif
     
     !
     !     diagonalize the reduced hamiltonian
     !
     call cdiaghg (nbase, 1, rc, sc, nvecx, ew, vcn )
     call ZCOPY (nvecx, vcn (1, 1), 1, vc, 1 )

     if (verb) then
        do n=1, nbase
           write(6,*) n,ew(n)
        enddo
        write(6,*)  
        do n=1,nbase
        enddo
        write(6,*) 'eigvec' 
        do n=1, nbase
           write(6,*) n, vc(n)
        enddo
     endif

     ! squared norm of current residual
     snorm = ew(1)
     !
     ! calculate new eigenvalues
     vcn(:,:)= (0.d0,0.d0)
     call ZGEMM ('n', 'n', nbase, 1, nbase, (1.d0 , 0.d0), hc, &
          nvecx, vc , nvecx, (0.d0, 0.d0), vcn (1, 1), nvecx)
     ec = DREAL( ZDOTC (nvecx, vc, 1, vcn (1, 1), 1) )
     call ZGEMM ('n', 'n', nbase, 1, nbase, (1.d0 , 0.d0), sc, &
          nvecx, vc , nvecx, (0.d0, 0.d0), vcn (1, 1), nvecx)
     ec = ec / DREAL( ZDOTC (nvecx, vc, 1, vcn (1, 1), 1) )

     if (verb) write(6,*) 'NORM RES=',snorm,'DELTA EIG=',ec-e(ib)

     !
     ! Convergence?
     ! Non occupied levels are converged with a lower precision than 
     ! occupied ones.
     !
     if (btype(ib) .eq. 0) then 
        if ( (snorm.lt.snorm0*0.3 .or. abs(ec-e(ib)).lt.ethr) &
             .and. ec.le.ew0) goto 20
     else
        if ( (snorm.lt.snorm0*0.3 .or. abs(ec-e(ib)).lt.max(ethr*50.0,1.D-4)) &
             .and. ec.le.ew0) goto 20
     endif
        
     if (minter.eq.1) then
        snorm0 = snorm
        ew0 = ec
     endif
     e(ib) = ec
     
     !
     ! Size of reduced basis is exceeded: refresh
     !
     if (nbase.ge.nvecx) then
        call ZGEMM ('n', 'n', ndim, 1, nbase, (1.d0, 0d0), psi, &
             ndmx, vc , nvecx, (0.d0, 0.d0), evc (1, ib), ndmx)
        if (verb) write(6,*) 'rotate band ',ib
        minter = kter + 1
        goto 10
     endif

     enddo ! iterate
20   continue

     e(ib) = ec
     diis_iter = diis_iter + kter

     if (kter .gt. maxter) then
        write (6, '("   WARNING: eigenvalue ",i5," not converged")') &
             ib
     else
        notcnv = notcnv - 1
        minter = 1
        if (verb) then
           write(6,*) 'BAND ',ib, ' CONVERGED'
           write(6,*)
        endif
     endif
     !
     ! calculate best approximated wavefunction and corresponding hpsi and spsi
     !
     call ZGEMM ('n', 'n', ndim, 1, nbase, (1.d0, 0d0), psi, &
          ndmx, vc , nvecx, (0.d0, 0.d0), evc (1, ib), ndmx)
     call ZGEMM ('n', 'n', ndim, 1, nbase, (1.d0, 0d0), hpsi, &
          ndmx, vc , nvecx, (0.d0, 0.d0), hevc (1, ib), ndmx)
     call ZGEMM ('n', 'n', ndim, 1, nbase, (1.d0, 0d0), spsi, &
          ndmx, vc , nvecx, (0.d0, 0.d0), sevc (1, ib), ndmx)

  enddo ! loop over bands     

  diis_iter = diis_iter / nvec

  deallocate( psi )
  deallocate(hpsi )
  deallocate(spsi )
  deallocate( res )
  deallocate( rc )
  deallocate( hc )
  deallocate( sc )
  deallocate( vc )
  deallocate(vcn )
  deallocate( ew )

  !
  ! orthonormalize bands
  !
!  if (mod(iter,6).eq.0) &
       call cgramg1 (ndmx, nvec, ndim, 1, nvec, evc, sevc, hevc)


  call stop_clock ('diis')
  return
end subroutine cdiisg


