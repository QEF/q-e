!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine regterg (ndim, ndmx, nvec, nvecx, evc, ethr, gstart, &
     e, notcnv, iter)
  !----------------------------------------------------------------------
  !
  !     iterative solution of the eigenvalue problem:
  !
  !     ( H - e S ) * evc = 0
  !
  !     where H is an hermitean operator, e is a real scalar,
  !     S is an overlap matrix, evc is a complex vector
  !     (real wavefunctions with only half plane waves stored)
  !
#include "machine.h"
  use parameters, only : DP
  use g_psi_mod
  implicit none
  ! on INPUT
  integer :: ndim, ndmx, nvec, nvecx, gstart
  ! dimension of the matrix to be diagonalized
  ! leading dimension of matrix evc, as declared in the calling pgm unit
  ! integer number of searched low-lying roots
  ! maximum dimension of the reduced basis set
  !    (the basis set is refreshed when its dimension would exceed nvecx)
  complex(kind=DP) :: evc (ndmx, nvec)
  real(kind=DP) :: ethr
  ! energy threshold for convergence: root improvement is stopped,
  ! when two consecutive estimates of the root differ by less than ethr.
  !
  ! on OUTPUT
  !  evc   contains the  refined estimates of the eigenvectors
  real(kind=DP) :: e (nvec)
  ! contains the estimated roots.

  integer :: iter, notcnv
  ! integer  number of iterations performed
  ! number of unconverged roots
  !
  ! LOCAL variables
  !
  integer, parameter :: maxter=20
  ! maximum number of iterations
  !
  integer :: kter, nbase, np, n, m
  ! counter on iterations
  ! dimension of the reduced basis
  ! counter on the reduced basis vectors
  ! do-loop counters
  real(kind=DP), allocatable :: hr (:,:),  sr (:,:), vr (:,:),  ew (:)
  ! Hamiltonian on the reduced basis
  ! S matrix on the reduced basis
  ! eigenvectors of the Hamiltonian
  ! eigenvalues of the reduced hamiltonian
  real(kind=DP), external :: DDOT
  complex(kind=DP), allocatable :: psi(:,:), hpsi (:,:),spsi (:,:)
  ! work space, contains psi
  ! the product of H and psi
  ! the product of S and psi
  logical, allocatable  :: conv (:)
  ! true if the root is converged
  !
  ! Called routines:
  external h_psi, g_psi
  ! h_psi(ndmx,ndim,nvec,psi,hpsi,spsi)
  !     calculates nvec H|psi> and S|psi> (if needed) products.
  !     Vectors psi,hpsi,spsi are dimensioned (ndmx,nvec)
  ! g_psi(ndmx,ndim,notcnv,psi,e)
  !    calculates (diag(h)-e)^-1 * psi, diagonal approx. to (h-e)^-1*psi
  !    the first nvec columns contain the trial eigenvectors
  !
  ! allocate the work arrays
  !

  test_new_preconditioning = .true.

  call start_clock ('cegterg')
  allocate ( psi( ndmx, nvecx))    
  allocate (hpsi( ndmx, nvecx))    
  allocate (spsi( ndmx, nvecx))    
  allocate (sr( nvecx, nvecx))    
  allocate (hr( nvecx, nvecx))    
  allocate (vr( nvecx, nvecx))    
  allocate (ew( nvecx))    
  allocate (conv( nvec))    

  if (nvec > nvecx / 2) call errore ('regter', 'nvecx is too small',1)
  !
  !     prepare the hamiltonian for the first iteration
  !
  notcnv = nvec
  nbase = nvec
  psi(:, 1:nvec) = evc(:, 1:nvec)
  !
  !     hpsi contains h times the basis vectors
  !
  call h_psi (ndmx, ndim, nvec, psi, hpsi, spsi)
  !
  !     hr contains the projection of the hamiltonian onto the reduced space
  !     vr contains the eigenvectors of hr
  !
  hr (:,:) = 0.d0
  vr (:,:) = 0.d0
  call DGEMM ('t', 'n', nbase, nbase, 2*ndim, 2.d0 , psi, &
       2*ndmx, hpsi, 2*ndmx, 0.d0, hr, nvecx)
  if (gstart.eq.2) call DGER (nbase, nbase, -1.d0, psi, 2*ndmx, &
       hpsi, 2*ndmx, hr, nvecx)
#ifdef __PARA
  call reduce (nbase * nvecx, hr)
#endif
  sr(:,:) = 0.d0
  call DGEMM ('t', 'n', nbase, nbase, 2*ndim, 2.d0, psi, &
       2*ndmx, spsi, 2*ndmx, 0.d0, sr, nvecx)
  if (gstart.eq.2) call DGER (nbase, nbase, -1.d0, psi, 2*ndmx, &
       spsi, 2*ndmx, sr, nvecx)
#ifdef __PARA
  call reduce (nbase * nvecx, sr)
#endif

  do n = 1, nbase
     e (n) = hr (n, n)
     conv (n) = .false.
     vr (n, n) = 1.d0
  enddo
  !
  !       iterate
  !
  do kter = 1, maxter
     iter = kter
     call start_clock ('update')
     np = 0
     do n = 1, nvec
        if ( .not.conv (n) ) then
           !
           !     this root not yet converged ... 
           !
           np = np + 1
           !
           ! reorder eigenvectors so that coefficients for unconverged
           ! roots come first. This allows to use quick matrix-matrix 
           ! multiplications to set a new basis vector (see below)
           !
           if (np .ne. n) vr(:,np) = vr(:,n)
           ! for use in g_psi
           ew(nbase+np) = e (n)
        endif
     enddo
     !
     !     expand the basis set with new basis vectors (h-es)psi ...
     !
     call DGEMM ('n', 'n', 2*ndim, notcnv, nbase, 1.d0, spsi, &
          2*ndmx, vr, nvecx, 0.d0, psi (1, nbase+1), 2*ndmx)
     do np = 1, notcnv
        psi (:,nbase+np) = - ew(nbase+np) * psi(:,nbase+np)
     end do
     call DGEMM ('n', 'n', 2*ndim, notcnv, nbase, 1.d0, hpsi, &
          2*ndmx, vr, nvecx, 1.d0, psi (1, nbase+1), 2*ndmx)

     call stop_clock ('update')
     !
     ! approximate inverse iteration
     !
     call g_psi (ndmx, ndim, notcnv, psi (1, nbase+1), ew(nbase+1) )
     !
     ! "normalize" correction vectors psi(*,nbase+1:nbase+notcnv) in order
     ! to improve numerical stability of subspace diagonalization rdiaghg
     ! ew is used as work array : ew = <psi_i|psi_i>, i=nbase+1,nbase+notcnv
     !
     do n = 1, notcnv
        ew (n) = 2.d0*DDOT(2*ndim, psi (1, nbase+n), 1, psi (1, nbase+n), 1)
        if (gstart.eq.2) ew (n) = ew(n) - psi (1, nbase+n)*psi (1, nbase+n)
     end do
#ifdef __PARA
     call reduce (notcnv, ew)
#endif
     do n = 1, notcnv
        call DSCAL (2 * ndim, 1.d0 / sqrt (ew (n) ), psi (1, nbase+n), 1)
     enddo
     !
     !   here compute the hpsi and spsi of the new functions
     !
     call h_psi (ndmx, ndim, notcnv, psi (1, nbase+1), &
                 hpsi (1, nbase+1), spsi (1, nbase+1) )
     !
     !     update the reduced hamiltonian
     !
     call start_clock ('overlap')
     !
     call DGEMM ('t', 'n', nbase+notcnv, notcnv, 2*ndim, 2.d0, &
          psi, 2*ndmx, hpsi (1, nbase+1), 2*ndmx, 0.d0, &
          hr (1, nbase+1) , nvecx)
     if (gstart.eq.2) call DGER (nbase+notcnv, notcnv, -1.d0, psi, 2*ndmx, &
          hpsi(1,nbase+1), 2*ndmx, hr (1, nbase+1), nvecx)
     !
#ifdef __PARA
     call reduce (nvecx * notcnv, hr (1, nbase+1) )
#endif
     call DGEMM ('t', 'n', nbase+notcnv, notcnv, 2*ndim, 2.d0, &
          psi, 2*ndmx, spsi (1, nbase+1) , 2*ndmx, 0.d0, &
          sr (1, nbase+1) , nvecx)
     if (gstart.eq.2) call DGER (nbase+notcnv, notcnv, -1.d0, psi, 2*ndmx, &
          spsi(1,nbase+1), 2*ndmx, sr (1, nbase+1), nvecx)
#ifdef __PARA
     call reduce (nvecx * notcnv, sr (1, nbase+1) )
#endif
     call stop_clock ('overlap')
     nbase = nbase+notcnv
     do n = 1, nbase
        do m = n + 1, nbase
           hr (m, n) = hr (n, m)
           sr (m, n) = sr (n, m)
        enddo
     enddo
     !
     !     diagonalize the reduced hamiltonian
     !
     call rdiaghg (nbase, nvec, hr, sr, nvecx, ew, vr)
     !
     !     test for convergence
     !
     notcnv = 0
     do n = 1, nvec
!!!        conv (n) = conv(n) .or. ( abs (ew (n) - e (n) ) <= ethr )
        conv (n) = ( abs (ew (n) - e (n) ) <= ethr )
        if ( .not. conv(n) ) notcnv = notcnv + 1
        e (n) = ew (n)
     enddo
     !
     !     if overall convergence has been achieved, OR
     !     the dimension of the reduced basis set is becoming too large, OR
     !     in any case if we are at the last iteration
     !     refresh the basis set. i.e. replace the first nvec elements
     !     with the current estimate of the eigenvectors;
     !     set the basis dimension to nvec.
     !
     if (notcnv == 0 .or. nbase+notcnv > nvecx .or. iter == maxter) then
        call start_clock ('last')
        call DGEMM ('n', 'n', 2*ndim, nvec, nbase, 1.d0, psi, &
             2*ndmx, vr, nvecx, 0.d0, evc, 2*ndmx)
        if (notcnv == 0) then
        !
        !     all roots converged: return
        !
           call stop_clock ('last')
           goto 10
        else if (iter == maxter) then
        !
        !     last iteration, some roots not converged: return
        !
#ifdef DEBUG_DAVIDSON
           do n = 1, nvec
              if ( .not.conv (n) ) write (6, '("   WARNING: e(",i3,") =",&
                   f10.5," is not converged to within ",1pe8.1)') n, e(n), ethr
           enddo
#else
           write (6, '("   WARNING: ",i5," eigenvalues not converged")') &
                notcnv
#endif
           call stop_clock ('last')
           goto 10
        end if
        !
        !     refresh psi, H*psi and S*psi
        !
        psi(:, 1:nvec) = evc(:, 1:nvec)

        call DGEMM ('n', 'n', 2*ndim, nvec, nbase, 1.d0, spsi, &
             2*ndmx, vr, nvecx, 0.d0, psi(1, nvec + 1), 2*ndmx)
        spsi(:, 1:nvec) = psi(:, nvec+1:2*nvec)

        call DGEMM ('n', 'n', 2*ndim, nvec, nbase, 1.d0, hpsi, &
             2*ndmx, vr, nvecx, 0.d0, psi(1, nvec + 1), 2*ndmx)
        hpsi(:, 1:nvec) = psi(:, nvec+1:2*nvec)
        !
        !     refresh the reduced hamiltonian
        !
        nbase = nvec
        hr (:, 1:nbase) = 0.d0
        sr (:, 1:nbase) = 0.d0
        vr (:, 1:nbase) = 0.d0
        do n = 1, nbase
           hr (n, n) = e(n)
           sr (n, n) = 1.d0
           vr (n, n) = 1.d0
        end do
        call stop_clock ('last')
     endif
  enddo

10 continue

  deallocate (conv)
  deallocate (ew)
  deallocate (vr)
  deallocate (hr)
  deallocate (sr)
  deallocate (spsi)
  deallocate (hpsi)
  deallocate ( psi)

  call stop_clock ('cegterg')
  return
end subroutine regterg

