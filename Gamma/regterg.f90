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
  !     ( H - e S ) * psi = 0
  !
  !     where H is an hermitean operator, e is a real scalar,
  !     S is an overlap matrix, psi is a complex vector
  !     (real wavefunctions with only half plane waves stored)
  !
#include "machine.h"
  use parameters, only : DP
  implicit none
  ! on INPUT
  integer :: ndim, ndmx, nvec, nvecx, gstart
  ! dimension of the matrix to be diagonalized
  ! leading dimension of matrix psi, as declared in the calling pgm unit
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
  real(kind=DP) :: DDOT
  complex(kind=DP), allocatable :: psi(:,:), hpsi (:,:),spsi (:,:)
  ! work space, contains psi
  ! the product of H and psi
  ! the product of S and psi
  integer, allocatable  :: conv (:)
  ! if 1 the root is converged if 0 no
  !
  ! Called routines:
  external h_psi, g_psi, DDOT
  ! h_psi(ndmx,ndim,nvec,psi,hpsi,spsi)
  !     calculates nvec H|psi> and S|psi> (if needed) products.
  !     Vectors psi,hpsi,spsi are dimensioned (ndmx,nvec)
  ! g_psi(ndmx,ndim,notcnv,psi,e)
  !    calculates (diag(h)-e)^-1 * psi, diagonal approx. to (h-e)^-1*psi
  !    the first nvec columns contain the trial eigenvectors
  !
  ! allocate the work arrays
  !

  call start_clock ('cegterg')
  allocate ( psi( ndmx, nvecx))    
  allocate (hpsi( ndmx, nvecx))    
  allocate (spsi( ndmx, nvecx))    
  allocate (sr( nvecx, nvecx))    
  allocate (hr( nvecx, nvecx))    
  allocate (vr( nvecx, nvecx))    
  allocate (ew( nvecx))    
  allocate (conv( nvec))    

  if (nvec.gt.nvecx / 2) call error ('regter', 'nvecx is too small',1)
  !
  !     prepare the hamiltonian for the first iteration
  !
  notcnv = nvec
  nbase = nvec
  call DCOPY (2*ndmx*nvec, evc, 1, psi, 1)
  !
  !     hpsi contains h times the basis vectors
  !
  call h_psi (ndmx, ndim, nvec, psi, hpsi, spsi)
  !
  !     hr contains the projection of the hamiltonian onto the reduced space
  !     vr contains the eigenvectors of hr
  !
  call setv (nvecx * nvecx, 0.d0, hr, 1)
  call setv (nvecx * nvecx, 0.d0, vr, 1)
  call DGEMM ('t', 'n', nbase, nbase, 2*ndim, 2.d0 , psi, &
       2*ndmx, hpsi, 2*ndmx, 0.d0, hr, nvecx)
  if (gstart.eq.2) call DGER (nbase, nbase, -1.d0, psi, 2*ndmx, &
       hpsi, 2*ndmx, hr, nvecx)
#ifdef PARA
  call reduce (nbase * nvecx, hr)
#endif
  call setv (nvecx * nvecx, 0.d0, sr, 1)
  call DGEMM ('t', 'n', nbase, nbase, 2*ndim, 2.d0, psi, &
       2*ndmx, spsi, 2*ndmx, 0.d0, sr, nvecx)
  if (gstart.eq.2) call DGER (nbase, nbase, -1.d0, psi, 2*ndmx, &
       spsi, 2*ndmx, sr, nvecx)
#ifdef PARA
  call reduce (nbase * nvecx, sr)
#endif

  do n = 1, nbase
     e (n) = hr (n, n)
     conv (n) = 0
     vr (n, n) = 1.d0
  enddo
  !
  !       iterate
  !
  do kter = 1, maxter
     iter = kter+1
     np = nbase
     call start_clock ('update')
     do n = 1, nvec
        if (conv (n) .eq.0) then
           !
           !     this root not yet converged ... set a new basis vector
           !     (position np) to (h-es)psi ...
           !
           np = np + 1
           ew(np) = e (n)
           call setv (2 * ndim, 0.d0, psi (1, np), 1)
           call DGEMV ('n', 2*ndim, nbase, 1.d0, hpsi, 2*ndmx, &
                vr (1, n) , 1, 0.d0, psi (1, np) , 1)
           call DGEMV ('n', 2*ndim, nbase, -e(n),spsi, 2*ndmx, &
                vr (1, n) , 1, 1.d0, psi (1, np) , 1)
        endif
     enddo
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
#ifdef PARA
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
#ifdef PARA
     call reduce (nvecx * notcnv, hr (1, nbase+1) )
#endif
     call DGEMM ('t', 'n', nbase+notcnv, notcnv, 2*ndim, 2.d0, &
          psi, 2*ndmx, spsi (1, nbase+1) , 2*ndmx, 0.d0, &
          sr (1, nbase+1) , nvecx)
     if (gstart.eq.2) call DGER (nbase+notcnv, notcnv, -1.d0, psi, 2*ndmx, &
          spsi(1,nbase+1), 2*ndmx, sr (1, nbase+1), nvecx)
#ifdef PARA
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
        if (conv (n) .eq.0) then
           if (.not.abs (ew (n) - e (n) ) .le.ethr) then
              notcnv = notcnv + 1
           else
              ! root converged
              conv (n) = 1
           endif
        endif
        e (n) = ew (n)
     end do
     !
     !     if overall convergence has been achieved or the dimension of the
     !     reduced basis set is becoming too large, refresh the basis set
     !     i.e. replace the first nvec elements with the current estimate
     !     of the eigenvectors and set the basis dimension to nvec.
     !
     if (notcnv.eq.0 .or. iter.gt.maxter .or. nbase+notcnv.gt.nvecx) then
        call start_clock ('last')
        call DGEMM ('n', 'n', 2*ndim, nvec, nbase, 1.d0, psi, &
             2*ndmx, vr, nvecx, 0.d0, evc, 2*ndmx)
        if (notcnv.eq.0.or. iter.gt.maxter) then
        !
        !     all roots converged: return
        !
           call stop_clock ('last')
           goto 10
        end if
        call DCOPY (2*ndmx*nvec, evc, 1, psi, 1)
        call DGEMM ('n', 'n', 2*ndim, nvec, nbase, 1.d0, spsi, &
             2*ndmx, vr, nvecx, 0.d0, psi(1,nvec+1), 2*ndmx)
        call DCOPY (2 * ndmx * nvec, psi(1,nvec+1), 1, spsi, 1)

        call DGEMM ('n', 'n', 2*ndim, nvec, nbase, 1.d0, hpsi, &
             2*ndmx, vr, nvecx, 0.d0, psi(1,nvec+1), 2*ndmx)
        call DCOPY (2 * ndmx * nvec, psi(1,nvec+1), 1, hpsi, 1)
        !
        !     modify the reduced hamiltonian accordingly
        !
        call setv (nvecx * nvecx, 0.d0, hr, 1)
        call DCOPY (nbase, e, 1, hr, nvecx + 1 )
        call setv (nvecx * nvecx, 0.d0, sr, 1)
        call setv (nvecx, 1.d0, sr, nvecx + 1 )
        call stop_clock ('last')
        nbase = nvec
        call setv (nvecx * nbase, 0.d0, vr, 1)
        do n = 1, nbase
           vr (n, n) = 1.d0
        enddo
     endif
  enddo

10 continue
  do n = 1, nvec
     if (conv (n) .eq.0) write (6, '("   WARNING: e(",i3,") =",&
          & f10.5," is not converged to within ",1pe8.1)') n, e(n), ethr
  enddo

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

