!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine cegterg (ndim, ndmx, nvec, nvecx, evc, ethr, overlap, &
     e, notcnv, iter)
  !----------------------------------------------------------------------
  !
  !     iterative solution of the eigenvalue problem:
  !
  !     ( H - e S ) * evc = 0
  !
  !     where H is an hermitean operator, e is a real scalar,
  !     S is an overlap matrix, evc is a complex vector
  !
#include "machine.h"
  USE io_global,  ONLY : stdout
  USE kinds, only : DP
  use g_psi_mod

  implicit none
  ! on INPUT
  integer :: ndim, ndmx, nvec, nvecx
  ! dimension of the matrix to be diagonalized
  ! leading dimension of matrix evc, as declared in the calling pgm unit
  ! integer number of searched low-lying roots
  ! maximum dimension of the reduced basis set
  !    (the basis set is refreshed when its dimension would exceed nvecx)
  complex(kind=DP) :: evc (ndmx, nvec)
  real(kind=DP) :: ethr
  ! energy threshold for convergence
  !   root improvement is stopped, when two consecutive estimates of the root
  !   differ by less than ethr.
  logical :: overlap
  ! if .false. : do not calculate S|psi>
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
  complex(kind=DP), allocatable :: hc (:,:),  sc (:,:), vc (:,:)
  ! Hamiltonian on the reduced basis
  ! S matrix on the reduced basis
  ! the eigenvectors of the Hamiltonian
  complex(kind=DP), allocatable :: psi(:,:), hpsi (:,:), spsi (:,:)
  ! work space, contains psi
  ! the product of H and psi
  ! the product of S and psi
  complex(kind=DP), external ::  ZDOTC
  ! scalar product routine
  complex(kind=DP) ::  eau
  ! auxiliary complex variable
  real(kind=DP), allocatable :: ew (:)
  ! eigenvalues of the reduced hamiltonian
  logical, allocatable  :: conv (:)
  ! true if the root is converged
  !
  ! Called routines:
  external h_psi, s_psi, g_psi
  ! h_psi(ndmx,ndim,nvec,psi,hpsi)
  !     calculates H|psi>
  ! s_psi(ndmx,ndim,nvec,spsi)
  !     calculates S|psi> (if needed)
  !     Vectors psi,hpsi,spsi are dimensioned (ndmx,nvec)
  ! g_psi(ndmx,ndim,notcnv,psi,e)
  !    calculates (diag(h)-e)^-1 * psi, diagonal approx. to (h-e)^-1*psi
  !    the first nvec columns contain the trial eigenvectors
  !
  ! allocate the work arrays
  !

  test_new_preconditioning = .true.

  call start_clock ('cegterg')
  allocate( psi (ndmx,  nvecx))
  allocate(hpsi (ndmx,  nvecx))
  if (overlap) allocate(spsi (ndmx,  nvecx))
  allocate(sc(nvecx, nvecx))
  allocate(hc (nvecx,  nvecx))
  allocate(vc (nvecx,  nvecx))
  allocate(ew (nvecx))
  allocate(conv (nvec))

  !      WRITE( stdout,*) 'eneter cegter',hc,vc,hpsi
  if (nvec > nvecx / 2) call errore ('cegter', 'nvecx is too small',1)
  !
  !     prepare the hamiltonian for the first iteration
  !
  notcnv = nvec
  nbase = nvec
  if (overlap) spsi = (0.d0, 0.d0)
  psi  = (0.d0, 0.d0)
  hpsi = (0.d0, 0.d0)
  psi(:, 1:nvec) = evc(:, 1:nvec)
  !
  !     hpsi contains h times the basis vectors
  !
  call h_psi (ndmx, ndim, nvec, psi, hpsi)
  if (overlap) call s_psi (ndmx, ndim, nvec, psi, spsi)
  !
  !   hc contains the projection of the hamiltonian onto the reduced space
  !   vc contains the eigenvectors of hc
  !
  hc(:,:) = (0.d0, 0.d0)
  sc(:,:) = (0.d0, 0.d0)
  vc(:,:) = (0.d0, 0.d0)
  call ZGEMM ('c', 'n', nbase, nbase, ndim, (1.d0, 0.d0) , psi, &
       ndmx, hpsi, ndmx, (0.d0, 0.d0) , hc, nvecx)
#ifdef __PARA
  call reduce (2 * nbase * nvecx, hc)
#endif
  if (overlap) then
     call ZGEMM ('c', 'n', nbase, nbase, ndim, (1.d0, 0.d0) , psi, &
          ndmx, spsi, ndmx, (0.d0, 0.d0) , sc, nvecx)
  else
     call ZGEMM ('c', 'n', nbase, nbase, ndim, (1.d0, 0.d0) , psi, &
          ndmx,  psi, ndmx, (0.d0, 0.d0) , sc, nvecx)
  endif
#ifdef __PARA
  call reduce (2 * nbase * nvecx, sc)
#endif

  do n = 1, nbase
     e (n) = hc (n, n)
     conv (n) = .false.
     vc (n, n) = (1.d0, 0.d0)
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
           if (np .ne. n) vc(:,np) = vc(:,n)
           ! for use in g_psi
           ew (nbase+np) = e (n)
        end if
     end do
     !
     !     expand the basis set with new basis vectors (h-es)psi ...
     !
     if (overlap) then
        call ZGEMM ('n', 'n', ndim, notcnv, nbase, (1.d0, 0.d0), spsi, &
             ndmx, vc, nvecx, (0.d0, 0.d0), psi (1, nbase+1), ndmx)
     else
        call ZGEMM ('n', 'n', ndim, notcnv, nbase, (1.d0, 0.d0), psi, &
             ndmx, vc, nvecx, (0.d0, 0.d0), psi (1, nbase+1), ndmx)
     end if
     do np = 1, notcnv
        psi (:,nbase+np) = - ew(nbase+np) * psi(:,nbase+np)
     end do
     call ZGEMM ('n', 'n', ndim, notcnv, nbase, (1.d0, 0d0), hpsi, &
          ndmx, vc, nvecx, (1.d0, 0.d0), psi (1, nbase+1), ndmx)

     call stop_clock ('update')
     !
     ! approximate inverse iteration
     !
     call g_psi (ndmx, ndim, notcnv, psi (1, nbase+1), ew (nbase+1) )

#ifdef DEBUG_DAVIDSON
     np = 0
     ew (1:nvec) = -1.0d0
     do n = 1, nvec
        if (.not.conv(n)) then
           np = np+1
           ew (n) = ZDOTC (ndim, psi (1, nbase+np), 1, psi (1, nbase+np), 1)
        endif
     end do
#ifdef __PARA
     call reduce (nvec, ew)
#endif
     WRITE( stdout,'(a,18f10.6)') 'NRM=',(ew(n),n=1,nvec)
#endif
     !
     ! "normalize" correction vectors psi(*,nbase+1:nbase+notcnv) in order
     ! to improve numerical stability of subspace diagonalization rdiaghg
     ! ew is used as work array : ew = <psi_i|psi_i>, i=nbase+1,nbase+notcnv
     !
     do n = 1, notcnv
        ew (n) = ZDOTC (ndim, psi (1, nbase+n), 1, psi (1, nbase+n), 1)
     enddo
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
                 hpsi (1, nbase+1) )
     if (overlap) call s_psi (ndmx, ndim, notcnv, psi (1, nbase+1), &
                 spsi (1, nbase+1) )
     !
     !     update the reduced hamiltonian
     !
     call start_clock ('overlap')
     call ZGEMM ('c', 'n', nbase+notcnv, notcnv, ndim, (1.d0, 0.d0) , &
          psi, ndmx, hpsi (1, nbase+1) , ndmx, (0.d0, 0.d0) , &
          hc (1, nbase+1) , nvecx)
#ifdef __PARA
     call reduce (2 * nvecx * notcnv, hc (1, nbase+1) )
#endif
     if (overlap) then
        call ZGEMM ('c', 'n', nbase+notcnv, notcnv, ndim, (1.d0, 0.d0), &
             psi, ndmx, spsi (1, nbase+1) , ndmx, (0.d0, 0.d0) , &
             sc (1, nbase+1) , nvecx)
     else
        call ZGEMM ('c', 'n', nbase+notcnv, notcnv, ndim, (1.d0, 0.d0), &
             psi, ndmx,  psi (1, nbase+1) , ndmx, (0.d0, 0.d0) , &
             sc (1, nbase+1) , nvecx)
     end if
#ifdef __PARA
     call reduce (2 * nvecx * notcnv, sc (1, nbase+1) )
#endif
     call stop_clock ('overlap')

     nbase = nbase+notcnv
     do n = 1, nbase
        !  the diagonal of hc must be strictly real 
        hc (n, n) = DCMPLX (DREAL (hc (n, n) ), 0.d0)
        do m = n + 1, nbase
           hc (m, n) = conjg (hc (n, m) )
        enddo
     enddo
     !
     !     diagonalize the reduced hamiltonian
     !
     do n = 1, nbase
        sc (n, n) = DCMPLX (DREAL (sc (n, n) ), 0.d0)
        do m = n + 1, nbase
           sc (m, n) = conjg (sc (n, m) )
        enddo
     enddo
     call cdiaghg (nbase, nvec, hc, sc, nvecx, ew, vc)
#ifdef DEBUG_DAVIDSON
     WRITE( stdout,'(a,18f10.6)') 'EIG=',(e(n),n=1,nvec)
     WRITE( stdout,'(a,18f10.6)') 'EIG=',(ew(n),n=1,nvec)
     WRITE( stdout,*) 
#endif
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
     if ( notcnv == 0 .or. nbase+notcnv > nvecx .or. iter == maxter) then
        call start_clock ('last')
        call ZGEMM ('n', 'n', ndim, nvec, nbase, (1.d0, 0.d0) , psi, &
             ndmx, vc, nvecx, (0.d0, 0.d0) , evc, ndmx)
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
              if ( .not.conv (n) ) WRITE( stdout, '("   WARNING: e(",i3,") =",&
                   f10.5," is not converged to within ",1pe8.1)') n, e(n), ethr
           enddo
#else
           WRITE( stdout, '("   WARNING: ",i5," eigenvalues not converged")') &
                notcnv
#endif
           call stop_clock ('last')
           goto 10
        end if
        !
        !     refresh psi, H*psi and S*psi
        !
        psi(:, 1:nvec) = evc(:, 1:nvec)

        if (overlap) then
           call ZGEMM ('n', 'n', ndim, nvec, nbase, (1.d0, 0.d0) , spsi, &
                ndmx, vc, nvecx, (0.d0, 0.d0) , psi(1, nvec + 1), ndmx)
           spsi(:, 1:nvec) = psi(:, nvec+1:2*nvec)
        endif

        call ZGEMM ('n', 'n', ndim, nvec, nbase, (1.d0, 0.d0) , hpsi, &
             ndmx, vc, nvecx, (0.d0, 0.d0) , psi (1, nvec + 1) , ndmx)
        hpsi(:, 1:nvec) = psi(:, nvec+1:2*nvec)
        !
        !     refresh the reduced hamiltonian 
        !
        nbase = nvec
        hc (:, 1:nbase) = (0.d0, 0.d0)
        sc (:, 1:nbase) = (0.d0, 0.d0)
        vc (:, 1:nbase) = (0.d0, 0.d0)
        do n = 1, nbase
           hc (n, n) = e(n)
           sc (n, n) = (1.d0, 0.d0)
           vc (n, n) = (1.d0, 0.d0)
        enddo
        call stop_clock ('last')

     endif
  enddo

10 continue
  deallocate (conv)
  deallocate (ew)
  deallocate (vc)
  deallocate (hc)
  deallocate (sc)
  if (overlap) deallocate (spsi)
  deallocate (hpsi)
  deallocate ( psi)

  call stop_clock ('cegterg')
  return
end subroutine cegterg

