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
  use parameters, only : DP
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
  ! if .true.  : use overlap matrix
  ! if .false. : use orthogonalization
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
  !  one parameter
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

  call start_clock ('cegterg')
  allocate( psi (ndmx,  nvecx))
  allocate(hpsi (ndmx,  nvecx))
  allocate(spsi (ndmx,  nvecx))
  if (overlap) allocate(sc(nvecx, nvecx))
  allocate(hc (nvecx,  nvecx))
  allocate(vc (nvecx,  nvecx))
  allocate(ew (nvecx))
  allocate(conv (nvec))

  !      write(6,*) 'eneter cegter',hc,vc,hpsi
  if (nvec > nvecx / 2) call errore ('cegter', 'nvecx is too small',1)
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
  !   hc contains the projection of the hamiltonian onto the reduced space
  !   vc contains the eigenvectors of hc
  !
  hc(:,:) = (0.d0, 0.d0)
  vc(:,:) = (0.d0, 0.d0)
  call ZGEMM ('c', 'n', nbase, nbase, ndim, (1.d0, 0.d0) , psi, &
       ndmx, hpsi, ndmx, (0.d0, 0.d0) , hc, nvecx)
#ifdef __PARA
  call reduce (2 * nbase * nvecx, hc)
#endif
  if (overlap) then
     sc(:,:) = (0.d0, 0.d0)
     call ZGEMM ('c', 'n', nbase, nbase, ndim, (1.d0, 0.d0) , psi, &
          ndmx, spsi, ndmx, (0.d0, 0.d0) , sc, nvecx)
#ifdef __PARA
     call reduce (2 * nbase * nvecx, sc)
#endif
  endif

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
     if (notcnv < nvec) then
        np = nbase
        do n = 1, nvec
           if ( .not.conv (n) ) then
              !
              !     this root not yet converged ... set a new basis vector
              !     (position np) to (h-es)psi ...
              !
              np = np + 1
              ! for use in g_psi
              ew (np) = e (n)
              call ZGEMV ('n', ndim, nbase, (1.d0, 0.d0) , hpsi, ndmx, &
                   vc (1, n) , 1, (0.d0, 0.d0) , psi (1, np) , 1)
              eau = DCMPLX ( - e (n), 0.d0)
              call ZGEMV ('n', ndim, nbase, eau, spsi, ndmx, vc (1, n) , 1, &
                   (1.d0, 0.d0) , psi (1, np) , 1)
           endif
        enddo
     else
        !
        !     expand the basis set with new basis vectors (h-es)psi ...
        !
        call ZGEMM ('n', 'n', ndim, nvec, nbase, (1.d0, 0.d0), spsi, &
             ndmx, vc, nvecx, (0.d0, 0.d0), psi (1, nbase+1), ndmx)
        do n = 1, nvec
           ! for use in g_psi
           ew (nbase+n) = e (n)
           psi (:,nbase+n) = - e(n) * psi(:,nbase+n)
        end do
        call ZGEMM ('n', 'n', ndim, nvec, nbase, (1.d0, 0d0), hpsi, &
             ndmx, vc, nvecx, (1.d0, 0.d0), psi (1, nbase+1), ndmx)
     end if

     call stop_clock ('update')
     !
     ! approximate inverse iteration
     !
     call g_psi (ndmx, ndim, notcnv, psi (1, nbase+1), ew (nbase+1) )
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
#ifdef DEBUG_DAVIDSON
     write (6,'(a,18f10.6)') 'NRM=',(ew(n),n=1,notcnv)
#endif
     !
     !   here compute the hpsi and spsi of the new functions
     !
     call h_psi (ndmx, ndim, notcnv, psi (1, nbase+1), &
                 hpsi (1, nbase+1), spsi (1, nbase+1) )
     !
     !  orthonormalize (if not using diagonalization with overlap matrix)
     !
     if (.not.overlap) call cgramg1 (ndmx, nvecx, ndim, nbase+1, &
                                     nbase+notcnv, psi, spsi, hpsi)
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
#ifdef __PARA
        call reduce (2 * nvecx * notcnv, sc (1, nbase+1) )
#endif
     endif

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
     if (overlap) then
        do n = 1, nbase
           sc (n, n) = DCMPLX (DREAL (sc (n, n) ), 0.d0)
           do m = n + 1, nbase
              sc (m, n) = conjg (sc (n, m) )
           enddo
        enddo
        call cdiaghg (nbase, nvec, hc, sc, nvecx, ew, vc)
#ifdef DEBUG_DAVIDSON
        write (6,'(a,18f10.6)') 'EIG=',(e(n),n=1,nvec)
        write (6,'(a,18f10.6)') 'EIG=',(ew(n),n=1,nvec)
        write (6,*) 
#endif
     else
        call cdiagh (nbase, hc, nvecx, ew, vc)
     endif
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
!!! TEST : update all eigenvectors if more than 1/4 are converged
     if (notcnv > nvec/4) then
        notcnv = nvec
        conv(:) = .false.
     end if
!!! END OF TEST
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

        call ZGEMM ('n', 'n', ndim, nvec, nbase, (1.d0, 0.d0) , spsi, &
             ndmx, vc, nvecx, (0.d0, 0.d0) , psi(1, nvec + 1), ndmx)
        spsi(:, 1:nvec) = psi(:, nvec+1:2*nvec)

        call ZGEMM ('n', 'n', ndim, nvec, nbase, (1.d0, 0.d0) , hpsi, &
             ndmx, vc, nvecx, (0.d0, 0.d0) , psi (1, nvec + 1) , ndmx)
        hpsi(:, 1:nvec) = psi(:, nvec+1:2*nvec)
        !
        !     refresh the reduced hamiltonian 
        !
        nbase = nvec
        hc (:, 1:nbase) = (0.d0, 0.d0)
        vc (:, 1:nbase) = (0.d0, 0.d0)
        do n = 1, nbase
           hc (n, n) = e(n)
           vc (n, n) = (1.d0, 0.d0)
        enddo
        if (overlap) then
           sc (:, 1:nbase) = (0.d0, 0.d0)
           do n = 1, nbase
              sc (n, n) = (1.d0, 0.d0)
           enddo
        endif
        call stop_clock ('last')

     endif
  enddo

10 continue
  deallocate (spsi)
  deallocate (hpsi)
  deallocate ( psi)
  if (overlap) deallocate (sc)
  deallocate (hc)
  deallocate (vc)
  deallocate (ew)
  deallocate (conv)

  call stop_clock ('cegterg')
  return
end subroutine cegterg

