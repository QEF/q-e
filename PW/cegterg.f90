!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine cegterg (ndim, ndmx, nvec, nvecx, psi, ethr, overlap, &
     e, notcnv, iter)
  !----------------------------------------------------------------------
  !
  !     iterative solution of the eigenvalue problem:
  !
  !     ( H - e S ) * psi = 0
  !
  !     where H is an hermitean operator, e is a real scalar,
  !     S is an overlap matrix, psi is a complex vector
  !
#include "machine.h"
  use parameters, only : DP
  implicit none
  ! on INPUT
  integer :: ndim, ndmx, nvec, nvecx
  ! dimension of the matrix to be diagonalized
  ! leading dimension of matrix psi, as declared in the calling pgm unit
  ! integer number of searched low-lying roots
  ! maximum dimension of the reduced basis set
  !    (the basis set is refreshed when its dimension would exceed nvecx)
  complex(kind=DP) :: psi (ndmx, nvecx)
  real(kind=DP) :: ethr
  ! energy threshold for convergence
  !   root improvement is stopped, when two consecutive estimates of the root
  !   differ by less than ethr.
  logical :: overlap
  ! if .true.  : use overlap matrix
  ! if .false. : use orthogonalization
  ! on OUTPUT
  !     complex*16 psi   the first nvec columns contain the
  !                      refined estimates of the eigenvectors
  real(kind=DP) :: e (nvecx)
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
  complex(kind=DP), allocatable :: hc (:,:),  smat (:,:), vc (:,:), aux (:)
  ! Hamiltonian on the reduced basis
  ! S matrix on the reduced basis
  ! the eigenvectors of the Hamiltonian
  ! work space
  complex(kind=DP), allocatable :: hpsi (:,:), spsi (:,:)
  ! the product of H and psi
  ! the product of S and psi
  complex(kind=DP) ::  ZDOTC, eau
  ! scalar product routine
  ! auxiliary complex variable
  real(kind=DP), allocatable :: ew (:)
  ! eigenvalues of the reduced hamiltonian
  integer, allocatable  :: conv (:)
  ! if 1 the root is converged if 0 no
  !
  ! Called routines:
  external h_psi, g_psi
  ! h_psi(ndmx,ndim,nvec,psi,hpsi,spsi)
  !     calculates nvec H|psi> and S|psi> (if needed) products.
  !     Vectors psi,hpsi,spsi are dimensioned (ndmx,nvec)
  ! g_psi(ndmx,ndim,notcnv,psi,e)
  !    calculates (diag(h)-e)^-1 * psi, diagonal approx. to (h-e)^-1*psi
  !    the first nvec columns contain the trial eigenvectors
  external ZDOTC
  !
  ! allocate the work arrays
  !

  call start_clock ('cegterg')
  if (overlap) allocate(smat(nvecx, nvecx))
  allocate(hc (nvecx,  nvecx))
  allocate(vc (nvecx,  nvecx))
  allocate(hpsi (ndmx,  nvecx))
  allocate(spsi (ndmx,  nvecx))
  allocate(aux (ndmx * nvec))
  allocate(ew (nvecx))
  allocate(conv (nvec))

  !      write(6,*) 'eneter cegter',hc,vc,hpsi
  if (nvec.gt.nvecx / 2) call error ('cegter', 'nvecx is too small',1)
  !
  !     prepare the hamiltonian for the first iteration
  !
  notcnv = nvec
  nbase = nvec
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
#ifdef PARA
  call reduce (2 * nbase * nvecx, hc)
#endif
  if (overlap) then
     smat(:,:) = (0.d0, 0.d0)
     call ZGEMM ('c', 'n', nbase, nbase, ndim, (1.d0, 0.d0) , psi, &
          ndmx, spsi, ndmx, (0.d0, 0.d0) , smat, nvecx)
#ifdef PARA
     call reduce (2 * nbase * nvecx, smat)
#endif
  endif

  do n = 1, nbase
     e (n) = hc (n, n)
     conv (n) = 0
     vc (n, n) = (1.d0, 0.d0)
  enddo
  !
  !       iterate
  !
  do kter = 1, maxter
     iter = kter
     np = nbase
     call start_clock ('update')

     do n = 1, nvec
        if (conv (n) .eq.0) then
           !
           !     this root not yet converged ... set a new basis vector
           !     (position np) to (h-es)psi ...
           !
           np = np + 1
           e (np) = e (n)
           psi(:,np) = (0.d0,0.d0)
           call ZGEMV ('n', ndim, nbase, (1.d0, 0.d0) , hpsi, ndmx, vc (1, &
                n) , 1, (0.d0, 0.d0) , psi (1, np) , 1)
           eau = DCMPLX ( - e (n), 0.d0)
           call ZGEMV ('n', ndim, nbase, eau, spsi, ndmx, vc (1, n) , 1, &
                (1.d0, 0.d0) , psi (1, np) , 1)
        endif
     enddo

     call stop_clock ('update')
     !
     ! approximate inverse iteration
     !
     call g_psi (ndmx, ndim, notcnv, psi (1, nbase+1), e (nbase+1) )
     !
     ! "normalize" correction vectors psi(*,nbase+1:nbase+notcnv) in order to
     ! the numerical stability of cdiaghg (use ew as work array for norm of p
     !
     do n = 1, notcnv
        ew (n) = ZDOTC (ndim, psi (1, nbase+n), 1, psi (1, nbase+n), 1)
     enddo
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
#ifdef PARA
     call reduce (2 * nvecx * notcnv, hc (1, nbase+1) )
#endif
     if (overlap) then
        call ZGEMM ('c', 'n', nbase+notcnv, notcnv, ndim, (1.d0, 0.d0) &
             , psi, ndmx, spsi (1, nbase+1) , ndmx, (0.d0, 0.d0) , smat (1, &
             nbase+1) , nvecx)
#ifdef PARA
        call reduce (2 * nvecx * notcnv, smat (1, nbase+1) )
#endif
     endif

     call stop_clock ('overlap')
     nbase = nbase+notcnv
     do n = 1, nbase
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
           smat (n, n) = DCMPLX (DREAL (smat (n, n) ), 0.d0)
           do m = n + 1, nbase
              smat (m, n) = conjg (smat (n, m) )
           enddo
        enddo
        call cdiaghg (nbase, nvec, hc, smat, nvecx, ew, vc)
     else
        call cdiagh (nbase, hc, nvecx, ew, vc)
     endif
     !
     !     test for convergence
     !
     notcnv = 0
     do n = 1, nvec
        if (conv (n) .eq.0) then
           if (.not.abs (ew (n) - e (n) ) .le.ethr) then
              notcnv = notcnv + 1
              e (n) = ew (n)
           else
              e (n) = ew (n)
              ! root converged
              conv (n) = 1
           endif
        endif
     enddo
     !     if overall convergence has been achieved or the dimension of the
     !     reduced basis set is becoming too large, then refresh the basis
     !     set. i.e. replace the first nvec elements with the current estimat
     !     of the eigenvectors and set the basis dimension to nvec.

     if (notcnv.eq.0.or.nbase+notcnv.gt.nvecx) then
        call start_clock ('last')
        aux(:) = (0.d0,0.d0)
        call ZGEMM ('n', 'n', ndim, nvec, nbase, (1.d0, 0.d0) , psi, &
             ndmx, vc, nvecx, (0.d0, 0.d0) , aux, ndmx)
        call ZCOPY (ndmx * nvec, aux, 1, psi, 1)
        if (notcnv.eq.0) then
        !
        !     all roots converged: return
        !
           call stop_clock ('last')
           goto 10
        end if

        call ZGEMM ('n', 'n', ndim, nvec, nbase, (1.d0, 0.d0) , spsi, &
             ndmx, vc, nvecx, (0.d0, 0.d0) , aux, ndmx)
        call ZCOPY (ndmx * nvec, aux, 1, spsi, 1)

        call ZGEMM ('n', 'n', ndim, nvec, nbase, (1.d0, 0.d0) , hpsi, &
             ndmx, vc, nvecx, (0.d0, 0.d0) , psi (1, nvec + 1) , ndmx)
        call ZCOPY (ndmx * nvec, psi (1, nvec + 1), 1, hpsi, 1)
        !
        !     modify the reduced hamiltonian accordingly
        !
        hc(:,:) = (0.d0, 0.d0)
        call DCOPY (nbase, e, 1, hc, 2 * (nvecx + 1) )
        if (overlap) then
           smat = (0.d0, 0.d0)
           do n = 1, nvecx
              smat (n, n) = (1.d0, 0.d0)
           enddo
        endif
        call stop_clock ('last')
        nbase = nvec
        vc(:,1:nbase) = (0.d0, 0.d0)
        do n = 1, nbase
           vc (n, n) = (1.d0, 0.d0)
        enddo
     endif
  enddo

  do n = 1, nvec
     if (conv (n) .eq.0) write (6, '("   WARNING: e(",i3,") =",&
          & f10.5," is not converged to within ",1pe8.1)') n, e(n), ethr
  enddo

10 continue
  deallocate (conv)
  deallocate (ew)
  deallocate (aux)
  deallocate (spsi)
  deallocate (hpsi)
  deallocate (vc)
  deallocate (hc)
  if (overlap) deallocate (smat)

  call stop_clock ('cegterg')
  return
end subroutine cegterg

