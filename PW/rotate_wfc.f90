!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine rotate_wfc (ndim, ndmx, nvec, nbnd, en, psi)
  !----------------------------------------------------------------------
  !
  !   Hamiltonian diagonalization in the subspace spanned
  !   by nvec states psi.
  !   Produces on output nbnd eigenvectors
  !   (nbnd <= nvec). Calls h_psi to calculate H|psi> ans S|psi>
  !
#include "machine.h"
  use parameters
  implicit none
  !

  integer :: ndim, ndmx, nvec, nbnd
  ! dimension of the matrix to be diagonalized
  ! leading dimension of matrix psi, as declared in the calling pgm unit
  ! input number of states
  ! output number of states

  real(kind=DP) :: en (nvec)
  ! eigenvalues

  complex(kind=DP) :: psi (ndmx,nbnd)
  ! eigenvectors

  ! auxiliary variables:
  complex(kind=DP), allocatable :: psip (:,:), psis (:,:), &
       hc (:,:), sc (:,:), vc (:,:)
  !
  !
  call start_clock ('wfcrot')
  allocate (psip(  ndmx , nvec))    
  allocate (psis(  ndmx , nvec))    
  allocate (hc( nvec , nvec))    
  allocate (sc( nvec , nvec))    
  allocate (vc( nvec , nvec))    

  call h_psi (ndmx, ndim, nvec, psi, psip, psis)

  hc(:,:) = (0.d0, 0.d0)
  sc(:,:) = (0.d0, 0.d0)
  call ZGEMM ('c', 'n', nvec, nvec, ndim, (1.d0, 0.d0) , psi, ndmx, &
       psip, ndmx, (0.d0, 0.d0) , hc, nvec)
#ifdef __PARA
  call reduce (2 * nvec * nvec, hc)
#endif
  call ZGEMM ('c', 'n', nvec, nvec, ndim, (1.d0, 0.d0) , psi, ndmx, &
       psis, ndmx, (0.d0, 0.d0) , sc, nvec)
#ifdef __PARA
  call reduce (2 * nvec * nvec, sc)
#endif

  call cdiaghg (nvec, nbnd, hc, sc, nvec, en, vc)

  call ZGEMM ('n', 'n', ndim, nbnd, nvec, (1.d0, 0.d0) , psi, ndmx, &
       vc, nvec, (0.d0, 0.d0) , psip, ndmx)
  psi(:,:) = psip(:,1:nbnd)

  deallocate (vc)
  deallocate (sc)
  deallocate (hc)
  deallocate (psis)
  deallocate (psip)

  call stop_clock ('wfcrot')
  return
end subroutine rotate_wfc

