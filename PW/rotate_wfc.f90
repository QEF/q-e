!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine rotate_wfc (npwx, npw, nstart, nbnd, psi, overlap, evc, e)
  !----------------------------------------------------------------------
  !
  !   Hamiltonian diagonalization in the subspace spanned
  !   by nstart states psi (atomic or random wavefunctions).
  !   Produces on output nbnd eigenvectors (nbnd <= nstart) in evc.
  !   Calls h_psi, s_psi to calculate H|psi> ans S|psi>
  !
#include "machine.h"
  USE kinds
  implicit none
  !

  integer :: npw, npwx, nstart, nbnd
  ! dimension of the matrix to be diagonalized
  ! leading dimension of matrix psi, as declared in the calling pgm unit
  ! input number of states
  ! output number of states
  logical :: overlap
  ! if .false. : S|psi> not needed

  complex(kind=DP) :: psi (npwx, nstart), evc (npwx,nbnd)
  ! input and output eigenvectors (may overlap)

  real(kind=DP) :: e (nbnd)
  ! eigenvalues

  ! auxiliary variables:
  complex(kind=DP), allocatable :: hpsi (:,:), spsi (:,:), &
       hc (:,:), sc (:,:), vc (:,:)
  real(kind=DP), allocatable :: en (:)
  !
  !
  call start_clock ('wfcrot')
  allocate (hpsi(  npwx , nstart))    
  if (overlap) allocate (spsi(  npwx , nstart))    
  allocate (hc( nstart , nstart))    
  allocate (sc( nstart , nstart))    
  allocate (vc( nstart , nstart))    
  allocate (en( nstart ))
  !
  ! Set up the Hamiltonian and Overlap matrix
  !
  call h_psi (npwx, npw, nstart, psi, hpsi)
  if (overlap) call s_psi (npwx, npw, nstart, psi, spsi)

  call ZGEMM ('c', 'n', nstart, nstart, npw, (1.d0, 0.d0) , psi, npwx, &
       hpsi, npwx, (0.d0, 0.d0) , hc, nstart)
#ifdef __PARA
  call reduce (2 * nstart * nstart, hc)
#endif
  if (overlap) then
     call ZGEMM ('c', 'n', nstart, nstart, npw, (1.d0, 0.d0) , psi, npwx, &
       spsi, npwx, (0.d0, 0.d0) , sc, nstart)
  else
     call ZGEMM ('c', 'n', nstart, nstart, npw, (1.d0, 0.d0) , psi, npwx, &
       psi, npwx, (0.d0, 0.d0) , sc, nstart)
  end if
#ifdef __PARA
  call reduce (2 * nstart * nstart, sc)
#endif
  !
  ! Diagonalize
  !
  call cdiaghg (nstart, nbnd, hc, sc, nstart, en, vc)
  !
  e (:) = en(1:nbnd)
  !
  !   update the basis set
  !  
  call ZGEMM ('n', 'n', npw, nbnd, nstart, (1.d0, 0.d0) , psi, npwx, &
       vc, nstart, (0.d0, 0.d0) , hpsi, npwx)
  evc(:, :) = hpsi(:, 1:nbnd)

  deallocate (en)
  deallocate (vc)
  deallocate (sc)
  deallocate (hc)
  if (overlap) deallocate (spsi)
  deallocate (hpsi)

  call stop_clock ('wfcrot')
  return
end subroutine rotate_wfc

