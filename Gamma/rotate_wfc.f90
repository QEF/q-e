!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine rotate_wfc &
     (npwx, npw, nstart, gstart, nbnd, psi, overlap, evc, e)
  !----------------------------------------------------------------------
  !
  !   Hamiltonian diagonalization in the subspace spanned
  !   by nstart states psi (atomic or random wavefunctions).
  !   Produces on output nbnd eigenvectors (nbnd <= nstart) in evc.
  !   (nbnd <= nvec). Calls h_psi, s_psi to calculate H|psi> and S|psi>.
  !   This version assumes real wavefunctions (k=0) with only
  !   half plane waves stored: psi(-G)=psi*(G), except G=0
  !
#include "machine.h"
  use parameters
  implicit none
  !
  integer :: npw, npwx, nstart, nbnd, gstart
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
  complex(kind=DP), allocatable :: hpsi (:,:), spsi (:,:)
  real(kind=8), allocatable :: hr (:,:), sr (:,:), vr (:,:), en(:)
  !
  call start_clock ('wfcrot')
  allocate (hpsi(  npwx , nstart))    
  if (overlap) allocate (spsi(  npwx , nstart))    
  allocate (hr( nstart , nstart))    
  allocate (sr( nstart , nstart))    
  allocate (vr( nstart , nstart))    
  allocate (en( nstart ))
  !
  ! Set up the Hamiltonian and Overlap matrix
  !
  call h_psi (npwx, npw, nstart, psi, hpsi)
  if (overlap) call s_psi (npwx, npw, nstart, psi, spsi)

  call DGEMM ('t', 'n', nstart, nstart, 2*npw, 2.d0 , psi, 2*npwx, &
       hpsi, 2*npwx, 0.d0, hr, nstart)
  if (gstart==2) call DGER (nstart, nstart,-1.d0, psi, 2*npwx, &
       hpsi, 2*npwx, hr, nstart)
#ifdef __PARA
  call reduce (nstart * nstart, hr)
#endif
  if (overlap) then
     call DGEMM ('t', 'n', nstart, nstart, 2*npw, 2.d0 , psi, 2*npwx, &
          spsi, 2*npwx, 0.d0, sr, nstart)
     if (gstart==2) call DGER (nstart, nstart,-1.d0, psi, 2*npwx, &
          spsi, 2*npwx, sr, nstart)
  else
     call DGEMM ('t', 'n', nstart, nstart, 2*npw, 2.d0 , psi, 2*npwx, &
          psi, 2*npwx, 0.d0, sr, nstart)
     if (gstart==2) call DGER (nstart, nstart,-1.d0, psi, 2*npwx, &
          psi, 2*npwx, sr, nstart)
  end if
#ifdef __PARA
  call reduce (nstart * nstart, sr)
#endif
  !
  ! Diagonalize
  !
  call rdiaghg (nstart, nbnd, hr, sr, nstart, en, vr)
  !
  e (:) = en(1:nbnd)
  !
  !   update the basis set
  !
  call DGEMM ('n', 'n', 2*npw, nbnd, nstart, 1.d0, psi, 2*npwx, &
       vr, nstart, 0.d0, hpsi, 2*npwx)
  evc(:, :) = hpsi(:, 1:nbnd)
  !
  deallocate (en)
  deallocate (vr)
  deallocate (sr)
  deallocate (hr)
  if (overlap) deallocate (spsi)
  deallocate (hpsi)

  call stop_clock ('wfcrot')
  return
end subroutine rotate_wfc

