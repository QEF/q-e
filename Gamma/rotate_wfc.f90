!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine rotate_wfc (ndim, ndmx, nvec, gstart, nbnd, en, psi)  
  !----------------------------------------------------------------------
  !
  !   Hamiltonian diagonalization in the subspace spanned
  !   by nvec states psi. Produces on output nbnd eigenvectors
  !   (nbnd <= nvec). Calls h_psi to calculate H|psi> and S|psi>.
  !   This version assumes real wavefunctions (k=0) with only
  !   half plane waves stored: psi(-G)=psi*(G), except G=0
  !
#include "machine.h"
  use parameters
  use allocate 
  implicit none
  !
  integer :: ndim, ndmx, nvec, nbnd, gstart
  ! dimension of the matrix to be diagonalized
  ! leading dimension of matrix psi, as declared in the calling pgm unit
  ! input number of states
  ! output number of states
  real(kind=DP) :: en (nvec)  
  ! energy eigenvalues
  complex(kind=DP) :: psi (ndmx,nbnd)  
  ! eigenvectors
  ! auxiliary variables:
  complex(kind=DP), pointer :: hpsi (:,:), spsi (:,:)
  real(kind=8), pointer :: hr (:,:), sr (:,:), vr (:,:)
  !
  call start_clock ('wfcrot')  
  call mallocate(hpsi,  ndmx , nvec)  
  call mallocate(spsi,  ndmx , nvec)  
  call mallocate(hr, nvec , nvec)  
  call mallocate(sr, nvec , nvec)  
  call mallocate(vr, nvec , nvec)  

  call h_psi (ndmx, ndim, nvec, psi, hpsi, spsi)  

  call DGEMM ('t', 'n', nvec, nvec, 2*ndim, 2.d0 , psi, 2*ndmx, &
       hpsi, 2*ndmx, 0.d0, hr, nvec)
  if (gstart==2) call DGER (nvec, nvec,-1.d0, psi, 2*ndmx, &
       hpsi, 2*ndmx, hr, nvec)
#ifdef PARA
  call reduce (nvec * nvec, hr)  
#endif
  call DGEMM ('t', 'n', nvec, nvec, 2*ndim, 2.d0 , psi, 2*ndmx, &
       spsi, 2*ndmx, 0.d0, sr, nvec)
  if (gstart==2) call DGER (nvec, nvec,-1.d0, psi, 2*ndmx, &
       spsi, 2*ndmx, sr, nvec)
#ifdef PARA
  call reduce (nvec * nvec, sr)  
#endif
  !
  call rdiaghg (nvec, nbnd, hr, sr, nvec, en, vr)  
  !
  hpsi(:,:) = (0.d0, 0.d0)
  call DGEMM ('n', 'n', 2*ndim, nbnd, nvec, 1.d0, psi, 2*ndmx, &
       vr, nvec, 0.d0, hpsi, 2*ndmx)
  psi(:,1:nbnd) = hpsi(:,1:nbnd)
  !
  call mfree (vr)  
  call mfree (sr)  
  call mfree (hr)  
  call mfree (spsi)  
  call mfree (hpsi)  

  call stop_clock ('wfcrot')  
  return  
end subroutine rotate_wfc

