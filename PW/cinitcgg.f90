!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine cinitcgg (npwx, npw, nstart, nbnd, psi, e)
  !-----------------------------------------------------------------------
  !
  !   Hamiltonian diagonalization in the subspace spanned
  !   by nstart states psi.
  !   Produces on output nbnd eigenvectors (nbnd <= nstart).
  !   Minimal memory use - calls h_1psi to calculate H|psi>, S|psi
  !
#include "machine.h"
  use parameters
  implicit none
  !
  integer :: npwx, npw, nstart, nbnd
  ! dimension of the matrix to be diagonalized
  ! leading dimension of matrix psi, as declared in the calling pgm unit
  ! input number of states
  ! output number of states

  complex(kind=DP) :: psi (npwx, nstart)
  ! eigenvectors

  real(kind=DP) :: e (nstart)
  ! eigenvalues
  !
  !   local variables
  !
  integer :: m, ibnd, i, j
  complex(kind=DP), allocatable ::  hpsi (:), spsi (:), haux (:,:,:), saux (:,:)
  complex(kind=DP) :: ZDOTC, ZDOTU
  real(kind=DP) :: DDOT
  !
  !
  allocate (spsi(  npwx))    
  allocate (hpsi(  npwx))    
  allocate (haux(  nstart , nstart , 2))    
  allocate (saux(  nstart , nstart))    
  !
  ! diagonalize on the reduced basis set formed by trial eigenvectors
  ! Set up the matrix in the reduced basis
  !
  do m = 1, nstart
     call h_1psi (npwx, npw, psi (1, m), hpsi, spsi)
     haux (m, m, 1) = DDOT (2 * npw, psi (1, m), 1, hpsi, 1)
     saux (m, m) = DDOT (2 * npw, psi (1, m), 1, spsi, 1)
     do j = m + 1, nstart
        haux (j, m, 1) = ZDOTC (npw, psi (1, j), 1, hpsi, 1)
        haux (m, j, 1) = conjg (haux (j, m, 1) )
        saux (j, m) = ZDOTC (npw, psi (1, j), 1, spsi, 1)
        saux (m, j) = conjg (saux (j, m) )
     enddo

  enddo
#ifdef PARA
  call reduce (2 * nstart * nstart, haux (1, 1, 1) )
  call reduce (2 * nstart * nstart, saux)
#endif
  !
  ! diagonalize
  !
  call cdiaghg (nstart, nbnd, haux, saux, nstart, e, haux (1, 1, 2))
  !
  !   update the reduced basis set
  !
  do i = 1, npw
     do ibnd = 1, nbnd
        haux (ibnd, 1, 1) = ZDOTU (nstart,haux(1,ibnd,2),1,psi(i,1),npwx)
     enddo
     do ibnd = 1, nbnd
        psi (i, ibnd) = haux (ibnd, 1, 1)
     enddo
  enddo
  !
  deallocate (saux)
  deallocate (haux)
  deallocate (hpsi)
  deallocate (spsi)

  return
end subroutine cinitcgg

