!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine cinitcgg (npwx, npw, nstart, nbnd, psi, evc, e)
  !-----------------------------------------------------------------------
  !
  !   Hamiltonian diagonalization in the subspace spanned
  !   by nstart states psi (atomic or random wavefunctions).
  !   Produces on output nbnd eigenvectors (nbnd <= nstart) in evc.
  !   Minimal memory use - evc and psi may overlap
  !   Calls h_1psi to calculate H|psi>, S|psi
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

  complex(kind=DP) :: psi (npwx, nstart), evc(npwx, nbnd)
  ! input and output eigenvectors (may overlap) 

  real(kind=DP) :: e (nbnd)
  ! eigenvalues
  !
  !   local variables
  !
  integer :: m, ibnd, i, j
  complex(kind=DP), allocatable ::  hpsi (:), spsi (:), hc (:,:,:), sc (:,:)
  complex(kind=DP) :: ZDOTC, ZDOTU
  real(kind=DP), allocatable :: en(:)
  real(kind=DP) :: DDOT
  !
  !
  call start_clock ('wfcrot1')
  allocate (spsi( npwx))    
  allocate (hpsi( npwx))    
  allocate (hc( nstart , nstart , 2))    
  allocate (sc( nstart , nstart))   
  allocate (en( nstart ))
  !
  ! Set up the Hamiltonian and Overlap matrix
  !
  do m = 1, nstart
     call h_1psi (npwx, npw, psi (1, m), hpsi, spsi)
     hc (m, m, 1) = DDOT (2 * npw, psi (1, m), 1, hpsi, 1)
     sc (m, m) = DDOT (2 * npw, psi (1, m), 1, spsi, 1)
     do j = m + 1, nstart
        hc (j, m, 1) = ZDOTC (npw, psi (1, j), 1, hpsi, 1)
        hc (m, j, 1) = conjg (hc (j, m, 1) )
        sc (j, m) = ZDOTC (npw, psi (1, j), 1, spsi, 1)
        sc (m, j) = conjg (sc (j, m) )
     enddo

  enddo
#ifdef __PARA
  call reduce (2 * nstart * nstart, hc (1, 1, 1) )
  call reduce (2 * nstart * nstart, sc)
#endif
  !
  ! diagonalize
  !
  call cdiaghg (nstart, nbnd, hc, sc, nstart, en, hc (1, 1, 2))
  !
  e (1:nbnd) = en(1:nbnd)
  !
  !   update the basis set
  !  
  do i = 1, npw
     do ibnd = 1, nbnd
        hc (ibnd, 1, 1) = ZDOTU (nstart, hc(1,ibnd,2), 1, psi(i,1), npwx)
     enddo
     do ibnd = 1, nbnd
        evc (i, ibnd) = hc (ibnd, 1, 1)
     enddo
  enddo
  !
  deallocate (en)
  deallocate (sc)
  deallocate (hc)
  deallocate (hpsi)
  deallocate (spsi)
  call stop_clock ('wfcrot1')

  return
end subroutine cinitcgg

