!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine cinitcgg_nc (npwx, npw, nstart, nbnd, psi, evc, e, overlap, npol)
  !-----------------------------------------------------------------------
  !
  !   Hamiltonian diagonalization in the subspace spanned
  !   by nstart states psi (atomic or random wavefunctions).
  !   Produces on output nbnd eigenvectors (nbnd <= nstart) in evc.
  !   Minimal memory use - evc and psi may overlap
  !   Calls h_1psi to calculate H|psi>, S|psi
  !
#include "f_defs.h"
  USE kinds
  implicit none
  !
  integer :: npw, npwx, nstart, nbnd, npol
  ! dimension of the matrix to be diagonalized
  ! leading dimension of matrix psi, as declared in the calling pgm unit
  ! input number of states
  ! output number of states

  complex(kind=DP) :: psi (npwx , npol, nstart), evc(npwx , npol, nbnd)
  ! input and output eigenvectors (may overlap) 

  real(kind=DP) :: e (nbnd)
  ! eigenvalues
  logical :: overlap
  !
  !   local variables
  !
  integer :: m, ibnd, i, j, ipol
  complex(kind=DP), allocatable ::  hpsi (:,:), spsi (:,:), hc (:,:,:), sc (:,:)
  complex(kind=DP) :: ZDOTC, ZDOTU
  real(kind=DP), allocatable :: en(:)
  real(kind=DP) :: DDOT
  !
  !
  if (nbnd.gt.nstart) call errore('cinitcgg_nc','nbnd larger than nstart',1)

  call start_clock ('wfcrot1')
  allocate (spsi( npwx,npol))    
  allocate (hpsi( npwx,npol))    
  allocate (hc( nstart , nstart , 2))    
  allocate (sc( nstart , nstart))   
  allocate (en( nstart ))

  spsi = (0.d0,0.d0)
  hpsi = (0.d0,0.d0)
  hc = (0.d0,0.d0)
  sc = (0.d0,0.d0)
  en = 0.d0
  !
  ! Set up the Hamiltonian and Overlap matrix
  !
  do m = 1, nstart
     call h_1psi_nc (npwx, npw, npol, psi(1, 1, m), hpsi, spsi)
     do ipol=1,npol
        hc (m, m, 1) = hc(m, m, 1)  + &
                 DDOT (2 * npw, psi (1, ipol, m), 1, hpsi(1,ipol), 1)
        sc (m, m) = sc(m,m) +         &
                 DDOT (2 * npw, psi (1, ipol, m), 1, spsi(1,ipol), 1) 
        do j = m + 1, nstart
           hc (j,m,1) = hc(j,m,1)+  &
                 ZDOTC (npw, psi (1, ipol, j), 1, hpsi(1,ipol), 1) 
           sc (j,m) = sc(j,m)    +  &
                 ZDOTC (npw, psi (1, ipol, j), 1, spsi(1,ipol), 1)
        enddo
     enddo
     do j = m + 1, nstart
        hc (m, j, 1) = conjg (hc (j, m, 1) )
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
  e (:) = en(1:nbnd)
  !
  !   update the basis set
  !  
  do ipol=1,npol
     do i = 1, npw
        do ibnd = 1, nbnd
           hc(ibnd,1,1)=ZDOTU(nstart,hc(1,ibnd,2),1,psi(i,ipol,1),npwx*npol)
        enddo
        do ibnd = 1, nbnd
           evc (i, ipol, ibnd) = hc (ibnd, 1, 1)
        enddo
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
end subroutine cinitcgg_nc
