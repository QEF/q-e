!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine rotate_wfc_nc &
     (npwx, npw, nstart, nbnd, psi, npol, overlap, evc, e)
  !----------------------------------------------------------------------
  !
  !   Hamiltonian diagonalization in the subspace spanned
  !   by nstart states psi (atomic or random wavefunctions).
  !   Produces on output nbnd eigenvectors (nbnd <= nstart) in evc.
  !   Calls h_psi to calculate H|psi> ans S|psi>
  !
#include "f_defs.h"
  USE kinds
  implicit none
  !
  integer :: npw, npwx, nstart, nbnd, npol, i, j, idx
  ! dimension of the matrix to be diagonalized
  ! leading dimension of matrix psi, as declared in the calling pgm unit
  ! input number of states
  ! output number of states
  ! number of coordinates of wfc
  ! counters
  logical :: overlap
  ! if .false. : S|psi> not needed

  complex(DP) :: psi (npwx, npol, nstart), evc (npwx, npol, nbnd),ZDOTC
  ! input and output eigenvectors (may overlap)

  real(DP) :: e (nbnd)
  ! eigenvalues

  ! auxiliary variables:
  complex(DP), allocatable :: hpsi (:,:,:), spsi (:,:,:), &
       hc (:,:), sc (:,:), vc (:,:) 
  real(DP), allocatable :: en (:)
  external ZDOTC
  !
  !
  call start_clock ('wfcrot')
  allocate (hpsi(  npwx ,npol, nstart))    
  if (overlap) allocate (spsi(  npwx , npol, nstart))    
  allocate (hc( nstart , nstart))    
  allocate (sc( nstart , nstart))    
  allocate (vc( nstart , nstart))    
  allocate (en( nstart ))
  !!!!!!!!!!!!!!!!! 
  hpsi=(0.d0,0.d0)
  if (overlap) spsi=(0.d0,0.d0)
  hc=(0.d0,0.d0)
  sc=(0.d0,0.d0)
  vc=(0.d0,0.d0)
  en=0.d0
  !!!!!!!!!!!!!!!!
  !
  ! Set up the Hamiltonian and Overlap matrix
  !
  call h_psi_nc (npwx, npw, nstart, psi(1,1,1), hpsi(1,1,1))
  if (overlap) call s_psi_nc (npwx, npw, nstart, psi(1,1,1), spsi(1,1,1))

  if (npol.eq.1) then
     call ZGEMM ('c','n',nstart,nstart,npw,(1.d0,0.d0),psi(1,1,1),npwx, &
       hpsi(1,1,1), npwx, (0.d0, 0.d0) , hc, nstart)
  else
     do j=1,nstart
        do i=1,nstart
           hc(i,j) = ZDOTC(npw,psi(1,1,i),1,hpsi(1,1,j),1) + &
                     ZDOTC(npw,psi(1,2,i),1,hpsi(1,2,j),1)
        enddo
     enddo
  endif
#ifdef __PARA
  call reduce (2 * nstart * nstart, hc)
#endif
  if (overlap) then
     if (npol.eq.1) then
        call ZGEMM ('c','n',nstart,nstart,npw,(1.d0,0.d0),psi(1,1,1),npwx, &
             spsi(1,1,1),npwx,(0.d0,0.d0),sc,nstart)
     else
        do j=1,nstart
           do i=1,nstart
              sc(i,j) = ZDOTC(npw,psi(:,1,i),1,spsi(:,1,j),1) + &
                        ZDOTC(npw,psi(:,2,i),1,spsi(:,2,j),1)
              ! WRITE( stdout,*) ' sc ', i, j, sc(i,j)
           enddo
        enddo
     end if
  else
     if (npol.eq.1) then
        call ZGEMM ('c','n',nstart,nstart,npw,(1.d0,0.d0),psi(1,1,1),npwx, &
             psi(1,1,1),npwx,(0.d0,0.d0),sc,nstart)
     else
        do j=1,nstart
           do i=1,nstart
              sc(i,j) = ZDOTC(npw,psi(:,1,i),1,psi(:,1,j),1) + &
                   ZDOTC(npw,psi(:,2,i),1,psi(:,2,j),1)
              ! WRITE( stdout,*) ' sc ', i, j, sc(i,j)
           enddo
        enddo
     endif
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
  if (npol.eq.1) then
     call ZGEMM ('n', 'n', npw, nbnd, nstart, (1.d0, 0.d0) , psi(1,1,1), npwx, &
                 vc, nstart, (0.d0, 0.d0) , hpsi(1,1,1), npwx)
  else
     do j=1,nbnd
        do i=1,npw
           hpsi(i,1,j) = (0.d0,0.d0)
           hpsi(i,2,j) = (0.d0,0.d0)
           do idx=1,nstart
              hpsi(i,1,j) = hpsi(i,1,j) + psi(i,1,idx)*vc(idx,j)
              hpsi(i,2,j) = hpsi(i,2,j) + psi(i,2,idx)*vc(idx,j)
           enddo
        enddo
     enddo
  endif
  evc(:, :,  :) = hpsi(:,:, 1:nbnd)

  deallocate (en)
  deallocate (vc)
  deallocate (sc)
  deallocate (hc)
  if (overlap) deallocate (spsi)
  deallocate (hpsi)

  call stop_clock ('wfcrot')
  return
end subroutine rotate_wfc_nc

