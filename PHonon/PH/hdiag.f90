!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine hdiag( npw, max_iter, avg_iter, et_ )
  !
  !! Diagonalizes the unperturbed Hamiltonian in a non-selfconsistent way
  !! by Conjugate Gradient (band-by-band).
  !
  USE kinds,     ONLY : DP
  USE gvect,     ONLY: g, gstart
  USE wvfct,     ONLY: g2kin, nbnd, npwx
  USE uspp,      ONLY: vkb, okvan
  USE noncollin_module,    ONLY: npol
  USE wavefunctions,ONLY: evc
  USE ramanm,    ONLY: eth_ns
  implicit none
  !
  integer :: npw
  !! number of plane waves
  integer :: max_iter
  !! maximum number of iterations
  real(DP) :: avg_iter
  !! iteration number in the diagonalization
  real(DP) :: et_(nbnd)
  !! eigenvalues of the diagonalization
  !
  ! ... local variables:
  !
  REAL(DP) :: cg_iter
  ! number of iteration in CG
  INTEGER  :: ig, ntry, notconv
  ! counter on G vectors
  ! number or repeated call to diagonalization in case of non convergence
  ! number of notconverged elements
  INTEGER, ALLOCATABLE :: btype(:)
    ! type of band: valence (1) or conduction (0)
  REAl(DP), ALLOCATABLE :: h_prec(:)
    ! preconditioning matrix (diagonal)
! CG diagonalization uses these external routines on a single band
   external hs_1psi, s_1psi
!  subroutine hs_1psi(npwx,npw,psi,hpsi,spsi)  computes H*psi and S*psi
!  subroutine s_1psi(npwx,npw,psi,spsi)  computes S*psi (if needed)

  call start_clock ('hdiag')

  allocate (h_prec( npwx), btype(nbnd))
  !
  !   various initializations
  !
  btype(:) = 1
  !
  ! Conjugate-Gradient diagonalization
  !
  h_prec=1.0_DP
  FORALL( ig = 1 : npwx )
      h_prec(ig) = 1.D0 + g2kin(ig) + SQRT( 1.D0 + ( g2kin(ig) - 1.D0 )**2 )
  END FORALL

  DO ntry = 1, 5
    if (ntry > 1) then
       call rotate_wfc &
         ( npwx, npw, nbnd, gstart, nbnd, evc, npol, okvan, evc, et_ )
       avg_iter = avg_iter + 1.d0
    endif
    CALL ccgdiagg( hs_1psi, s_1psi, h_prec, &
         npwx, npw, nbnd, npol, evc, et_, btype, eth_ns, &
         max_iter, .true., notconv, cg_iter)
    avg_iter = avg_iter + cg_iter
 
    IF(notconv == 0 ) EXIT
  ENDDO

  deallocate (btype, h_prec)
  call stop_clock ('hdiag')

  return
end subroutine hdiag

