!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine dynmat_nlcc (imode0, drhop, npe)
  !
  !! This routine adds a contribution to the dynamical matrix due
  !! to the NLCC.
  !

  USE kinds, only : DP
  USE mp_bands,  ONLY: intra_bgrp_comm
  USE mp,        ONLY: mp_sum
  USE ions_base, ONLY : nat
  USE cell_base, ONLY : omega
  USE fft_base, ONLY : dfftp
  USE noncollin_module, ONLY : nspin_lsda, nspin_mag
  USE dynmat, ONLY : dyn, dyn_rec
  USE modes,  ONLY : nirr, npert, u
  USE uspp,   ONLY : nlcc_any
  USE dv_of_drho_lr,    ONLY : dv_of_drho_xc
  USE control_lr, ONLY: lmultipole

  implicit none

  integer, intent(in) :: imode0
  !! input: the starting mode
  integer, intent(in) :: npe
  !! input: the number of perturbations
  complex(DP), intent(in) :: drhop (dfftp%nnr, nspin_mag, npe)
  !! input: the change of density due to perturbation

  integer :: nrtot, ipert, jpert, is, irr, mode, mode1
  ! the total number of points
  ! counter on perturbations
  ! counter on spin
  ! counter on representations
  ! counter on modes

  complex(DP) :: dyn1 (3 * nat, 3 * nat)
  ! auxiliary dynamical matrix
  complex(DP), allocatable :: drhoc (:), dvaux (:,:)
  ! the change of the core
  ! the change of the potential

  real(DP) :: fac
  ! auxiliary factor

  if (.not.nlcc_any) return

  allocate (drhoc(  dfftp%nnr))
  allocate (dvaux(  dfftp%nnr, nspin_mag))

  dyn1 (:,:) = (0.d0, 0.d0)
!
!  compute the exchange and correlation potential for this mode
!
  nrtot = dfftp%nr1 * dfftp%nr2 * dfftp%nr3
  fac = 1.d0 / DBLE (nspin_lsda)
!
!  Compute the change of xc potential due to the perturbation
!
  do ipert = 1, npe
     mode = imode0 + ipert
     call addcore(u(1, mode), drhoc)
     !
     dvaux(:, :) = (0.d0, 0.d0)
     IF (.NOT. lmultipole) THEN
       CALL dv_of_drho_xc(dvaux, drho = drhop(1, 1, ipert), drhoc = drhoc)
     ELSE
       CALL dv_of_drho_xc(dvaux, drho = drhop(1, 1, ipert))
     ENDIF
     !
     mode1 = 0
     do irr = 1, nirr
        do jpert = 1, npert (irr)
           mode1 = mode1 + 1
           call addcore(u(1, mode1), drhoc)
           do is = 1, nspin_lsda
              dyn1 (mode, mode1) = dyn1 (mode, mode1) + &
                   dot_product (dvaux (:,is), drhoc) * &
                   omega * fac / DBLE (nrtot)
           enddo
        enddo
     enddo
  enddo
  !
  ! collect contributions from all r/G points.
  !
  call mp_sum ( dyn1, intra_bgrp_comm )

  dyn (:,:) = dyn(:,:) + dyn1(:,:)
  dyn_rec(:,:)=dyn_rec(:,:)+dyn1(:,:)
  deallocate (dvaux)
  deallocate (drhoc)
  return
end subroutine dynmat_nlcc
