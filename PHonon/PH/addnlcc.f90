!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine addnlcc (imode0, drhoscf, npe)
  !
  !! This routine adds a contribution to the dynamical matrix due
  !! to the NLCC.
  !

  USE kinds, only : DP
  USE ions_base, ONLY : nat
  USE funct,  only : dft_is_nonlocc
  USE xc_lib, only : xclib_dft_is
  USE cell_base, ONLY : omega
  use scf, only : rho, rho_core
  USE gvect, ONLY : g, ngm
  USE fft_base, ONLY : dfftp
  USE noncollin_module, ONLY : nspin_lsda, nspin_gga, nspin_mag
  USE dynmat, ONLY : dyn, dyn_rec
  USE modes,  ONLY : nirr, npert, u
  USE uspp,   ONLY : nlcc_any

  USE mp_bands,  ONLY: intra_bgrp_comm
  USE mp,        ONLY: mp_sum

  USE qpoint,  ONLY : xq
  USE eqv,     ONLY : dmuxc
  USE gc_lr,   ONLY: grho,  dvxc_rr,  dvxc_sr,  dvxc_ss, dvxc_s
  USE dv_of_drho_lr,    ONLY : dv_of_drho_xc

  implicit none

  integer :: imode0
  !! input: the starting mode
  integer :: npe
  !! input: the number of perturbations
  complex(DP) :: drhoscf (dfftp%nnr, nspin_mag, npe)
  !! input: the change of density due to perturbation

  integer :: nrtot, ipert, jpert, is, is1, irr, ir, mode, mode1
  ! the total number of points
  ! counter on perturbations
  ! counter on spin
  ! counter on representations
  ! counter on real space points
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
  dvaux(:, :) = (0.d0, 0.d0)
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
     CALL dv_of_drho_xc(dvaux, drho = drhoscf(1, 1, ipert), drhoc = drhoc)
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
end subroutine addnlcc
