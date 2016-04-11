!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine addnlcc (imode0, drhoscf, npe)
  !
  !     This routine adds a contribution to the dynamical matrix due
  !     to the NLCC
  !

  USE kinds, only : DP
  USE ions_base, ONLY : nat
  use funct, only : dft_is_gradient, dft_is_nonlocc
  USE cell_base, ONLY : omega, alat
  use scf, only : rho, rho_core
  USE gvect, ONLY : g, ngm, nl
  USE fft_base, ONLY : dfftp
  USE noncollin_module, ONLY : nspin_lsda, nspin_gga, nspin_mag
  USE dynmat, ONLY : dyn, dyn_rec
  USE modes,  ONLY : nirr, npert
  USE uspp,   ONLY : nlcc_any

  USE mp_bands,  ONLY: intra_bgrp_comm
  USE mp,        ONLY: mp_sum

  USE qpoint,  ONLY : xq
  USE eqv,     ONLY : dmuxc
  USE gc_lr,   ONLY: grho,  dvxc_rr,  dvxc_sr,  dvxc_ss, dvxc_s

  implicit none

  integer :: imode0, npe
  ! input: the starting mode
  ! input: the number of perturbations
  ! input: the change of density due to perturbation

  complex(DP) :: drhoscf (dfftp%nnr, nspin_mag, npe)

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
  complex(DP), external :: zdotc
  ! the scalar product function


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
! add core charge to the density
!
  DO is=1,nspin_lsda
     rho%of_r(:,is) = rho%of_r(:,is) + fac * rho_core(:)
  ENDDO
!
!  Compute the change of xc potential due to the perturbation
!
  do ipert = 1, npe
     mode = imode0 + ipert
     dvaux (:,:) = (0.d0, 0.d0)
     call addcore (mode, drhoc)
     do is = 1, nspin_lsda
        call daxpy (2 * dfftp%nnr, fac, drhoc, 1, drhoscf (1, is, ipert), 1)
     enddo
     do is = 1, nspin_lsda
        do is1 = 1, nspin_mag
           do ir = 1, dfftp%nnr
              dvaux (ir, is) = dvaux (ir, is) + dmuxc (ir, is, is1) * &
                                                drhoscf ( ir, is1, ipert)
           enddo
        enddo
     enddo
     !
     ! add gradient correction to xc, NB: if nlcc is true we need to add here
     ! its contribution. grho contains already the core charge
     !
     if ( dft_is_gradient() ) &
       call dgradcorr (rho%of_r, grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s, xq, &
          drhoscf (1, 1, ipert), dfftp%nnr, nspin_mag, nspin_gga, nl, ngm, g, alat,&
          dvaux)
     if (dft_is_nonlocc()) &
       call dnonloccorr(rho%of_r, drhoscf (1, 1, ipert), xq, dvaux)

     do is = 1, nspin_lsda
        call daxpy (2 * dfftp%nnr, - fac, drhoc, 1, drhoscf (1, is, ipert), 1)
     enddo
     mode1 = 0
     do irr = 1, nirr
        do jpert = 1, npert (irr)
           mode1 = mode1 + 1
           call addcore (mode1, drhoc)
           do is = 1, nspin_lsda
              dyn1 (mode, mode1) = dyn1 (mode, mode1) + &
                   zdotc (dfftp%nnr, dvaux (1, is), 1, drhoc, 1) * &
                   omega * fac / DBLE (nrtot)
           enddo
        enddo
     enddo
  enddo
  DO is=1,nspin_lsda
     rho%of_r(:,is) = rho%of_r(:,is) - fac * rho_core(:)
  ENDDO
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
