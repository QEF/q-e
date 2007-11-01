!
! Copyright (C) 2001 PWSCF group
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
#include "f_defs.h"

  USE ions_base, ONLY : nat
  use funct, only : dft_is_gradient
  use pwcom
  use scf, only : rho, rho_core
  USE kinds, only : DP
  use phcom
  implicit none

  integer :: imode0, npe
  ! input: the starting mode
  ! input: the number of perturbations
  ! input: the change of density due to perturbatio

  complex(DP) :: drhoscf (nrxx, nspin, npertx)

  integer :: nrtot, ipert, jpert, is, is1, irr, ir, mode, mode1, nspin0, nspin1
  ! the total number of points
  ! counter on perturbations
  ! counter on spin
  ! counter on representations
  ! counter on real space points
  ! counter on modes

  complex(DP) :: ZDOTC, dyn1 (3 * nat, 3 * nat)
  ! the scalar product function
  ! auxiliary dynamical matrix
  complex(DP), allocatable :: drhoc (:), dvaux (:,:)
  ! the change of the core
  ! the change of the potential

  real(DP) :: fac
  ! auxiliary factor


  if (.not.nlcc_any) return

  nspin0=nspin
  nspin1=nspin
  if (nspin==4) then
     nspin0=1
     nspin1=1
     if (domag) nspin1=2
  endif

  allocate (drhoc(  nrxx))    
  allocate (dvaux(  nrxx , nspin))    

  dyn1 (:,:) = (0.d0, 0.d0)
  !
  !  compute the exchange and correlation potential for this mode
  !
  nrtot = nr1 * nr2 * nr3
  fac = 1.d0 / DBLE (nspin0)

  do ipert = 1, npe
     mode = imode0 + ipert
     dvaux (:,:) = (0.d0, 0.d0)
     call addcore (mode, drhoc)
     do is = 1, nspin0
        call DAXPY (nrxx, fac, rho_core, 1, rho%of_r(1, is), 1)
        call DAXPY (2 * nrxx, fac, drhoc, 1, drhoscf (1, is, ipert), 1)
     enddo
     do is = 1, nspin
        do is1 = 1, nspin
           do ir = 1, nrxx
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
          drhoscf (1, 1, ipert), nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, &
          nspin, nspin1, nl, ngm, g, alat, omega, dvaux)
     do is = 1, nspin0
        call DAXPY (nrxx, - fac, rho_core, 1, rho%of_r(1, is), 1)
        call DAXPY (2 * nrxx, - fac, drhoc, 1, drhoscf (1, is, ipert), 1)
     enddo
     mode1 = 0
     do irr = 1, nirr
        do jpert = 1, npert (irr)
           mode1 = mode1 + 1
           call addcore (mode1, drhoc)
           do is = 1, nspin0
              dyn1 (mode, mode1) = dyn1 (mode, mode1) + &
                   ZDOTC (nrxx, dvaux (1, is), 1, drhoc, 1) * &
                   omega * fac / DBLE (nrtot)
           enddo
        enddo
     enddo
  enddo
#ifdef __PARA
  !
  ! collect contributions from all pools (sum over k-points)
  !
  call reduce (18 * nat * nat, dyn1)
#endif
  dyn (:,:) = dyn(:,:) + dyn1(:,:) 
  deallocate (dvaux)
  deallocate (drhoc)
  return
end subroutine addnlcc
