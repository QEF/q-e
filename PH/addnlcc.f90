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
#include "machine.h"

  use funct
  use pwcom
  use parameters, only : DP
  use phcom
  implicit none

  integer :: imode0, npe
  ! input: the starting mode
  ! input: the number of perturbations
  ! input: the change of density due to perturbatio

  complex(kind=DP) :: drhoscf (nrxx, nspin, npertx)

  integer :: nrtot, ipert, jpert, is, is1, irr, ir, mode, mode1
  ! the total number of points
  ! counter on perturbations
  ! counter on spin
  ! counter on representations
  ! counter on real space points
  ! counter on modes

  complex(kind=DP) :: ZDOTC, dyn1 (3 * nat, 3 * nat)
  complex(kind=DP), allocatable :: drhoc (:),&
       dvaux (:,:)
  ! the scalar product function
  ! auxiliary dynamical matrix
  ! the change of the core
  ! the change of the potential

  real(kind=DP) :: fac
  ! auxiliary factor


  if (.not.nlcc_any) return

  allocate (drhoc(  nrxx))    
  allocate (dvaux(  nrxx , nspin))    

  call setv (2 * 3 * nat * 3 * nat, 0.d0, dyn1, 1)
  !
  !  compute the exchange and correlation potential for this mode
  !
  nrtot = nr1 * nr2 * nr3
  fac = 1.d0 / float (nspin)
  do ipert = 1, npe

     mode = imode0 + ipert
     call setv (2 * nrxx * nspin, 0.d0, dvaux, 1)
     call addcore (mode, drhoc)
     do is = 1, nspin
        call DAXPY (nrxx, fac, rho_core, 1, rho (1, is), 1)
        call DAXPY (2 * nrxx, fac, drhoc, 1, drhoscf (1, is, ipert), &
             1)

     enddo
     do is = 1, nspin
        do is1 = 1, nspin
           do ir = 1, nrxx
              dvaux (ir, is) = dvaux (ir, is) + dmuxc (ir, is, is1) * drhoscf ( &
                   ir, is1, ipert)
           enddo
        enddo

     enddo
     !
     ! add gradient correction to xc, NB: if nlcc is true we need to add here
     ! its contribution. grho contains already the core charge
     !

     if (igcx.ne.0.or.igcc.ne.0) &
       call dgradcorr (rho, grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s, xq, &
          drhoscf (1, 1, ipert), nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, &
          nspin, nl, ngm, g, alat, omega, dvaux)
     do is = 1, nspin
        call DAXPY (nrxx, - fac, rho_core, 1, rho (1, is), 1)
        call DAXPY (2 * nrxx, - fac, drhoc, 1, drhoscf (1, is, ipert), &
             1)

     enddo
     mode1 = 0
     do irr = 1, nirr
        do jpert = 1, npert (irr)
           mode1 = mode1 + 1
           call addcore (mode1, drhoc)
           do is = 1, nspin
              dyn1 (mode, mode1) = dyn1 (mode, mode1) + ZDOTC (nrxx, dvaux (1, &
                   is), 1, drhoc, 1) * omega * fac / float (nrtot)
           enddo
        enddo
     enddo
  enddo
#ifdef PARA
  !
  ! collect contributions from all pools (sum over k-points)
  !
  call reduce (18 * nat * nat, dyn1)
#endif

  call ZAXPY (9 * nat * nat, (1.d0, 0.d0), dyn1, 1, dyn, 1)
  deallocate (dvaux)
  deallocate (drhoc)
  return
end subroutine addnlcc
