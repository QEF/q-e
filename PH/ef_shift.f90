!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine ef_shift (drhoscf, ldos, ldoss, dos_ef, irr, npe, flag)
!-----------------------------------------------------------------------
!    This routine takes care of the effects of a shift of Ef, due to the
!    perturbation, that can take place in a metal at q=0
!
#include "machine.h"


use pwcom
use parameters, only : DP
use phcom
implicit none
!
! input/output variables
!

integer :: npe
                             ! input: the number of perturbation

complex(kind=DP) :: drhoscf (nrxx, nspin, npe), ldos (nrxx, nspin), &
 ldoss (nrxxs, nspin)
                                      ! inp/out:the change of the charge
                                    ! inp: local DOS at Ef
                                    ! inp: local DOS at Ef without augme

real(kind=DP) :: dos_ef
                              ! inp: density of states at Ef

integer :: irr
                       ! inp: index of the current irr. rep.
logical :: flag
                      ! inp: if true the eigenfunctions are updated
!
! local variables
!
                 !--> these quantities may be complex since perturbation

complex(kind=DP) :: delta_n, def (3), wfshift
                       ! the change in electron number
                       ! the change of the Fermi energy for each pert.
                       ! the shift coefficient for the wavefunction

real(kind=DP) :: w0gauss
                       ! the smeared delta function

integer :: ibnd, ik, is, ipert, nrec, ikrec
                       ! counter on occupied bands
                       ! counter on k-point
                       ! counter on spin polarizations
                       ! counter on perturbations
                       ! record number
                       ! record position of wfc at k
save def

external w0gauss
!
! determines Fermi energy shift (such that each pertubation is neutral)
!
call start_clock ('ef_shift')
if (.not.flag) then
   write (6, * )
   do ipert = 1, npert (irr)
   delta_n = (0.d0, 0.d0)
   do is = 1, nspin
   call cft3 (drhoscf (1, is, ipert), nr1, nr2, nr3, nrx1, nrx2, &
    nrx3, - 1)
   if (gg (1) .lt.1.0d-8) delta_n = delta_n + omega * drhoscf (nl &
    (1), is, ipert)
   call cft3 (drhoscf (1, is, ipert), nr1, nr2, nr3, nrx1, nrx2, &
    nrx3, + 1)
   enddo
   call reduce (2, delta_n)
   def (ipert) = - delta_n / dos_ef
!check
!         write (6,*) 'delta_n , dos_ef, def(ipert)'
!scal         write (6,*) 'delta_n , dos_ef, def(ipert)'
!         write (6,*)  delta_n , dos_ef, def(ipert)
!check
   enddo
!
! symmetrizes the Fermi energy shift
!
   call sym_def (def, irr)
write (6, '(5x,"Pert. #",i3, &
&                   ": Fermi energy shift (Ryd) =", &
&         2f10.4)')  (ipert, def (ipert) , ipert = 1, npert (irr) )
!
! corrects the density response accordingly...
!
   do ipert = 1, npert (irr)
   call ZAXPY (nrxx * nspin, def (ipert), ldos, 1, drhoscf (1, 1, &
    ipert), 1)
!check
!         delta_n= (0.d0,0.d0)
!         do is=1,nspin
!            call cft3(drhoscf(1,is,ipert),nr1,nr2,nr3,nrx1,nrx2,nrx3,-1
!            if (gg(1).lt.1.0d-8)
!     +         delta_n = delta_n + omega*drhoscf(nl(1),is,ipert)
!            call cft3(drhoscf(1,is,ipert),nr1,nr2,nr3,nrx1,nrx2,nrx3,+1
!         end do
!         call reduce(2,delta_n)
!         write (6,*) 'new delta_n , dos_ef, def(ipert)',nd_nmbr
!         write (6,*)      delta_n , dos_ef, def(ipert)
!check
   enddo
else
!
! does the same for perturbed wfc
!
   do ik = 1, nksq
   npw = ngk (ik)
!
! reads unperturbed wavefuctions psi_k in G_space, for all bands
!
   ikrec = ik
   if (nksq.gt.1) call davcio (evc, lrwfc, iuwfc, ikrec, - 1)
!
! reads delta_psi from iunit iudwf, k=kpoint
!
   do ipert = 1, npert (irr)
   nrec = (ipert - 1) * nksq + ik
   if (nksq.gt.1.or.npert (irr) .gt.1) call davcio (dpsi, lrdwf, &
    iudwf, nrec, - 1)
   do ibnd = 1, nbnd_occ (ik)
   wfshift = 0.5d0 * def (ipert) * w0gauss ( (ef - et (ibnd, ik) ) &
    / degauss, ngauss) / degauss
   call ZAXPY (npw, wfshift, evc (1, ibnd), 1, dpsi (1, ibnd), &
    1)
   enddo
!
! writes corrected delta_psi to iunit iudwf, k=kpoint,
!
   if (nksq.gt.1.or.npert (irr) .gt.1) call davcio (dpsi, lrdwf, &
    iudwf, nrec, + 1)
   enddo

   enddo
   do ipert = 1, npert (irr)
   do is = 1, nspin
   call ZAXPY (nrxxs, def (ipert), ldoss (1, is), 1, drhoscf (1, &
    is, ipert), 1)
   enddo
   enddo
endif
call stop_clock ('ef_shift')
return
end subroutine ef_shift
