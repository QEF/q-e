!------------------------------------------------
subroutine addnlcc_zstar_eu_us (drhoscf) 
!----------===================-------------------

#include "machine.h"

  use funct
  use pwcom
  use parameters, only : DP
  use phcom
  implicit none
  
  complex(kind=DP) :: drhoscf (nrxx,nspin,3)


  integer :: nrtot, ipert, jpert, is, is1, irr, ir, mode, mode1
  integer :: imode0, npe, ipol

  real(kind=DP) :: fac
  
  complex(kind=DP), dimension(nrxx) :: drhoc
  complex(kind=DP), dimension(nrxx,nspin) :: dvaux

  if (.not.nlcc_any) return
  
  do ipol = 1, 3
     imode0 = 0
     do irr = 1, nirr
        npe = npert(irr)
        !
        !  compute the exchange and correlation potential for this mode
        !
        nrtot = nr1 * nr2 * nr3
        fac = 1.d0 / float (nspin)
        do ipert = 1, npe
           mode = imode0 + ipert
           
           dvaux = (0.0_dp,0.0_dp)
           call addcore (mode, drhoc)
           
           do is = 1, nspin
              rho(:,is) = rho(:,is) + fac * rho_core
           end do

           do is = 1, nspin
              do is1 = 1, nspin
                 do ir = 1, nrxx
                    dvaux (ir, is) = dvaux (ir, is) +     &
                         dmuxc (ir, is, is1) *            &
                         drhoscf (ir, is1, ipol)
                 enddo
              enddo
           end do
           !
           ! add gradient correction to xc, NB: if nlcc is true we need to add here
           ! its contribution. grho contains already the core charge
           !

           if (igcx.ne.0.or.igcc.ne.0) &
                call dgradcorr (rho, grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s, xq, &
                drhoscf (1, 1, ipert), nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, &
                nspin, nl, ngm, g, alat, omega, dvaux)
        
           do is = 1, nspin
              rho(:,is) = rho(:,is) - fac * rho_core
           end do
           
           do is = 1, nspin
              zstareu0(ipol,mode) = zstareu0(ipol,mode) -                  &
                   omega * fac / real(nrtot, kind = dp) *         &
                   dot_product(dvaux(1:nrxx,is),drhoc(1:nrxx)) 
           end do
        end do
        imode0 = imode0 + npe
     end do
  end do

  return
end subroutine addnlcc_zstar_eu_us
