!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine dielec
!-----------------------------------------------------------------------
!
!      calculates the dielectric tensor
!
#include "machine.h"


use pwcom
use parameters, only : DP
use phcom
implicit none

integer :: ibnd, ipol, jpol, nrec, ik
                        ! counter on polarizations
                       ! counter on records
                       ! counter on k points
real(kind=DP) :: w, weight

complex(kind=DP) :: ZDOTC
call setv (9, 0.0d0, epsilon, 1)
if (nksq.gt.1) rewind (unit = iunigk)
do ik = 1, nksq
   if (nksq.gt.1) read (iunigk) npw, igk
   weight = wk (ik)
   w = fpi * weight / omega
   do ipol = 1, 3
      nrec = (ipol - 1) * nksq + ik
      call davcio (dvpsi, lrbar, iubar, nrec, - 1)
      do jpol = 1, 3
         nrec = (jpol - 1) * nksq + ik
         call davcio (dpsi, lrdwf, iudwf, nrec, - 1)
         do ibnd = 1, nbnd_occ (ik)
!
!  this is the real part of <DeltaV*psi(E)|DeltaPsi(E)>
!
            epsilon(ipol,jpol)=epsilon(ipol,jpol)-4.d0*w*DREAL( &
                    ZDOTC (npw, dvpsi (1, ibnd), 1, dpsi (1, ibnd), 1) )
         enddo
      enddo
   enddo
enddo
#ifdef PARA
call reduce (9, epsilon)

call poolreduce (9, epsilon)
#endif
!
!      symmetrize
!
!       write(6,'(/,10x,"Unsymmetrized in crystal axis ",/)')
!       write(6,'(10x,"(",3f15.5," )")') ((epsilon(ipol,jpol),
!     +                                ipol=1,3),jpol=1,3)


call symtns (epsilon, nsym, s)
!
!    pass to cartesian axis
!
!      write(6,'(/,10x,"Symmetrized in crystal axis ",/)')
!      write(6,'(10x,"(",3f15.5," )")') ((epsilon(ipol,jpol),
!     +                                ipol=1,3),jpol=1,3)
call trntns (epsilon, at, bg, 1)
!
! add the diagonal part
!
do ipol = 1, 3
   epsilon (ipol, ipol) = epsilon (ipol, ipol) + 1.d0
enddo
!
!  and print the result
!
write (6, '(/,10x,"Dielectric constant in cartesian axis ",/)')

write (6, '(10x,"(",3f18.9," )")') ((epsilon(ipol,jpol), ipol=1,3), jpol=1,3)
return
end subroutine dielec
