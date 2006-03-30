!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine polariz ( iw )
  !-----------------------------------------------------------------------
  !
  !      calculates the frequency depwndent polarizability 
  !
#include "f_defs.h"

  USE io_global,  ONLY : stdout
  USE io_files, ONLY: iunigk
  use pwcom
  USE kinds, only : DP
  use phcom
  use cell_base,  only : omega
  implicit none
  !
  ! I/O variables
  !
  real(kind=DP) :: iw
  !
  ! local variables
  !
  integer :: ibnd, ipol, jpol, nrec, ik
  ! counter on polarizations
  ! counter on records
  ! counter on k points
  real(kind=DP) :: w, weight

  complex(kind=DP) :: ZDOTC

  call start_clock ('polariz')
  epsilon(:,:) = 0.d0
  if (nksq > 1) rewind (unit = iunigk)
  do ik = 1, nksq
     if (nksq > 1) read (iunigk) npw, igk
     weight = wk (ik)
     w = fpi * weight / omega
     do ipol = 1, 3
        nrec = (ipol - 1) * nksq + ik
        call davcio (dvpsi, lrebar, iuebar, nrec, - 1)
        do jpol = 1, 3
           nrec = (jpol - 1) * nksq + ik
           call davcio (dpsi, lrdwf, iudwf, nrec, - 1)
           do ibnd = 1, nbnd_occ (ik)
              !
              !  this is the real part of <DeltaV*psi(E)|DeltaPsi(E)>
              !
              epsilon(ipol,jpol)=epsilon(ipol,jpol)-4.d0*w*REAL( &
                   ZDOTC (npw, dvpsi (1, ibnd), 1, dpsi (1, ibnd), 1) )
           enddo
        enddo
     enddo
  enddo
#ifdef __PARA
  call reduce (9, epsilon)
  call poolreduce (9, epsilon)
#endif
  !
  !      symmetrize
  !
  !       WRITE( stdout,'(/,10x,"Unsymmetrized in crystal axis ",/)')
  !       WRITE( stdout,'(10x,"(",3f15.5," )")') ((epsilon(ipol,jpol),
  !     +                                ipol=1,3),jpol=1,3)

  call symtns (epsilon, nsym, s)
  !
  !    pass to cartesian axis
  !
  !      WRITE( stdout,'(/,10x,"Symmetrized in crystal axis ",/)')
  !      WRITE( stdout,'(10x,"(",3f15.5," )")') ((epsilon(ipol,jpol),
  !     +                                ipol=1,3),jpol=1,3)
  call trntns (epsilon, at, bg, 1)
  !
  ! add the diagonal part
  !
  do ipol = 1, 3
     epsilon (ipol, ipol) = epsilon (ipol, ipol) + 1.d0
  enddo
  !
  ! compute the polarization
  !
  do ipol = 1, 3
     do jpol = 1, 3
        if ( epsilon (ipol, jpol) .gt. 1.d-4 ) &
        epsilon (ipol, jpol) = (3.d0*omega/fpi) * ( epsilon (ipol, jpol) - 1.d0 ) / &
                                                  ( epsilon (ipol, jpol) + 2.d0 )
     enddo
  enddo
  !
  !  and print the result
  !
  WRITE( stdout, '(/,10x,"Polarizability in cartesian axis at frequency ",f5.2,/)') iw

  WRITE( stdout, '(10x,"(",3f18.9," )")') ((epsilon(ipol,jpol), ipol=1,3), jpol=1,3)

  call stop_clock ('polariz')

  return
end subroutine polariz
