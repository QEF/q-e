!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine zstar_eu
  !-----------------------------------------------------------------------
  ! calculate the effective charges Z(E,Us) (E=scf,Us=bare)
  !
  ! epsil =.true. is needed for this calculation to be meaningful
  !
#include "machine.h"

  USE io_global,      ONLY : stdout
  use pwcom
  USE wavefunctions,  ONLY: evc
  use parameters, only : DP
  use phcom
  implicit none

  integer :: ibnd, ipol, jpol, icart, na, nu, mu, imode0, irr, &
       imode, nrec, mode, ik
  ! counter on bands
  ! counter on polarization
  ! counter on cartesian coordinates
  ! counter on atoms and modes
  ! counter on modes
  ! counter on modes
  ! counter on records
  ! counter on modes
  ! counter on k points

  real(kind=DP) :: work (3, 3, nat), weight
  !  auxiliary space
  !

  complex(kind=DP) :: ZDOTC
  !  scalar product
  call start_clock ('zstar_eu')

  call setv (2 * 9 * nat, 0.d0, zstareu0, 1)
  call setv (9 * nat, 0.d0, zstareu, 1)

  if (nksq.gt.1) rewind (iunigk)
  do ik = 1, nksq
     if (nksq.gt.1) read (iunigk) npw, igk
     npwq = npw
     weight = wk (ik)
     if (nksq.gt.1) call davcio (evc, lrwfc, iuwfc, ik, - 1)
     call init_us_2 (npw, igk, xk (1, ik), vkb)
     imode0 = 0
     do irr = 1, nirr
        do imode = 1, npert (irr)
           mode = imode+imode0
           dvpsi(:,:) = (0.d0, 0.d0)
           !
           ! recalculate  DeltaV*psi(ion) for mode nu
           !
           call dvqpsi_us (ik, mode, u (1, mode), .not.okvan)
           do jpol = 1, 3
              nrec = (jpol - 1) * nksq + ik
              !
              ! read DeltaV*psi(scf) for electric field in jpol direction
              !
              call davcio (dpsi, lrdwf, iudwf, nrec, - 1)
              do ibnd = 1, nbnd
                 zstareu0(jpol,mode)=zstareu0(jpol, mode)-2.d0*weight*&
                      ZDOTC(npw,dpsi(1,ibnd),1,dvpsi(1,ibnd),1)
              enddo
           enddo
        enddo
        imode0 = imode0 + npert (irr)
     enddo
  enddo
  !
  ! Now we add the terms which are due to the USPP
  !
  if (okvan) call zstar_eu_us

#ifdef __PARA
  call reduce (18 * nat, zstareu0)
  call poolreduce (18 * nat, zstareu0)
#endif
  !
  ! bring the mode index to cartesian coordinates
  !
  do jpol = 1, 3
     do mu = 1, 3 * nat
        na = (mu - 1) / 3 + 1
        icart = mu - 3 * (na - 1)
        do nu = 1, 3 * nat
           zstareu (jpol, icart, na) = zstareu (jpol, icart, na) + conjg (u ( &
                mu, nu) ) * zstareu0 (jpol, nu)
        enddo
     enddo

  enddo
  call setv (9 * nat, 0.d0, work, 1)
  !
  ! bring to crystal axis for symmetrization
  ! NOTA BENE: the electric fields are already in crystal axis
  !
  do na = 1, nat
     do ipol = 1, 3
        do jpol = 1, 3
           do icart = 1, 3
              work (jpol, ipol, na) = work (jpol, ipol, na) + zstareu (jpol, &
                   icart, na) * at (icart, ipol)
           enddo
        enddo
     enddo
  enddo

  !      WRITE( stdout,'(/,10x,"Effective charges E-U in crystal axis ",/)')
  !      do na=1,nat
  !         WRITE( stdout,'(10x," atom ",i6)') na
  !         WRITE( stdout,'(10x,"(",3f15.5," )")') ((work(jpol,ipol,na),
  !     +                                ipol=1,3),jpol=1,3)
  !      enddo

  call symz (work, nsym, s, nat, irt)
  do na = 1, nat
     call trntns (work (1, 1, na), at, bg, 1)
  enddo

  call DCOPY (9 * nat, work, 1, zstareu, 1)
  !
  ! add the diagonal part
  !
  do ipol = 1, 3
     do na = 1, nat
        zstareu (ipol, ipol, na) = zstareu (ipol, ipol, na) + zv (ityp ( na) )
     enddo
  enddo

  WRITE( stdout, '(/,10x,"Effective charges E-U in cartesian axis ",/)')
  do na = 1, nat
     WRITE( stdout, '(10x," atom ",i6)') na
     WRITE( stdout, '(10x,"(",3f15.5," )")') ( (zstareu (ipol, jpol, na) &
          , ipol = 1, 3) , jpol = 1, 3)

  enddo
  call stop_clock ('zstar_eu')
  return
end subroutine zstar_eu
