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
#include "f_defs.h"
  !
  USE ions_base, ONLY : nat, zv, ityp
  USE io_global,      ONLY : stdout
  USE io_files, ONLY: iunigk
  use pwcom
  USE wavefunctions_module,  ONLY: evc
  USE kinds, only : DP
  use phcom
  implicit none

  integer :: ibnd, ipol, jpol, icart, na, nu, mu, imode0, irr, &
       imode, nrec, mode, ik
  ! counters
  real(DP) :: work (3, 3, nat), weight
  !  auxiliary space
  complex(DP), external :: ZDOTC
  !  scalar product
  !
  call start_clock ('zstar_eu')

  zstareu0(:,:) = (0.d0,0.d0)
  zstareu (:,:,:) = 0.d0

  if (nksq > 1) rewind (iunigk)
  do ik = 1, nksq
     if (nksq > 1) read (iunigk) npw, igk
     npwq = npw
     weight = wk (ik)
     if (nksq > 1) call davcio (evc, lrwfc, iuwfc, ik, - 1)
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
              do ibnd = 1, nbnd_occ(ik)
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
           zstareu (jpol, icart, na) = zstareu (jpol, icart, na) + &
                CONJG(u (mu, nu) ) * zstareu0 (jpol, nu)
        enddo
     enddo
  enddo
  !
  work(:,:,:) = 0.d0
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
  zstareu(:,:,:) = work(:,:,:)
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
          , jpol = 1, 3) , ipol = 1, 3)
  enddo
  call stop_clock ('zstar_eu')
  return
end subroutine zstar_eu
