!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine sym_and_write_zue
  !-----------------------------------------------------------------------
  !
  !  symmetrize the effective charges in the U-E case (Us=scf,E=bare)
  !  and write them on iudyn and standard output
  !
#include "f_defs.h"
  !
  USE ions_base, ONLY : nat, zv, ityp
  USE io_global,  ONLY : stdout
  use pwcom
  USE kinds, only : DP
  use phcom
  USE mp_global,  ONLY : inter_pool_comm, intra_pool_comm
  USE mp,         ONLY : mp_sum

  implicit none

  integer :: ipol, jpol, icart, jcart, na, nu, mu
  ! counter on polarization
  ! counter on cartesian coordinates
  ! counter on atoms and modes
  ! counter on modes

  real(DP) :: work (3, 3, nat)
  ! auxiliary space (note the order of indices)
  !
#ifdef __PARA
  call mp_sum ( zstarue0, intra_pool_comm )
  call mp_sum ( zstarue0, inter_pool_comm )
#endif
  !
  zstarue(:,:,:) = 0.d0
  do jcart = 1, 3
     do mu = 1, 3 * nat
        na = (mu - 1) / 3 + 1
        icart = mu - 3 * (na - 1)
        do nu = 1, 3 * nat
           zstarue (icart, na, jcart) = zstarue (icart, na, jcart) + &
                u (mu, nu) * zstarue0 (nu, jcart)
        enddo
     enddo

  enddo
  !
  ! copy to work (a vector with E-U index order) and transform to crystal
  ! NOTA BENE: the E index is already in crystal axis
  !
  work(:,:,:) = 0.d0
  do na = 1, nat
     do icart = 1, 3
        do jcart = 1, 3
           work (jcart, icart, na) = at (1, icart) * zstarue (1, na, jcart) &
                                   + at (2, icart) * zstarue (2, na, jcart) &
                                   + at (3, icart) * zstarue (3, na, jcart)
        enddo
     enddo
  enddo
  !
  ! symmetrize
  !
  call symz (work, nsym, s, nat, irt)
  !
  ! back to cartesian axis and U-E ordering
  !
  do na = 1, nat
     call trntns (work (1, 1, na), at, bg, 1)
     do icart = 1, 3
        do jcart = 1, 3
           zstarue (icart, na, jcart) = work (jcart, icart, na)
        enddo
     enddo
  enddo
  !
  ! add the diagonal part
  !
  do ipol = 1, 3
     do na = 1, nat
        zstarue (ipol, na, ipol) = zstarue (ipol, na, ipol) + zv (ityp (na) )
     enddo
  enddo
  !
  ! write Z_{s,alpha}{beta} on iudyn
  !
  write (iudyn, '(/5x, &
       &               "Effective Charges U-E: Z_{s,alpha}{beta}",/)')
  do na = 1, nat
     write (iudyn, '(5x,"atom # ",i4)') na
     write (iudyn, '(3e24.12)') ( (zstarue (ipol, na, jpol) , jpol = 1, &
          3) , ipol = 1, 3)
  enddo
  !
  ! write Z_{s,alpha}{beta} on standard output
  !
  WRITE( stdout, '(/,10x,"Effective charges U-E in cartesian axis ",/)' &
       &)
  ! WRITE( stdout, '(10x,  "          Z_{s,alpha}{beta} ",/)')
  do na = 1, nat
     WRITE( stdout, '(10x," atom ",i6)') na
     WRITE( stdout, '(10x,"(",3f15.5," )")') ( (zstarue (ipol, na, jpol) &
          , jpol = 1, 3) , ipol = 1, 3)

  enddo
  return
end subroutine sym_and_write_zue
