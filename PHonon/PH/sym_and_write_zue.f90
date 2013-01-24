!
! Copyright (C) 2001-2010 Quantum ESPRESSO group
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
  !
  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : nat, zv, atm, ityp
  USE io_global,  ONLY : stdout
  USE cell_base,  ONLY : at, bg
  USE symme,      ONLY : symtensor
  USE efield_mod, ONLY : zstarue, zstarue0
  USE modes,      ONLY : u
  USE ph_restart, ONLY : ph_writefile
  USE control_ph, ONLY : zue, done_zue, xmldyn
  USE units_ph,   ONLY : iudyn

  implicit none

  integer :: ipol, jpol, icart, jcart, na, nu, mu, ierr
  ! counter on polarization
  ! counter on cartesian coordinates
  ! counter on atoms and modes
  ! counter on modes

  real(DP) :: work (3, 3, nat)
  ! auxiliary space (note the order of indices)
  !
  IF (.NOT.zue.OR.done_zue) RETURN

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
  ! copy to work (a vector with E-U index order) and transform to
  ! cartesian axis (NOTA BENE: the E index is in crystal axis)
  !
  work(:,:,:) = 0.d0
  do jcart = 1, 3
     do icart = 1, 3
        work (jcart,icart,:) = zstarue(icart,:,1) * bg(jcart,1) + &
                               zstarue(icart,:,2) * bg(jcart,2) + &
                               zstarue(icart,:,3) * bg(jcart,3)
     enddo
  enddo
  !
  ! symmetrize
  !
  call symtensor (nat, work)
  !
  ! back to U-E ordering
  !
  do icart = 1, 3
     do jcart = 1, 3
        zstarue (icart, :, jcart) = work (jcart, icart, :)
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
  IF (.NOT. xmldyn) THEN
     write (iudyn, '(/5x, &
       &               "Effective Charges U-E: Z_{s,alpha}{beta}",/)')
     do na = 1, nat
        write (iudyn, '(5x,"atom # ",i4)') na
        write (iudyn, '(3e24.12)') ( (zstarue (ipol, na, jpol) , jpol = 1, &
          3) , ipol = 1, 3)
     enddo
  ENDIF
  !
  ! write Z_{s,alpha}{beta} on standard output
  !
  done_zue=.true.
  CALL summarize_zue()
  CALL ph_writefile('tensors', 0, 0, ierr)

  return
end subroutine sym_and_write_zue
