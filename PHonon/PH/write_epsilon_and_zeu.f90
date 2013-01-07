!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine write_epsilon_and_zeu (zstareu, epsilon, nat, iudyn)
  !-----------------------------------------------------------------------
  USE kinds, only : DP
  USE control_ph, ONLY : xmldyn
  USE io_global, ONLY : ionode
  implicit none
  ! input variables
  integer :: iudyn, nat
  ! unit number
  ! number of atom in the unit cell

  real(DP) :: zstareu (3, 3, nat), epsilon (3, 3)
  !  the effective charges
  !  the dielectric tensor
  ! local variables
  integer :: na, icar, jcar
  ! counter on atoms
  ! cartesian coordinate counters
  !
  ! write dielectric tensor and Z(E,Us) effective charges on iudyn
  !
  IF (.NOT.xmldyn.AND.ionode) THEN
     write (iudyn, '(/,5x,"Dielectric Tensor:",/)')
     write (iudyn, '(3f24.12)') ((epsilon(icar,jcar), jcar=1,3), icar=1,3)
     write (iudyn, '(/5x, "Effective Charges E-U: Z_{alpha}{s,beta}",/)')
     do na = 1, nat
        write (iudyn, '(5x,"atom # ",i4)') na
        write (iudyn, '(3f24.12)') ((zstareu(icar,jcar,na), jcar=1,3), icar=1,3)
     enddo
  ENDIF
  !
  ! write dielectric tensor and Z(E,Us) effective charges on standard output
  !
  CALL summarize_epsilon()

  CALL summarize_zeu()

  return
end subroutine write_epsilon_and_zeu
