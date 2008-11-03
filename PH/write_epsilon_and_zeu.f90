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
  USE io_global,  ONLY : stdout
  USE constants,  ONLY : fpi, bohr_radius_angs
  USE kinds, only : DP
  USE cell_base, only : omega
  USE ions_base, only : ityp, atm
  USE control_ph, ONLY : lgamma_gamma
  implicit none
  ! input variables
  integer :: iudyn, nat
  ! unit number
  ! number of atom in the unit cell

  real(DP) :: zstareu (3, 3, nat), epsilon (3, 3), chi(3, 3)
  !  the effective charges
  !  the dielectric tensor
  ! local variables
  integer :: na, icar, jcar
  ! counter on atoms
  ! cartesian coordinate counters
  !
  ! write dielectric tensor and Z(E,Us) effective charges on iudyn
  !
  write (iudyn, '(/,5x,"Dielectric Tensor:",/)')
  write (iudyn, '(3f24.12)') ((epsilon(icar,jcar), jcar=1,3), icar=1,3)
  write (iudyn, '(/5x, "Effective Charges E-U: Z_{alpha}{s,beta}",/)')
  do na = 1, nat
     write (iudyn, '(5x,"atom # ",i4)') na
     write (iudyn, '(3f24.12)') ((zstareu(icar,jcar,na), jcar=1,3), icar=1,3)
  enddo
  !
  ! write dielectric tensor and Z(E,Us) effective charges on standard output
  !
  WRITE( stdout, '(/,10x,"Dielectric constant in cartesian axis ",/)')

  WRITE( stdout, '(10x,"(",3f15.5," )")') &
                                   ((epsilon(icar,jcar), jcar=1,3), icar=1,3)

  IF (lgamma_gamma) THEN
!
! The system is probably a molecule. Try to estimate the polarizability
!
     DO icar=1,3
        DO jcar=1,3
           IF (icar == jcar) THEN
              chi(icar,jcar) = (epsilon(icar,jcar)-1.0_DP)*omega/fpi
           ELSE
              chi(icar,jcar) = epsilon(icar,jcar)*omega/fpi
           END IF
        END DO
     END DO
     WRITE(stdout,'(/5x,"Polarizability (a.u.)^3",20x,"polarizability (A^3)")')
     WRITE(stdout,'(3f10.2,5x,3f14.4)') ( (chi(icar,jcar), jcar=1,3), &
                   (chi(icar,jcar)*BOHR_RADIUS_ANGS**3, jcar=1,3), icar=1,3)
  ENDIF



  WRITE( stdout, '(/,10x,"Effective charges (d Force / dE) in cartesian axis",/)')
  ! WRITE( stdout, '(10x,  "          Z_{alpha}{s,beta} ",/)')
  do na = 1, nat
     WRITE( stdout, '(10x," atom ",i6,a6)') na, atm(ityp(na))
     WRITE( stdout, '(6x,"Ex  (",3f15.5," )")') &
                            (zstareu(1,jcar,na), jcar=1,3)
     WRITE( stdout, '(6x,"Ey  (",3f15.5," )")') &
                            (zstareu(2,jcar,na), jcar=1,3)
     WRITE( stdout, '(6x,"Ez  (",3f15.5," )")') &
                            (zstareu(3,jcar,na), jcar=1,3)
  enddo
  return
end subroutine write_epsilon_and_zeu
