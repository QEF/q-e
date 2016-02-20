!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE summarize_epsilon()
  !-----------------------------------------------------------------------
  !
  !      write the dielectric tensor on output
  !

  USE kinds, only : DP
  USE io_global,  ONLY : stdout
  USE constants, ONLY: fpi, bohr_radius_angs
  USE cell_base, ONLY: omega
  USE noncollin_module, ONLY : npol
  USE efield_mod, ONLY : epsilon
  USE control_ph, ONLY : lgamma_gamma, lnoloc, done_epsil
  USE control_lr, ONLY : lrpa

  IMPLICIT NONE

  INTEGER :: ipol, jpol
  ! counter on polarizations
  ! counter on records
  ! counter on k points
  REAL(DP) :: chi(3,3)
  !
  IF (.NOT. done_epsil) RETURN
  !
  IF (lnoloc) THEN
      WRITE( stdout, '(/,10x,"Dielectric constant in cartesian axis (DV_Hxc=0)",/)')
  ELSE IF (lrpa) THEN
     WRITE( stdout, '(/,10x,"RPA dielectric constant in cartesian axis (DV_xc=0)",/)')
  ELSE
      WRITE( stdout, '(/,10x,"Dielectric constant in cartesian axis ",/)')
  ENDIF

  WRITE( stdout, '(10x,"(",3f18.9," )")') ((epsilon(ipol,jpol), ipol=1,3), jpol=1,3)

  IF (lgamma_gamma) THEN
!
! The system is probably a molecule. Try to estimate the polarizability
!
     DO ipol=1,3
        DO jpol=1,3
           IF (ipol == jpol) THEN
              chi(ipol,jpol) = (epsilon(ipol,jpol)-1.0_DP)*3.0_DP*omega/fpi &
                               /(epsilon(ipol,jpol)+2.0_DP)
           ELSE
              chi(ipol,jpol) = epsilon(ipol,jpol)*omega/fpi
           END IF
        END DO
     END DO

     WRITE(stdout,'(/5x,"Polarizability (a.u.)^3",20x,"Polarizability (A^3)")')
     WRITE(stdout,'(3f10.2,5x,3f14.4)') ( (chi(ipol,jpol), jpol=1,3), &
                   (chi(ipol,jpol)*bohr_radius_angs**3, jpol=1,3), ipol=1,3)
  ENDIF

  RETURN
END SUBROUTINE summarize_epsilon
!
!-----------------------------------------------------------------------
SUBROUTINE summarize_zeu()
  !-----------------------------------------------------------------------
  !
  !  write the zue effective charges on output
  !
  USE kinds,     ONLY : DP
  USE ions_base, ONLY : nat, ityp, atm
  USE io_global, ONLY : stdout

  USE efield_mod,   ONLY : zstareu
  USE control_ph,   ONLY : done_zeu

  IMPLICIT NONE

  INTEGER :: jpol, na
  ! counters
  !

  IF (.NOT. done_zeu) RETURN

  WRITE( stdout, '(/,10x,"Effective charges (d Force / dE) in cartesian axis",/)')
  DO na = 1, nat
     WRITE( stdout, '(10x," atom ",i6, a6)') na, atm(ityp(na))
     WRITE( stdout, '(6x,"Ex  (",3f15.5," )")')  (zstareu (1, jpol, na), &
            jpol = 1, 3)
     WRITE( stdout, '(6x,"Ey  (",3f15.5," )")')  (zstareu (2, jpol, na), &
            jpol = 1, 3)
     WRITE( stdout, '(6x,"Ez  (",3f15.5," )")')  (zstareu (3, jpol, na), &
            jpol = 1, 3)
  ENDDO

  RETURN
END SUBROUTINE summarize_zeu

!-----------------------------------------------------------------------
SUBROUTINE summarize_zue
!-----------------------------------------------------------------------
!
!  Write the zue effective charges on output
!
  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : nat, atm, ityp
  USE io_global,  ONLY : stdout
  USE efield_mod, ONLY : zstarue
  USE control_ph, ONLY : done_zue

  IMPLICIT NONE

  INTEGER :: ipol, na
  ! counter on polarization
  ! counter on atoms
  !
  IF (.NOT. done_zue) RETURN

  WRITE( stdout, '(/,10x,"Effective charges (d P / du) in cartesian axis ",/)')
  !
  DO na = 1, nat
     WRITE( stdout, '(10x," atom ",i6,a6)') na, atm(ityp(na))
     WRITE( stdout, '(6x,"Px  (",3f15.5," )")') (zstarue (ipol, na, 1), &
            ipol = 1, 3)
     WRITE( stdout, '(6x,"Py  (",3f15.5," )")') (zstarue (ipol, na, 2), &
            ipol = 1, 3)
     WRITE( stdout, '(6x,"Pz  (",3f15.5," )")') (zstarue (ipol, na, 3), &
            ipol = 1, 3)
  ENDDO
  !
  RETURN
END SUBROUTINE summarize_zue
!
!-----------------------------------------------------------------------
SUBROUTINE summarize_elopt()
  !-----------------------------------------------------------------------
  !
  ! write the electro-optic tensor on output
  !
  USE io_global,  ONLY : stdout
  USE ramanm,     ONLY : eloptns, done_elop
IMPLICIT NONE

INTEGER :: ipa, ipb, ipc

IF (.NOT. done_elop) RETURN

WRITE(stdout, '(/,10x,''    Electro-optic tensor is defined as '')' )
WRITE(stdout, '(10x  ,''  the derivative of the dielectric tensor '')' )
WRITE(stdout, '(10x  ,''    with respect to one electric field '')' )
WRITE(stdout, '(10x  ,''       units are Rydberg a.u. '',/)' )
WRITE(stdout, '(10x  ,''  to obtain the static chi^2 multiply by 1/2  '',/)' )
WRITE(stdout, '(10x  ,''  to convert to pm/Volt multiply per 2.7502  '',/)' )
WRITE(stdout, '(/,10x,''Electro-optic tensor in cartesian axis: '',/)' )

DO ipc = 1, 3
   DO ipb = 1, 3
      WRITE(stdout,'(10x,''('',3f18.9,'' )'')')      &
                   (eloptns (ipa, ipb, ipc), ipa = 1, 3)
   ENDDO
   WRITE(6,'(10x)')
ENDDO

RETURN
END SUBROUTINE summarize_elopt

SUBROUTINE summarize_fpol ()
  !-----------------------------------------------------------------------

  USE kinds, only : DP
  USE freq_ph,    ONLY : nfs, done_iu

  IMPLICIT NONE
  ! input variables

  INTEGER :: iu
  ! counter on frequencies
  !
  DO iu=nfs, 1, -1
     IF ( done_iu(iu) ) CALL write_polariz(iu)
  END DO

  RETURN
END SUBROUTINE summarize_fpol
!
