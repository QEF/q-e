!
! Copyright (C) 2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE check_if_partial_dyn(u, nirr, npert, comp_irr)
!
!  This routine decides which irreducible representation have to
!  be computed for each q point on the basis of start_irr and last_irr
!  or on the basis of nat_todo. It sets the array comp_irr. 
!  If this is part of a dispersion calculation the routine has to be
!  called separately for each q and in that case different displacement
!  patterns are given as input. Note that this routine is called before
!  distributing the irrep among the images and only the irrep selected
!  here are distributed.
!
USE kinds, ONLY : DP
USE ions_base, ONLY : nat
USE symm_base, ONLY : irt
USE partial, ONLY : nat_todo, atomo
USE control_ph, ONLY : start_irr, last_irr, ldiag
USE control_flags, ONLY : modenum

USE lr_symm_base, ONLY : nsymq

IMPLICIT NONE

COMPLEX(DP), INTENT(IN) :: u( 3*nat, 3*nat )
INTEGER, INTENT(IN) :: nirr, npert( 3*nat )
LOGICAL, INTENT(OUT) :: comp_irr(0:3*nat)

INTEGER, ALLOCATABLE :: ifat(:)

INTEGER :: na, isym, mu, nu, ipol, imode0, irr, ipert
INTEGER :: last_irr_eff
!
!  If modenum is specified, only this mode is calculated and nat_todo is
!  ignored
!
  comp_irr = .FALSE.
  comp_irr(0) = .TRUE.
  IF (modenum /= 0) THEN
     comp_irr(modenum)=.TRUE.
     RETURN
  ENDIF

  ALLOCATE ( ifat(nat) )
  IF (nat_todo > 0) THEN
     !
     !   Sets the atoms which must be computed: the requested
     !   atoms and all the symmetry related atoms
     !
     ifat = 0
     DO na = 1, nat_todo
        IF ( atomo(na)>nat .OR. atomo(na)<1 ) &
           CALL errore('phq_setup', &
                       'one of atoms to do (nat_todo) is < 0 or > nat', 1)
        ifat (atomo (na) ) = 1
        DO isym = 1, nsymq
           ifat (irt (isym, atomo (na) ) ) = 1
        ENDDO
     ENDDO
     !
     !  Find the irreducible representations where the required atoms moves
     !
     imode0 = 0
     do irr = 1, nirr
        do ipert = 1, npert (irr)
           mu = imode0 + ipert
           do na = 1, nat
              if (ifat (na) == 1 .and. .NOT.comp_irr (irr) ) then
                 do ipol = 1, 3
                    nu = 3 * (na - 1) + ipol
                    if (abs (u (nu, mu) ) > 1.d-6)  comp_irr (irr) = .TRUE.
                 enddo
              endif
           enddo
        enddo
        imode0 = imode0 + npert (irr)
     enddo
  ELSE
     !
     !  nat_todo is not given. In principle all representation have
     !  to be calculated if not found somewhere else or limited below
     !
     comp_irr=.TRUE.
  ENDIF
!
! The representations that are smaller than start_irr or larger than
! last_irr_eff must be removed
!
  last_irr_eff=last_irr
  IF (last_irr > nirr.or.last_irr<0) last_irr_eff=nirr
  IF (start_irr > 1) comp_irr(0:start_irr-1) = .FALSE.
  IF (last_irr_eff < nirr ) comp_irr(last_irr_eff+1:nirr) = .FALSE.

!
! When ldiag=.true. the partial dynamical matrix is diagonalized.
! It must contain also the part calculated by dynmat0
!  
  IF (ldiag) comp_irr(0)=.TRUE.

  DEALLOCATE(ifat)

  RETURN
END SUBROUTINE check_if_partial_dyn
