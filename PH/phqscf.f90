!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
SUBROUTINE phqscf
  !-----------------------------------------------------------------------
  !
  !     This subroutine is the main driver of the self consistent cycle
  !     which gives as output the change of the wavefunctions and the
  !     change of the self-consistent potential due to a phonon of
  !     a fixed q or to an electric field.
  !

  USE io_global,  ONLY : stdout, ionode
  USE check_stop, ONLY: max_seconds
  USE pwcom
  USE kinds, ONLY : DP
  USE phcom

  IMPLICIT NONE

  INTEGER :: irr, irr1, irrc, imode0
  ! counter on the representations
  ! counter on the representations
  ! number of representation computed
  ! counter on the modes

  REAL(DP) :: tcpu, get_clock
  ! timing variables

  LOGICAL :: exst
  ! used to test the recover file

  EXTERNAL get_clock
  ! the change of density due to perturbations

  COMPLEX(DP), ALLOCATABLE :: drhoscf (:,:,:)

  CALL start_clock ('phqscf')
  !
  !    For each irreducible representation we compute the change
  !    of the wavefunctions
  !
  irrc = 0

  ALLOCATE (drhoscf( nrxx , nspin, npertx))    
  DO irr = 1, nirr
     IF ( (comp_irr (irr) == 1) .AND. (done_irr (irr) == 0) ) THEN
        irrc = irrc + 1
        imode0 = 0
        DO irr1 = 1, irr - 1
           imode0 = imode0 + npert (irr1)
        ENDDO
        IF (npert (irr) == 1) THEN
           WRITE( stdout, '(//,5x,"Representation #", i3," mode # ",i3)') &
                              irr, imode0 + 1
        ELSE
           WRITE( stdout, '(//,5x,"Representation #", i3," modes # ",4i3)') &
                              irr, (imode0+irr1, irr1=1,npert(irr))
        ENDIF
        !
        !    then for this irreducible representation we solve the linear system
        !
        WRITE( stdout, '(/,5x,"Self-consistent Calculation")')
        CALL solve_linter (irr, imode0, npert (irr), drhoscf)
        WRITE( stdout, '(/,5x,"End of self-consistent calculation")')
        !
        !   Add the contribution of this mode to the dynamical matrix
        !
        IF (convt) THEN
           CALL drhodv (imode0, npert (irr), drhoscf)
           !
           !   add the contribution of the modes imode0+1 -> imode+npe
           !   to the effective charges Z(Us,E) (Us=scf,E=bare)
           !
           IF (zue) CALL add_zstar_ue (imode0, npert (irr) )
           IF (zue.AND. okvan) CALL add_zstar_ue_us(imode0, npert (irr) )

        ENDIF
        IF (convt) THEN
           WRITE( stdout, '(/,5x,"Convergence has been achieved ")')
           done_irr (irr) = 1
        ELSE
           WRITE( stdout, '(/,5x,"No convergence has been achieved ")')
           CALL stop_ph (.FALSE.)
        ENDIF
        !
        tcpu = get_clock ('PHONON')
        ! if (tcpu > max_second) then
        ! temporary disabled: recover does not work if program stop here
           !
           ! WRITE( stdout, '(/,5x,"Stopping for time limit ",2f10.0)') &
           !      tcpu, max_second
           ! CALL stop_ph (.FALSE.)
        ! ENDIF
        !
        !   We test here if we have done the appropriate number of
        !   representation
        !
        IF (irrc >= maxirr) THEN
           WRITE( stdout, '(/,5x,"Stopping at Representation #",i6)') irr
#ifdef DEBUG
           IF ( ionode ) CLOSE (6)
#endif
           CALL stop_ph (.FALSE.)
        ENDIF
     ENDIF

  ENDDO
  DEALLOCATE (drhoscf)

  CALL stop_clock ('phqscf')
  RETURN
END SUBROUTINE phqscf
