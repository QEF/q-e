!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE gen_dwfc (isw_sl)
  !-----------------------------------------------------------------------
  !
  !    Calculates and writes  | d/du(0) psi(k+q) >
  !
  !    Several cases are possible:
  !  isw_sl = 1   : it calculates | d/du(q) psi_k >
  !  isw_sl = 2   : it calculates | d/du(0) psi_k+q >
  !  isw_sl = 3,4 : it calculates | d/du(0) psi_k >
  !
  USE io_global,  ONLY : stdout, ionode
  USE pwcom
  USE phcom
  USE d3com
  !
  IMPLICIT NONE
  !
  INTEGER isw_sl, nirr_x, irr, irr1, imode0
  ! switch
  ! the number of irreducible representation
  ! counter on the representations
  ! counter on the representations
  ! counter on the modes
  INTEGER, POINTER ::  npert_x (:)
  ! the number of perturbations per IR

  IF (isw_sl.EQ.1) THEN
     nirr_x = nirr
     npert_x => npert
  ELSE
     nirr_x = nirrg0
     npert_x => npertg0
  ENDIF
  !
  !    For each irreducible representation we compute the change
  !    of the wavefunctions
  !
  DO irr = 1, nirr_x
     imode0 = 0
     DO irr1 = 1, irr - 1
        imode0 = imode0 + npert_x (irr1)
     ENDDO
     IF (npert_x (irr) .EQ.1) THEN
        WRITE( stdout, '(//,5x,"Representation #", i3, &
             &                        " mode # ",i3)') irr, imode0 + 1
     ELSE
        WRITE( stdout, '(//,5x,"Representation #", i3, &
             &                 " modes # ",3i3)') irr,  (imode0 + irr1, irr1 = &
             & 1, npert_x (irr) )

     ENDIF
     CALL solve_linter_d3 (irr, imode0, npert_x (irr), isw_sl)
  ENDDO
  !
  ! Writes FermiEnergy shift on a file
  !
  IF ( ionode ) THEN
     !
     IF (isw_sl.EQ.3.AND.degauss.NE.0.d0) THEN
        REWIND (unit = iuef)
        WRITE (iuef) ef_sh
     ENDIF
     !
  END IF
  !
  ! closes and opens some units --useful in case of interrupted run--
  !

  CALL close_open (isw_sl)
  RETURN
END SUBROUTINE gen_dwfc
