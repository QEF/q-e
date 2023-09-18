!
! Copyright (C) 2004-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE beef_interface
    !----------------------------------------------------------------------
    !! This module contains fortran wrapper to BEEF library functions.
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    !
    PUBLIC :: beefx, beeflocalcorr, beeflocalcorrspin, beefsetmode, &
        beefrandinit, beefrandinitdef, beefensemble, beef_set_type
    !
#if !defined(__NOBEEF)
    INTERFACE
    !
    SUBROUTINE beefx( r, g, e, dr, dg, addlda ) BIND(C, NAME="beefx_")
        USE iso_c_binding
        !$acc routine seq
        REAL (C_DOUBLE)            :: r, g, e, dr, dg
        INTEGER(C_INT), INTENT(IN) :: addlda
    END SUBROUTINE beefx
    !
    SUBROUTINE beeflocalcorr( r, g, e, dr, dg, addlda) BIND(C, NAME="beeflocalcorr_")
        USE iso_c_binding
        !$acc routine seq
        REAL (C_DOUBLE), INTENT(INOUT) :: r, g, e, dr, dg
        INTEGER(C_INT), INTENT(IN) :: addlda
    END SUBROUTINE beeflocalcorr
    !
    SUBROUTINE beeflocalcorrspin(r, z, g, e, drup, drdown, dg, addlda) BIND(C, NAME="beeflocalcorrspin_")
        USE iso_c_binding
        !$acc routine seq
        REAL (C_DOUBLE), INTENT(INOUT) :: r, z, g, e, drup, drdown, dg
        INTEGER(C_INT), INTENT(IN) :: addlda
    END SUBROUTINE beeflocalcorrspin
    !
    SUBROUTINE beefsetmode(mode) BIND(C, NAME="beefsetmode_")
        USE iso_c_binding
        INTEGER(C_INT), INTENT(IN) :: mode
    END SUBROUTINE beefsetmode
    !
    SUBROUTINE beefrandinit(seed) BIND(C, NAME="beefrandinit_")
        USE iso_c_binding
        INTEGER(C_INT), INTENT(IN) :: seed
    END SUBROUTINE beefrandinit
    !
    SUBROUTINE beefrandinitdef() BIND(C, NAME="beefrandinitdef_")
    END SUBROUTINE beefrandinitdef
    !
    SUBROUTINE beefensemble(beefxc, ensemble) BIND(C, NAME="beefensemble_")
        USE iso_c_binding
        REAL (C_DOUBLE), INTENT(INOUT) :: beefxc(*), ensemble(*)
    END SUBROUTINE beefensemble
    !
    FUNCTION beef_set_type_interface(tbeef, ionode) &
            BIND(C,name="beef_set_type_") RESULT(r)
        USE iso_c_binding
        INTEGER(C_INT), INTENT(IN) :: tbeef, ionode
        INTEGER(C_INT)             :: r
    END FUNCTION beef_set_type_interface
    !
    END INTERFACE
    !
#endif
    !
    CONTAINS
    ! ====================================================================
    !
#if !defined(__NOBEEF)
    FUNCTION beef_set_type(tbeef, ionode) RESULT(r)
        INTEGER, INTENT(IN) :: tbeef
        LOGICAL, INTENT(IN) :: ionode
        LOGICAL             :: r
        ! ... local variables ...
        INTEGER             :: ionode_ = 0
        INTEGER             :: r_
        !
        IF ( ionode ) ionode_ = 1
        !
        r_ = beef_set_type_interface(tbeef, ionode_)
        !
        IF ( r_ /= 0 ) THEN
            r = .TRUE.
        ELSE
            r = .FALSE.
        END IF
        !
    END FUNCTION beef_set_type
    !
#else
    !
    ! Empty routines to prevent compilation errors.
    ! For simplicity the comments for Ford doc generator have been put here
    ! and not in the above interfaces.
    !
    SUBROUTINE beefx( r, g, e, dr, dg, addlda )
      !! Evaluate bee exchange energy and its derivatives \(d\epsilon/d\rho\) and 
      !! \( (d\epsilon/d|\nabla \rho| ) / |\nabla\rho|\).
      USE kind_l, ONLY : dp
      !$acc routine seq
      REAL (dp) :: r, g, e, dr, dg
      INTEGER :: addlda
    END SUBROUTINE beefx
    !
    SUBROUTINE beeflocalcorr( r, g, e, dr, dg, addlda)
      !! Evaluate local part of bee correlation and its derivatives \(d\epsilon/drho\)
      !! and \( (d\epsilon/d|\nabla\rho|) / |\nabla\rho| \).
      USE kind_l, ONLY : dp
      !$acc routine seq
      REAL (dp), INTENT(INOUT) :: r, g, e, dr, dg
      INTEGER :: addlda
    END SUBROUTINE beeflocalcorr
    !
    SUBROUTINE beeflocalcorrspin(r, z, g, e, drup, drdown, dg, addlda)
      !! Evaluate local part of bee correlation for spin polarized system.
      USE kind_l, ONLY : dp
      !$acc routine seq
      REAL (dp), INTENT(INOUT) :: r, z, g, e, drup, drdown, dg
      INTEGER :: addlda
    END SUBROUTINE beeflocalcorrspin
    !
    SUBROUTINE beefsetmode(mode)
      !! \(\text{mode} \geq 0\): for perturbed parameters --- calculate Legendre
      !! order mode only:  
      !! -1:        standard beefxc expansion coefficients;  
      !! -2:        PBE correlation only;  
      !! -3:        LDA correlation only;  
      !! else:      no correlation either.
      INTEGER :: mode
    END SUBROUTINE beefsetmode
    !
    SUBROUTINE beefrandinit(seed)
      !! Initialize pseudo random number generator
      INTEGER :: seed
    END SUBROUTINE beefrandinit
    !
    SUBROUTINE beefrandinitdef()
      !! Initialize pseudo random number generator with default seed
    END SUBROUTINE beefrandinitdef
    !
    SUBROUTINE beefensemble(beefxc, ensemble)
      !! Calculate ensemble energies
      USE kind_l, ONLY : dp
      REAL (dp) :: beefxc(:), ensemble(:)
    END SUBROUTINE beefensemble
    !
    LOGICAL FUNCTION beef_set_type(tbeef, ionode)
      !! Set type of beef functional to be used, returns true on success.
      !! 0: BEEF-vdW
      INTEGER :: tbeef
      LOGICAL :: ionode
      CALL xclib_error('beef_set_type','no beef! support for BEEF not compiled',1)
    END FUNCTION beef_set_type
#endif
    !
END MODULE beef_interface
!
