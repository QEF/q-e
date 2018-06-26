!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE plugin_tddfpt_potential(drho,dv)
!
! This routine is used for printing plugins clocks
! DO NOT REMOVE THE TAGS ! ***ADDSON_NAME KIND_OF_PATCH***
!
USE io_global,        ONLY : stdout, ionode
USE kinds,            ONLY : DP
USE fft_base,         ONLY : dfftp
USE lsda_mod,         ONLY : nspin
USE lr_variables,     ONLY : davidson
USE plugin_flags
!
! ***Environ MODULES BEGIN***
USE scf,              ONLY : rho
USE environ_main,     ONLY : calc_dvenviron
! ***Environ MODULES END***
!
IMPLICIT NONE
!
! ***Environ VARIABLES BEGIN***
! ***Environ VARIABLES END***
!
! ***Environ CALLS BEGIN***
IF ( use_environ ) THEN
   !
   IF (.not.davidson) THEN
      WRITE( stdout, '(5x,"ENVIRON: Calculate the response &
           & polarization and dielectric potentials")' )
   ENDIF
   !
   CALL calc_dvenviron( dfftp%nnr, nspin, rho%of_r, rho_1, dv )
   !
ENDIF
! ***Environ CALLS END***
!
RETURN
!
END SUBROUTINE plugin_tddfpt_potential
