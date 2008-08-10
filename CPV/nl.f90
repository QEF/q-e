!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "f_defs.h"

!

   SUBROUTINE nlrh_x( c0, tstress, bec, eigr, enl, denl )

      !  this routine computes:
      !  Kleinman-Bylander pseudopotential terms (see nlsm1)
      !  enl: nonlocal potential contribution to total energy (see ene_nl)
      !
      ! ... include modules

      USE kinds,                   ONLY: DP
      USE read_pseudo_module_fpmd, ONLY: nspnl
      USE electrons_base,          ONLY: iupdwn, nupdwn, nspin, f
      USE gvecw,                   ONLY: ngw
      USE uspp,                    ONLY: becsum, nkb
      USE cp_interfaces,           ONLY: stress_nl

      IMPLICIT NONE

      ! ... declare subroutine arguments

      COMPLEX(DP)                 :: eigr(:,:)     ! exp(i G dot r)
      COMPLEX(DP)                 :: c0(:,:)       ! wave functions
      LOGICAL,     INTENT(IN)     :: tstress       ! if .TRUE. compute nl Q-M stress
      REAL(DP)                    :: bec(:,:)
      REAL(DP),    INTENT(OUT)    :: enl
      REAL(DP),    INTENT(OUT)    :: denl( 6 )

      REAL(DP)    :: ennl
      EXTERNAL    :: ennl

      ! ... declare other variables
      !
      INTEGER     :: iss, i, j
      REAL(DP), ALLOCATABLE :: btmp( :, :, : )

      ! ... end of declarations
      !

      DO iss = 1, nspin
         !
         CALL nlsm1 ( nupdwn( iss ), 1, nspnl, eigr(1,1),    &
                      c0( 1, iupdwn( iss ) ), bec(1, iupdwn( iss ) ) )
         !
      END DO
      
      enl = ennl( becsum, bec )
      !
      IF( tstress ) THEN
         !
         ! ... compute enl (non-local) contribution
         !
         CALL stress_nl( denl, c0, f, eigr, bec, enl )
         !
      END IF
      !
      RETURN
   END SUBROUTINE nlrh_x

