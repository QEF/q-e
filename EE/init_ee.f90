!
! Copyright (C) 2007-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! original version by I. Dabo and N. Marzari (MIT)
!
! contributions by E. Lamas and S. de Gironcoli (SISSA/DEMOCRITOS)
!
!--------------------------------------------------------------------
      SUBROUTINE init_ee(nrxx)
!--------------------------------------------------------------------
 
#include "f_defs.h"
      
      ! ... Declares modules
      USE kinds,             ONLY: DP
      USE uspp_param,       ONLY : nhm
      USE ions_base,        ONLY : nat, ntyp => nsp
      USE ee_mod
      !
      IMPLICIT NONE      
      !
      !
      INTEGER, INTENT (IN) :: nrxx
      !
      INTEGER :: i
      INTEGER :: j
      INTEGER :: k
      !
!      n_cycle = 0
      !
      ! ... Allocates self-interaction variables
      !
      ! ... Allocates charge compensation variables
      !
      IF( do_comp ) THEN
        ALLOCATE( vcomp( nrxx ) )
        ALLOCATE( vloccoul( mr1 * mr2 * mr3 ) )
        ALLOCATE( rhoion( mr1 * mr2 * mr3 ) )
        ALLOCATE( vcoul( mr1 * mr2 * mr3 ) )
        vcomp = 0.D0
        rhoion = 0.D0
        ecomp = 0.D0
        vcoul = 0.D0
        omegafact = 1.D0 / ( cellmax( 1 ) - cellmin( 1 ) )             &
                  / ( cellmax( 2 ) - cellmin( 2 ) )                    &
                  / ( cellmax( 3 ) - cellmin( 3 ) )
      END IF
      !
      RETURN

!--------------------------------------------------------------------
      END SUBROUTINE init_ee
!--------------------------------------------------------------------
