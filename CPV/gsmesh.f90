!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"

!  AB INITIO COSTANT PRESSURE MOLECULAR DYNAMICS
!  ----------------------------------------------
!  Car-Parrinello Parallel Program
!  Carlo Cavazzoni - Gerardo Ballabio
!  SISSA, Trieste, Italy - 1997-99
!  Last modified: Fri Oct  8 15:08:56 MDT; 1999
!  ----------------------------------------------
!  BEGIN manual

!=----------------------------------------------------------------------------=!
      MODULE reciprocal_space_mesh
!=----------------------------------------------------------------------------=!

!  (describe briefly what this module does...)
!  ----------------------------------------------
!  routines in this module:
!  SUBROUTINE gmeshset(ngw,ng)
!  INTEGER FUNCTION owner_of_gvec(ig)
!  SUBROUTINE newg(gv,htm1)
!  ----------------------------------------------
!  END manual

! ...   declare included modules

        USE kinds

        USE cell_base, ONLY: tpiba, tpiba2

        IMPLICIT NONE
        SAVE

        PRIVATE

        PUBLIC :: recvecs_units 


!=----------------------------------------------------------------------------=!
      CONTAINS
!=----------------------------------------------------------------------------=!


   SUBROUTINE recvecs_units( alat )
     USE constants, ONLY: tpi
     REAL(DP) :: alat
     IF( alat == 0.0d0 ) &
       CALL errore(' recvecs_units ', ' alat is zero ! ', 1 )
     tpiba = tpi / alat
     tpiba2 = tpiba ** 2
     RETURN
   END SUBROUTINE recvecs_units


!  ----------------------------------------------


   INTEGER FUNCTION owner_of_gvec(mill)

     USE stick_base, ONLY: stown => sticks_owner

     IMPLICIT NONE

     INTEGER, INTENT(IN) :: mill(:)

     owner_of_gvec = stown( mill(1), mill(2) )

     RETURN
   END FUNCTION owner_of_gvec


!=----------------------------------------------------------------------------=!
      END MODULE reciprocal_space_mesh
!=----------------------------------------------------------------------------=!

