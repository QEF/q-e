!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


!  AB INITIO COSTANT PRESSURE MOLECULAR DYNAMICS
!  ----------------------------------------------
!  Car-Parrinello Parallel Program
!  Carlo Cavazzoni - Gerardo Ballabio
!  SISSA, Trieste, Italy - 1997-99
!  Last modified: Fri Oct  8 13:34:20 MDT; 1999
!  ----------------------------------------------
      MODULE cutoffs

!  (describe briefly what this module does...)
!  ----------------------------------------------
!  routines in this module:
!  SUBROUTINE cutoffs_print_info(unit)
!  SUBROUTINE cutoffs_setup(alat,ecut_inp,ecfix_inp, gcutz_inp, &
!                           esig_inp,tk_inp,nk_inp,kpoints_inp)
!  ----------------------------------------------

        USE kinds

        IMPLICIT NONE
        SAVE

        PRIVATE

! ...   declare module-scope variables

        PUBLIC :: cutoffs_print_info

!  end of module-scope declarations
!  ----------------------------------------------

      CONTAINS

!  subroutines
!  ----------------------------------------------
!  ----------------------------------------------
        SUBROUTINE cutoffs_print_info(unit)

!  (describe briefly what this routine does...)
!  ----------------------------------------------

        USE gvecw, ONLY: ecutwfc => ecutw,  gcutw
        USE gvecp, ONLY: ecutrho => ecutp,  gcutp
        USE gvecw, ONLY: ecfix, gcfix, ecutz, gcutz, ecsig, gcsig, tecfix
        USE gvecw, ONLY: ekcut, gkcut
        USE gvecs, ONLY: ecuts, gcuts
        USE gvecb, ONLY: ecutb, gcutb

! ...     declare subroutine argument
          INTEGER, INTENT(IN) :: unit

          WRITE(unit,100) ecutwfc, ecutrho, sqrt(gcutw), sqrt(gcutp) 
100       FORMAT(/,3X,'Energy Cut-offs',/ &
                  ,3X,'---------------',/ &
                  ,3X,'Ecutwfc = ',F6.1,'Ryd., ', 3X,'Ecutrho = ',F6.1,'Ryd., ',/ &
                  ,3X,'Gcutwfc = ',F6.1,'Ryd., ', 3X,'Gcutrho = ',F6.1,'Ryd.')
          RETURN
        END SUBROUTINE cutoffs_print_info

!  ----------------------------------------------
!  ----------------------------------------------
 
      END MODULE cutoffs

