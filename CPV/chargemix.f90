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
!  Last modified: Sun Nov 21 13:29:38 MET 1999
!  ----------------------------------------------
!  BEGIN manual

!=----------------------------------------------------------------------------=!
      MODULE charge_mix
!=----------------------------------------------------------------------------=!

!  (describe briefly what this module does...)
!  ----------------------------------------------
!  routines in this module:
!  SUBROUTINE allocate_charge_mix(ng_l)
!  SUBROUTINE deallocate_charge_mix
!  SUBROUTINE charge_mix_print_info(unit)
!  SUBROUTINE newrho(rhor,nfi,tcel,drho)
!  SUBROUTINE invgen(aa,dimaa)
!  ----------------------------------------------
!  END manual

        USE kinds
        IMPLICIT NONE
        SAVE

        PRIVATE

! ...   declare module-scope variables
        REAL(DP)  :: achmix
        REAL(DP)  :: g1met2
        REAL(DP)  :: g0chmix2
        INTEGER    :: daamax

        REAL(DP), ALLOCATABLE :: aa_save(:,:)
        COMPLEX(DP), ALLOCATABLE :: rho(:,:)
        COMPLEX(DP), ALLOCATABLE :: rr(:,:)
        COMPLEX(DP), ALLOCATABLE :: chmix(:)
        COMPLEX(DP), ALLOCATABLE :: metric(:)

! ...   end of module-scope declarations
!  ----------------------------------------------

        PUBLIC :: charge_mix_setup
        PUBLIC :: allocate_charge_mix, deallocate_charge_mix
        PUBLIC :: charge_mix_print_info
        PUBLIC :: chmix, metric, rho, rr, aa_save
        PUBLIC :: achmix, g1met2, g0chmix2, daamax

!=----------------------------------------------------------------------------=!
      CONTAINS
!=----------------------------------------------------------------------------=!

!  subroutines
!  ----------------------------------------------
!  ----------------------------------------------

        SUBROUTINE charge_mix_setup(achmix_inp, g0chmix2_inp, daamax_inp, g1met2_inp)
          REAL(DP), INTENT(IN) ::  achmix_inp, g0chmix2_inp
          REAL(DP), INTENT(IN) ::  g1met2_inp
          INTEGER, INTENT(IN) :: daamax_inp
          achmix = achmix_inp
          g0chmix2 = g0chmix2_inp
          daamax = daamax_inp
          g1met2 = g1met2_inp
          RETURN
        END SUBROUTINE charge_mix_setup

!  ----------------------------------------------
!  ----------------------------------------------

        SUBROUTINE allocate_charge_mix(ng)
          INTEGER, INTENT(IN) :: ng
          INTEGER :: ierr
          ALLOCATE( rho(ng, daamax), STAT=ierr )
          IF( ierr /= 0 ) CALL errore(' allocate_charge_mix ', ' allocating rho ', ierr)
          ALLOCATE( rr(ng, daamax), STAT=ierr )
          IF( ierr /= 0 ) CALL errore(' allocate_charge_mix ', ' allocating rr ', ierr)
          ALLOCATE( aa_save(daamax, daamax), STAT=ierr )
          IF( ierr /= 0 ) CALL errore(' allocate_charge_mix ', ' allocating aa_save ', ierr)
          ALLOCATE( chmix(ng), STAT=ierr )
          IF( ierr /= 0 ) CALL errore(' allocate_charge_mix ', ' allocating chmix ', ierr)
          ALLOCATE( metric(ng), STAT=ierr )
          IF( ierr /= 0 ) CALL errore(' allocate_charge_mix ', ' allocating metric ', ierr)
          RETURN
        END SUBROUTINE allocate_charge_mix

!  ----------------------------------------------
!  ----------------------------------------------

        SUBROUTINE deallocate_charge_mix
          INTEGER :: ierr
          IF( ALLOCATED(rho) ) THEN
            DEALLOCATE(rho, STAT=ierr)
            IF( ierr /= 0 ) CALL errore(' deallocate_charge_mix ', ' deallocating rho ', ierr)
          END IF
          IF( ALLOCATED(rr) ) THEN
            DEALLOCATE(rr, STAT=ierr)
            IF( ierr /= 0 ) CALL errore(' deallocate_charge_mix ', ' deallocating rr ', ierr)
          END IF
          IF( ALLOCATED(aa_save) ) THEN
            DEALLOCATE(aa_save, STAT=ierr)
            IF( ierr /= 0 ) CALL errore(' deallocate_charge_mix ', ' deallocating aa_save ', ierr)
          END IF
          IF( ALLOCATED(chmix) ) THEN
            DEALLOCATE(chmix, STAT=ierr)
            IF( ierr /= 0 ) CALL errore(' deallocate_charge_mix ', ' deallocating chmix ', ierr)
          END IF
          IF( ALLOCATED(metric) ) THEN
            DEALLOCATE(metric, STAT=ierr)
            IF( ierr /= 0 ) CALL errore(' deallocate_charge_mix ', ' deallocating metric ', ierr)
          END IF
          RETURN
        END SUBROUTINE deallocate_charge_mix

!  ----------------------------------------------
!  ----------------------------------------------

        SUBROUTINE charge_mix_print_info(unit)
          INTEGER, INTENT(IN) :: unit
            WRITE(unit,300)
            WRITE(unit,310) achmix, g0chmix2, g1met2
            WRITE(unit,320) daamax
300         FORMAT(/,3X,'Charge mixing parameters:')
310         FORMAT(  3X,'A = ', D14.6, ' G0^2 = ', D14.6,' G1^2 = ',D14.6)
320         FORMAT(  3X,'charge mixing matrix maximum size = ', I5)
          RETURN
        END SUBROUTINE charge_mix_print_info

!  ----------------------------------------------
!  ----------------------------------------------

!=----------------------------------------------------------------------------=!
      END MODULE charge_mix
!=----------------------------------------------------------------------------=!

