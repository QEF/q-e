!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!  ----------------------------------------------
!  BEGIN manual

      MODULE stick

        USE kinds
        USE fft_types, ONLY: fft_dlay_descriptor, fft_dlay_allocate, fft_dlay_deallocate, &
          fft_dlay_set
        USE fft_base, ONLY: dfftp, dffts

        IMPLICIT NONE
        PRIVATE
        SAVE

        ! dfftp  ! fft descriptor for potentials
        ! dffts  ! fft descriptor for potentials

        PUBLIC :: deallocate_stick, fft_dlay_descriptor 
        PUBLIC :: dfftp, dffts

!  ----------------------------------------------
!  END manual


!=----------------------------------------------------------------------=
      CONTAINS
!=----------------------------------------------------------------------=

        SUBROUTINE deallocate_stick
          CALL fft_dlay_deallocate( dfftp )
          CALL fft_dlay_deallocate( dffts )
          RETURN
        END SUBROUTINE deallocate_stick

!=----------------------------------------------------------------------=


      END MODULE STICK
