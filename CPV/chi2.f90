!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
      module chi2

        USE kinds
        IMPLICIT NONE
        SAVE

        PRIVATE
 
        COMPLEX(DP), allocatable :: VLOCAL(:)
        COMPLEX(DP), allocatable :: RHOCHI(:)
        logical :: tchi2

        PUBLIC :: allocate_chi2, deallocate_chi2
        PUBLIC :: rhochi

      contains

        subroutine allocate_chi2(ng)
          integer, intent(in) :: ng
          allocate(VLOCAL(ng))
          allocate(RHOCHI(ng))
        end subroutine allocate_chi2

        subroutine deallocate_chi2
          IF( ALLOCATED( vlocal ) ) DEALLOCATE(vlocal)
          IF( ALLOCATED( rhochi ) ) DEALLOCATE(rhochi)
        end subroutine deallocate_chi2

      SUBROUTINE PRINTCHI2(box)

      USE cell_base, ONLY: tpiba
      USE cell_module, only: boxdimensions
      use mp, ONLY: mp_sum
      USE mp_global, ONLY: intra_image_comm
      USE io_global, ONLY: ionode
      USE io_files, ONLY: chiunit, chifile
      USE reciprocal_vectors, ONLY: gstart, gx
      USE gvecp, ONLY: ngm

      IMPLICIT NONE

!---------------------------------------------------ARGUMENT
      type (boxdimensions), intent(in) :: box
      REAL(DP) GXR,GYR,GZR
!---------------------------------------------------COMMON 


      COMPLEX(DP) :: CHI(3,3,3)
      COMPLEX(DP) :: ctmp
      REAL(DP)    :: omega
      INTEGER      :: idum, ig, i, j, k, l, m, s
      INTEGER      :: NFI, ierr
 
      IF( .NOT. tchi2 ) RETURN

      omega      = box%deth

      idum = 0

!=======================================================================
!==           FFT: RHO(R) --> RHO(G)  (IN ARRAY V)                    ==
!=======================================================================

      CHI = 0.0d0
      DO IG = gstart, ngm 
        ctmp = CONJG(RHOCHI(IG))*VLOCAL(IG)
        do i=1,3
          GXR       = gx(i,IG)*TPIBA
          do j=1,3
            GYR       = gx(j,IG)*TPIBA
            do k=1,3
              GZR       = gx(k,IG)*TPIBA
              CHI(i,j,k) = CHI(i,j,k) - GXR*GYR*GZR*ctmp
            end do
          end do
        end do
      END DO
      DO IG = gstart, ngm 
        ctmp = RHOCHI(IG)*CONJG(VLOCAL(IG))
        do i=1,3
          GXR       = -gx(i,IG)*TPIBA
          do j=1,3
            GYR       =  -gx(j,IG)*TPIBA
            do k=1,3
              GZR       =  -gx(k,IG)*TPIBA
              CHI(i,j,k) = CHI(i,j,k) - GXR*GYR*GZR*ctmp
            end do
          end do
        end do
      END DO

      CALL mp_sum(CHI,intra_image_comm)
!
      CHI = CHI * OMEGA * CMPLX(0.0d0,1.0d0)

      ierr = 0
      IF( ionode ) THEN

         OPEN(UNIT=chiunit, FILE=chifile, STATUS='unknown', POSITION='append', ERR=300)

         WRITE(chiunit,*) ' == CHI2 ==' 
         do i=1,3
           do j=1,3
             do k=1,3
                WRITE(chiunit,'(2X,3I3,2F12.4)')  i, j, k, chi( i, j, k )
             end do
           end do
         end do

         CLOSE(UNIT=chiunit, ERR=300)

         GO TO 310
 300     ierr = 1
 310     CONTINUE

      END IF

      CALL mp_sum( ierr, intra_image_comm )
      IF( ierr > 0 ) &
         CALL errore( ' printchi2 ', ' writing to file '//TRIM( chifile ), 1 )

      RETURN
      END SUBROUTINE PRINTCHI2

      end module chi2
