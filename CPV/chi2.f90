!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

      module chi2

        USE kinds
        IMPLICIT NONE
        SAVE

        PRIVATE
 
        COMPLEX(dbl), allocatable :: VLOCAL(:)
        COMPLEX(dbl), allocatable :: RHOCHI(:)
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

      SUBROUTINE PRINTCHI2(gv,box)

      USE cell_base, ONLY: tpiba
      use cp_types, ONLY: recvecs
      USE cell_module, only: boxdimensions
      use mp, ONLY: mp_sum
      USE mp_global, ONLY: nproc, mpime, group
      USE io_files, ONLY: chiunit, chifile

      IMPLICIT NONE

!---------------------------------------------------ARGUMENT
      type (recvecs), intent(in) :: gv
      type (boxdimensions), intent(in) :: box
      REAL(dbl) GX,GY,GZ
!---------------------------------------------------COMMON 


      COMPLEX(dbl) :: CHI(3,3,3)
      COMPLEX(dbl) :: ctmp
      REAL(dbl)    :: omega
      INTEGER      :: idum, ig, i, j, k, l, m, s
      INTEGER      :: NFI, ierr
 
      IF( .NOT. tchi2 ) RETURN

      omega      = box%deth

      idum = 0

!=======================================================================
!==           FFT: RHO(R) --> RHO(G)  (IN ARRAY V)                    ==
!=======================================================================

      CHI = 0.0d0
      DO IG=gv%gstart,gv%ng_l 
        ctmp = CONJG(RHOCHI(IG))*VLOCAL(IG)
        do i=1,3
          GX        = gv%gx_l(i,IG)*TPIBA
          do j=1,3
            GY        = gv%gx_l(j,IG)*TPIBA
            do k=1,3
              GZ        = gv%gx_l(k,IG)*TPIBA
              CHI(i,j,k) = CHI(i,j,k) - GX*GY*GZ*ctmp
            end do
          end do
        end do
      END DO
      DO IG=gv%gstart,gv%ng_l 
        ctmp = RHOCHI(IG)*CONJG(VLOCAL(IG))
        do i=1,3
          GX        = -gv%gx_l(i,IG)*TPIBA
          do j=1,3
            GY        =  -gv%gx_l(j,IG)*TPIBA
            do k=1,3
              GZ        =  -gv%gx_l(k,IG)*TPIBA
              CHI(i,j,k) = CHI(i,j,k) - GX*GY*GZ*ctmp
            end do
          end do
        end do
      END DO

      CALL mp_sum(CHI,group)
!
      CHI = CHI * OMEGA * CMPLX(0.0d0,1.0d0)

      ierr = 0
      IF( mpime == 0 ) THEN

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

      CALL mp_sum( ierr, group )
      IF( ierr > 0 ) &
         CALL errore( ' printchi2 ', ' writing to file '//TRIM( chifile ), 1 )

      RETURN
      END SUBROUTINE PRINTCHI2

      end module chi2
