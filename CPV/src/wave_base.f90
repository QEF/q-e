!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!  BEGIN manual

!==----------------------------------------------==!
        MODULE wave_base
!==----------------------------------------------==!


!  (describe briefly what this module does...)
!  ----------------------------------------------

!  END manual

          USE kinds

          IMPLICIT NONE
          SAVE
          PRIVATE

          REAL(DP) :: frice  = 0.0_DP
          !! friction parameter for electronic damped dynamics
          REAL(DP) :: grease = 0.0_DP
          !! friction parameter for electronic damped dynamics

          PUBLIC :: rande_base

          PUBLIC :: wave_steepest
          PUBLIC :: wave_verlet
          PUBLIC :: wave_speed2

          PUBLIC :: frice, grease

          PUBLIC :: print_norm_square_difference

!==----------------------------------------------==!
        CONTAINS
!==----------------------------------------------==!

      SUBROUTINE rande_base(wf,ampre)
      !! Randomize wave functions coefficients.
!  ----------------------------------------------
      USE random_numbers, ONLY : randy
      IMPLICIT NONE
! ... declare subroutine arguments
      COMPLEX(DP)          :: wf(:,:)
      REAL(DP), INTENT(IN) :: ampre

! ... declare other variables
      INTEGER i, j
      REAL(DP)  rranf1, rranf2
! ... end of declarations
!  ----------------------------------------------
      DO i = 1, SIZE(wf, 2)
        DO j = 1, SIZE( wf, 1)
          rranf1 = 0.5_DP - randy()
          rranf2 = 0.5_DP - randy()
          wf(j,i) = wf(j,i) + ampre * CMPLX(rranf1, rranf2, KIND=DP)
        END DO
      END DO
      RETURN
      END SUBROUTINE rande_base

!==----------------------------------------------==!
!==----------------------------------------------==!

   SUBROUTINE wave_steepest( CP, C0, dt2m, grad, ngw, idx )
      IMPLICIT NONE
      COMPLEX(DP), INTENT(OUT) :: CP(:)
      COMPLEX(DP), INTENT(IN) :: C0(:)
      COMPLEX(DP), INTENT(IN) :: grad(:)
      REAL(DP), INTENT(IN) ::  dt2m(:)
      INTEGER, OPTIONAL, INTENT(IN) :: ngw, idx
      !
      IF( PRESENT( ngw ) .AND. PRESENT( idx ) ) THEN
         CP( : )  = C0( : )  + dt2m(:) * grad( (idx-1)*ngw+1 : idx*ngw )
      ELSE
         CP( : )  = C0( : )  + dt2m(:) * grad(:)
      END IF
      !
      RETURN
   END SUBROUTINE wave_steepest

!==----------------------------------------------==!
!==----------------------------------------------==!

   SUBROUTINE wave_verlet( cm, c0, ver1, ver2, ver3, grad, ngw, idx )
      IMPLICIT NONE
      COMPLEX(DP), INTENT(INOUT) :: cm(:)
      COMPLEX(DP), INTENT(IN) :: c0(:)
      COMPLEX(DP), INTENT(IN) :: grad(:)
      REAL(DP), INTENT(IN) ::  ver1, ver2, ver3(:)
      INTEGER, OPTIONAL, INTENT(IN) :: ngw, idx
      !
      IF( PRESENT( ngw ) .AND. PRESENT( idx ) ) THEN
         cm( : )  = ver1 * c0( : ) + ver2 * cm( : ) + ver3( : ) * grad( (idx-1)*ngw+1:idx*ngw)
      ELSE
         cm( : )  = ver1 * c0( : ) + ver2 * cm( : ) + ver3( : ) * grad( : )
      END IF
      !
      RETURN
   END SUBROUTINE wave_verlet

!==----------------------------------------------==!
!==----------------------------------------------==!

   FUNCTION wave_speed2( cp, cm, wmss, fact )
     IMPLICIT NONE
     COMPLEX(DP), INTENT(IN) :: cp(:)
     COMPLEX(DP), INTENT(IN) :: cm(:)
     REAL(DP) :: wmss(:), fact
     REAL(DP) :: wave_speed2
     REAL(DP) :: ekinc
     COMPLEX(DP) :: speed
     INTEGER :: j
     speed  = ( cp(1) - cm(1) )
     ekinc  = fact * wmss(1) * CONJG( speed ) * speed
     DO j = 2, SIZE( cp )
       speed  = ( cp(j) - cm(j) )
       ekinc  = ekinc + wmss(j) * CONJG( speed ) * speed
     END DO
     wave_speed2 = ekinc
     RETURN
   END FUNCTION wave_speed2

      subroutine print_norm_square_difference(c1,c2,ngw,nbnd, nam,&
         ionode,comm )
      use mp, only: mp_sum
      implicit none
      complex(dp), intent(in) :: c1(ngw,nbnd), c2(ngw,nbnd)
      real(dp) :: x
      logical, intent(in) :: ionode
      integer, intent(in) :: ngw,nbnd,comm
      character(len=*), intent(in) :: nam
      integer :: i,ig
      if (ionode) &
           write (*,*) ' CHECKING NORM SQUARE DIFFERENCE OF ', trim(nam)
       x=0.0_dp
       do i =1,nbnd
          do ig = 1, ngw
                 x=x+( aimag(c1(ig,i))-aimag(c2(ig,i)))**2 &
                   + ( real(c1(ig,i))-real(c2(ig,i)))**2
          enddo
       enddo
       call mp_sum(x,comm)
       if (ionode) &
          write (*,*) ' :',  x

      end subroutine


!==----------------------------------------------==!
       END MODULE wave_base
!==----------------------------------------------==!
