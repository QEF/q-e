!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!  ----------------------------------------------
!  BEGIN manual

      MODULE pseudotab_base

!  (describe briefly what this module does...)
!  ----------------------------------------------
!  routines in this module:
!  ----------------------------------------------
!  END manual

! ... declare modules

      USE kinds
      USE constants, ONLY: gsmall
      USE cell_base, ONLY: tpiba

      IMPLICIT NONE

      SAVE

      PRIVATE

      PUBLIC ::  nlintab_base, formftab_base, corecortab_base
      PUBLIC ::  chkpstab

      CONTAINS

!  subroutines
!  ----------------------------------------------
!  ----------------------------------------------
      LOGICAL FUNCTION chkpstab(hg, xgtabmax)
        USE mp, ONLY: mp_max
        USE io_global, ONLY: stdout

        IMPLICIT none
        REAL(dbl), INTENT(IN) :: hg(:)
        REAL(dbl), INTENT(IN) :: xgtabmax
        REAL(dbl) :: xgmax

        chkpstab = .FALSE.
        xgmax = tpiba * SQRT( MAXVAL( hg ) )
        CALL mp_max(xgmax)
        IF( xgmax > xgtabmax ) THEN
          chkpstab = .TRUE.
          WRITE( stdout, fmt='(  "CHKPSTAB: recalculate pseudopotential table" )' )
        END IF
        RETURN
      END FUNCTION chkpstab


      SUBROUTINE nlintab_base(hg, wnl, lnl, wnl_sp, xgmax)

        USE splines, ONLY: spline_data, spline

        IMPLICIT NONE

! ...   subroutine arguments
        REAL(dbl), INTENT(OUT) :: wnl(:,:)
        TYPE(spline_data)      :: wnl_sp(:)
        REAL(dbl), INTENT(IN)  :: xgmax
        REAL(dbl), INTENT(IN) :: hg(:)
        INTEGER, INTENT(IN) :: lnl

! ...   other variables
        REAL(dbl) :: xg
        INTEGER :: ig, l

        DO ig = 1, SIZE( wnl, 1 )
          IF( hg(ig) < gsmall ) THEN
            DO l = 1,lnl
              wnl(ig,l) = wnl_sp(l)%y(1)
            END DO
          ELSE
            xg = SQRT( hg(ig) ) * tpiba
            DO l = 1, lnl
              wnl(ig,l) = spline( wnl_sp(l), xg ) 
            END DO
          END IF
        END DO
        RETURN
      END SUBROUTINE nlintab_base


!  ----------------------------------------------
!  ----------------------------------------------

      SUBROUTINE formftab_base(hg, vps, dvps, vps_sp, dvps_sp, xgmax, omega)

!  this routine computes:
!  pseudopotentials are given as numerical tables
!  ----------------------------------------------

        USE splines, ONLY: spline_data, spline

        IMPLICIT NONE

! ...   declare subroutine arguments
        REAL(dbl), INTENT(OUT) :: vps(:), dvps(:)
        TYPE(spline_data)      ::  vps_sp
        TYPE(spline_data)      :: dvps_sp
        REAL(dbl), INTENT(IN)  :: xgmax, omega
        REAL(dbl), INTENT(IN)  :: hg(:)

! ...   declare other variables
        REAL(dbl) :: xg, cost1
        INTEGER :: ig
! ...   end of declarations
!  ----------------------------------------------
          cost1 = 1.0d0/omega
          DO ig = 1, SIZE( vps )
            IF( hg(ig) < gsmall ) THEN
              vps(ig)  = vps_sp%y(1)  * cost1
              dvps(ig) = dvps_sp%y(1) * cost1
            ELSE
              xg = SQRT( hg(ig) ) * tpiba
              vps(ig)  = spline( vps_sp, xg ) 
              dvps(ig) = spline( dvps_sp, xg ) 
              vps(ig)  = vps(ig)  * cost1
              dvps(ig) = dvps(ig) * cost1
            END IF
          END DO
        RETURN
      END SUBROUTINE formftab_base

!  ----------------------------------------------

      SUBROUTINE corecortab_base(hg, rhoc1, rhocp, rhoc1_sp, rhocp_sp, xgmax, omega  )

!  this routine computes:
!  ----------------------------------------------

        USE splines, ONLY: spline_data, spline

        IMPLICIT NONE

! ...   declare subroutine arguments
        REAL(dbl), INTENT(OUT) :: rhoc1(:)     , rhocp(:)
        TYPE(spline_data)      :: rhoc1_sp  , rhocp_sp
        REAL(dbl), INTENT(IN)  :: xgmax, omega
        REAL(dbl), INTENT(IN)  :: hg(:)
! ...   declare other variables
        REAL(dbl) :: xg, cost1
        INTEGER :: ig
! ...   end of declarations
!  ----------------------------------------------
          cost1 = 1.0d0/omega
          DO ig = 1, SIZE(rhoc1)
            IF( hg(ig) < gsmall ) THEN
              rhoc1(ig) = rhoc1_sp%y(1) * cost1
              rhocp(ig) = 0.0d0
            ELSE
              xg = SQRT(hg(ig)) * tpiba
              rhoc1(ig) = spline( rhoc1_sp, xg ) 
              rhocp(ig) = spline( rhocp_sp, xg ) 

              rhoc1(ig) = rhoc1(ig) * cost1
              rhocp(ig) = rhocp(ig) * cost1
            END IF
          END DO
        RETURN
      END SUBROUTINE corecortab_base

  END MODULE pseudotab_base
