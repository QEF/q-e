!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------

subroutine sort_gvec( ng, g2, mill )
  !---------------------------------------------------------------------
  !
  !     first the input variables
  !
  use kinds, ONLY: dbl
  use constants, ONLY: eps8
  implicit none
  INTEGER, INTENT(IN) :: ng
  REAL(dbl) :: g2( ng )
  INTEGER   :: mill( 3, ng )

  REAL(dbl), ALLOCATABLE :: gsort( : )
  INTEGER, ALLOCATABLE :: index( : )
  INTEGER :: ig, icurr, it, im

  ALLOCATE( gsort( ng ) )
  ALLOCATE( index( ng ) )

  DO ig = 1, ng
    IF ( g2(ig) > eps8 ) THEN
      gsort(ig) = g2(ig)
    ELSE
      gsort(ig) = 0.d0
    END IF
  END DO

  index(1) = 0
  CALL hpsort_eps( ng, gsort( 1 ), index( 1 ), eps8 )

  ! ... sort indices accordingly
  DO ig = 1, ng-1
    icurr = ig
30  IF( index(icurr) /= ig ) THEN
      ! ...     swap g-vec indices
      im = mill(1,icurr); mill(1,icurr) = mill(1,index(icurr)); mill(1,index(icurr)) = im
      im = mill(2,icurr); mill(2,icurr) = mill(2,index(icurr)); mill(2,index(icurr)) = im
      im = mill(3,icurr); mill(3,icurr) = mill(3,index(icurr)); mill(3,index(icurr)) = im
      ! ...     swap modules
      gsq = g2( icurr ); g2( icurr ) = g2( index(icurr) ); g2( index(icurr) ) = gsq
      ! ...     swap indices
      it = icurr; icurr = index(icurr); index(it) = it
      IF( index(icurr) == ig ) THEN
        index(icurr) = icurr
      ELSE
        GOTO 30
      END IF
    END IF
  END DO

  DEALLOCATE( gsort )
  DEALLOCATE( index )

  return
end subroutine
