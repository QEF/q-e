!
! Copyright (C) 2010 Quant ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------------=!
   MODULE gvecb
!=----------------------------------------------------------------------------=!
     USE kinds, ONLY: DP

     IMPLICIT NONE
     SAVE

     ! ...   G vectors less than the box grid cut-off ( ? )
     INTEGER :: ngb  = 0  ! local number of G vectors
     INTEGER :: ngbt = 0  ! in parallel execution global number of G vectors,
                       ! in serial execution this is equal to ngw
     INTEGER :: ngbl = 0  ! number of G-vector shells up to ngw
     INTEGER :: ngbx = 0  ! maximum local number of G vectors

     REAL(DP), ALLOCATABLE :: gb(:), gxb(:,:), glb(:)
     INTEGER, ALLOCATABLE :: npb(:), nmb(:)
     INTEGER, ALLOCATABLE :: mill_b(:,:)

     REAL(DP) :: ecutb = 0.0_DP
     REAL(DP) :: gcutb = 0.0_DP

   CONTAINS

     SUBROUTINE gvecb_set( ecut, tpibab )
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: ecut, tpibab
         ecutb = ecut
         gcutb = ecut / tpibab / tpibab
       RETURN
     END SUBROUTINE gvecb_set

     SUBROUTINE deallocate_gvecb()
       IF( ALLOCATED( gb ) ) DEALLOCATE( gb )
       IF( ALLOCATED( gxb ) ) DEALLOCATE( gxb )
       IF( ALLOCATED( glb ) ) DEALLOCATE( glb )
       IF( ALLOCATED( npb ) ) DEALLOCATE( npb )
       IF( ALLOCATED( nmb ) ) DEALLOCATE( nmb )
       IF( ALLOCATED( mill_b ) ) DEALLOCATE( mill_b )
     END SUBROUTINE deallocate_gvecb

!=----------------------------------------------------------------------------=!
   END MODULE gvecb
!=----------------------------------------------------------------------------=!
