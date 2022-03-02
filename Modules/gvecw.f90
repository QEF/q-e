!
! Copyright (C) 2010-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------------=!
   MODULE gvecw
!=----------------------------------------------------------------------------=!
     !! G vectors module.
     
     USE kinds, ONLY: DP

     IMPLICIT NONE
     SAVE

     PRIVATE
     PUBLIC :: ngw, ngw_g, ngwx, ecutwfc, gcutw, ekcut, gkcut
     PUBLIC :: gg2kin, ecfixed, qcutz, q2sigma
     PUBLIC :: gvecw_init, gg2kin_init, deallocate_gvecw

     ! ...   G vectors less than the wave function cut-off ( ecutwfc )
     INTEGER :: ngw  = 0
     !! local number of G vectors
     INTEGER :: ngw_g= 0
     !! in parallel execution global number of G vectors, in serial execution
     !! this is equal to \(\text{ngw}\).
     INTEGER :: ngwx = 0
     !! maximum local number of G vectors

     REAL(DP) :: ecutwfc = 0.0_DP
     !! wave function cut-off
     REAL(DP) :: gcutw = 0.0_DP

     ! values for costant cut-off computations

     REAL(DP) :: ecfixed=0.0_DP
     !! value of the constant cut-off
     REAL(DP) :: qcutz = 0.0_DP
     !! height of the penalty function (above ecfix)
     REAL(DP) :: q2sigma=0.0_DP
     !! spread of the penalty function around ecfix
     ! augmented cut-off for k-point calculation

     REAL(DP) :: ekcut = 0.0_DP
     REAL(DP) :: gkcut = 0.0_DP
    
     ! array of G vectors module plus penalty function for constant cut-off 
     ! simulation.
     
     REAL(DP), ALLOCATABLE :: gg2kin(:)
     !! \(\text{gg2kin} = g + (\text{agg} / \text{tpiba}^2)\cdot(1+\text{erf}
     !! ((\text{tpiba2}\cdot g-\text{e0gg})/\text{sgg}))\)

   CONTAINS

     SUBROUTINE gvecw_init( ngw_ , comm )
       !
       USE mp, ONLY: mp_max, mp_sum
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: ngw_
       INTEGER, INTENT(IN) :: comm
       !
       ngw = ngw_
       !
       !  calculate maximum over all processors
       !
       ngwx = ngw
       CALL mp_max( ngwx, comm )
       !
       !  calculate sum over all processors
       !
       ngw_g = ngw
       CALL mp_sum( ngw_g, comm )
       !
       !  allocate kinetic energy
       !
       ALLOCATE( gg2kin(ngw) )
       !$acc enter data create(gg2kin)
       !
       RETURN 
       !
     END SUBROUTINE gvecw_init
     !
     SUBROUTINE gg2kin_init( ggg, tpiba2 )
       !! Initialize kinetic energy
       USE gvect, ONLY : gg
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: ggg(:), tpiba2
       REAL(DP) :: gcutz
       INTEGER :: ig
       !
       gcutz  = qcutz / tpiba2
       IF( gcutz > 0.0d0 ) THEN
          DO ig=1,ngw
             gg2kin(ig) = gg(ig) + gcutz * &
                     ( 1.0d0 + ERF( ( tpiba2 *gg(ig) - ecfixed )/q2sigma ) )
          ENDDO
       ELSE
          gg2kin( 1 : ngw ) = gg( 1 : ngw )
       END IF
       !
       !$acc update device(gg2kin)
       !
       RETURN 
       !
     END SUBROUTINE gg2kin_init
     !
     SUBROUTINE deallocate_gvecw
       !$acc exit data delete(gg2kin)
       IF( ALLOCATED( gg2kin ) ) DEALLOCATE( gg2kin )
     END SUBROUTINE deallocate_gvecw
     !=----------------------------------------------------------------------------=!
   END MODULE gvecw
!=----------------------------------------------------------------------------=!
