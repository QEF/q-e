!
! Copyright (C) 2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!=----------------------------------------------------------------------------=!
   MODULE gvecw
!=----------------------------------------------------------------------------=!
     USE kinds, ONLY: DP

     IMPLICIT NONE
     SAVE

     ! ...   G vectors less than the wave function cut-off ( ecutwfc )
     INTEGER :: ngw  = 0  ! local number of G vectors
     INTEGER :: ngw_g= 0  ! in parallel execution global number of G vectors,
                       ! in serial execution this is equal to ngw
     INTEGER :: ngwx = 0  ! maximum local number of G vectors
     INTEGER :: ng0  = 0  ! first G-vector with nonzero modulus
                       ! needed in the parallel case (G=0 is on one node only!)

     REAL(DP) :: ecutwfc = 0.0_DP
     REAL(DP) :: gcutw = 0.0_DP

     !   values for costant cut-off computations

     REAL(DP) :: ecfixed=0.0_DP     ! value of the constant cut-off
     REAL(DP) :: qcutz = 0.0_DP     ! height of the penalty function (above ecfix)
     REAL(DP) :: q2sigma=0.0_DP     ! spread of the penalty function around ecfix
     LOGICAL   :: tecfix = .FALSE.  ! .TRUE. if constant cut-off is in use

     ! augmented cut-off for k-point calculation

     REAL(DP) :: ekcut = 0.0_DP
     REAL(DP) :: gkcut = 0.0_DP
    
     ! array of G vectors module plus penalty function for constant cut-off 
     ! simulation.
     !
     ! ggp = g + ( agg / tpiba**2 ) * ( 1 + erf( ( tpiba2 * g - e0gg ) / sgg ) )

     REAL(DP), ALLOCATABLE, TARGET :: ggp(:)

   CONTAINS

     SUBROUTINE deallocate_gvecw
       IF( ALLOCATED( ggp ) ) DEALLOCATE( ggp )
     END SUBROUTINE deallocate_gvecw

!=----------------------------------------------------------------------------=!
   END MODULE gvecw
!=----------------------------------------------------------------------------=!

!=----------------------------------------------------------------------------=!
   MODULE recvecs_subroutines
!=----------------------------------------------------------------------------=!

     IMPLICIT NONE
     SAVE

   CONTAINS

     SUBROUTINE recvecs_init( ngm_ , ngw_ , ngs_ )
       USE mp_global, ONLY: intra_image_comm
       USE mp, ONLY: mp_max, mp_sum
       USE gvecw, ONLY: ngw, ngwx, ngw_g
       USE gvect, ONLY: ngm, ngmx, ngm_g
       USE gvecs, ONLY: ngms, ngsx, ngms_g

       IMPLICIT NONE
       INTEGER, INTENT(IN) :: ngm_ , ngw_ , ngs_

       ngm = ngm_
       ngw = ngw_
       ngms= ngs_

       !
       !  calculate maxima over all processors
       !
       ngwx = ngw
       ngmx = ngm
       ngsx = ngms
       CALL mp_max( ngwx, intra_image_comm )
       CALL mp_max( ngmx, intra_image_comm )
       CALL mp_max( ngsx, intra_image_comm )
       !
       !  calculate SUM over all processors
       !
       ngw_g = ngw
       ngm_g= ngm
       ngms_g=ngms
       CALL mp_sum( ngw_g, intra_image_comm )
       CALL mp_sum( ngm_g, intra_image_comm )
       CALL mp_sum( ngms_g,intra_image_comm )

       RETURN 
     END SUBROUTINE recvecs_init

!=----------------------------------------------------------------------------=!
   END MODULE recvecs_subroutines
!=----------------------------------------------------------------------------=!
