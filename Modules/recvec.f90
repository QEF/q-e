!
! Copyright (C) 2002-2010 Quantum ESPRESSO group
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

     REAL(DP) :: ecfix = 0.0_DP     ! value of the constant cut-off
     REAL(DP) :: ecutz = 0.0_DP     ! height of the penalty function (above ecfix)
     REAL(DP) :: ecsig = 0.0_DP     ! spread of the penalty function around ecfix
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
   MODULE gvecp
!=----------------------------------------------------------------------------=!
     USE kinds, ONLY: DP

     IMPLICIT NONE
     SAVE

     ! ...   G vectors less than the potential cut-off ( ecutrho )
     INTEGER :: ngm  = 0  ! local number of G vectors
     INTEGER :: ngm_g= 0  ! in parallel execution global number of G vectors,
                          ! in serial execution this is equal to ngm
     INTEGER :: ngl = 0   ! number of G-vector shells up to ngw
     INTEGER :: ngmx = 0  ! maximum local number of G vectors

     REAL(DP) :: ecutrho = 0.0_DP
     REAL(DP) :: gcutm = 0.0_DP

     !     nl      = fft index for G>
     !     nlm     = fft index for G<

     INTEGER, ALLOCATABLE :: nl(:), nlm(:)

   CONTAINS

     SUBROUTINE deallocate_gvecp()
       IF( ALLOCATED( nl ) ) DEALLOCATE( nl )
       IF( ALLOCATED( nlm ) ) DEALLOCATE( nlm )
     END SUBROUTINE deallocate_gvecp

!=----------------------------------------------------------------------------=!
   END MODULE gvecp
!=----------------------------------------------------------------------------=!

!=----------------------------------------------------------------------------=!
   MODULE gvecs
!=----------------------------------------------------------------------------=!
     USE kinds, ONLY: DP

     IMPLICIT NONE
     SAVE

     ! ...   G vectors less than the smooth grid cut-off ( ? )
     INTEGER :: ngms = 0  ! local number of G vectors
     INTEGER :: ngms_g=0  ! in parallel execution global number of G vectors,
                       ! in serial execution this is equal to ngw
     INTEGER :: ngsx = 0  ! maximum local number of G vectors

     INTEGER, ALLOCATABLE :: nls(:), nlsm(:)

     REAL(DP) :: ecuts = 0.0_DP
     REAL(DP) :: gcutms= 0.0_DP

     REAL(DP) :: dual = 0.0_DP
     LOGICAL   :: doublegrid = .FALSE.

   CONTAINS

     SUBROUTINE deallocate_gvecs()
       IF( ALLOCATED( nls ) ) DEALLOCATE( nls )
       IF( ALLOCATED( nlsm ) ) DEALLOCATE( nlsm )
     END SUBROUTINE deallocate_gvecs

!=----------------------------------------------------------------------------=!
   END MODULE gvecs
!=----------------------------------------------------------------------------=!

!=----------------------------------------------------------------------------=!
   MODULE reciprocal_vectors
!=----------------------------------------------------------------------------=!

     USE kinds, ONLY: DP
     USE gvecp
     USE gvecs
     USE gvecw

     IMPLICIT NONE
     SAVE

     ! ...   declare module-scope variables

     INTEGER :: gstart = 2 ! index of the first G vector whose module is > 0
                           ! Needed in parallel execution: gstart=2 for the
                           ! proc that holds G=0, gstart=1 for all others

     !     G^2 in increasing order (in units of tpiba2=(2pi/a)^2)
     !
     REAL(DP), ALLOCATABLE, TARGET :: gg(:) 

     !     shells of G^2
     !
     REAL(DP), ALLOCATABLE, TARGET :: gl(:) 

     !     G-vectors cartesian components ( units tpiba =(2pi/a)  )
     !
     REAL(DP), ALLOCATABLE, TARGET :: g(:,:) 

     !     g2_g    = all G^2 in increasing order, replicated on all procs
     !
     REAL(DP), ALLOCATABLE, TARGET :: g2_g(:)

     !     mill_g  = miller index of G vecs (increasing order), replicated on all procs
     !
     INTEGER, ALLOCATABLE, TARGET :: mill_g(:,:)

     !     mill    = miller index of G vectors (local to each processor)
     !
     INTEGER, ALLOCATABLE, TARGET :: mill(:,:)

     !     ig_l2g  = "l2g" means local to global, this array convert a local
     !               G-vector index into the global index, in other words
     !               the index of the G-v. in the overall array of G-vectors
     !
     INTEGER, ALLOCATABLE, TARGET :: ig_l2g(:)

     !     sortedig_l2g = array obtained by sorting ig_l2g
     !
     !
     INTEGER, ALLOCATABLE, TARGET :: sortedig_l2g(:)

   CONTAINS

     SUBROUTINE deallocate_recvecs
       IF( ALLOCATED( gg ) ) DEALLOCATE( gg )
       IF( ALLOCATED( gl ) ) DEALLOCATE( gl )
       IF( ALLOCATED( g ) )  DEALLOCATE( g )
       IF( ALLOCATED( g2_g ) ) DEALLOCATE( g2_g )
       IF( ALLOCATED( mill_g ) ) DEALLOCATE( mill_g )
       IF( ALLOCATED( mill ) ) DEALLOCATE( mill )
       IF( ALLOCATED( ig_l2g ) ) DEALLOCATE( ig_l2g )
       IF( ALLOCATED( sortedig_l2g ) ) DEALLOCATE( sortedig_l2g )
       CALL deallocate_gvecw( )
       CALL deallocate_gvecs( )
     END SUBROUTINE deallocate_recvecs

!=----------------------------------------------------------------------------=!
   END MODULE reciprocal_vectors
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
       USE gvecp, ONLY: ngm, ngmx, ngm_g
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
