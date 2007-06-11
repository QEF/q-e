!
! Copyright (C) 2002 FPMD group
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
     INTEGER :: ngwt = 0  ! in parallel execution global number of G vectors,
                       ! in serial execution this is equal to ngw
     INTEGER :: ngwl = 0  ! number of G-vector shells up to ngw
     INTEGER :: ngwx = 0  ! maximum local number of G vectors
     INTEGER :: ng0  = 0  ! first G-vector with nonzero modulus
                       ! needed in the parallel case (G=0 is on one node only!)

     REAL(DP) :: ecutw = 0.0_DP
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
     INTEGER :: ngmt = 0  ! in parallel execution global number of G vectors,
                       ! in serial execution this is equal to ngm
     INTEGER :: ngml = 0  ! number of G-vector shells up to ngw
     INTEGER :: ngmx = 0  ! maximum local number of G vectors

     REAL(DP) :: ecutp = 0.0_DP
     REAL(DP) :: gcutp = 0.0_DP

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
     INTEGER :: ngs  = 0  ! local number of G vectors
     INTEGER :: ngst = 0  ! in parallel execution global number of G vectors,
                       ! in serial execution this is equal to ngw
     INTEGER :: ngsl = 0  ! number of G-vector shells up to ngw
     INTEGER :: ngsx = 0  ! maximum local number of G vectors

     INTEGER, ALLOCATABLE :: nps(:), nms(:)

     REAL(DP) :: ecuts = 0.0_DP
     REAL(DP) :: gcuts = 0.0_DP

     REAL(DP) :: dual = 0.0_DP
     LOGICAL   :: doublegrid = .FALSE.

   CONTAINS

     SUBROUTINE deallocate_gvecs()
       IF( ALLOCATED( nps ) ) DEALLOCATE( nps )
       IF( ALLOCATED( nms ) ) DEALLOCATE( nms )
     END SUBROUTINE deallocate_gvecs

!=----------------------------------------------------------------------------=!
   END MODULE gvecs
!=----------------------------------------------------------------------------=!

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
     INTEGER, ALLOCATABLE :: npb(:), nmb(:), iglb(:)
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
       IF( ALLOCATED( iglb ) ) DEALLOCATE( iglb )
       IF( ALLOCATED( mill_b ) ) DEALLOCATE( mill_b )
     END SUBROUTINE deallocate_gvecb

!=----------------------------------------------------------------------------=!
   END MODULE gvecb
!=----------------------------------------------------------------------------=!


!=----------------------------------------------------------------------------=!
   MODULE reciprocal_vectors
!=----------------------------------------------------------------------------=!

     USE kinds, ONLY: DP
     USE gvecp
     USE gvecb
     USE gvecs
     USE gvecw

     IMPLICIT NONE
     SAVE

     ! ...   declare module-scope variables

     LOGICAL :: gzero  = .TRUE.   ! .TRUE. if the first G vectors on this processor is
                                  ! the null G vector ( i.e. |G| == 0 )
     INTEGER :: gstart = 2        ! index of the first G vectors whose module is greather
                                  ! than 0 . 
                                  ! gstart = 2 when gzero == .TRUE., gstart = 1 otherwise 

     !     G^2 in increasing order (in units of tpiba2=(2pi/a)^2)
     !
     REAL(DP), ALLOCATABLE, TARGET :: g(:) 

     !     shells of G^2
     !
     REAL(DP), ALLOCATABLE, TARGET :: gl(:) 

     !     G-vectors cartesian components ( units tpiba =(2pi/a)  )
     !
     REAL(DP), ALLOCATABLE, TARGET :: gx(:,:) 

     !     g2_g    = all G^2 in increasing order, replicated on all procs
     !
     REAL(DP), ALLOCATABLE, TARGET :: g2_g(:)

     !     mill_g  = miller index of G vecs (increasing order), replicated on all procs
     !
     INTEGER, ALLOCATABLE, TARGET :: mill_g(:,:)

     !     mill_l  = miller index of G vecs local to the processors
     !
     INTEGER, ALLOCATABLE, TARGET :: mill_l(:,:)

     !     ig_l2g  = "l2g" means local to global, this array convert a local
     !               G-vector index into the global index, in other words
     !               the index of the G-v. in the overall array of G-vectors
     !
     INTEGER, ALLOCATABLE, TARGET :: ig_l2g(:)

     !     sortedig_l2g = array obtained by sorting ig_l2g
     !
     !
     INTEGER, ALLOCATABLE, TARGET :: sortedig_l2g(:)

     !     igl = index of the g-vector shells
     !
     INTEGER, ALLOCATABLE, TARGET :: igl(:)

     !     bi  = base vector used to generate the reciprocal space
     !
     REAL(DP) :: bi1(3) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     REAL(DP) :: bi2(3) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     REAL(DP) :: bi3(3) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)

   CONTAINS

     SUBROUTINE deallocate_recvecs
       IF( ALLOCATED( g ) ) DEALLOCATE( g )
       IF( ALLOCATED( gl ) ) DEALLOCATE( gl )
       IF( ALLOCATED( gx ) ) DEALLOCATE( gx )
       IF( ALLOCATED( g2_g ) ) DEALLOCATE( g2_g )
       IF( ALLOCATED( mill_g ) ) DEALLOCATE( mill_g )
       IF( ALLOCATED( mill_l ) ) DEALLOCATE( mill_l )
       IF( ALLOCATED( ig_l2g ) ) DEALLOCATE( ig_l2g )
       IF( ALLOCATED( sortedig_l2g ) ) DEALLOCATE( sortedig_l2g )
       IF( ALLOCATED( igl ) ) DEALLOCATE( igl )
       CALL deallocate_gvecw( )
       CALL deallocate_gvecs( )
       CALL deallocate_gvecb( )
     END SUBROUTINE deallocate_recvecs

!=----------------------------------------------------------------------------=!
   END MODULE reciprocal_vectors
!=----------------------------------------------------------------------------=!


!=----------------------------------------------------------------------------=!
   MODULE recvecs_indexes
!=----------------------------------------------------------------------------=!

     IMPLICIT NONE
     SAVE

     !     np      = fft index for G>
     !     nm      = fft index for G<
     !     in1p,in2p,in3p = G components in crystal axis


     INTEGER, ALLOCATABLE :: np(:), nm(:), in1p(:), in2p(:), in3p(:)

   CONTAINS

     SUBROUTINE deallocate_recvecs_indexes
       IF( ALLOCATED( np ) ) DEALLOCATE( np )
       IF( ALLOCATED( nm ) ) DEALLOCATE( nm )
       IF( ALLOCATED( in1p ) ) DEALLOCATE( in1p )
       IF( ALLOCATED( in2p ) ) DEALLOCATE( in2p )
       IF( ALLOCATED( in3p ) ) DEALLOCATE( in3p )
      END SUBROUTINE deallocate_recvecs_indexes

!=----------------------------------------------------------------------------=!
   END MODULE recvecs_indexes
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
       USE gvecw, ONLY: ngw, ngwx, ngwt
       USE gvecp, ONLY: ngm, ngmx, ngmt
       USE gvecs, ONLY: ngs, ngsx, ngst

       IMPLICIT NONE
       INTEGER, INTENT(IN) :: ngm_ , ngw_ , ngs_

       ngm = ngm_
       ngw = ngw_
       ngs = ngs_

       !
       !  calculate maxima over all processors
       !
       ngwx = ngw
       ngmx = ngm
       ngsx = ngs
       CALL mp_max( ngwx, intra_image_comm )
       CALL mp_max( ngmx, intra_image_comm )
       CALL mp_max( ngsx, intra_image_comm )
       !
       !  calculate SUM over all processors
       !
       ngwt = ngw
       ngmt = ngm
       ngst = ngs
       CALL mp_sum( ngwt, intra_image_comm )
       CALL mp_sum( ngmt, intra_image_comm )
       CALL mp_sum( ngst, intra_image_comm )

       RETURN 
     END SUBROUTINE recvecs_init

!=----------------------------------------------------------------------------=!
   END MODULE recvecs_subroutines
!=----------------------------------------------------------------------------=!
