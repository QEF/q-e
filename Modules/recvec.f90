!
! Copyright (C) 2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------------=!
   MODULE gvect
!=----------------------------------------------------------------------------=!

     ! ... variables describing the reciprocal lattice vectors
     ! ... G vectors with |G|^2 < ecutrho, cut-off for charge density

     USE kinds, ONLY: DP

     IMPLICIT NONE
     SAVE

     INTEGER :: ngm  = 0  ! local  number of G vectors (on this processor)
                          ! with gamma tricks, only half sphere (G>)
     INTEGER :: ngm_g= 0  ! global number of G vectors (summed on all procs)
                          ! in serial execution, ngm_g = ngm
     INTEGER :: ngl = 0   ! number of G-vector shells
     INTEGER :: ngmx = 0  ! local number of G vectors, maximum across all procs

     REAL(DP) :: ecutrho = 0.0_DP ! energy cut-off for charge density 
     REAL(DP) :: gcutm = 0.0_DP   ! ecutrho/(2 pi/a)^2, cut-off for G|^2

     ! nl  = fft index for G-vectors (with gamma tricks, only G>)
     ! nlm = as above, for -G (with gamma tricks)

     INTEGER, ALLOCATABLE :: nl(:), nlm(:)

     INTEGER :: gstart = 2 ! index of the first G vector whose module is > 0
                           ! Needed in parallel execution: gstart=2 for the
                           ! proc that holds G=0, gstart=1 for all others

     !     G^2 in increasing order (in units of tpiba2=(2pi/a)^2)
     !
     REAL(DP), ALLOCATABLE, TARGET :: gg(:) 

     !     gl(i) = i-th shell of G^2
     !     igtongl(n) = shell index for n-th G vector
     !
     REAL(DP), POINTER :: gl(:) 
     INTEGER, ALLOCATABLE, TARGET :: igtongl(:) 
     !
     !     G-vectors cartesian components ( units tpiba =(2pi/a)  )
     !
     REAL(DP), ALLOCATABLE, TARGET :: g(:,:) 

     !     mill    = miller index of G vectors (local to each processor)
     !
     INTEGER, ALLOCATABLE, TARGET :: mill(:,:)
     
     !     g2_g    = all G^2 in increasing order, replicated on all procs
     !
     REAL(DP), ALLOCATABLE, TARGET :: g2_g(:)

     !     mill_g  = miller index of G vecs (increasing order), replicated on all procs
     !
     INTEGER, ALLOCATABLE, TARGET :: mill_g(:,:)

     !     ig_l2g  = "l2g" means local to global, this array convert a local
     !               G-vector index into the global index, in other words
     !               the index of the G-v. in the overall array of G-vectors
     !
     INTEGER, ALLOCATABLE, TARGET :: ig_l2g(:)

     !     sortedig_l2g = array obtained by sorting ig_l2g
     !
     INTEGER, ALLOCATABLE, TARGET :: sortedig_l2g(:)
     !
     ! the phases e^{-iG*tau_s}
     !
     COMPLEX(DP), ALLOCATABLE :: eigts1(:,:), eigts2(:,:), eigts3(:,:)
     !
   CONTAINS

     SUBROUTINE deallocate_gvect()
       IF( ALLOCATED( gg ) ) DEALLOCATE( gg )
       IF( ASSOCIATED( gl ) ) DEALLOCATE( gl )
       IF( ALLOCATED( g ) )  DEALLOCATE( g )
       IF( ALLOCATED( g2_g ) ) DEALLOCATE( g2_g )
       IF( ALLOCATED( mill_g ) ) DEALLOCATE( mill_g )
       IF( ALLOCATED( mill ) ) DEALLOCATE( mill )
       IF( ALLOCATED( igtongl ) ) DEALLOCATE( igtongl )
       IF( ALLOCATED( ig_l2g ) ) DEALLOCATE( ig_l2g )
       IF( ALLOCATED( sortedig_l2g ) ) DEALLOCATE( sortedig_l2g )
       IF( ALLOCATED( eigts1 ) ) DEALLOCATE( eigts1 )
       IF( ALLOCATED( eigts2 ) ) DEALLOCATE( eigts2 )
       IF( ALLOCATED( eigts3 ) ) DEALLOCATE( eigts3 )
       IF( ALLOCATED( nl ) ) DEALLOCATE( nl )
       IF( ALLOCATED( nlm ) ) DEALLOCATE( nlm )
     END SUBROUTINE deallocate_gvect

!=----------------------------------------------------------------------------=!
   END MODULE gvect
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
