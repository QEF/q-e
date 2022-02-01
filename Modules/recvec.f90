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
     !! Variables describing the reciprocal lattice vectors
     !! G vectors with |G|^2 < ecutrho, cut-off for charge density
     !! With gamma tricks, G-vectors are divided into two half-spheres,
     !! G> and G<, containing G and -G (G=0 is in G>)
     !! This is referred to as the "dense" (or "hard", or "thick") grid

     USE kinds, ONLY: DP

     IMPLICIT NONE
     SAVE

     INTEGER :: ngm = 0
     !! local  number of G vectors (on this processor) with gamma tricks, only
     !! vectors in G>
     INTEGER :: ngm_g = 0
     !! global number of G vectors (summed on all procs) in serial execution,
     !! \(\text{ngm}_g = \text{ngm}\).
     INTEGER :: ngl = 0
     !! number of G-vector shells
     INTEGER :: ngmx = 0
     !! local number of G vectors, maximum across all procs

     REAL(DP) :: ecutrho = 0.0_DP
     !! energy cut-off for charge density 
     REAL(DP) :: gcutm = 0.0_DP
     !! \(\text{ecutrho}/(2 \pi/a)^2\), cut-off for \(\|G\|^2\)

     INTEGER :: gstart = 2
     !! index of the first G vector whose module is > 0. Needed in parallel
     !execution: gstart=2 for the proc that holds G=0, gstart=1 for all others.

     REAL(DP), ALLOCATABLE, TARGET :: gg(:) 
     !! \(G^2\) in increasing order (in units of \(\text{tpiba2}=(2\pi/a)^2\) )
     REAL(DP), ALLOCATABLE         :: gg_d(:)
     !! device

     REAL(DP), POINTER, PROTECTED :: gl(:)
     !! gl(i) = i-th shell of G^2 (in units of tpiba2)
     INTEGER, ALLOCATABLE, TARGET, PROTECTED :: igtongl(:)
     !! igtongl(n) = shell index for n-th G-vector
     
     ! Duplicate of the above variables (new style duplication).
     REAL(DP), ALLOCATABLE                   :: gl_d(:)      ! device
     INTEGER, ALLOCATABLE, TARGET            :: igtongl_d(:) ! device
     !
     REAL(DP), ALLOCATABLE, TARGET :: g(:,:) 
     !! G-vectors cartesian components ( in units \(\text{tpiba} =(2\pi/a)\) )
     REAL(DP), ALLOCATABLE         :: g_d(:,:)   ! device
     !
     INTEGER, ALLOCATABLE, TARGET :: mill(:,:)
     !! miller index of G vectors (local to each processor)
     !! G(:) = mill(1)*bg(:,1)+mill(2)*bg(:,2)+mill(3)*bg(:,3) 
     !! where bg are the reciprocal lattice basis vectors.
     INTEGER, ALLOCATABLE         :: mill_d(:,:) ! device
     !
     INTEGER, ALLOCATABLE, TARGET :: ig_l2g(:)
     !! converts a local G-vector index into the global index
     !! ("l2g" means local to global): ig\_l2g(i) = index of i-th
     !! local G-vector in the global array of G-vectors
     !
     INTEGER, ALLOCATABLE, TARGET :: mill_g(:,:)
     !! Miller index of all G vectors
     !
     COMPLEX(DP), ALLOCATABLE :: eigts1(:,:)
     !! the phases \(e^{-iG\text{tau}_s}\) used to calculate structure factors.
     COMPLEX(DP), ALLOCATABLE :: eigts2(:,:), eigts3(:,:)
     COMPLEX(DP), ALLOCATABLE :: eigts1_d(:,:)    ! device
     COMPLEX(DP), ALLOCATABLE :: eigts2_d(:,:)
     COMPLEX(DP), ALLOCATABLE :: eigts3_d(:,:)
     !   
#if defined(__CUDA)
     attributes(DEVICE) :: gl_d, igtongl_d
     attributes(DEVICE) :: gg_d, g_d, mill_d, eigts1_d, eigts2_d, eigts3_d
#endif
   CONTAINS

     SUBROUTINE gvect_init( ngm_ , comm )
       !
       !! Set local and global dimensions, allocate arrays
       !
       USE mp, ONLY: mp_max, mp_sum
       USE control_flags, ONLY : use_gpu
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: ngm_
       INTEGER, INTENT(IN) :: comm
       !! communicator of the group on which g-vecs are distributed
       !
       ngm = ngm_
       !
       !  calculate maximum over all processors
       !
       ngmx = ngm
       CALL mp_max( ngmx, comm )
       !
       !  calculate sum over all processors
       !
       ngm_g = ngm
       CALL mp_sum( ngm_g, comm )
       !
       !  allocate arrays - only those that are always kept until the end
       !
       ALLOCATE( gg(ngm) )
       ALLOCATE( g(3, ngm) )
       ALLOCATE( mill(3, ngm) )
       ALLOCATE( ig_l2g(ngm) )
       ALLOCATE( igtongl(ngm) )
       !
       IF (use_gpu) THEN
          ALLOCATE( mill_d(3, ngm) )
          ALLOCATE( igtongl_d(ngm) )
          ALLOCATE( gl_d(ngm) )
       ENDIF  
       !$acc enter data create( mill(1:3, 1:ngm), g(1:3, 1:ngm) ) 
       !
       RETURN 
       !
     END SUBROUTINE gvect_init

     SUBROUTINE deallocate_gvect(vc)
       !! Deallocate G-vector related arrays.
       !
       USE control_flags, ONLY : use_gpu
       IMPLICIT NONE
       !
       LOGICAL, OPTIONAL, INTENT(IN) :: vc
       LOGICAL :: vc_
       !
       vc_ = .false.
       IF (PRESENT(vc)) vc_ = vc
       IF ( .NOT. vc_ ) THEN
          IF ( ASSOCIATED( gl ) ) DEALLOCATE ( gl )
       END IF
       !
       !$acc exit data delete(eigts1,  eigts2,  eigts3)
       !$acc exit data delete(mill, g)
       !  
       IF( ALLOCATED( gg ) )     DEALLOCATE( gg )
       IF( ALLOCATED( g ) )      DEALLOCATE( g )
       IF( ALLOCATED( mill_g ) ) DEALLOCATE( mill_g )
       IF( ALLOCATED( mill ) )   DEALLOCATE( mill )
       IF( ALLOCATED( igtongl )) DEALLOCATE( igtongl )
       IF( ALLOCATED( ig_l2g ) ) DEALLOCATE( ig_l2g )
       IF( ALLOCATED( eigts1 ) ) DEALLOCATE( eigts1 )
       IF( ALLOCATED( eigts2 ) ) DEALLOCATE( eigts2 )
       IF( ALLOCATED( eigts3 ) ) DEALLOCATE( eigts3 )
       !
       ! GPU vars
       IF (use_gpu) THEN
          IF (ALLOCATED( igtongl_d )) DEALLOCATE( igtongl_d )
          IF (ALLOCATED( gl_d ) )     DEALLOCATE( gl_d )
          IF (ALLOCATED( gg_d ) )     DEALLOCATE( gg_d )
          IF (ALLOCATED( g_d ) )      DEALLOCATE( g_d )
          IF (ALLOCATED( mill_d ) )   DEALLOCATE( mill_d )
          IF (ALLOCATED( eigts1_d ) ) DEALLOCATE( eigts1_d )
          IF (ALLOCATED( eigts2_d ) ) DEALLOCATE( eigts2_d )
          IF (ALLOCATED( eigts3_d ) ) DEALLOCATE( eigts3_d )
       ENDIF
       ! 
     END SUBROUTINE deallocate_gvect

     SUBROUTINE deallocate_gvect_exx()
       IF( ALLOCATED( gg ) )      DEALLOCATE( gg )
       IF( ALLOCATED( g ) )       DEALLOCATE( g )
       IF( ALLOCATED( mill ) )    DEALLOCATE( mill )
       IF( ALLOCATED( igtongl ) ) DEALLOCATE( igtongl )
       IF( ALLOCATED( ig_l2g ) )  DEALLOCATE( ig_l2g )
     END SUBROUTINE deallocate_gvect_exx
     !
     !-----------------------------------------------------------------------
     SUBROUTINE gshells ( vc )
        !----------------------------------------------------------------------
        !! Calculate number of G shells: ngl, and the index ng = igtongl(ig)
        !! that gives the shell index ng for (local) G-vector of index ig.
        !
        USE kinds,              ONLY : DP
        USE constants,          ONLY : eps8
        USE control_flags,      ONLY : use_gpu
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: vc
        !
        INTEGER :: ng, igl
        !
        IF ( vc ) THEN
           !
           ! in case of a variable cell run each G vector has its shell
           !
           ngl = ngm
           gl => gg
           DO ng = 1, ngm
              igtongl (ng) = ng
           ENDDO
        ELSE
           !
           ! G vectors are grouped in shells with the same norm
           !
           ngl = 1
           igtongl (1) = 1
           DO ng = 2, ngm
              IF (gg (ng) > gg (ng - 1) + eps8) THEN
                 ngl = ngl + 1
              ENDIF
              igtongl (ng) = ngl
           ENDDO

           ALLOCATE (gl( ngl))
           gl (1) = gg (1)
           igl = 1
           DO ng = 2, ngm
              IF (gg (ng) > gg (ng - 1) + eps8) THEN
                 igl = igl + 1
                 gl (igl) = gg (ng)
              ENDIF
           ENDDO

           IF (igl /= ngl) CALL errore ('gshells', 'igl <> ngl', ngl)

        ENDIF
        IF (use_gpu)      gl_d = gl
        IF (use_gpu) igtongl_d = igtongl
     END SUBROUTINE gshells
!=----------------------------------------------------------------------------=!
   END MODULE gvect
!=----------------------------------------------------------------------------=!

!=----------------------------------------------------------------------------=!
   MODULE gvecs
!=----------------------------------------------------------------------------=!
     !! G vectors with \(\|G\|^2 < 4\cdot\text{ecutwfc}\), cut-off for wavefunctions
     !! ("smooth" grid). Gamma tricks and units as for the "dense" grid.
     
     USE kinds, ONLY: DP

     IMPLICIT NONE
     SAVE

     ! ... G vectors with |G|^2 < 4*ecutwfc, cut-off for wavefunctions
     ! ... ("smooth" grid). Gamma tricks and units as for the "dense" grid
     !
     INTEGER :: ngms = 0
     !! local  number of smooth vectors (on this processor)
     INTEGER :: ngms_g=0
     !! global number of smooth vectors (summed on procs) 
     !! in serial execution this is equal to \(\text{ngms}\)
     INTEGER :: ngsx = 0
     !! local number of smooth vectors, max across procs

     REAL(DP) :: ecuts = 0.0_DP
     !! energy cut-off = \(4\cdot\text{ecutwfc}\)

     REAL(DP) :: gcutms= 0.0_DP
     !! ecuts/(2 \pi/a)^2, cut-off for |G|^2

     REAL(DP) :: dual = 0.0_DP
     !! \(\text{ecutrho}=\text{dual}\cdot\text{ecutwfc}
     LOGICAL  :: doublegrid = .FALSE.
     !! TRUE if smooth and dense grid differ. \(\text{doublegrid}=(\text{dual}>4)\)


   CONTAINS

     SUBROUTINE gvecs_init( ngs_ , comm )
       !! G-vector initialization
       USE mp, ONLY: mp_max, mp_sum
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: ngs_
       INTEGER, INTENT(IN) :: comm
       !! communicator of the group on which g-vecs are distributed
       !
       ngms = ngs_
       !
       !  calculate maximum over all processors
       !
       ngsx = ngms
       CALL mp_max( ngsx, comm )
       !
       !  calculate sum over all processors
       !
       ngms_g = ngms
       CALL mp_sum( ngms_g, comm )
       !
       !  allocate arrays 
       !
       ! ALLOCATE( nls (ngms) )
       ! ALLOCATE( nlsm(ngms) )
       !
       RETURN 
       !
     END SUBROUTINE gvecs_init

!=----------------------------------------------------------------------------=!
   END MODULE gvecs
!=----------------------------------------------------------------------------=!

