!
! Copyright (C) 2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!=----------------------------------------------------------------------=
MODULE recvec_subs
  !=----------------------------------------------------------------------=

  !  ... subroutines generating G-vectors and variables nl* needed to map
  !  ... G-vector components onto the FFT grid(s) in reciprocal space
  !
  USE kinds, ONLY : dp
  USE fft_types, ONLY: fft_stick_index, fft_type_descriptor
  USE fft_ggen, ONLY : fft_set_nl
  !  
  PRIVATE
  SAVE
  
  PUBLIC :: ggen, ggens

  !=----------------------------------------------------------------------=
CONTAINS
  !=----------------------------------------------------------------------=
  !
  !-----------------------------------------------------------------------
  SUBROUTINE ggen ( dfftp, gamma_only, at, bg,  gcutm, ngm_g, ngm, &
       g, gg, mill, ig_l2g, gstart, no_global_sort )
    !----------------------------------------------------------------------
    !
    !     This routine generates all the reciprocal lattice vectors
    !     contained in the sphere of radius gcutm. Furthermore it
    !     computes the indices nl which give the correspondence
    !     between the fft mesh points and the array of g vectors.
    !
    USE mp, ONLY: mp_rank, mp_size, mp_sum
    USE constants, ONLY : eps8
    !
    IMPLICIT NONE
    !
    TYPE(fft_type_descriptor),INTENT(INOUT) :: dfftp
    LOGICAL,  INTENT(IN) :: gamma_only
    REAL(DP), INTENT(IN) :: at(3,3), bg(3,3), gcutm
    INTEGER, INTENT(IN) :: ngm_g
    INTEGER, INTENT(INOUT) :: ngm
    REAL(DP), INTENT(OUT) :: g(:,:), gg(:)
    INTEGER, INTENT(OUT) :: mill(:,:), ig_l2g(:), gstart
    !  if no_global_sort is present (and it is true) G vectors are sorted only
    !  locally and not globally. In this case no global array needs to be
    !  allocated and sorted: saves memory and a lot of time for large systems.
    !
    LOGICAL,  OPTIONAL, INTENT(IN) :: no_global_sort
    !
    !     here a few local variables
    !
    REAL(DP) ::  t (3), tt
    INTEGER :: ngm_save, n1, n2, n3, ngm_offset, ngm_max
    !
    REAL(DP), ALLOCATABLE :: g2sort_g(:)
    ! array containing all g vectors, on all processors: replicated data
    ! when no_global_sort is present (and it is true) only g vectors for
    ! the current processor are stored
    INTEGER, ALLOCATABLE :: mill_g(:,:), mill_unsorted(:,:)
    ! array containing all g vectors generators, on all processors
    ! (replicated data). When no_global_sort is present and .true.,
    ! only g-vectors for the current processor are stored
    INTEGER, ALLOCATABLE :: igsrt(:)
    !
    INTEGER :: ni, nj, nk, i, j, k, ipol, ng, igl, indsw
    INTEGER :: mype, npe
    LOGICAL :: global_sort
    INTEGER, ALLOCATABLE :: ngmpe(:)
    !
    global_sort = .TRUE.
    IF( PRESENT( no_global_sort ) ) THEN
       global_sort = .NOT. no_global_sort
    END IF
    !
    IF( .NOT. global_sort ) THEN
       ngm_max = ngm
    ELSE
       ngm_max = ngm_g
    END IF
    !
    ! save current value of ngm
    !
    ngm_save  = ngm
    !
    ngm = 0
    !
    !    set the total number of fft mesh points and and initial value of gg
    !    The choice of gcutm is due to the fact that we have to order the
    !    vectors after computing them.
    !
    gg(:) = gcutm + 1.d0
    !
    !    and computes all the g vectors inside a sphere
    !
    ALLOCATE( mill_g( 3, ngm_max ),mill_unsorted( 3, ngm_max ) )
    ALLOCATE( igsrt( ngm_max ) )
    ALLOCATE( g2sort_g( ngm_max ) )
    !
    g2sort_g(:) = 1.0d20
    !
    ! max miller indices (same convention as in module stick_set)
    !
    ni = (dfftp%nr1-1)/2
    nj = (dfftp%nr2-1)/2
    nk = (dfftp%nr3-1)/2
    !
    iloop: DO i = -ni, ni
       !
       ! gamma-only: exclude space with x < 0
       !
       IF ( gamma_only .and. i < 0) CYCLE iloop
       jloop: DO j = -nj, nj
          !
          ! gamma-only: exclude plane with x = 0, y < 0
          !
          IF ( gamma_only .and. i == 0 .and. j < 0) CYCLE jloop
          
          IF( .NOT. global_sort ) THEN
             IF ( fft_stick_index( dfftp, i, j ) == 0) CYCLE jloop
          END IF
          
          kloop: DO k = -nk, nk
             !
             ! gamma-only: exclude line with x = 0, y = 0, z < 0
             !
             IF ( gamma_only .and. i == 0 .and. j == 0 .and. k < 0) CYCLE kloop
             t(:) = i * bg (:,1) + j * bg (:,2) + k * bg (:,3)
             tt = t(1)**2+t(2)**2+t(3)**2
             IF (tt <= gcutm) THEN
                ngm = ngm + 1
                IF (ngm > ngm_max) CALL errore ('ggen 1', 'too many g-vectors', ngm)
                mill_unsorted( :, ngm ) = (/ i,j,k /)
                IF ( tt > eps8 ) THEN
                   g2sort_g(ngm) = tt
                ELSE
                   g2sort_g(ngm) = 0.d0
                ENDIF
             ENDIF
          ENDDO kloop
       ENDDO jloop
    ENDDO iloop
    
    IF (ngm  /= ngm_max) &
         CALL errore ('ggen', 'g-vectors missing !', abs(ngm - ngm_max))
    !
    igsrt(1) = 0
    IF( .NOT. global_sort ) THEN
       CALL hpsort_eps( ngm, g2sort_g, igsrt, eps8 )
    ELSE
       CALL hpsort_eps( ngm_g, g2sort_g, igsrt, eps8 )
    END IF
    mill_g(1,:) = mill_unsorted(1,igsrt(:))
    mill_g(2,:) = mill_unsorted(2,igsrt(:))
    mill_g(3,:) = mill_unsorted(3,igsrt(:))
    DEALLOCATE( g2sort_g, igsrt, mill_unsorted )
    
    IF( .NOT. global_sort ) THEN
       !
       ! compute adeguate offsets in order to avoid overlap between
       ! g vectors once they are gathered on a single (global) array
       !
       mype = mp_rank( dfftp%comm )
       npe  = mp_size( dfftp%comm )
       ALLOCATE( ngmpe( npe ) )
       ngmpe = 0
       ngmpe( mype + 1 ) = ngm
       CALL mp_sum( ngmpe, dfftp%comm )
       ngm_offset = 0
       DO ng = 1, mype
          ngm_offset = ngm_offset + ngmpe( ng )
       END DO
       DEALLOCATE( ngmpe )
       !
    END IF
    
    ngm = 0
    !
    ngloop: DO ng = 1, ngm_max
       
       i = mill_g(1, ng)
       j = mill_g(2, ng)
       k = mill_g(3, ng)
       
       IF( dfftp%lpara .AND. global_sort ) THEN
          IF ( fft_stick_index( dfftp, i, j ) == 0) CYCLE ngloop
       END IF
       
       ngm = ngm + 1
       
       !  Here map local and global g index !!! N.B: :
       !  the global G vectors arrangement depends on the number of processors
       !
       IF( .NOT. global_sort ) THEN
          ig_l2g( ngm ) = ng + ngm_offset
       ELSE
          ig_l2g( ngm ) = ng
       END IF
       
       g (1:3, ngm) = i * bg (:, 1) + j * bg (:, 2) + k * bg (:, 3)
       gg (ngm) = sum(g (1:3, ngm)**2)
       
    ENDDO ngloop
    
    DEALLOCATE( mill_g )
    
    IF (ngm /= ngm_save) &
         CALL errore ('ggen', 'g-vectors (ngm) missing !', abs(ngm - ngm_save))
    !
    !     determine first nonzero g vector
    !
    IF (gg(1).le.eps8) THEN
       gstart=2
    ELSE
       gstart=1
    ENDIF
    !
    !     Now set nl and nls with the correct fft correspondence
    !
    CALL fft_set_nl( dfftp, at, g, mill )
    !
  END SUBROUTINE ggen
  !
  !-----------------------------------------------------------------------
  SUBROUTINE ggens( dffts, gamma_only, at, g, gg, mill, gcutms, ngms, &
       gs, ggs )
    !-----------------------------------------------------------------------
    !
    ! Initialize number and indices of  g-vectors for a subgrid,
    ! for exactly the same ordering as for the original FFT grid
    !
    !--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: gamma_only
    TYPE (fft_type_descriptor), INTENT(INOUT) :: dffts
    ! primitive lattice vectors
    REAL(dp), INTENT(IN) :: at(3,3)
    ! G-vectors in FFT grid
    REAL(dp), INTENT(IN) :: g(:,:), gg(:)    
    ! Miller indices for G-vectors of FFT grid
    INTEGER, INTENT(IN) :: mill(:,:)
    ! cutoff for subgrid
    REAL(DP), INTENT(IN):: gcutms
    ! Local number of G-vectors in subgrid
    INTEGER, INTENT(OUT):: ngms
    ! Optionally: G-vectors and modules
    REAL(DP), INTENT(OUT), POINTER, OPTIONAL:: gs(:,:), ggs(:)
    !
    INTEGER :: i, ng, ngm
    !
    ngm  = SIZE(gg)
    ngms = dffts%ngm
    IF ( ngms > ngm  ) CALL errore ('ggens','wrong  number of G-vectors',1)
    !
    IF ( PRESENT(gs) ) ALLOCATE ( gs(3,ngms) )
    IF ( PRESENT(ggs)) ALLOCATE ( ggs(ngms) )
    ng = 0
    DO i = 1, ngm
       IF ( gg(i) > gcutms ) exit
       IF ( PRESENT(gs) ) gs (:,i) = g(:,i)
       IF ( PRESENT(ggs)) ggs(i) = gg(i)
       ng = i
    END DO
    IF ( ng /= ngms ) CALL errore ('ggens','mismatch in number of G-vectors',2)
    !
    CALL fft_set_nl ( dffts, at, g )
    !
  END SUBROUTINE ggens
  !
  !=----------------------------------------------------------------------=
END MODULE recvec_subs
!=----------------------------------------------------------------------=
!
!-----------------------------------------------------------------------
SUBROUTINE gshells ( vc )
   !----------------------------------------------------------------------
   !
   ! calculate number of G shells: ngl, and the index ng = igtongl(ig)
   ! that gives the shell index ng for (local) G-vector of index ig
   !
   USE kinds,              ONLY : DP
   USE gvect,              ONLY : gg, ngm, gl, ngl, igtongl
   USE constants,          ONLY : eps8
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

   END SUBROUTINE gshells
