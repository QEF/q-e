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

!  ... Most important dependencies: next three modules
   USE gvect,              ONLY : ig_l2g, g, gg, ngm, ngm_g, gcutm, &
                                  mill,  nl, gstart
   USE gvecs,              ONLY : ngms, gcutms, ngms_g, nls
   USE fft_base,           ONLY : dfftp, dffts
!
   USE kinds,              ONLY : DP
   USE constants,          ONLY : eps8

   PRIVATE
   SAVE

   PUBLIC :: ggen

!=----------------------------------------------------------------------=
CONTAINS
!=----------------------------------------------------------------------=
!
   !-----------------------------------------------------------------------
   SUBROUTINE ggen ( gamma_only, at, bg, comm, no_global_sort )
   !----------------------------------------------------------------------
   !
   !     This routine generates all the reciprocal lattice vectors
   !     contained in the sphere of radius gcutm. Furthermore it
   !     computes the indices nl which give the correspondence
   !     between the fft mesh points and the array of g vectors.
   !
   USE mp, ONLY: mp_rank, mp_size, mp_sum
   !
   IMPLICIT NONE
   !
   LOGICAL,  INTENT(IN) :: gamma_only
   REAL(DP), INTENT(IN) :: at(3,3), bg(3,3)
   INTEGER,  OPTIONAL, INTENT(IN) :: comm
   LOGICAL,  OPTIONAL, INTENT(IN) :: no_global_sort
   !  if no_global_sort is present (and it is true) G vectors are sorted only
   !  locally and not globally. In this case no global array needs to be
   !  allocated and sorted: saves memory and a lot of time for large systems.
   !
   !     here a few local variables
   !
   REAL(DP) ::  t (3), tt
   INTEGER :: ngm_save, ngms_save, n1, n2, n3, n1s, n2s, n3s, ngm_offset, ngm_max, ngms_max
   !
   REAL(DP), ALLOCATABLE :: g2sort_g(:)
   ! array containing all g vectors, on all processors: replicated data
   ! when no_global_sort is present (and it is true) only g vectors for the current processor are stored
   INTEGER, ALLOCATABLE :: mill_g(:,:), mill_unsorted(:,:)
   ! array containing all g vectors generators, on all processors: replicated data
   ! when no_global_sort is present (and it is true) only g vectors for the current processor are stored
   INTEGER, ALLOCATABLE :: igsrt(:)
   !
   INTEGER :: m1, m2, mc
   INTEGER :: ni, nj, nk, i, j, k, ipol, ng, igl, indsw
   INTEGER :: mype, npe
   LOGICAL :: global_sort
   INTEGER, ALLOCATABLE :: ngmpe(:)
   !
   IF( PRESENT( no_global_sort ) .AND. .NOT. PRESENT( comm ) ) THEN
      CALL errore ('ggen', ' wrong subroutine arguments, communicator is missing ', 1)
   END IF
   IF( .NOT. PRESENT( no_global_sort ) .AND. PRESENT( comm ) ) THEN
      CALL errore ('ggen', ' wrong subroutine arguments, parameter no_global_sort is missing ', 1)
   END IF
   !
   global_sort = .TRUE.
   !
   IF( PRESENT( no_global_sort ) ) THEN
      global_sort = .NOT. no_global_sort
   END IF
   !
   IF( .NOT. global_sort ) THEN
      mype = mp_rank( comm )
      npe  = mp_size( comm )
      ALLOCATE( ngmpe( npe ) )
      ngmpe = 0
      ngm_max = ngm
      ngms_max = ngms
   ELSE
      ngm_max = ngm_g
      ngms_max = ngms_g
   END IF
   !
   ! save current value of ngm and ngms
   !
   ngm_save  = ngm
   ngms_save = ngms
   !
   ngm = 0
   ngms = 0
   !
   ! counters
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
            m1 = mod (i, dfftp%nr1) + 1
            IF (m1 < 1) m1 = m1 + dfftp%nr1
            m2 = mod (j, dfftp%nr2) + 1
            IF (m2 < 1) m2 = m2 + dfftp%nr2
            mc = m1 + (m2 - 1) * dfftp%nr1x
            IF ( dfftp%isind ( mc ) == 0) CYCLE jloop
         END IF

         kloop: DO k = -nk, nk
            !
            ! gamma-only: exclude line with x = 0, y = 0, z < 0
            !
            IF ( gamma_only .and. i == 0 .and. j == 0 .and. k < 0) CYCLE kloop
            t(:) = i * bg (:,1) + j * bg (:,2) + k * bg (:,3)
            !tt = sum(t(:)**2)
            tt = t(1)**2+t(2)**2+t(3)**2
            IF (tt <= gcutm) THEN
               ngm = ngm + 1
               IF (tt <= gcutms) ngms = ngms + 1
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

   IF( .NOT. global_sort ) THEN
      ngmpe( mype + 1 ) = ngm
      CALL mp_sum( ngmpe, comm )
   END IF
   IF (ngm  /= ngm_max) &
         CALL errore ('ggen', 'g-vectors missing !', abs(ngm - ngm_max))
   IF (ngms /= ngms_max) &
         CALL errore ('ggen', 'smooth g-vectors missing !', abs(ngms - ngms_max))

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
      ! compute adeguate offsets in order to avoid overlap between
      ! g vectors once they are gathered on a single (global) array
      !
      ngm_offset = 0
      DO ng = 1, mype
         ngm_offset = ngm_offset + ngmpe( ng )
      END DO
   END IF

   ngm = 0
   ngms = 0
   !
   ngloop: DO ng = 1, ngm_max

      i = mill_g(1, ng)
      j = mill_g(2, ng)
      k = mill_g(3, ng)

#if defined(__MPI)
      IF( global_sort ) THEN
         m1 = mod (i, dfftp%nr1) + 1
         IF (m1 < 1) m1 = m1 + dfftp%nr1
         m2 = mod (j, dfftp%nr2) + 1
         IF (m2 < 1) m2 = m2 + dfftp%nr2
         mc = m1 + (m2 - 1) * dfftp%nr1x
         IF ( dfftp%isind ( mc ) == 0) CYCLE ngloop
      END IF
#endif

      ngm = ngm + 1

      !  Here map local and global g index !!!
      !  N.B. the global G vectors arrangement depends on the number of processors
      !
      IF( .NOT. global_sort ) THEN
         ig_l2g( ngm ) = ng + ngm_offset
      ELSE
         ig_l2g( ngm ) = ng
      END IF

      g (1:3, ngm) = i * bg (:, 1) + j * bg (:, 2) + k * bg (:, 3)
      gg (ngm) = sum(g (1:3, ngm)**2)

      IF (gg (ngm) <= gcutms) ngms = ngms + 1
      IF (ngm > ngm_save) CALL errore ('ggen 2', 'too many g-vectors', ngm)
   ENDDO ngloop

   IF (ngm /= ngm_save) &
      CALL errore ('ggen', 'g-vectors (ngm) missing !', abs(ngm - ngm_save))
   IF (ngms /= ngms_save) &
      CALL errore ('ggen', 'g-vectors (ngms) missing !', abs(ngm - ngms_save))
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
   DO ng = 1, ngm
      n1 = nint (sum(g (:, ng) * at (:, 1))) + 1
      mill (1,ng) = n1 - 1
      n1s = n1
      IF (n1<1) n1 = n1 + dfftp%nr1
      IF (n1s<1) n1s = n1s + dffts%nr1

      n2 = nint (sum(g (:, ng) * at (:, 2))) + 1
      mill (2,ng) = n2 - 1
      n2s = n2
      IF (n2<1) n2 = n2 + dfftp%nr2
      IF (n2s<1) n2s = n2s + dffts%nr2

      n3 = nint (sum(g (:, ng) * at (:, 3))) + 1
      mill (3,ng) = n3 - 1
      n3s = n3
      IF (n3<1) n3 = n3 + dfftp%nr3
      IF (n3s<1) n3s = n3s + dffts%nr3

      IF (n1>dfftp%nr1 .or. n2>dfftp%nr2 .or. n3>dfftp%nr3) &
         CALL errore('ggen','Mesh too small?',ng)

#if defined (__MPI) && !defined (__USE_3D_FFT)
      nl (ng) = n3 + ( dfftp%isind (n1 + (n2 - 1) * dfftp%nr1x) - 1) * dfftp%nr3x
      IF (ng <= ngms) &
         nls (ng) = n3s + ( dffts%isind (n1s+(n2s-1)*dffts%nr1x) - 1 ) * dffts%nr3x
#else
      nl (ng) = n1 + (n2 - 1) * dfftp%nr1x + (n3 - 1) * dfftp%nr1x * dfftp%nr2x
      IF (ng <= ngms) &
         nls (ng) = n1s + (n2s - 1) * dffts%nr1x + (n3s - 1) * dffts%nr1x * dffts%nr2x
#endif
   ENDDO
   !
   DEALLOCATE( mill_g )

   IF ( gamma_only) CALL index_minusg()

   IF( ALLOCATED( ngmpe ) ) DEALLOCATE( ngmpe )

   END SUBROUTINE ggen
   !
   !-----------------------------------------------------------------------
   SUBROUTINE index_minusg()
   !----------------------------------------------------------------------
   !
   !     compute indices nlm and nlms giving the correspondence
   !     between the fft mesh points and -G (for gamma-only calculations)
   !
   USE gvect,    ONLY : ngm, nlm, mill
   USE gvecs,    ONLY : nlsm, ngms
   USE fft_base, ONLY : dfftp, dffts
   !
   IMPLICIT NONE
   !
   INTEGER :: n1, n2, n3, n1s, n2s, n3s, ng
   !
   DO ng = 1, ngm
      n1 = -mill (1,ng) + 1
      n1s = n1
      IF (n1 < 1) THEN
         n1 = n1 + dfftp%nr1
         n1s = n1s + dffts%nr1
      END IF

      n2 = -mill (2,ng) + 1
      n2s = n2
      IF (n2 < 1) THEN
         n2 = n2 + dfftp%nr2
         n2s = n2s + dffts%nr2
      END IF
      n3 = -mill (3,ng) + 1
      n3s = n3
      IF (n3 < 1) THEN
         n3 = n3 + dfftp%nr3
         n3s = n3s + dffts%nr3
      END IF

      IF (n1>dfftp%nr1 .or. n2>dfftp%nr2 .or. n3>dfftp%nr3) THEN
         CALL errore('index_minusg','Mesh too small?',ng)
      ENDIF

#if defined (__MPI) && !defined (__USE_3D_FFT)
      nlm(ng) = n3 + (dfftp%isind (n1 + (n2 - 1) * dfftp%nr1x) - 1) * dfftp%nr3x
      IF (ng<=ngms) &
         nlsm(ng) = n3s + (dffts%isind (n1s+(n2s-1) * dffts%nr1x) - 1) * dffts%nr3x
#else
      nlm(ng) = n1 + (n2 - 1) * dfftp%nr1x + (n3 - 1) * dfftp%nr1x * dfftp%nr2x
      IF (ng<=ngms) &
         nlsm(ng) = n1s + (n2s - 1) * dffts%nr1x + (n3s-1) * dffts%nr1x * dffts%nr2x
#endif
   ENDDO

   END SUBROUTINE index_minusg
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
   ! that gives the shell index ng for (lacal) G-vector of index ig
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
