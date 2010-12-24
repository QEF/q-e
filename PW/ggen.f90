!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE ggen()
   !----------------------------------------------------------------------
   !
   !     This routine generates all the reciprocal lattice vectors
   !     contained in the sphere of radius gcutm. Furthermore it
   !     computes the indices nl which give the correspondence
   !     between the fft mesh points and the array of g vectors.
   !
   USE kinds,              ONLY : DP
   USE cell_base,          ONLY : at, bg
   USE gvect, ONLY : ig_l2g
   USE gvect,              ONLY : g, gg, ngm, ngm_g, gcutm, &
                                  mill,  nl, gstart, gl, ngl, igtongl
   USE gvecs,            ONLY : ngms, gcutms, ngms_g, nls
   USE control_flags,      ONLY : gamma_only
   USE cellmd,             ONLY : lmovecell
   USE constants,          ONLY : eps8
   USE fft_base,           ONLY : dfftp, dffts

   IMPLICIT NONE
   !
   !     here a few local variables
   !
   REAL(DP) ::  t (3), tt, swap
   !
   INTEGER :: ngmx, n1, n2, n3, n1s, n2s, n3s
   !
   REAL(DP), ALLOCATABLE :: g2sort_g(:)
   ! array containing all g vectors, on all processors: replicated data
   INTEGER, ALLOCATABLE :: mill_g(:,:), mill_unsorted(:,:)
   ! array containing all g vectors generators, on all processors:
   !     replicated data
   INTEGER, ALLOCATABLE :: igsrt(:)
   !
#ifdef __PARA
   INTEGER :: m1, m2, mc
   !
#endif
   INTEGER :: i, j, k, ipol, ng, igl, iswap, indsw
   !
   ! counters
   !
   !    set the total number of fft mesh points and and initial value of gg
   !    The choice of gcutm is due to the fact that we have to order the
   !    vectors after computing them.
   !
   gg(:) = gcutm + 1.d0
   !
   !     set d vector for unique ordering
   !
   !    and computes all the g vectors inside a sphere
   !
   ALLOCATE( ig_l2g( ngm ) )
   ALLOCATE( mill_g( 3, ngm_g ),mill_unsorted( 3, ngm_g ) )
   ALLOCATE( igsrt( ngm_g ) )
   ALLOCATE( g2sort_g( ngm_g ) )
   g2sort_g(:) = 1.0d20
   !
   ! save present value of ngm in ngmx variable
   !
   ngmx = ngm
   !
   ngm = 0
   ngms = 0
   iloop: DO i = -dfftp%nr1-1, dfftp%nr1+1
      !
      ! gamma-only: exclude space with x < 0
      !
      IF ( gamma_only .and. i < 0) CYCLE iloop
      jloop: DO j = -dfftp%nr2-1, dfftp%nr2+1
         !
         ! gamma-only: exclude plane with x = 0, y < 0
         !
         IF ( gamma_only .and. i == 0 .and. j < 0) CYCLE jloop
         kloop: DO k = -dfftp%nr3-1, dfftp%nr3+1
            !
            ! gamma-only: exclude line with x = 0, y = 0, z < 0
            !
            IF ( gamma_only .and. i == 0 .and. j == 0 .and. k < 0) CYCLE kloop
            t(:) = i * bg (:,1) + j * bg (:,2) + k * bg (:,3)
            tt = sum(t(:)**2)
            IF (tt <= gcutm) THEN
               ngm = ngm + 1
               IF (tt <= gcutms) ngms = ngms + 1
               IF (ngm > ngm_g) CALL errore ('ggen', 'too many g-vectors', ngm)
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

   IF (ngm  /= ngm_g ) &
         CALL errore ('ggen', 'g-vectors missing !', abs(ngm - ngm_g))
   IF (ngms /= ngms_g) &
         CALL errore ('ggen', 'smooth g-vectors missing !', abs(ngms - ngms_g))

   igsrt(1) = 0
   CALL hpsort_eps( ngm_g, g2sort_g, igsrt, eps8 )
   mill_g(1,:) = mill_unsorted(1,igsrt(:))
   mill_g(2,:) = mill_unsorted(2,igsrt(:))
   mill_g(3,:) = mill_unsorted(3,igsrt(:))
   DEALLOCATE( g2sort_g, igsrt, mill_unsorted )

   ngm = 0
   ngms = 0
   ngloop: DO ng = 1, ngm_g
      i = mill_g(1, ng)
      j = mill_g(2, ng)
      k = mill_g(3, ng)

#ifdef __PARA
      m1 = mod (i, dfftp%nr1) + 1
      IF (m1 < 1) m1 = m1 + dfftp%nr1
      m2 = mod (j, dfftp%nr2) + 1
      IF (m2 < 1) m2 = m2 + dfftp%nr2
      mc = m1 + (m2 - 1) * dfftp%nr1x
      IF ( dfftp%isind ( mc ) == 0) CYCLE ngloop
#endif

      ngm = ngm + 1

      !  Here map local and global g index !!!
      ig_l2g( ngm ) = ng

      g (1:3, ngm) = i * bg (:, 1) + j * bg (:, 2) + k * bg (:, 3)
      gg (ngm) = sum(g (1:3, ngm)**2)

      IF (gg (ngm) <= gcutms) ngms = ngms + 1
      IF (ngm > ngmx) CALL errore ('ggen', 'too many g-vectors', ngm)
   ENDDO ngloop

   IF (ngm /= ngmx) &
      CALL errore ('ggen', 'g-vectors missing !', abs(ngm - ngmx))

   !
   !  here to initialize berry_phase
   !  CALL berry_setup(ngm, ngm_g, nr1, nr2, nr3, mill_g)
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

#if defined (__PARA) && !defined (__USE_3D_FFT)
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
   !
   ! calculate number of G shells: ngl
   !
   IF (lmovecell) THEN
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

      IF (igl /= ngl) CALL errore ('setup', 'igl <> ngl', ngl)

   ENDIF

   IF ( gamma_only) CALL index_minusg()

END SUBROUTINE ggen

!
!-----------------------------------------------------------------------
SUBROUTINE index_minusg()
   !----------------------------------------------------------------------
   !
   !     compute indices nlm and nlms giving the correspondence
   !     between the fft mesh points and -G (for gamma-only calculations)
   !
   USE gvect,   ONLY : ngm, nlm, mill
   USE gvecs, ONLY : nlsm, ngms
   USE fft_base,  ONLY : dfftp, dffts
   IMPLICIT NONE
   !
   INTEGER :: n1, n2, n3, n1s, n2s, n3s, ng
   !
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

#if defined (__PARA) && !defined (__USE_3D_FFT)
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
