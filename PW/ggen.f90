!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
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
  USE reciprocal_vectors, ONLY : ig_l2g
  USE gvect,              ONLY : g, gg, ngm, ngm_g, ngm_l, nr1, nr2, nr3, &
                                 gcutm, nrx1, nrx2, nrx3, ig1, ig2, ig3,  &
                                 nl, gstart, gl, ngl, igtongl
  USE gsmooth,            ONLY : ngms, gcutms, ngms_g, nr1s, nr2s, nr3s, &
                                 nrx1s, nrx3s, nls
  USE control_flags,      ONLY : gamma_only
  USE cellmd,             ONLY : lmovecell
  USE constants,          ONLY : eps8
  USE sticks,             ONLY : dfftp, dffts

  IMPLICIT NONE
  !
  !     here a few local variables
  !
  REAL(DP) ::  t (3), tt, swap
  REAL(DP), ALLOCATABLE ::  esort (:)
  !
  INTEGER :: ngmx, n1, n2, n3, n1s, n2s, n3s
  !
  REAL(DP), ALLOCATABLE :: g2sort_g(:)
  ! array containing all g vectors, on all processors: replicated data
  INTEGER, ALLOCATABLE :: mill_g(:,:)
  ! array containing all g vectors generators, on all processors:
  !     replicated data
  INTEGER, ALLOCATABLE :: igsrt(:)
  !
#ifdef __PARA
  INTEGER :: m1, m2, m3, mc
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
  ALLOCATE( ig_l2g( ngm_l ) )
  ALLOCATE( mill_g( 3, ngm_g ) )
  ALLOCATE( igsrt( ngm_g ) )
  ALLOCATE( g2sort_g( ngm_g ) )
  g2sort_g(:) = 1.0d20
  !
  n1 = nr1 + 1
  n2 = nr2 + 1
  n3 = nr3 + 1
  !
  ! save present value of ngm in ngmx variable
  !
  ngmx = ngm
  !
  ngm = 0
  ngms = 0
  DO i = - n1, n1
     !
     ! Gamma-only: exclude space with x < 0
     !
     IF ( gamma_only .AND. i < 0) go to 10
     DO j = - n2, n2
        !
        ! exclude plane with x = 0, y < 0
        !
        IF ( gamma_only .AND. i == 0 .AND. j < 0) go to 11
        DO k = - n3, n3
           !
           ! exclude line with x = 0, y = 0, z < 0
           !
           IF ( gamma_only .AND. i == 0 .AND. j == 0 .AND. k < 0) go to 12
           tt = 0.d0
           DO ipol = 1, 3
              t (ipol) = i * bg (ipol, 1) + j * bg (ipol, 2) + k * bg (ipol, 3)
              tt = tt + t (ipol) * t (ipol)
           ENDDO
           IF (tt <= gcutm) THEN
              ngm = ngm + 1
              IF (tt <= gcutms) ngms = ngms + 1
              IF (ngm > ngm_g) CALL errore ('ggen', 'too many g-vectors', ngm)
              mill_g( 1, ngm ) = i
              mill_g( 2, ngm ) = j
              mill_g( 3, ngm ) = k
              IF ( tt > eps8 ) THEN
                 g2sort_g(ngm) = tt
              ELSE
                 g2sort_g(ngm) = 0.d0
              ENDIF
           END IF
12         CONTINUE
        ENDDO
11      CONTINUE
     ENDDO
10   CONTINUE
  ENDDO

  IF (ngm  /= ngm_g ) &
       CALL errore ('ggen', 'g-vectors missing !', ABS(ngm - ngm_g))
  IF (ngms /= ngms_g) &
       CALL errore ('ggen', 'smooth g-vectors missing !', ABS(ngms - ngms_g))

  igsrt(1) = 0
  CALL hpsort_eps( ngm_g, g2sort_g, igsrt, eps8 )
  DEALLOCATE( g2sort_g )
  DO ng = 1, ngm_g-1
    indsw = ng
7   IF(igsrt(indsw) /= ng) THEN
! ..  swap indices
      DO i = 1, 3
        iswap = mill_g(i,indsw)
        mill_g(i,indsw) = mill_g(i,igsrt(indsw))
        mill_g(i,igsrt(indsw)) = iswap
      END DO
! ..  swap indices
      iswap = indsw; indsw = igsrt(indsw); igsrt(iswap) = iswap
      IF(igsrt(indsw) == ng) THEN
        igsrt(indsw)=indsw
      ELSE
        GOTO 7
      END IF
    END IF
  END DO

  DEALLOCATE( igsrt )

  ! WRITE( stdout, fmt="(//,' --- Executing new GGEN Loop ---',//)" )

  ALLOCATE(esort(ngm) )
  esort(:) = 1.0d20
  ngm = 0
  ngms = 0
  DO ng = 1, ngm_g
    i = mill_g(1, ng)
    j = mill_g(2, ng)
    k = mill_g(3, ng)

#ifdef __PARA
    m1 = MOD (i, nr1) + 1
    IF (m1.LT.1) m1 = m1 + nr1
    m2 = MOD (j, nr2) + 1
    IF (m2.LT.1) m2 = m2 + nr2
    mc = m1 + (m2 - 1) * nrx1
    IF ( dfftp%isind ( mc ) .EQ.0) GOTO 1
#endif

    tt = 0.d0
    DO ipol = 1, 3
      t (ipol) = i * bg (ipol, 1) + j * bg (ipol, 2) + k * bg (ipol, 3)
      tt = tt + t (ipol) * t (ipol)
    ENDDO

    ngm = ngm + 1
    IF (tt <= gcutms) ngms = ngms + 1
    IF (ngm > ngmx) CALL errore ('ggen', 'too many g-vectors', ngm)
    !
    !  Here map local and global g index !!!
    !
    ig_l2g( ngm ) = ng
    !
    g (1:3, ngm) = t (1:3)
    gg (ngm) = tt

    IF (tt > eps8) THEN
      esort (ngm) = tt 
    ELSE
      esort (ngm) = 0.d0
    ENDIF

1   CONTINUE
  ENDDO

     IF (ngm.NE.ngmx) &
          CALL errore ('ggen', 'g-vectors missing !', ABS(ngm - ngmx))
     !
     !   reorder the g's in order of increasing magnitude. On exit
     !   from hpsort esort is ordered, and nl contains the new order.
     !
     !   initialize the index inside sorting routine

     nl (1) = 0
     CALL hpsort_eps ( ngm, esort, nl, eps8 )
     !
     DEALLOCATE( esort  )
     !
     !   reorder also the g vectors, and nl
     !
     DO ng = 1, ngm - 1
20      indsw = nl (ng)
        IF (indsw.NE.ng) THEN
           DO ipol = 1, 3
              swap = g (ipol, indsw)
              g (ipol, indsw) = g (ipol, nl (indsw) )
              g (ipol, nl (indsw) ) = swap
           ENDDO
           swap = gg (indsw)
           gg (indsw) = gg (nl (indsw) )
           gg (nl (indsw) ) = swap

          !
          !  Remember: ig_l2g is the index of a given G vectors in the
          !  sorted global array containing all G vectors, it is used to
          !  collect all wave function components
          !
          iswap = ig_l2g( indsw )
          ig_l2g( indsw ) = ig_l2g( nl(indsw) )
          ig_l2g( nl(indsw) ) = iswap

           iswap = nl (ng)
           nl (ng) = nl (indsw)
           nl (indsw) = iswap

           GOTO 20
        ENDIF

     ENDDO
     !
     !  here to initialize berry_phase
     !  work in progress ...
     !  CALL berry_setup(ngm, ngm_g, nr1, nr2, nr3, mill_g)
     !
     !     determine first nonzero g vector
     !
     IF (gg(1).LE.eps8) THEN
        gstart=2
     ELSE
        gstart=1
     END IF
     !
     !     Now set nl and nls with the correct fft correspondence
     !
     DO ng = 1, ngm
        n1 = NINT (g (1, ng) * at (1, 1) + g (2, ng) * at (2, 1) + g (3, &
             ng) * at (3, 1) ) + 1
        ig1 (ng) = n1 - 1
        n1s = n1
        IF (n1.LT.1) n1 = n1 + nr1
        IF (n1s.LT.1) n1s = n1s + nr1s
        n2 = NINT (g (1, ng) * at (1, 2) + g (2, ng) * at (2, 2) + g (3, &
             ng) * at (3, 2) ) + 1
        ig2 (ng) = n2 - 1
        n2s = n2
        IF (n2.LT.1) n2 = n2 + nr2
        IF (n2s.LT.1) n2s = n2s + nr2s
        n3 = NINT (g (1, ng) * at (1, 3) + g (2, ng) * at (2, 3) + g (3, &
             ng) * at (3, 3) ) + 1
        ig3 (ng) = n3 - 1
        n3s = n3
        IF (n3.LT.1) n3 = n3 + nr3
        IF (n3s.LT.1) n3s = n3s + nr3s
        IF (n1.LE.nr1.AND.n2.LE.nr2.AND.n3.LE.nr3) THEN
#if defined (__PARA) && !defined (__USE_3D_FFT)
           nl (ng) = n3 + ( dfftp%isind (n1 + (n2 - 1) * nrx1) - 1) * nrx3
           IF (ng.LE.ngms) nls (ng) = n3s + ( dffts%isind (n1s + (n2s - 1) &
                * nrx1s) - 1) * nrx3s
#else
           nl (ng) = n1 + (n2 - 1) * nrx1 + (n3 - 1) * nrx1 * nrx2
           IF (ng.LE.ngms) nls (ng) = n1s + (n2s - 1) * nrx1s + (n3s - 1) &
                * nrx1s * nr2s
#endif
        ELSE
           CALL errore('ggen','Mesh too small?',ng)
        ENDIF
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

        IF (igl.NE.ngl) CALL errore ('setup', 'igl <> ngl', ngl)

     ENDIF

     CALL index_minusg()

     RETURN
   END SUBROUTINE ggen

!
!-----------------------------------------------------------------------
SUBROUTINE index_minusg()
  !----------------------------------------------------------------------
  !
  !     compute indices nlm and nlms giving the correspondence
  !     between the fft mesh points and -G (for gamma-only calculations)
  !
  USE gvect,   ONLY : ngm, nr1, nr2, nr3, &
                      nrx1, nrx2, nrx3, nlM, ig1, ig2, ig3
  USE gsmooth, ONLY : nr1s, nr2s, nr3s, nrx1s, nrx3s, nlsm, ngms
  USE sticks,  ONLY : dfftp, dffts
  USE control_flags, ONLY : gamma_only
  IMPLICIT NONE
  !
  INTEGER :: n1, n2, n3, n1s, n2s, n3s, ng
  !
  !
  IF (gamma_only) THEN
     DO ng = 1, ngm
        n1 = -ig1 (ng) + 1
        n1s = n1
        IF (n1 < 1) n1 = n1 + nr1
        IF (n1s < 1) n1s = n1s + nr1s
        n2 = -ig2 (ng) + 1
        n2s = n2
        IF (n2 < 1) n2 = n2 + nr2
        IF (n2s < 1) n2s = n2s + nr2s
        n3 = -ig3 (ng) + 1
        n3s = n3
        IF (n3 < 1) n3 = n3 + nr3
        IF (n3s < 1) n3s = n3s + nr3s
        IF (n1.LE.nr1 .AND. n2.LE.nr2 .AND. n3.LE.nr3) THEN
#if defined (__PARA) && !defined (__USE_3D_FFT)
           nlm(ng) = n3 + (dfftp%isind (n1 + (n2 - 1) * nrx1) - 1) * nrx3
           IF (ng.LE.ngms) nlsm(ng) = n3s + (dffts%isind (n1s + (n2s - 1) &
                * nrx1s) - 1) * nrx3s
#else
           nlm(ng) = n1 + (n2 - 1) * nrx1 + (n3 - 1) * nrx1 * nrx2
           IF (ng.LE.ngms) nlsm(ng) = n1s + (n2s - 1) * nrx1s + (n3s - 1) &
                * nrx1s * nr2s
#endif
        ELSE
           CALL errore('index_minusg','Mesh too small?',ng)
        ENDIF
     ENDDO
  END IF
  RETURN
END SUBROUTINE index_minusg

