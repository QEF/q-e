!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine ggen
  !----------------------------------------------------------------------
  !
  !     This routine generates all the reciprocal lattice vectors
  !     contained in the sphere of radius gcutm. Furthermore it
  !     computes the indices nl which give the correspondence
  !     between the fft mesh points and the array of g vectors.
  !
#include "machine.h"
  use pwcom
#ifdef __PARA
  use para
#endif
  implicit none
  !
  real(kind=DP), parameter :: eps = 1.d-8
  !
  !     here a few local variables
  !
  real(kind=DP) ::  t (3), tt, swap, d (3), dnorm
  real(kind=DP), allocatable ::  esort (:)
  !
  integer :: ngmx, n1, n2, n3, n1s, n2s, n3s
  !
  real(kind=DP), allocatable :: g2sort_g(:)
  ! array containing all g vectors, on all processors: replicated data
  integer, allocatable :: mill_g(:,:)
  ! array containing all g vectors generators, on all processors: replicated data
  integer, allocatable :: igsrt(:)
  !
#ifdef __PARA
  integer :: m1, m2, m3, mc
  !
#endif
  integer :: i, j, k, ipol, ng, igl, iswap, indsw
  !
  ! counters
  !
  !    set the total number of fft mesh points and and initial value of gg
  !    The choice of gcutm is due to the fact that we have to order the
  !    vectors after computing them.
  !
  gg(:) = gcutm + 1.d0
  allocate(esort(ngm) )
  esort(:) = 1.0d20
  !
  !     set d vector for unique ordering
  !
  d (1) = 0.25657642786d0
  d (2) = 0.35342818974d0
  d (3) = 0.56421652427d0
  dnorm = sqrt (d (1) * d (1) + d (2) * d (2) + d (3) * d (3) )
  d (1) = d (1) / dnorm
  d (2) = d (2) / dnorm
  d (3) = d (3) / dnorm
  !
  !    and computes all the g vectors inside a sphere
  !
  allocate( igsrt( ngm_g ) )
  allocate( g2sort_g( ngm_g ) )
  g2sort_g(:) = 1.0d20
  allocate( mill_g( 3, ngm_g ) )
  allocate( ig_l2g( ngm_l ) )
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
  do i = - n1, n1
    do j = - n2, n2
      do k = - n3, n3
        tt = 0.d0
        do ipol = 1, 3
          t (ipol) = i * bg (ipol, 1) + j * bg (ipol, 2) + k * bg (ipol, 3)
          tt = tt + t (ipol) * t (ipol)
        enddo
        if (tt <= gcutm) then
          ngm = ngm + 1
          if (tt <= gcutms) ngms = ngms + 1
          if (ngm > ngm_g) call errore ('ggen', 'too many g-vectors', ngm)
          mill_g( 1, ngm ) = i
          mill_g( 2, ngm ) = j
          mill_g( 3, ngm ) = k
          if ( tt > eps ) then
             g2sort_g(ngm) = 1.d4 * tt + &
               (t (1) * d (1) + t (2) * d (2) + t (3) * d (3) ) / sqrt (tt)
          else
             g2sort_g(ngm) = 0.d0
          endif
        end if
      enddo
    enddo
  enddo

  if (ngm  /= ngm_g ) call errore ('ggen', 'g-vectors missing !', abs(ngm - ngm_g))
  if (ngms /= ngms_g) call errore ('ggen', 'smooth g-vectors missing !', abs(ngms - ngms_g))

  igsrt(1) = 0
  call hpsort(ngm_g, g2sort_g, igsrt)
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

#ifndef __OLD_GGEN_LOOP

  write(6, fmt="(//,' --- Executing new GGEN Loop ---',//)" )

  ngm = 0
  ngms = 0
  do ng = 1, ngm_g
    i = mill_g(1, ng)
    j = mill_g(2, ng)
    k = mill_g(3, ng)

#ifdef __PARA
    m1 = mod (i, nr1) + 1
    if (m1.lt.1) m1 = m1 + nr1
    m2 = mod (j, nr2) + 1
    if (m2.lt.1) m2 = m2 + nr2
    mc = m1 + (m2 - 1) * nrx1
    if (ipc (mc) .eq.0) goto 1
#endif

    tt = 0.d0
    do ipol = 1, 3
      t (ipol) = i * bg (ipol, 1) + j * bg (ipol, 2) + k * bg (ipol, 3)
      tt = tt + t (ipol) * t (ipol)
    enddo

    ngm = ngm + 1
    if (tt.le.gcutms) ngms = ngms + 1
    if (ngm > ngmx) call errore ('ggen', 'too many g-vectors', ngm)
    !
    !  Here map local and global g index !!!
    !
    ig_l2g( ngm ) = ng
    !
    g (1:3, ngm) = t (1:3)
    gg (ngm) = tt

    if (tt.gt.eps) then
      esort (ngm) = 1.d4 * tt + (t (1) * d (1) + t (2) * d (2) &
      + t (3) * d (3) ) / sqrt (tt)
    else
      esort (ngm) = 0.d0
    endif

1   continue
  enddo

! ... Uncomment to make tests and comparisons with other codes
!  DO ng=1,ngm_g
!    WRITE( 202, fmt="( I6, 3I4, 2D25.16 )" ) &
!      ng, mill_g(1,ng), mill_g(2,ng), mill_g(3,ng), gg( ng ), g2sort_g( ng )
!  END DO
!  CLOSE( 202 )

#else
  !
  !    and computes all the g vectors inside a sphere
  !
  n1 = nr1 + 1
  n2 = nr2 + 1
  n3 = nr3 + 1
  ngmx = ngm
  ngm = 0
  ngms = 0

  do i = - n1, n1
#ifdef __PARA
     m1 = mod (i, nr1) + 1
     if (m1.lt.1) m1 = m1 + nr1
     do j = - n2, n2
        m2 = mod (j, nr2) + 1
        if (m2.lt.1) m2 = m2 + nr2
        mc = m1 + (m2 - 1) * nrx1
        if (ipc (mc) .eq.0) goto 1
#else
        do j = - n2, n2
#endif
           do k = - n3, n3
              tt = 0.d0
              do ipol = 1, 3
                 t (ipol) = i * bg (ipol, 1) + j * bg (ipol, 2) + k * bg (ipol, 3)
                 tt = tt + t (ipol) * t (ipol)
              enddo
              if (tt.le.gcutm) then
                 ngm = ngm + 1
                 if (tt.le.gcutms) ngms = ngms + 1
                 if (ngm.gt.ngmx) call errore ('ggen', 'too many g-vectors', ngm)
                 do ipol = 1, 3
                    g (ipol, ngm) = t (ipol)
                 enddo
                 gg (ngm) = tt
                 if (tt.gt.eps) then
                    esort (ngm) = 1.d4 * tt + (t (1) * d (1) + t (2) * d (2) &
                         + t (3) * d (3) ) / sqrt (tt)
                 else
                    esort (ngm) = 0.d0
                 endif
              endif
           enddo
1          continue
        enddo

     enddo

#endif

     if (ngm.ne.ngmx) call errore ('ggen', 'g-vectors missing !', abs(ngm - ngmx))
     !
     !   reorder the g's in order of increasing magnitude. On exit
     !   from hpsort esort is ordered, and nl contains the new order.
     !
     !   initialize the index inside sorting routine

     nl (1) = 0
     call hpsort (ngm, esort, nl)
     !
     !   reorder also the g vectors, and nl
     !
     do ng = 1, ngm - 1
20      indsw = nl (ng)
        if (indsw.ne.ng) then
           do ipol = 1, 3
              swap = g (ipol, indsw)
              g (ipol, indsw) = g (ipol, nl (indsw) )
              g (ipol, nl (indsw) ) = swap
           enddo
           swap = gg (indsw)
           gg (indsw) = gg (nl (indsw) )
           gg (nl (indsw) ) = swap

#ifndef __OLD_GGEN_LOOP
          !
          !  Remember: ig_l2g is the index of a given G vectors in the
          !  sorted global array containing all G vectors, it is used to
          !  collect all wave function components
          !
          iswap = ig_l2g( indsw )
          ig_l2g( indsw ) = ig_l2g( nl(indsw) )
          ig_l2g( nl(indsw) ) = iswap
#endif

           iswap = nl (ng)
           nl (ng) = nl (indsw)
           nl (indsw) = iswap

           goto 20
        endif

     enddo
     !
     !     determine first nonzero g vector
     !
     if (gg(1).le.eps) then
        gstart=2
     else
        gstart=1
     end if
     !
     !     Now set nl and nls with the correct fft correspondence
     !
     do ng = 1, ngm
        n1 = nint (g (1, ng) * at (1, 1) + g (2, ng) * at (2, 1) + g (3, &
             ng) * at (3, 1) ) + 1
        ig1 (ng) = n1 - 1
        n1s = n1
        if (n1.lt.1) n1 = n1 + nr1
        if (n1s.lt.1) n1s = n1s + nr1s
        n2 = nint (g (1, ng) * at (1, 2) + g (2, ng) * at (2, 2) + g (3, &
             ng) * at (3, 2) ) + 1
        ig2 (ng) = n2 - 1
        n2s = n2
        if (n2.lt.1) n2 = n2 + nr2
        if (n2s.lt.1) n2s = n2s + nr2s
        n3 = nint (g (1, ng) * at (1, 3) + g (2, ng) * at (2, 3) + g (3, &
             ng) * at (3, 3) ) + 1
        ig3 (ng) = n3 - 1
        n3s = n3
        if (n3.lt.1) n3 = n3 + nr3
        if (n3s.lt.1) n3s = n3s + nr3s
        if (n1.le.nr1.and.n2.le.nr2.and.n3.le.nr3) then
#ifdef __PARA
           nl (ng) = n3 + (ipc (n1 + (n2 - 1) * nrx1) - 1) * nrx3
           if (ng.le.ngms) nls (ng) = n3s + (ipcs (n1s + (n2s - 1) &
                * nrx1s) - 1) * nrx3s
#else
           nl (ng) = n1 + (n2 - 1) * nrx1 + (n3 - 1) * nrx1 * nrx2
           if (ng.le.ngms) nls (ng) = n1s + (n2s - 1) * nrx1s + (n3s - 1) &
                * nrx1s * nr2s
#endif
        else
           call errore('ggen','Mesh too small?',ng)
        endif
     enddo
     !
     ! calculate number of G shells: ngl
     !
     if (lmovecell) then
        !
        ! in case of a variable cell run each G vector has its shell
        !
        ngl = ngm
        gl => gg
        do ng = 1, ngm
           igtongl (ng) = ng

        enddo

     else
        !
        ! G vectors are grouped in shells with the same norm
        !
        ngl = 1
        igtongl (1) = 1
        do ng = 2, ngm
           if (gg (ng) .gt.gg (ng - 1) + eps) then
              ngl = ngl + 1
           endif
           igtongl (ng) = ngl

        enddo

        allocate (gl( ngl))    
        gl (1) = gg (1)
        igl = 1
        do ng = 2, ngm
           if (gg (ng) .gt.gg (ng - 1) + eps) then
              igl = igl + 1
              gl (igl) = gg (ng)
           endif

        enddo

        if (igl.ne.ngl) call errore ('setup', 'igl <> ngl', ngl)

     endif


     deallocate( esort  )
     deallocate( g2sort_g )
     deallocate( mill_g )
     deallocate( igsrt )

     return
   end subroutine ggen

