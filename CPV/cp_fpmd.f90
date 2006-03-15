!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine ggenb (b1b, b2b, b3b, nr1b ,nr2b, nr3b, nr1bx ,nr2bx, nr3bx, gcutb )
!-----------------------------------------------------------------------
   !
   ! As ggen, for the box grid. A "b" is appended to box variables.
   ! The documentation for ggen applies
   !
   use gvecb, only: ngb, ngbt, ngbl, ngbx, gb, gxb, glb, npb, nmb
   use gvecb, only: iglb, mill_b
   use io_global, only: stdout, ionode
   use mp_global, only: nproc
   use control_flags, only: iprsta
!
   implicit none
!
   integer nr1b, nr2b, nr3b, nr1bx, nr2bx, nr3bx
   real(8) b1b(3), b2b(3), b3b(3), gcutb
!
   integer, allocatable:: index(:)
   integer n1pb, n2pb, n3pb, n1mb, n2mb, n3mb
   integer it, icurr, nr1m1, nr2m1, nr3m1, ir, ig, i,j,k, itv(3), idum, ip
   real(8) t(3), g2
!
      nr1m1=nr1b-1
      nr2m1=nr2b-1
      nr3m1=nr3b-1
      ngb=0
!
!     first step : count the number of vectors with g2 < gcutb
!
!     exclude space with x<0
!
      do i= 0,nr1m1
         do j=-nr2m1,nr2m1
!
!     exclude plane with x=0, y<0
!
            if(i.eq.0.and.j.lt.0) go to 10
!
            do k=-nr3m1,nr3m1
!
!     exclude line with x=0, y=0, z<0
!
               if(i.eq.0.and.j.eq.0.and.k.lt.0) go to 20
               g2=0.d0
               do ir=1,3
                  t(ir) = DBLE(i)*b1b(ir) + DBLE(j)*b2b(ir) + DBLE(k)*b3b(ir)
                  g2=g2+t(ir)*t(ir)
               end do
               if(g2.gt.gcutb) go to 20
               ngb=ngb+1
 20            continue
            end do
 10         continue
         end do
      end do
!
!     second step: allocate space
!
      allocate(gxb(3,ngb))
      allocate(gb(ngb))
      allocate(npb(ngb))
      allocate(nmb(ngb))
      allocate(iglb(ngb))
      allocate(mill_b(3,ngb))
      allocate(index(ngb))
!
!     third step : find the vectors with g2 < gcutb
!
      ngb=0
!
!     exclude space with x<0
!
      do i= 0,nr1m1
         do j=-nr2m1,nr2m1
!
!     exclude plane with x=0, y<0
!
            if(i.eq.0.and.j.lt.0) go to 15
!
            do k=-nr3m1,nr3m1
!
!     exclude line with x=0, y=0, z<0
!
               if(i.eq.0.and.j.eq.0.and.k.lt.0) go to 25
               g2=0.d0
               do ir=1,3
                  t(ir) = DBLE(i)*b1b(ir) + DBLE(j)*b2b(ir) + DBLE(k)*b3b(ir)
                  g2=g2+t(ir)*t(ir)
               end do
               if(g2.gt.gcutb) go to 25
               ngb=ngb+1
               gb(ngb)=g2
               mill_b(1,ngb)=i
               mill_b(2,ngb)=j
               mill_b(3,ngb)=k
 25            continue
            end do
 15         continue
         end do
      end do

      IF( iprsta > 3 ) THEN
        WRITE( stdout,*)
        WRITE( stdout,170) ngb
 170    format(' ggenb: # of gb vectors < gcutb ngb = ',i6)
      END IF

      call kb07ad_cp90 (gb,ngb,index)

      do ig=1,ngb-1
         icurr=ig
 30      if(index(icurr).ne.ig) then
            itv=mill_b(:,icurr)
            mill_b(:,icurr)=mill_b(:,index(icurr))
            mill_b(:,index(icurr))=itv

            it=icurr
            icurr=index(icurr)
            index(it)=it
            if(index(icurr).eq.ig) then
               index(icurr)=icurr
               goto 35
            endif
            goto 30
         endif
 35      continue
      end do
!
      deallocate(index)
!
! costruct fft indexes (n1b,n2b,n3b) for the box grid
!
      do ig=1,ngb
         i=mill_b(1,ig)
         j=mill_b(2,ig)
         k=mill_b(3,ig)
         n1pb=i+1
         n2pb=j+1
         n3pb=k+1
!
! n1pb,n2pb,n3pb: indexes of G
! negative indexes are refolded (note that by construction i.ge.0)
!
         if(i.lt.0) n1pb=n1pb+nr1b
         if(j.lt.0) n2pb=n2pb+nr2b
         if(k.lt.0) n3pb=n3pb+nr3b
!
! n1mb,n2mb,n3mb: indexes of -G
!
         if(i.eq.0) then
            n1mb=1
         else
            n1mb=nr1b-n1pb+2
         end if
         if(j.eq.0) then
            n2mb=1
         else
            n2mb=nr2b-n2pb+2
         end if
         if(k.eq.0) then
            n3mb=1
         else
            n3mb=nr3b-n3pb+2
         end if
!
! conversion from (i,j,k) index to combined 1-d ijk index:
! ijk = 1 + (i-1)+(j-1)*ix+(k-1)*ix*jx
! where the (i,j,k) array is assumed to be dimensioned (ix,jx,kx)
!
         npb(ig) = n1pb+(n2pb-1)*nr1bx+(n3pb-1)*nr1bx*nr2bx
         nmb(ig) = n1mb+(n2mb-1)*nr1bx+(n3mb-1)*nr1bx*nr2bx
      end do
!
! shells of G - first calculate their number and position
!

      CALL gshcount( ngbl, idum, idum, iglb, ngb, gb, -1.0d0, -1.0d0 )

      IF( iprsta > 3 ) THEN
        WRITE( stdout,180) ngbl
 180    format(' ggenb: # of gb shells  < gcutb ngbl= ',i6)
      END IF
!
! then allocate the array glb
!
      allocate(glb(ngbl))
!
! and finally fill glb with the values of the shells
!
      glb(iglb(1))=gb(1)
      do ig=2,ngb
         if(iglb(ig).ne.iglb(ig-1)) glb(iglb(ig))=gb(ig)
      end do
!
! calculation of G-vectors
!
      do ig=1,ngb
         i=mill_b(1,ig)
         j=mill_b(2,ig)
         k=mill_b(3,ig)
         gxb(1,ig)=i*b1b(1)+j*b2b(1)+k*b3b(1)
         gxb(2,ig)=i*b1b(2)+j*b2b(2)+k*b3b(2)
         gxb(3,ig)=i*b1b(3)+j*b2b(3)+k*b3b(3)
      end do
!
      return
end subroutine ggenb



!-----------------------------------------------------------------------
      subroutine gcalb( alatb, b1b_ , b2b_ , b3b_  )
!-----------------------------------------------------------------------
!
      use control_flags, only: iprint
      use gvecb
!
      implicit none
      real(8), intent(in) :: alatb, b1b_ (3), b2b_ (3), b3b_ (3)
      real(8) :: b1b(3), b2b(3), b3b(3)
!
      integer i, i1,i2,i3,ig

      b1b = b1b_ * alatb
      b2b = b2b_ * alatb
      b3b = b3b_ * alatb
!
!     calculation of gxb(3,ngbx)
!
      do ig=1,ngb
         i1=mill_b(1,ig)
         i2=mill_b(2,ig)
         i3=mill_b(3,ig)
         gxb(1,ig)=i1*b1b(1)+i2*b2b(1)+i3*b3b(1)
         gxb(2,ig)=i1*b1b(2)+i2*b2b(2)+i3*b3b(2)
         gxb(3,ig)=i1*b1b(3)+i2*b2b(3)+i3*b3b(3)
         gb(ig)=gxb(1,ig)**2 + gxb(2,ig)**2 + gxb(3,ig)**2
      enddo
!
      return
      end subroutine gcalb


!-------------------------------------------------------------------------
      subroutine ggencp ( b1, b2, b3, nr1, nr2, nr3, nr1s, nr2s, nr3s,               &
     &      gcut, gcuts, gcutw, lgam )
!-----------------------------------------------------------------------
!   generates the reciprocal lattice vectors (g>) with length squared
!   less than gcut and returns them in order of increasing length.
!      g=i*b1+j*b2+k*b3,
!   where b1,b2,b3 are the vectors defining the reciprocal lattice
!
!   Only half of the g vectors (g>) are stored:
!   if g is present, -g is not (with the exception of g=0)
!   The set g> is defined by
!          g> = line(i=j=0,k>0)+plane(i=0,j>0)+space(i>0)
!
!   n1p,n2p, and n3p are the fast-fourier transform indexes of g> :
!      n1p=i+1      if i.ge.0
!      n1p=i+1+nr1  if i.lt.0
!   and the similar definitions for n2p and n3p.
!
!   n1m,n2m, and n3m are the fft indexes for g<, that is, the set
!   of vectors g=-i*b1-j*b2-k*b3 . These can be shown to be:
!      n1m=1          if i.eq.0 (or n1p.eq.1)
!      n1m=nr1-n1p+2  if i.ne.0
!   and the similar definitions for n2m and n3m.
!
!   the indexes (n1p,n2p,n3p) are collapsed into a one-dimensional
!   index np, and the same applies to negative vectors indexes
!
!   The fft indices are displaced by one unit so that g=0 corresponds
!   to element (1,1,1) (and not (0,0,0)) for purely historical reasons.
!   Negative coefficients are refolded to positive coefficients,
!   introducing a factor  exp(m*2pi*i)=1 in the fourier transform.
!
!   For a transform length n and for a single axis, if n odd:
!   -n-1                       n-1 n+1
!   ----, ..., -1, 0, 1, ...., ---,---,....,n-1  is the "true" index i,
!     2         |  |  |         2   2
!     |         |  |  |         |
!     |         |  V  V         V
!     |         |              n+1 n+3
!     |         |  1, 2, ...., ---,---,....,n    is the fft index n1 of G
!     |         |               2   2
!     | folding  \_________________ | ______|
!     |_____________________________|
!
! so: if (n1.le.(n+1)/2) i=n1-1 , otherwise, i=n1-n-1
!
! If n is even:
!     n                         n  n
!    -- , ..., -1, 0, 1, ....,  - ,-+1,....,n-1  is the "real" index i,
!     2         |  |  |         2  2
!     |         |  |  |         |
!     |         |  V  V         V
!     |         |              n   n
!     |         |  1, 2, ...., -+1,-+2,....,n    is the fft index n1 of G
!     |         |              2   2
!     | folding  \_____________ | __________|
!     |_________________________|
!
! so: if (n1.le.n/2+1) i=n1-1 ; if(n1.gt.n/2+1) i=n1-n-1 ;
!     if (n1.eq.n/2+1) i=n1-1 or i=n1-n-1, depending on how
!     the G vectors are refolded
!
!   The indices mill_l and mill_g  are the i,j,k values.
!   They are used to quickly calculate the structure factors
!      eigt=exp(-i*g*tau)        (i=imaginary unit!)
!   by decomposing eigt into products of exponentials:
!      eigt=ei1(i)*ei2(j)*ei3(k) (i=index, see above!).
!
!   ng is the total number of vectors with length squared less than gcut.
!
!   The smooth grid of g with length squared less than gcuts
!   (gcuts.le.gcut) is calculated in this routine.
!   Smooth grid variables have an "s" appended.
!
!   ngw is the total number of vectors with length squared less than gcutw
!   (gcutw.le.gcut).
!
!   the g's are in units of 2pi/a.
!
      use reciprocal_vectors, only: g, gx, igl, mill_g, g2_g, gl
      use reciprocal_vectors, only: mill_l, ig_l2g
      use reciprocal_vectors, only: gzero, gstart, sortedig_l2g
      use recvecs_indexes,    only: nm, np
      use gvecs,              only: ngs, nms, ngsl, nps
      use gvecw,              only: ngw, ngwl, ngwt, ggp
      use gvecp,              only: ng => ngm, ngl => ngml, ng_g => ngmt
      use io_global,          only: stdout
      USE fft_base,           ONLY: dfftp, dffts, fft_dlay_descriptor
      use mp,                 ONLY: mp_sum, mp_max
      use io_global,          only: ionode
      use constants,          only: eps8
      use control_flags,      only: iprsta
      !
      implicit none
      !
      real(8) :: b1(3), b2(3), b3(3), gcut, gcuts, gcutw
      real(8) :: t(3), g2
      logical      :: lgam
      integer      :: nr1,nr2,nr3, nr1s,nr2s,nr3s
      integer      :: n1p, n2p, n3p, n1m, n2m, n3m
      integer      :: n1ps, n2ps, n3ps, n1ms, n2ms, n3ms
      integer      :: it, icurr, nr1m1, nr2m1, nr3m1, nrefold, ir, ig, i,j,k
      integer      :: ichk
      integer      :: mill(3)
      !
      !  First of all count the number of G vectors according with the FFT mesh 
      ! 
      CALL gcount( ng, ngs, ngw, b1, b2, b3, nr1, nr2, nr3, gcut, gcuts, gcutw,&
                   dfftp%isind, dfftp%nr1x, lgam )
      !
      !     Second step. Compute and sort all G vectors, and build non
      !     distributed reciprocal space vectors arrays (ng_g = global
      !     number og Gs )
      !
      ng_g = ng
      ngwt = ngw

      CALL mp_sum( ng_g )
      CALL mp_sum( ngwt )

      !
      !     Temporary global and replicated arrays, used for sorting
      !
      allocate( g2_g( ng_g ) )
      allocate( mill_g( 3, ng_g ) )

      CALL gglobal( ng_g, g2_g, mill_g, b1, b2, b3, nr1, nr2, nr3, gcut, lgam )

      !
      !     third step: allocate space
      !     ng is the number of Gs local to this processor
      !
      allocate( gx ( 3, ng ) )
      allocate( g  ( ng ) )
      allocate( ggp( ngw ) )
      allocate( np ( ng ) )
      allocate( nm ( ng ) )
      allocate( igl( ng ) )

      allocate( ig_l2g( ng ) )
      allocate( mill_l( 3, ng ) )
      allocate( sortedig_l2g( ng ) )

      !
      !     fourth step : find the vectors with g2 < gcut
      !     local to each processor
      !
      CALL glocal( ng, g, ig_l2g, mill_l, ng_g, g2_g, mill_g, nr1, nr2, nr3, dfftp%isind, dfftp%nr1x  )

      IF( iprsta > 3 ) THEN
        WRITE( stdout,*)
        WRITE( stdout,150) ng
 150    format(' ggen:  # of g vectors < gcut   ng= ',i6)
        WRITE( stdout,160) ngs
 160    format(' ggen:  # of g vectors < gcuts ngs= ',i6)
        WRITE( stdout,170) ngw
 170    format(' ggen:  # of g vectors < gcutw ngw= ',i6)
      END IF

      !
      !     check for the presence of refolded G-vectors (dense grid)
      !

      CALL gchkrefold( ng, mill_l, nr1, nr2, nr3 )

      !
      !     costruct fft indexes (n1,n2,n3) for the dense grid
      !
      CALL gfftindex( np, nm, ng, mill_l, nr1, nr2, nr3, &
                      dfftp%isind, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x )

! ... Uncomment to make tests and comparisons with other codes
!      IF ( ionode ) THEN
!        DO ig=1,ng
!          WRITE( 201, fmt="( 3I6 )" ) ig, &
!             ( np( ig ) - 1 ) / dfftp%nr3x + 1, &
!             MOD( ( np( ig ) - 1 ), dfftp%nr3x ) + 1
!        END DO
!        CLOSE( 201 )
!      END IF


      !
      ! check for the presence of refolded G-vectors (smooth  grid)
      !
      CALL gchkrefold( ngs, mill_l, nr1s, nr2s, nr3s )

      !
      ! costruct fft indexes (n1s,n2s,n3s) for the smooth grid
      !
      allocate(nps(ngs))
      allocate(nms(ngs))
!
      CALL gfftindex( nps, nms, ngs, mill_l, nr1s, nr2s, nr3s, &
                      dffts%isind, dffts%nr1x, dffts%nr2x, dffts%nr3x )


! ... Uncomment to make tests and comparisons with other codes
!      IF ( ionode ) THEN
!        DO ig=1,ngs
!          WRITE( 202, fmt="( I6, 2I6, 3I4 )" ) ig, nps(ig), nms(ig), mill_l(1,ig), mill_l(2,ig), mill_l(3,ig)
!        END DO
!        CLOSE( 202 )
!      END IF

      !  ... here igl is used as temporary storage area
      !  ... sortedig_l2g is used to find out local G index given the global G index
      !
      DO ig = 1, ng
        sortedig_l2g( ig ) = ig
      END DO
      DO ig = 1, ng
        igl( ig ) = ig_l2g( ig )
      END DO
      CALL  ihpsort( ng, igl, sortedig_l2g )
      igl = 0

!
! shells of G - first calculate their number and position
!

      CALL gshcount( ngl, ngsl, ngwl, igl, ng, g, gcuts, gcutw )

!
! then allocate the array gl
!
      allocate(gl(ngl))
!
! and finally fill gl with the values of the shells
!
      gl(igl(1))=g(1)
      do ig=2,ng
         if(igl(ig).ne.igl(ig-1)) gl(igl(ig))=g(ig)
      end do
!
! gstart is the index of the first nonzero G-vector
! needed in the parallel case (G=0 is found on one node only!)
!
      if ( g(1) < 1.e-6 ) then
         gstart = 2
         gzero  = .TRUE.
      else
         gstart = 1
         gzero  = .FALSE.
      end if

      ichk = gstart
      CALL mp_max( ichk )
      IF( ichk /= 2 ) &
        CALL errore( ' ggencp ', ' inconsistent value for gstart ', ichk )
!
      IF( iprsta > 3 ) THEN
        WRITE( stdout,180) ngl
 180    format(' ggen:  # of g shells  < gcut  ngl= ',i6)
        WRITE( stdout,*)
      END IF
!
! calculation of G-vectors
!
      do ig=1,ng
         i=mill_l(1,ig)
         j=mill_l(2,ig)
         k=mill_l(3,ig)
         gx(1,ig)=i*b1(1)+j*b2(1)+k*b3(1)
         gx(2,ig)=i*b1(2)+j*b2(2)+k*b3(2)
         gx(3,ig)=i*b1(3)+j*b2(3)+k*b3(3)
      end do

      return
      end subroutine ggencp



!-------------------------------------------------------------------------
SUBROUTINE gcount( ng, ngs, ngw, b1, b2, b3, nr1, nr2, nr3, gcut, gcuts, gcutw, isind, ldis, lgam )
!-------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER ng, ngs, ngw
  integer nr1,nr2,nr3
  real(8) b1(3), b2(3), b3(3), gcut, gcuts, gcutw
  INTEGER :: isind(*), ldis
  LOGICAL :: lgam

  INTEGER :: nr1m1, nr2m1, nr3m1
  INTEGER :: i, j, k, n1p, n2p, ir

  real(8) :: g2, t(3)

  if( gcut < gcuts ) call errore(' gcount ', ' gcut .lt. gcuts ', 1 )

  ng  = 0
  ngs = 0
  ngw = 0

!
! NOTA BENE: these limits are larger than those actually needed
! (-nr/2,..,+nr/2  for nr even; -(nr-1)/2,..,+(nr-1)/2  for nr odd).
! This allows to use a slightly undersized fft grid, with some degree
! of G-vector refolding, at your own risk
!
      nr1m1=nr1-1
      nr2m1=nr2-1
      nr3m1=nr3-1

!
!     first step : count the number of vectors with g2 < gcut
!
!     exclude space with x<0
!
      loop_x: do i= -nr1m1, nr1m1
         if( lgam .AND. ( i < 0 ) ) cycle loop_x
         loop_y: do j=-nr2m1,nr2m1
!
!     exclude plane with x=0, y<0
!
            if( lgam .AND. ( i.eq.0.and.j.lt.0 ) ) cycle loop_y
!
            loop_z: do k=-nr3m1,nr3m1
!
!     exclude line with x=0, y=0, z<0
!
               if( lgam .AND. ( i.eq.0.and.j.eq.0.and.k.lt.0 ) ) cycle loop_z
!
!     consider only columns that belong to this node
!

#if defined __PARA
               n1p = i + 1
               if (n1p.lt.1) n1p = n1p + nr1
               n2p = j + 1
               if (n2p.lt.1) n2p = n2p + nr2
               if ( isind( n1p + (n2p-1)*ldis ) .eq. 0 ) cycle loop_z
#endif

               g2=0.d0
               do ir=1,3
                  t(ir) = DBLE(i)*b1(ir) + DBLE(j)*b2(ir) + DBLE(k)*b3(ir)
                  g2=g2+t(ir)*t(ir)
               end do
               if(g2.gt.gcut) cycle loop_z
               ng=ng+1
               if(g2.lt.gcutw) ngw=ngw+1
               if(g2.lt.gcuts) ngs=ngs+1
            end do loop_z
         end do loop_y
      end do loop_x

  RETURN
END SUBROUTINE gcount



!-------------------------------------------------------------------------
SUBROUTINE gglobal( ng_g, g2_g, mill_g, b1, b2, b3, nr1, nr2, nr3, gcut, lgam )
!-------------------------------------------------------------------------

  use io_global, only: ionode

  IMPLICIT NONE

  INTEGER :: ng_g
  INTEGER :: mill_g(3,*)
  real(8) :: g2_g(*)
  integer :: nr1, nr2, nr3
  real(8) :: b1(3), b2(3), b3(3), gcut
  LOGICAL :: lgam

  INTEGER :: nr1m1, nr2m1, nr3m1
  INTEGER :: i, j, k, ir, ng, ig

  real(8) :: g2, t(3)


      nr1m1=nr1-1
      nr2m1=nr2-1
      nr3m1=nr3-1

      ng = 0
!
!     exclude space with x<0
!
      loopx: do i= -nr1m1,nr1m1
         if( lgam .AND. ( i < 0 ) ) cycle loopx
         loopy: do j=-nr2m1,nr2m1
! ...       exclude plane with x=0, y<0
            if( lgam .AND. ( i.eq.0.and.j.lt.0) ) cycle loopy
            loopz: do k=-nr3m1,nr3m1
! ...          exclude line with x=0, y=0, z<0
               if( lgam .AND. (i.eq.0.and.j.eq.0.and.k.lt.0)) cycle loopz
               g2=0.d0
               do ir=1,3
                  t(ir) = DBLE(i)*b1(ir)+DBLE(j)*b2(ir)+DBLE(k)*b3(ir)
                  g2=g2+t(ir)*t(ir)
               end do
               if(g2 <= gcut) then
                 ng=ng+1
                 if( ng > ng_g ) call errore( ' gglobal ', ' too many G vectors ', ng )
                 g2_g(ng)=g2
                 mill_g(1,ng)=i
                 mill_g(2,ng)=j
                 mill_g(3,ng)=k
               end if
            end do loopz
         end do loopy
      end do loopx

      if( ng /= ng_g ) call errore( ' gglobal ', ' inconsistent number of G vectors ', ng )

      CALL sort_gvec( ng, g2_g, mill_g )

! ... Uncomment to make tests and comparisons with other codes
!      IF ( ionode ) THEN
!        DO ig=1,ng_g
!          WRITE( 201, fmt="( I6, 3I4, 1D25.16 )" ) &
!            ig, mill_g(1,ig), mill_g(2,ig), mill_g(3,ig), g2_g( ig )
!        END DO
!        CLOSE( 201 )
!      END IF

  RETURN
END SUBROUTINE gglobal



!-------------------------------------------------------------------------
SUBROUTINE glocal( ng, g, ig_l2g, mill_l, ng_g, g2_g, mill_g, nr1, nr2, nr3, isind, ldis )
!-------------------------------------------------------------------------

  use io_global, only: ionode

  IMPLICIT NONE

  INTEGER :: ng_g, ng
  INTEGER :: mill_g(3,*), ig_l2g(*), mill_l(3,*)
  real(8) :: g2_g(*), g(*)
  integer :: nr1, nr2, nr3, isind(*), ldis

  INTEGER :: i, j, k, ig, n1p, n2p, ng_l
  INTEGER :: icurr, it
  INTEGER :: mill(3)
  integer, allocatable:: index(:)

      ng_l=0
      loop_allg: do ig = 1, ng_g
        i = mill_g(1,ig)
        j = mill_g(2,ig)
        k = mill_g(3,ig)

#if defined __PARA
        n1p = i + 1
        if (n1p.lt.1) n1p = n1p + nr1
        n2p = j + 1
        if (n2p.lt.1) n2p = n2p + nr2
        if (isind(n1p+(n2p-1)*ldis).eq.0) cycle loop_allg
#endif

        ng_l=ng_l+1
        g(ng_l)=g2_g(ig)
        ig_l2g(ng_l) = ig
        mill_l(1:3,ng_l) = mill_g(1:3,ig)
      end do loop_allg

      if( ng /= ng_l ) call errore( ' glocal ', ' inconsistent number of G vectors ', ng_l )

      allocate(index(ng))
!
!     reorder the local g's in order of increasing magnitude.
!
      call kb07ad_cp90(g,ng,index)
!
      do ig=1,ng-1
         icurr=ig
 30      if(index(icurr).ne.ig) then

            it=ig_l2g(icurr)
            ig_l2g(icurr)=ig_l2g(index(icurr))
            ig_l2g(index(icurr))=it

            mill=mill_l(:,icurr)
            mill_l(:,icurr)=mill_l(:,index(icurr))
            mill_l(:,index(icurr))=mill
!
            it=icurr
            icurr=index(icurr)
            index(it)=it
            if(index(icurr).eq.ig) then
               index(icurr)=icurr
               goto 35
            endif
            goto 30
         endif
 35      continue
      end do

! ... Uncomment to make tests and comparisons with other codes
!      IF ( ionode ) THEN
!        DO ig=1,ng
!          WRITE( 201, fmt="( I6, 3I4 )" ) &
!            ig, mill_l(1,ig), mill_l(2,ig), mill_l(3,ig)
!        END DO
!        CLOSE( 201 )
!      END IF


      deallocate( index )

  RETURN
END SUBROUTINE glocal



!-------------------------------------------------------------------------
SUBROUTINE gchkrefold( ng, mill_l, nr1, nr2, nr3 )
!-------------------------------------------------------------------------

  use io_global, only: stdout

  IMPLICIT NONE

  INTEGER :: ng
  INTEGER :: mill_l(3,*)
  integer :: nr1, nr2, nr3

  INTEGER :: nr1m1, nr2m1, nr3m1
  INTEGER :: nrefold, ig

      nrefold=0
      if (mod(nr1,2).eq.0) then
         nr1m1=nr1/2-1
      else
         nr1m1=(nr1-1)/2
      end if
      if (mod(nr2,2).eq.0) then
         nr2m1=nr2/2-1
      else
         nr2m1=(nr2-1)/2
      end if
      if (mod(nr3,2).eq.0) then
         nr3m1=nr3/2-1
      else
         nr3m1=(nr3-1)/2
      end if
      do ig=1,ng
         if ( mill_l(1,ig).lt.-nr1m1.or.mill_l(1,ig).gt.nr1m1 .or.              &
     &        mill_l(2,ig).lt.-nr2m1.or.mill_l(2,ig).gt.nr2m1 .or.              &
     &        mill_l(3,ig).lt.-nr3m1.or.mill_l(3,ig).gt.nr3m1      )            &
     &        nrefold=nrefold+1
      end do
      if (nrefold.ne.0) WRITE( stdout, '('' WARNING: '',i6,                   &
     &     '' G-vectors refolded into FFT grid (ng,nrefold)'')') ng, nrefold

  RETURN
END SUBROUTINE gchkrefold


!-------------------------------------------------------------------------

SUBROUTINE gfftindex( np, nm, ng, mill_l, nr1, nr2, nr3, isind, nr1x, nr2x, nr3x )
  !
  IMPLICIT NONE

  INTEGER :: ng
  INTEGER :: isind(*), nr1x, nr2x, nr3x
  INTEGER :: mill_l(3,*), np(*), nm(*)
  integer :: nr1, nr2, nr3

  INTEGER :: n1p, n2p, n3p
  INTEGER :: n1m, n2m, n3m
  INTEGER :: i, j, k, ig, isp, ism


      do ig = 1, ng

         i = mill_l(1,ig)
         j = mill_l(2,ig)
         k = mill_l(3,ig)

         !
         ! n1p,n2p,n3p: indexes of G
         ! negative indexes are refolded (note that by construction i.ge.0)
         !
         n1p=i+1
         n2p=j+1
         n3p=k+1
         if(i.lt.0) n1p=n1p+nr1
         if(j.lt.0) n2p=n2p+nr2
         if(k.lt.0) n3p=n3p+nr3

         !
         ! n1m,n2m,n3m: indexes of -G
         !
         if(i.eq.0) then
            n1m=1
         else
            n1m=nr1-n1p+2
         end if
         if(j.eq.0) then
            n2m=1
         else
            n2m=nr2-n2p+2
         end if
         if(k.eq.0) then
            n3m=1
         else
            n3m=nr3-n3p+2
         end if

         !
         ! conversion from (i,j,k) index to combined 1-d ijk index:
         ! ijk = 1 + (i-1)+(j-1)*ix+(k-1)*ix*jx
         ! where the (i,j,k) array is assumed to be dimensioned (ix,jx,kx)
         !
         ! for the parallel case: columns along z are stored contiguously
         !

#if defined __PARA && !defined __USE_3D_FFT

         isp = isind( n1p + ( n2p - 1 ) * nr1x )
         IF( isp <= 0 ) &
           CALL errore( ' gfftindex ', ' wrong index: isp', 1 )
         IF( n3p > nr3x ) &
           CALL errore( ' gfftindex ', ' wrong index: n3p ', 1 )

         ism = isind( n1m + ( n2m - 1 ) * nr1x )
         IF( ism <= 0 ) &
           CALL errore( ' gfftindex ', ' wrong index: ism ', 1 )
         IF( n3m > nr3x ) &
           CALL errore( ' gfftindex ', ' wrong index: n3m ', 1 )

         np(ig) = n3p + ( isp - 1 ) * nr3x
         nm(ig) = n3m + ( ism - 1 ) * nr3x

#else

         np(ig) = n1p + (n2p-1)*nr1x + (n3p-1)*nr1x*nr2x
         nm(ig) = n1m + (n2m-1)*nr1x + (n3m-1)*nr1x*nr2x

#endif

      end do

  RETURN
END SUBROUTINE gfftindex


!-------------------------------------------------------------------------
SUBROUTINE gshcount( ngl, ngsl, ngwl, igl, ng, g, gcuts, gcutw )
!-------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: ngl, ngsl, ngwl
  INTEGER :: igl(*)
  INTEGER :: ng
  REAL(8) :: g(*), gcuts, gcutw

  INTEGER :: ig

      ngl=1
      igl(1)=ngl
      do ig=2,ng
         if(abs(g(ig)-g(ig-1)).gt.1.e-6)then
            ngl=ngl+1
            if (g(ig).lt.gcuts) ngsl=ngl
            if (g(ig).lt.gcutw) ngwl=ngl
         endif
         igl(ig)=ngl
      end do

  RETURN
END SUBROUTINE gshcount


!-------------------------------------------------------------------------
      subroutine gcal( alat, b1_ , b2_ , b3_ , gmax )
!-----------------------------------------------------------------------
!   calculates the values of g-vectors to be assigned to the lattice
!   points generated in subroutine ggen. these values are derived
!   from the actual values of lattice parameters, with fixed number
!   of plane waves and a cut-off function to keep energy cut-off fixed.
!
!      g=i*b1+j*b2+k*b3,
!
!   where b1,b2,b3 are the vectors defining the reciprocal lattice,
!   i go from 1 to +(nr-1) and j,k go from -(nr-1) to +(nr-1).
!
!   the g's are in units of 2pi/a.
!
      use constants, only: tpi
      use control_flags, only: iprint
      use reciprocal_vectors, only: g, gx, mill_l
      use gvecp, only: ngm
      use gvecw, only: ngw
      use gvecw, only: ggp, ecutz, ecsig, ecfix
      implicit none
!
      real(8) :: alat, b1_(3),b2_(3),b3_(3), gmax
      real(8), external :: erf
      real(8) :: b1(3),b2(3),b3(3), tpiba2, gcutz
!
      integer i1,i2,i3,ig

      b1 = b1_ * alat
      b2 = b2_ * alat
      b3 = b3_ * alat
!
!     calculation of gx(3,ng)
!
      gmax=0.
      do ig=1,ngm
         i1=mill_l(1,ig)
         i2=mill_l(2,ig)
         i3=mill_l(3,ig)
         gx(1,ig)=i1*b1(1)+i2*b2(1)+i3*b3(1)
         gx(2,ig)=i1*b1(2)+i2*b2(2)+i3*b3(2)
         gx(3,ig)=i1*b1(3)+i2*b2(3)+i3*b3(3)
         g(ig)=gx(1,ig)**2 + gx(2,ig)**2 + gx(3,ig)**2
         if(g(ig).gt.gmax) gmax=g(ig)
      enddo
 
      tpiba2 = ( tpi / alat ) ** 2
      gcutz  = ecutz / tpiba2
!
      IF( gcutz > 0.0d0 ) THEN
        do ig=1,ngw
           ggp(ig) = g(ig) + gcutz * ( 1.0d0 + erf( ( tpiba2 * g(ig) - ecfix ) / ecsig ) )
        enddo
      ELSE
        ggp( 1 : ngw ) = g( 1 : ngw )
      END IF
!
      return
      end subroutine gcal


!=----------------------------------------------------------------------------=!

        SUBROUTINE newgb( a1, a2, a3, omega, alat )
!
!     re-generation of little box g-vectors
!
          USE kinds, ONLY: DP
          USE grid_dimensions, only: nr1, nr2, nr3
          USE smallbox_grid_dimensions, only: nr1b, nr2b, nr3b
          USE small_box, only: a1b, a2b, a3b, ainvb, omegab, tpibab
          USE constants, ONLY: pi

          IMPLICIT NONE
          REAL(DP) :: a1( 3 ), a2( 3 ), a3( 3 ), omega, alat

          INTEGER :: i
          REAL(DP) :: alatb, b1b(3),b2b(3),b3b(3)

          alatb  = alat / nr1*nr1b
          tpibab = 2.d0*pi / alatb
          do i=1,3
            a1b(i)=a1(i)/nr1*nr1b
            a2b(i)=a2(i)/nr2*nr2b
            a3b(i)=a3(i)/nr3*nr3b
          enddo

          omegab=omega/nr1*nr1b/nr2*nr2b/nr3*nr3b
!
          call recips( a1b, a2b, a3b, b1b, b2b, b3b )
          !
          call gcalb( alatb, b1b, b2b, b3b )
!
          do i=1,3
            ainvb(1,i)=b1b(i)
            ainvb(2,i)=b2b(i)
            ainvb(3,i)=b3b(i)
          end do

          RETURN
        END SUBROUTINE newgb

!------------------------------------------------------------------------------!
!
!
!------------------------------------------------------------------------------!

        SUBROUTINE ecutoffs_setup( ecutwfc, ecutrho, ecfixed, qcutz, q2sigma,  &
                                   refg_ )
 
          USE kinds, ONLY: DP
          USE constants, ONLY: eps8
          USE gvecw, ONLY: ecutw
          USE gvecw, ONLY: ecfix, ecutz, ecsig
          USE gvecp, ONLY: ecutp
          USE gvecs, ONLY: ecuts, dual, doublegrid
          use betax, only: mmx, refg

          IMPLICIT NONE
          REAL(DP), INTENT(IN) ::  ecutwfc, ecutrho, ecfixed, qcutz, q2sigma
          REAL(DP), INTENT(IN) ::  refg_

          ecutw = ecutwfc

          IF ( ecutrho <= 0.D0 ) THEN
             !
             dual = 4.D0
             !
          ELSE
             !
             dual = ecutrho / ecutwfc
             !
             IF ( dual <= 1.D0 ) &
                CALL errore( ' ecutoffs_setup ', ' invalid dual? ', 1 )
             !
          END IF

          ecutp = dual * ecutwfc

          doublegrid = ( dual > 4.D0 )
          !
          IF ( doublegrid ) THEN
             !
             ecuts = 4.D0 * ecutwfc
             !
          ELSE
             !
             ecuts = ecutp
             !
          END IF

          !
          ecfix = ecfixed
          ecutz = qcutz
          ecsig = q2sigma

          refg = refg_
          mmx  = NINT( 1.2d0 * ecutp / refg )

          RETURN
        END SUBROUTINE ecutoffs_setup


        SUBROUTINE gcutoffs_setup( alat, tk_inp, nk_inp, kpoints_inp )

!  (describe briefly what this routine does...)
!  ----------------------------------------------

          USE kinds, ONLY: DP
          USE gvecw, ONLY: ecutwfc => ecutw,  gcutw
          USE gvecp, ONLY: ecutrho => ecutp,  gcutp
          USE gvecs, ONLY: ecuts, gcuts
          USE gvecb, ONLY: ecutb, gcutb
          USE gvecw, ONLY: ecfix, ecutz, ecsig
          USE gvecw, ONLY: ekcut, gkcut
          USE constants, ONLY: eps8, pi

          IMPLICIT NONE

! ...     declare subroutine arguments
          REAL(DP), INTENT(IN) :: alat
          LOGICAL, INTENT(IN) :: tk_inp
          INTEGER, INTENT(IN) :: nk_inp
          REAL(DP), INTENT(IN) :: kpoints_inp(3,*)

! ...     declare other variables
          INTEGER   :: i
          REAL(DP) :: kcut, ksq
          REAL(DP) :: tpiba

!  end of declarations
!  ----------------------------------------------

! ...   Set Values for the cutoff


          IF( alat < eps8 ) THEN
            CALL errore(' cut-off setup ', ' alat too small ', 0)
          END IF

          tpiba = 2.0d0 * pi / alat 

          ! ...  Constant cutoff simulation parameters

          gcutw = ecutwfc / tpiba**2  ! wave function cut-off
          gcutp = ecutrho / tpiba**2  ! potential cut-off
          gcuts = ecuts   / tpiba**2  ! smooth mesh cut-off

          kcut = 0.0_DP
          IF ( tk_inp ) THEN
! ...       augment plane wave cutoff to include all k+G's
            DO i = 1, nk_inp
! ...         calculate modulus
              ksq = kpoints_inp( 1, i ) ** 2 + kpoints_inp( 2, i ) ** 2 + kpoints_inp( 3, i ) ** 2
              IF ( ksq > kcut ) kcut = ksq
            END DO
          END IF

          gkcut = ( sqrt( kcut ) + sqrt( gcutw ) ) ** 2

          ekcut = gkcut * tpiba ** 2

          RETURN
        END SUBROUTINE gcutoffs_setup

!  ----------------------------------------------

      SUBROUTINE cutoffs_print_info()

        !  Print out informations about different cut-offs

        USE gvecw, ONLY: ecutwfc => ecutw,  gcutw
        USE gvecp, ONLY: ecutrho => ecutp,  gcutp
        USE gvecw, ONLY: ecfix, ecutz, ecsig
        USE gvecw, ONLY: ekcut, gkcut
        USE gvecs, ONLY: ecuts, gcuts
        USE gvecb, ONLY: ecutb, gcutb
        use betax, only: mmx, refg
        USE io_global, ONLY: stdout

        WRITE( stdout, 100 ) ecutwfc, ecutrho, ecuts, sqrt(gcutw), sqrt(gcutp), sqrt(gcuts)
        IF( ecutz > 0.0d0 ) THEN
          WRITE( stdout, 150 ) ecutz, ecsig, ecfix
        END IF

        WRITE( stdout,200) refg, mmx

100     FORMAT(/,3X,'Energy Cut-offs',/ &
                ,3X,'---------------',/ &
                ,3X,'Ecutwfc = ',F6.1,' Ryd., ', 3X,'Ecutrho = ',F6.1,' Ryd., ', 3X,'Ecuts = ',F6.1,' Ryd.',/ &
                ,3X,'Gcutwfc = ',F6.1,'     , ', 3X,'Gcutrho = ',F6.1,'       ', 3X,'Gcuts = ',F6.1)
150     FORMAT(  3X,'modified kinetic energy functional, with parameters:',/,   &
                 3X,'ecutz = ',f8.4,'  ecsig = ', f7.4,'  ecfix = ',f6.2)
200     FORMAT(  3X,'NOTA BENE: refg, mmx = ', f10.6,I6 )

        RETURN
      END SUBROUTINE cutoffs_print_info

!  ----------------------------------------------

      SUBROUTINE orthogonalize_info( )
        USE control_flags, ONLY: ortho_eps, ortho_max
        USE io_global, ONLY: stdout
        IMPLICIT NONE
           WRITE(stdout, 585)
           WRITE(stdout, 511) ortho_eps, ortho_max
  511   FORMAT(   3X,'Orthog. with lagrange multipliers : eps = ',E10.2, ',  max = ',I3)
  585   FORMAT(   3X,'Eigenvalues calculated without the kinetic term contribution')
        RETURN
      END SUBROUTINE orthogonalize_info


!  ----------------------------------------------


      SUBROUTINE electrons_print_info( )

          USE electrons_base, ONLY: nbnd, nspin, nel, nelt, nupdwn, iupdwn, f
          USE io_global, ONLY: stdout

          IMPLICIT NONE
          INTEGER :: i

          IF( nspin == 1) THEN
            WRITE(stdout,6) nelt, nbnd
            WRITE(stdout,7) ( f( i ), i = 1, nbnd )
          ELSE
            WRITE(stdout,8) nelt
            WRITE(stdout,9) nel(1)
            WRITE(stdout,7) ( f( i ), i = 1, nupdwn(1))
            WRITE(stdout,10) nel(2)
            WRITE(stdout,7) ( f( i ), i = iupdwn(2), ( iupdwn(2) + nupdwn(2) - 1 ) )
          END IF
6         FORMAT(/,3X,'Electronic states',/  &
                  ,3X,'-----------------',/  &
                  ,3X,'Number of Electron = ',I5,', of States = ',I5,/ &
                  ,3X,'Occupation numbers :')
7         FORMAT(2X,10F5.2)
8         FORMAT(/,3X,'Electronic states',/  &
                  ,3X,'-----------------',/  &
                  ,3X,'Local Spin Density calculation',/ &
                  ,3X,'Number of Electron = ',I5)
9         FORMAT(  3X,'Spins up   = ', I5, ', occupations: ')
10        FORMAT(  3X,'Spins down = ', I5, ', occupations: ')

          RETURN
      END SUBROUTINE electrons_print_info


!  ----------------------------------------------


      SUBROUTINE exch_corr_print_info()

        USE funct, ONLY: get_iexch, get_icorr, get_igcx, get_igcc, write_dft_name
        USE io_global, ONLY: stdout

        IMPLICIT NONE

        CHARACTER(LEN = 60) :: exch_info
        CHARACTER(LEN = 60) :: corr_info
        CHARACTER(LEN = 60) :: exgc_info
        CHARACTER(LEN = 60) :: cogc_info

        WRITE(stdout,800)

          ! ...     iexch => Exchange functional form
          ! ...     icorr => Correlation functional form
          ! ...     igcx  => Gradient Correction to the Exchange potential
          ! ...     igcc  => Gradient Correction to the Correlation potential

          SELECT CASE ( get_iexch() )
            CASE (0)
              exch_info = 'NONE'
            CASE (1)
              exch_info = 'SLATER'
            CASE (2)
              exch_info = 'SLATER (alpha=1)'
            CASE DEFAULT
              exch_info = 'UNKNOWN'
          END SELECT
          SELECT CASE ( get_icorr() )
            CASE (0)
              corr_info = 'NONE'
            CASE (1)
              corr_info = 'PERDEW AND ZUNGER'
            CASE (2)
              corr_info = 'VOSKO, WILK AND NUSAIR'
            CASE (3)
              corr_info = 'LEE, YANG, AND PARR'
            CASE (4)
              corr_info = 'PERDEW AND WANG'
            CASE (9)
              corr_info = 'PADE APPROXIMATION'
            CASE DEFAULT
              corr_info = 'UNKNOWN'
          END SELECT
          SELECT CASE ( get_igcx() )
            CASE (0)
              exgc_info = 'NONE'
            CASE (1)
              exgc_info = 'BECKE'
            CASE (2)
              exgc_info = 'PERDEW'
            CASE (3)
              exgc_info = 'PERDEW BURKE ERNZERHOF'
            CASE (7)
              exgc_info = 'META-TPSS'
            CASE DEFAULT
              exgc_info = 'UNKNOWN'
          END SELECT
          SELECT CASE ( get_igcc() )
            CASE (0)
              cogc_info = 'NONE'
            CASE (1)
              cogc_info = 'PERDEW'
            CASE (2)
              cogc_info = 'LEE, YANG AND PARR'
            CASE (3)
              cogc_info = 'PERDEW AND WANG'
            CASE (4)
              cogc_info = 'PERDEW BURKE ERNZERHOF'
            CASE (6)
              cogc_info = 'META-TPSS'
            CASE DEFAULT
              cogc_info = 'UNKNOWN'
          END SELECT

          WRITE(stdout,910)
          WRITE(stdout,fmt='(5X,"Exchange functional: ",A)') exch_info
          WRITE(stdout,fmt='(5X,"Correlation functional: ",A)') corr_info
          IF( ( get_igcx() > 0 ) .OR. ( get_igcc() > 0 ) ) THEN
            WRITE(stdout,810)
            WRITE(stdout,fmt='(5X,"Exchange functional: ",A)') exgc_info
            WRITE(stdout,fmt='(5X,"Correlation functional: ",A)') cogc_info
          END IF

        call write_dft_name

800 FORMAT(//,3X,'Exchange and correlations functionals',/ &
             ,3X,'-------------------------------------')
810 FORMAT(   3X,'Using Generalized Gradient Corrections with')
910 FORMAT(   3X,'Using Local Density Approximation with')

        RETURN
      END SUBROUTINE exch_corr_print_info



!  ----------------------------------------------



       SUBROUTINE ions_print_info( )
            
         !  Print info about input parameter for ion dynamic

         USE io_global,     ONLY: ionode, stdout
         USE control_flags, ONLY: tranp, amprp, tnosep, tolp, tfor, tsdp, tzerop, &
                                  tv0rd, taurdr, nv0rd, nbeg, tcp, tcap, &
                                  program_name
         USE ions_base,     ONLY: tau_srt, tau_units, if_pos, ind_srt, nsp, na, &
                                  pmass, nat, fricp, greasp, rcmax
         USE ions_nose,     ONLY: tempw, ndega
         USE constants,     ONLY: scmass

         IMPLICIT NONE
              
         integer is, ia, k, ic, isa
         LOGICAL :: ismb( 3 ) 
                
         WRITE( stdout, 50 ) 

         IF( .NOT. tfor ) THEN
           WRITE( stdout, 518 )
         ELSE
           WRITE( stdout, 520 )
           IF( tsdp ) THEN
             WRITE( stdout, 521 )
           ELSE
             WRITE( stdout, 522 )
           END IF
           WRITE( stdout, 523 ) ndega
           WRITE( stdout, 524 ) fricp, greasp
           IF( tzerop ) then
             IF( tv0rd ) THEN
               WRITE( stdout, 850 ) nv0rd
             ELSE
               WRITE( stdout, 635 )
             ENDIF 
           ENDIF
         END IF 
              
         DO is = 1, nsp
           IF( tranp(is) ) THEN
             WRITE( stdout,510)
             WRITE( stdout,512) is, amprp(is)
           END IF
         END DO

         WRITE(stdout,660) 
         isa = 0
         DO IS = 1, nsp
           WRITE(stdout,1000) is, na(is), pmass(is), pmass(is) / scmass, rcmax(is)
           DO IA = 1, na(is)
             isa = isa + 1
             WRITE(stdout,1010) ( tau_srt(k,isa), K = 1,3 )
           END DO
         END DO    

         IF ( ( nbeg > -1 ) .AND. ( .NOT. taurdr ) ) THEN
            WRITE(stdout,661)
         ELSE
            WRITE(stdout,662)
         ENDIF

         IF( tfor ) THEN

            IF( ANY( ( if_pos( 1:3, 1:nat ) == 0 )  ) ) THEN

              WRITE(stdout,1020)
              WRITE(stdout,1022)

              DO isa = 1, nat
                ia = ind_srt( isa )
                ismb( 1 ) = ( if_pos(1,ia) /= 0 )
                ismb( 2 ) = ( if_pos(2,ia) /= 0 )
                ismb( 3 ) = ( if_pos(3,ia) /= 0 )
                IF( .NOT. ALL( ismb ) ) THEN
                  WRITE( stdout, 1023 ) isa, ( ismb(k), K = 1, 3 )
                END IF
              END DO

            ELSE

              WRITE(stdout,1021)

            END IF
         END IF

         IF( tfor ) THEN
           if( ( tcp .or. tcap .or. tnosep ) .and. tsdp ) then
             call errore(' ions_print_info',' t contr. for ions when tsdp=.t.',1)
           endif
           IF(.not. tcp .and. .not. tcap .and. .not. tnosep ) THEN
              WRITE( stdout,550)
           ELSE IF( tcp .and. tcap ) then
             call errore(' ions_print_info',' tcp and tcap both true',1)
           ELSE IF( tcp .and. tnosep ) then
             call errore(' ions_print_info',' tcp and tnosep both true',1)
           ELSE IF(tcap .and. tnosep ) then
             call errore(' ions_print_info',' tcap and tnosep both true',1)
           ELSE IF(tcp) THEN
             WRITE( stdout,555) tempw,tolp
           ELSE IF(tcap) THEN
             WRITE( stdout,560) tempw,tolp
           ELSE IF(tnosep) THEN
             WRITE( stdout,595)
           ELSE
             WRITE( stdout,550)
           END IF
         END IF

   50 FORMAT(//,3X,'Ions Simulation Parameters',/ &
               ,3X,'--------------------------')

  510 FORMAT(   3X,'Initial random displacement of ionic coordinates',/, & 
                3X,' specie  amplitude')
  512 FORMAT(   3X,I7,2X,F9.6)

  518 FORMAT(   3X,'Ions are not allowed to move')
  520 FORMAT(   3X,'Ions are allowed to move')
  521 FORMAT(   3X,'Ions dynamics with steepest descent')
  522 FORMAT(   3X,'Ions dynamics with newton equations')
  523 format(   3X,'the temperature is computed for ',i5,' degrees of freedom')
  524 format(   3X,'ion dynamics with fricp = ',f7.4,' and greasp = ',f7.4)
  550 FORMAT(   3X,'Ionic temperature is not controlled')
  555 FORMAT(   3X,'Ionic temperature control via ', &
                   'rescaling of velocities :',/ &
               ,3X,'temperature required = ',F10.5,'K, ', &
                   'tolerance = ',F10.5,'K')
  560 FORMAT(   3X,'Ionic temperature control via ', &
                   'canonical velocities rescaling :',/ &
               ,3X,'temperature required = ',F10.5,'K, ', &
                   'tolerance = ',F10.5,'K')
  595 FORMAT(   3X,'Ionic temperature control via nose thermostat')
  635 FORMAT(   3X,'Zero initial momentum for ions')

  660 FORMAT(   3X,'Ionic position (from input)', /, &
                3X,'sorted by specie, and converted to real a.u. coordinates')
  661 FORMAT(   3X,'Ionic position will be re-read from restart file')
  662 FORMAT(   3X,'Ionic position read from input file')

  850 FORMAT(   3X,'Initial ion velocities read from unit : ',I4)

 1000 FORMAT(3X,'Species ',I3,' atoms = ',I4,' mass = ',F12.2, ' (a.u.), ', &
               & F12.2, ' (amu)', ' rcmax = ', F6.2, ' (a.u.)' )
 1010 FORMAT(3X,3(1X,F12.6))
 1020 FORMAT(/,3X,'NOT all atoms are allowed to move ')
 1021 FORMAT(/,3X,'All atoms are allowed to move')
 1022 FORMAT(  3X,' indx  ..x.. ..y.. ..z..')
 1023 FORMAT(  3X,I4,3(1X,L5))



         RETURN
       END SUBROUTINE ions_print_info


!  ----------------------------------------------

        subroutine cell_print_info( )

          USE constants, ONLY: au_gpa
          USE control_flags, ONLY: thdyn, tsdc, tzeroc, tbeg, nbeg, tpre
          USE control_flags, ONLY: tnoseh
          USE io_global, ONLY: stdout
          USE cell_base, ONLY: press, frich, greash, wmass

          IMPLICIT NONE

          WRITE(stdout,545 )
          IF ( tpre ) WRITE( stdout, 600 )
          IF ( tbeg ) THEN
            WRITE(stdout,546)
          ELSE
            WRITE(stdout,547)
            IF( nbeg > -1 ) WRITE( stdout, 548 )
          END IF

          IF( .NOT. thdyn ) THEN
            WRITE( stdout,525)
            WRITE( stdout,606)
          ELSE
            IF( tsdc ) THEN
              WRITE( stdout,526)
            ELSE
              IF( frich /= 0.0d0 ) THEN
                WRITE( stdout,602) frich, greash
              ELSE
                WRITE( stdout,527)
              END IF
              IF( tnoseh ) then
                WRITE( stdout,604) 
              ELSE
                WRITE( stdout,565)
              END IF
              ! if( thdiag ) WRITE( stdout,608)
              IF( tzeroc ) THEN
                WRITE( stdout,563)
              ENDIF
            END IF
            WRITE( stdout,530) press * au_gpa, wmass
          END IF


 545     FORMAT(//,3X,'Cell Dynamics Parameters (from STDIN)',/ &
                  ,3X,'-------------------------------------')
 546     FORMAT(   3X,'Simulation cell read from STDIN')
 547     FORMAT(   3X,'Starting cell generated from CELLDM')
 548     FORMAT(   3X,'Cell parameters will be re-read from restart file')
 525     FORMAT(   3X,'Constant VOLUME Molecular dynamics')
 606     format(   3X,'cell parameters are not allowed to move')
 526     FORMAT(   3X,'Volume dynamics with steepest descent')
 527     FORMAT(   3X,'Volume dynamics with newton equations')
 530     FORMAT(   3X,'Constant PRESSURE Molecular dynamics:',/ &
                  ,3X,'External pressure (GPa) = ',F11.2,/ &
                  ,3X,'Volume mass             = ',F11.2)
 563     FORMAT(   3X,'Zero initial momentum for cell variables')
 565     FORMAT(   3X,'Volume dynamics: the temperature is not controlled')
 604     format(   3X,'cell parameters dynamics with nose` temp. control' )

 600  format( 3X, 'internal stress tensor calculated')
 602  format( 3X, 'cell parameters dynamics with frich = ',f7.4,            &
     &        3X, 'and greash = ',f7.4 )
 608  format( 3X, 'frozen off-diagonal cell parameters'//)

        return
      end subroutine cell_print_info


!----------------------------------------------
SUBROUTINE gmeshinfo( )
!----------------------------------------------
   !
   !   Print out the number of g vectors for the different mesh
   !
   USE mp_global, ONLY: nproc, mpime, group
   USE io_global, ONLY: ionode, ionode_id, stdout
   USE mp,        ONLY: mp_max, mp_gather
   use gvecb,     only: ngb
   USE reciprocal_vectors, only: ngst, ngs, ngsx,  &
              ngw_g  => ngwt,   &
              ngw_l  => ngw ,   &
              ngw_lx => ngwx,   &
              ng_g   => ngmt,   &
              ng_l   => ngm ,   &
              ng_lx  => ngmx

   IMPLICIT NONE

   INTEGER :: ip, ng_snd(3), ng_rcv(3,nproc)

   IF(ionode) THEN
      WRITE( stdout,*)
      WRITE( stdout,*) '  Reciprocal Space Mesh'
      WRITE( stdout,*) '  ---------------------'
   END IF

   ng_snd(1) = ng_g
   ng_snd(2) = ng_l
   ng_snd(3) = ng_lx
   CALL mp_gather(ng_snd, ng_rcv, ionode_id, group)
   !
   IF(ionode) THEN
      WRITE( stdout,1000)
      DO ip = 1, nproc
         WRITE( stdout,1010) ip, ng_rcv(1,ip), ng_rcv(2,ip), ng_rcv(3,ip)
      END DO
   END IF
   !
   ng_snd(1) = ngst
   ng_snd(2) = ngs
   ng_snd(3) = ngsx
   CALL mp_gather(ng_snd, ng_rcv, ionode_id, group)
   !
   IF(ionode) THEN
      WRITE( stdout,1001)
      DO ip = 1, nproc
         WRITE( stdout,1010) ip, ng_rcv(1,ip), ng_rcv(2,ip), ng_rcv(3,ip)
      END DO
   END IF
   !
   ng_snd(1) = ngw_g
   ng_snd(2) = ngw_l
   ng_snd(3) = ngw_lx
   CALL mp_gather(ng_snd, ng_rcv, ionode_id, group)
   !
   IF(ionode) THEN
      WRITE( stdout,1002)
      DO ip = 1, nproc
         WRITE( stdout,1010) ip, ng_rcv(1,ip), ng_rcv(2,ip), ng_rcv(3,ip)
      END DO
   END IF
   !
   IF(ionode) THEN
      WRITE( stdout,1050)
      WRITE( stdout,1060) ngb
   END IF

   1000    FORMAT(16X,'Large Mesh',/, &
           3X,'PE   Global(ngmt)     Local(ngm) MaxLocal(ngmx)') 
   1001    FORMAT(16X,'Smooth Mesh',/, &
           3X,'PE   Global(ngst)     Local(ngs) MaxLocal(ngsx)') 
   1002    FORMAT(16X,'Wave function Mesh',/, &
           3X,'PE   Global(ngwt)     Local(ngw) MaxLocal(ngwx)') 
   1010    FORMAT( I5,3I15 )
   1050    FORMAT(/,16X,'Small box Mesh')
   1060    FORMAT( 3X, 'ngb = ', I12, ' not distributed to processors' )

   RETURN

END SUBROUTINE gmeshinfo
