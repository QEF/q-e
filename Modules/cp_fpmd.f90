! this file contains common subroutines and modules between
! CP and FPMD


module cvan
  !     ionic pseudo-potential variables
  use parameters, only: nsx
  implicit none
  save
  integer nvb, ish(nsx), ipp(nsx)
  !     nvb    = number of species with Vanderbilt PPs
  !     ipp(is)= pseudopotential type - to be removed
  !     ish(is)= used for indexing the nonlocal projectors betae
  !              with contiguous indices inl=ish(is)+(iv-1)*na(is)+1
  !              where "is" is the species and iv=1,nh(is)
end module cvan


module qrl_mod

  use parameters, only: nsx, ndmx, nbrx, lqmax
  implicit none
  save
!
! qrl       q(r) functions (old format)
! cmesh     used only for Herman-Skillman mesh (old format)
!
  real(kind=8) :: qrl(ndmx,nbrx,nbrx,lqmax,nsx)
  real(kind=8) :: cmesh(nsx)

end module qrl_mod


!-----------------------------------------------------------------------
      subroutine fill_qrl(is)
!-----------------------------------------------------------------------
!
! for compatibility with old Vanderbilt formats
!
      use uspp_param, only: qfunc, nqf, qfcoef, rinner, lll, nbeta, &
                       kkbeta
      use qrl_mod, only: qrl
      use atom, only: r
      ! the above module variables has no dependency from iosys
!
      implicit none
      integer :: is
!
      integer :: iv, jv, lmin, lmax, l, ir, i
!
      do iv=1,nbeta(is)
         do jv=iv,nbeta(is)
            lmin=lll(jv,is)-lll(iv,is)+1
            lmax=lmin+2*lll(iv,is)
            do l=lmin,lmax
               do ir=1,kkbeta(is)
                  if (r(ir,is).ge.rinner(l,is)) then
                     qrl(ir,iv,jv,l,is)=qfunc(ir,iv,jv,is)
                  else
                     qrl(ir,iv,jv,l,is)=qfcoef(1,l,iv,jv,is)
                     do i = 2, nqf(is)
                        qrl(ir,iv,jv,l,is)=qrl(ir,iv,jv,l,is) +      &
                             qfcoef(i,l,iv,jv,is)*r(ir,is)**(2*i-2)
                     end do
                     qrl(ir,iv,jv,l,is) = qrl(ir,iv,jv,l,is) * r(ir,is)**(l+1)
                  end if
               end do
            end do
         end do
      end do
    end subroutine fill_qrl


!-----------------------------------------------------------------------
      subroutine ggenb (b1b, b2b, b3b, nr1b ,nr2b, nr3b, nr1bx ,nr2bx, nr3bx, gcutb )
!-----------------------------------------------------------------------
!
! As ggen, for the box grid. A "b" is appended to box variables.
! The documentation for ggen applies
!
      use gvecb, only: ngb, ngbt, ngbl, ngbx, gb, gxb, glb, npb, nmb
      use gvecb, only: iglb, mill_b
      use io_global, only: stdout

!
      implicit none
!
      integer nr1b, nr2b, nr3b, nr1bx, nr2bx, nr3bx
      real(kind=8) b1b(3), b2b(3), b3b(3), gcutb
!
      integer, allocatable:: index(:)
      integer n1pb, n2pb, n3pb, n1mb, n2mb, n3mb
      integer it, icurr, nr1m1, nr2m1, nr3m1, ir, ig, i,j,k, itv(3), idum
!       cray:
!      integer jwork(257)
      real(kind=8) t(3), g2
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
                  t(ir) = dble(i)*b1b(ir) + dble(j)*b2b(ir) + dble(k)*b3b(ir)
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
                  t(ir) = dble(i)*b1b(ir) + dble(j)*b2b(ir) + dble(k)*b3b(ir)
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
!
      WRITE( stdout,*)
      WRITE( stdout,170) ngb
 170  format(' ggenb: # of gb vectors < gcutb ngb = ',i6)
!
!   reorder the g's in order of increasing magnitude.
!       cray:
!       call orders (2,jwork,gb,index,ngb,1,8,1)
!       generic:
      call kb07ad_cp90 (gb,ngb,index)
      do ig=1,ngb-1
         icurr=ig
 30      if(index(icurr).ne.ig) then
!     comment if not using cray orders from here
!         g2=gb(icurr)
!         gb(icurr)=gb(index(icurr))
!         gb(index(icurr))=g2
!       to here.
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

      WRITE( stdout,180) ngbl
 180  format(' ggenb: # of gb shells  < gcutb ngbl= ',i6)
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
      end subroutine

!-----------------------------------------------------------------------
      subroutine gcalb(b1b,b2b,b3b)
!-----------------------------------------------------------------------
!
      use control_flags, only: iprint
      use gvecb
!
      implicit none
      integer nr1b,nr2b,nr3b
      real(kind=8) b1b(3),b2b(3),b3b(3)
!
      integer i, i1,i2,i3,ig
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
      end


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
      use reciprocal_vectors, only: gstart, sortedig_l2g
      use recvecs_indexes, only: nm, np
      use gvecs, only: ngs, nms, ngsl, nps
      use gvecw, only: ngw, ngwl, ngwt, ggp
      use gvecp, only: ng => ngm, ngl => ngml, ng_g => ngmt
      use io_global, only: stdout
      USE fft_base, ONLY: dfftp, dffts, fft_dlay_descriptor
      use mp, ONLY: mp_sum
      use io_global, only: ionode
      use constants, only: eps8

      implicit none
      integer nr1,nr2,nr3, nr1s,nr2s,nr3s
      real(kind=8) b1(3), b2(3), b3(3), gcut, gcuts, gcutw
      logical lgam
      integer n1p, n2p, n3p, n1m, n2m, n3m
      integer n1ps, n2ps, n3ps, n1ms, n2ms, n3ms
!       cray:
!      integer jwork(257)
      integer it, icurr, nr1m1, nr2m1, nr3m1, nrefold, ir, ig, i,j,k
      integer mill(3)
      real(kind=8) t(3), g2
!

      CALL gcount( ng, ngs, ngw, b1, b2, b3, nr1, nr2, nr3,  &
     &      gcut, gcuts, gcutw, dfftp%isind, dfftp%nr1x, lgam )

      ! WRITE(6,*) 'DEBUG gcount: ng, ngs, ngw = ', ng, ngs, ngw

!
!     Second step. Compute and sort all G vectors, and build non
!     distributed reciprocal space vectors arrays (ng_g = global
!     number og Gs )
!
      ng_g = ng
      ngwt = ngw

      CALL mp_sum( ng_g )
      CALL mp_sum( ngwt )

      allocate(g2_g(ng_g))
      allocate(mill_g(3,ng_g))

      CALL gglobal( ng_g, g2_g, mill_g, b1, b2, b3, nr1, nr2, nr3, gcut, lgam )

      ! WRITE(6,*) 'DEBUG gglobal: ng_g  = ', ng_g

!
!     third step: allocate space
!     ng is the number of Gs local to this processor
!
      allocate(gx(3,ng))
      allocate(g(ng))
      allocate(ggp(ngw))
      allocate(np(ng))
      allocate(nm(ng))
      allocate(igl(ng))

      allocate(ig_l2g(ng))
      allocate(mill_l(3,ng))
      allocate(sortedig_l2g(ng))

!
!     fourth step : find the vectors with g2 < gcut
!     local to each processor
!
      CALL glocal( ng, g, ig_l2g, mill_l, ng_g, g2_g, mill_g, nr1, nr2, nr3, dfftp%isind, dfftp%nr1x  )

      ! WRITE(6,*) 'DEBUG glocal: ng, ng_g  = ', ng, ng_g
      ! call hangup
      ! stop

      WRITE( stdout,*)
      WRITE( stdout,150) ng
 150  format(' ggen:  # of g vectors < gcut   ng= ',i6)
      WRITE( stdout,160) ngs
 160  format(' ggen:  # of g vectors < gcuts ngs= ',i6)
      WRITE( stdout,170) ngw
 170  format(' ggen:  # of g vectors < gcutw ngw= ',i6)

!
! check for the presence of refolded G-vectors (dense grid)
!

      CALL gchkrefold( ng, mill_l, nr1, nr2, nr3 )
!
! costruct fft indexes (n1,n2,n3) for the dense grid
!

      CALL gfftindex( np, nm, ng, mill_l, nr1, nr2, nr3, dfftp%isind, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x )

!
! check for the presence of refolded G-vectors (smooth  grid)
!

      CALL gchkrefold( ngs, mill_l, nr1s, nr2s, nr3s )

!
! costruct fft indexes (n1s,n2s,n3s) for the small grid
!
      allocate(nps(ngs))
      allocate(nms(ngs))
!
      CALL gfftindex( nps, nms, ngs, mill_l, nr1s, nr2s, nr3s, dffts%isind, dffts%nr1x, dffts%nr2x, dffts%nr3x )

      !  ... here igl is used as temporary storage area
      !  ... sortedig_l2g is used to find out local G index given the global G index
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
      if (g(1).lt.1.e-6) then
         gstart=2
      else
         gstart=1
      end if
!
      WRITE( stdout,180) ngl
 180  format(' ggen:  # of g shells  < gcut  ngl= ',i6)
      WRITE( stdout,*)
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
      end subroutine



!-------------------------------------------------------------------------
SUBROUTINE gcount( ng, ngs, ngw, b1, b2, b3, nr1, nr2, nr3, gcut, gcuts, gcutw, isind, ldis, lgam )
!-------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER ng, ngs, ngw
  integer nr1,nr2,nr3
  real(kind=8) b1(3), b2(3), b3(3), gcut, gcuts, gcutw
  INTEGER :: isind(*), ldis
  LOGICAL :: lgam

  INTEGER :: nr1m1, nr2m1, nr3m1
  INTEGER :: i, j, k, n1p, n2p, ir

  real(kind=8) :: g2, t(3)

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
               n1p = i + 1
               if (n1p.lt.1) n1p = n1p + nr1
               n2p = j + 1
               if (n2p.lt.1) n2p = n2p + nr2
               if ( isind( n1p + (n2p-1)*ldis ) .eq. 0 ) cycle loop_z

               g2=0.d0
               do ir=1,3
                  t(ir) = dble(i)*b1(ir) + dble(j)*b2(ir) + dble(k)*b3(ir)
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
END SUBROUTINE



!-------------------------------------------------------------------------
SUBROUTINE gglobal( ng_g, g2_g, mill_g, b1, b2, b3, nr1, nr2, nr3, gcut, lgam )
!-------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: ng_g
  INTEGER :: mill_g(3,*)
  real(kind=8) :: g2_g(*)
  integer :: nr1, nr2, nr3
  real(kind=8) :: b1(3), b2(3), b3(3), gcut
  LOGICAL :: lgam

  INTEGER :: nr1m1, nr2m1, nr3m1
  INTEGER :: i, j, k, ir, ng

  real(kind=8) :: g2, t(3)


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
                  t(ir) = dble(i)*b1(ir)+dble(j)*b2(ir)+dble(k)*b3(ir)
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
!          WRITE( 201, fmt="( I6, 3I4, 2D25.16 )" ) &
!            ig, mill_g(1,ig), mill_g(2,ig), mill_g(3,ig), g2_g( ig ), g2sort_g( ig )
!        END DO
!        CLOSE( 201 )
!      END IF

  RETURN
END SUBROUTINE



!-------------------------------------------------------------------------
SUBROUTINE glocal( ng, g, ig_l2g, mill_l, ng_g, g2_g, mill_g, nr1, nr2, nr3, isind, ldis )
!-------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: ng_g, ng
  INTEGER :: mill_g(3,*), ig_l2g(*), mill_l(3,*)
  real(kind=8) :: g2_g(*), g(*)
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
        n1p = i + 1
        if (n1p.lt.1) n1p = n1p + nr1
        n2p = j + 1
        if (n2p.lt.1) n2p = n2p + nr2
        if (isind(n1p+(n2p-1)*ldis).eq.0) cycle loop_allg
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

      deallocate( index )

  RETURN
END SUBROUTINE



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
END SUBROUTINE


!-------------------------------------------------------------------------
SUBROUTINE gfftindex( np, nm, ng, mill_l, nr1, nr2, nr3, isind, nr1x, nr2x, nr3x )
!-------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: ng
  INTEGER :: isind(*), nr1x, nr2x, nr3x
  INTEGER :: mill_l(3,*), np(*), nm(*)
  integer :: nr1, nr2, nr3

  INTEGER :: n1p, n2p, n3p
  INTEGER :: n1m, n2m, n3m
  INTEGER :: i, j, k, ig


      do ig=1,ng
         i=mill_l(1,ig)
         j=mill_l(2,ig)
         k=mill_l(3,ig)
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
! conversion from (i,j,k) index to combined 1-d ijk index
! for the parallel case: columns along z are stored contiguously
!
#ifdef __PARA
         np(ig) = n3p + (isind(n1p+(n2p-1)*nr1x)-1)*nr3x
         nm(ig) = n3m + (isind(n1m+(n2m-1)*nr1x)-1)*nr3x
#else
         np(ig) = n1p + (n2p-1)*nr1x + (n3p-1)*nr1x*nr2x
         nm(ig) = n1m + (n2m-1)*nr1x + (n3m-1)*nr1x*nr2x
#endif
!
! conversion from (i,j,k) index to combined 1-d ijk index:
! ijk = 1 + (i-1)+(j-1)*ix+(k-1)*ix*jx
! where the (i,j,k) array is assumed to be dimensioned (ix,jx,kx)
!
      end do

  RETURN
END SUBROUTINE


!-------------------------------------------------------------------------
SUBROUTINE gshcount( ngl, ngsl, ngwl, igl, ng, g, gcuts, gcutw )
!-------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: ngl, ngsl, ngwl
  INTEGER :: igl(*)
  INTEGER :: ng
  REAL(kind=8) :: g(*), gcuts, gcutw

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
END SUBROUTINE


!-------------------------------------------------------------------------
      subroutine gcal(b1,b2,b3,gmax)
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
      use control_flags, only: iprint
      use reciprocal_vectors, only: g, gx, mill_l
      use gvecp, only: ngm
      use gvecw, only: ngw
      use gvecw, only: ggp, agg => ecutz, sgg => ecsig, e0gg => ecfix
      use cell_base, only: tpiba2
      implicit none
!
      real(kind=8) b1(3),b2(3),b3(3), gmax
      real(kind=8), external :: erf
!
      integer i1,i2,i3,ig
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
!
      do ig=1,ngw
         ggp(ig) = g(ig) +                                              &
     &             (agg/tpiba2)*(1.0+erf((tpiba2*g(ig)-e0gg)/sgg))
      enddo
!
      return
      end subroutine


!=----------------------------------------------------------------------------=!

        SUBROUTINE newgb( a1, a2, a3, omega, alat )
!
!     re-generation of little box g-vectors
!
          USE kinds, ONLY: dbl
          USE grid_dimensions, only: nr1, nr2, nr3
          USE smallbox_grid_dimensions, only: nr1b, nr2b, nr3b
          USE small_box, only: a1b, a2b, a3b, ainvb, omegab, tpibab
          USE constants, ONLY: pi

          IMPLICIT NONE
          REAL(dbl) :: a1( 3 ), a2( 3 ), a3( 3 ), omega, alat

          INTEGER :: i
          REAL(dbl) :: alatb, b1b(3),b2b(3),b3b(3)

          alatb  = alat/nr1*nr1b
          tpibab = 2.d0*pi/alatb
          do i=1,3
            a1b(i)=a1(i)/nr1*nr1b
            a2b(i)=a2(i)/nr2*nr2b
            a3b(i)=a3(i)/nr3*nr3b
          enddo
          omegab=omega/nr1*nr1b/nr2*nr2b/nr3*nr3b
!
          call recips(a1b,a2b,a3b,b1b,b2b,b3b)
          b1b = b1b * alatb
          b2b = b2b * alatb
          b3b = b3b * alatb
          call gcalb(b1b,b2b,b3b,nr1b,nr2b,nr3b)
!
          do i=1,3
            ainvb(1,i)=b1b(i)/alatb
            ainvb(2,i)=b2b(i)/alatb
            ainvb(3,i)=b3b(i)/alatb
          end do

          RETURN
        END SUBROUTINE

!------------------------------------------------------------------------------!
!
!
!------------------------------------------------------------------------------!

        SUBROUTINE ecutoffs_setup( ecutwfc, ecutrho, ecfixed, qcutz, q2sigma )
 
          USE kinds, ONLY: dbl
          USE constants, ONLY: eps8
          USE gvecw, ONLY: ecutw
          USE gvecw, ONLY: ecfix, ecutz, ecsig, tecfix
          USE gvecp, ONLY: ecutp
          USE gvecs, ONLY: ecuts

          IMPLICIT NONE
          REAL(dbl), INTENT(IN) ::  ecutwfc, ecutrho, ecfixed, qcutz, q2sigma
          REAL(dbl) :: dual

          ecutw = ecutwfc
          ecutp = ecutrho
          IF( ecutp <= 0.d0 ) THEN
              ecutp = 4.0d0 * ecutw
          ELSE IF( ecutp < ecutw ) THEN
              CALL errore(' ecutoffs_setup ',' invalid ecutrho ', INT( ecutp ) )
          END IF
          dual  = 4.0d0
          ecuts = dual * ecutwfc
          ecfix = ecfixed
          ecutz = qcutz
          ecsig = q2sigma
          IF( ABS( ecutz ) < eps8 ) THEN
            tecfix = .FALSE.
          ELSE
            tecfix = .TRUE.
          ENDIF

          RETURN
        END SUBROUTINE


        SUBROUTINE gcutoffs_setup( alat, tk_inp, nk_inp, kpoints_inp )

!  (describe briefly what this routine does...)
!  ----------------------------------------------

          USE kinds, ONLY: dbl
          USE gvecw, ONLY: ecutwfc => ecutw,  gcutw
          USE gvecp, ONLY: ecutrho => ecutp,  gcutp
          USE gvecs, ONLY: ecuts, gcuts
          USE gvecb, ONLY: ecutb, gcutb
          USE gvecw, ONLY: ecfix, gcfix, ecutz, gcutz, esig => ecsig, g2sig => gcsig, tecfix
          USE gvecw, ONLY: ekcut, gkcut
          USE constants, ONLY: eps8, pi

          IMPLICIT NONE

! ...     declare subroutine arguments
          REAL(dbl), INTENT(IN) :: alat
          INTEGER, INTENT(IN) :: nk_inp
          LOGICAL, INTENT(IN) :: tk_inp
          REAL(dbl), INTENT(IN) :: kpoints_inp(3,*)

! ...     declare other variables
          INTEGER i
          REAL(dbl) kcut, ksq, dual
          REAL(dbl) tpiba

!  end of declarations
!  ----------------------------------------------

! ...   Set Values for the cutoff


          IF( alat < eps8 ) THEN
            CALL errore(' cut-off setup ', ' alat too small ', 0)
          END IF

          tpiba = 2.0d0 * pi / alat 

          ! ...  Constant cutoff simulation parameters

          gcfix = ecfix / tpiba**2
          gcutz = ecutz / tpiba**2
          g2sig = esig  / tpiba**2

          gcutw = ecutwfc / tpiba**2  ! wave function cut-off
          gcutp = ecutrho / tpiba**2  ! potential cut-off

          gcuts = ecuts / tpiba**2    ! smooth mesh cut-off

          kcut = 0.0_dbl
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

