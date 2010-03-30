MODULE pw2blip
   USE kinds, ONLY: DP

   USE io_global, ONLY: ionode, ionode_id
   USE mp_global, ONLY: me_pool,nproc_pool,intra_pool_comm
   USE mp, ONLY: mp_get
   USE control_flags, ONLY: gamma_only
   USE constants, ONLY: tpi
   USE cell_base, ONLY: at,alat
   USE fft_scalar, ONLY: allowed, good_fft_dimension

   PRIVATE
   PUBLIC pw2blip_init,pw2blip_cleanup,pw2blip_transform,pw2blip_transform2,&
    &blipgrid,cavc,avc1,avc2,pw2blip_get,pw2blip_stat,blipeval,blip3dk,g_int,phase1,phase2

   INTEGER,PUBLIC :: blipreal = 0
   ! blipreal == 0 -- complex wfn1
   ! blipreal == 1 -- one real wfn (.not.gamma_only)
   ! blipreal == 2 -- two real wfn (.not.gamma_only)
   ! blipreal == -1 -- one real wfn (gamma_only)
   ! blipreal == -2 -- two real wfn (gamma_only)

   INTEGER :: ngtot
   COMPLEX(dp),ALLOCATABLE :: psic(:),cavc_flat(:)
   INTEGER :: blipgrid(3),ld_bg(3),bg_vol

   REAL(dp),ALLOCATABLE :: gamma(:)

   INTEGER,PARAMETER :: gamma_approx = 1
   REAL(dp),PARAMETER :: pi = 3.14159265358979324d0

   INTEGER,ALLOCATABLE :: map_igk_to_fft(:)
   INTEGER,ALLOCATABLE :: map_minus_igk_to_fft(:) ! gamma_only
   INTEGER,ALLOCATABLE :: map_neg_igk(:)          ! blipreal.and..not.gamma_only
   LOGICAL,ALLOCATABLE :: unique_igk(:)           ! blipreal.and..not.gamma_only
   INTEGER,ALLOCATABLE :: do_fft_x(:),do_fft_y(:)

   REAL(dp) :: norm_real(2),norm_imag(2)
   COMPLEX(DP) :: phase1,phase2

   INTEGER :: nr(3)
   INTEGER,ALLOCATABLE :: g_int(:,:)
   REAL(dp) :: rnr(3),rnr2(3),bg(3,3),lvp(6)

CONTAINS

   SUBROUTINE pw2blip_init(ngtot_in,g_vec,multiplicity)
      INTEGER,INTENT(in) :: ngtot_in
      REAL(dp),INTENT(in) :: g_vec(3,ngtot_in)
      REAL(dp),INTENT(in) :: multiplicity
      REAL(dp) :: da(3),k,k2,k4,cosk
      INTEGER :: ig,ig2,d,g_idx(3),dummy(0)
      INTEGER,PARAMETER :: nmax = 5000

      ngtot = ngtot_in

      allocate(g_int(3,ngtot))
      do ig=1,ngtot
         g_int(1,ig) = nint (sum(g_vec(:,ig) * at (:,1)))
         g_int(2,ig) = nint (sum(g_vec(:,ig) * at (:,2)))
         g_int(3,ig) = nint (sum(g_vec(:,ig) * at (:,3)))
      enddo

      if(any(g_int(:,1)/=0))then
         CALL errore('pw2blip_init','first G vector is not zero',0)
      endif

      ! choose size of blip grid in real space
      do d=1,3
         blipgrid(d) = 2*ceiling(dble(maxval(abs(g_int(d,:))))*multiplicity)+2
         do while(.not.allowed(blipgrid(d)))
            blipgrid(d) = blipgrid(d) + 1
         enddo
         if (blipgrid(d)>nmax) &
            call errore ('pw2blip_init', 'blipgrid is unreasonably large', blipgrid(d))
      enddo

      nr(:) = blipgrid(:)
      rnr(:) = dble(nr(:))
      rnr2(:) = rnr(:)*rnr(:)

      call inve(at,bg)
      bg=transpose(bg)
      lvp(1)=bg(1,1)**2+bg(2,1)**2+bg(3,1)**2
      lvp(2)=bg(1,2)**2+bg(2,2)**2+bg(3,2)**2
      lvp(3)=bg(1,3)**2+bg(2,3)**2+bg(3,3)**2
      lvp(4)=2.d0*(bg(1,1)*bg(1,2)+bg(2,1)*bg(2,2)+bg(3,1)*bg(3,2))
      lvp(5)=2.d0*(bg(1,2)*bg(1,3)+bg(2,2)*bg(2,3)+bg(3,2)*bg(3,3))
      lvp(6)=2.d0*(bg(1,3)*bg(1,1)+bg(2,3)*bg(2,1)+bg(3,3)*bg(3,1))

      ! set up leading dimensions of fft data array
      ld_bg(1) = good_fft_dimension(blipgrid(1))
      ld_bg(2) = blipgrid(2)
      ld_bg(3) = blipgrid(3)
      bg_vol = ld_bg(1)*ld_bg(2)*ld_bg(3)

! Set up indices to fft grid: map_igk_to_fft
      allocate(map_igk_to_fft(ngtot))
!      map_igk_to_fft(1) = 1
      if(blipreal<0)then ! gamma_only
         allocate(map_minus_igk_to_fft(ngtot))
!         map_minus_igk_to_fft(1) = 1
      elseif(blipreal>0)then
         allocate(map_neg_igk(ngtot),unique_igk(ngtot))
         map_neg_igk(:)=0
!         map_neg_igk(1)=1
         unique_igk(:)=.true.
      endif
      allocate(do_fft_x(blipgrid(3)*ld_bg(2)),do_fft_y(blipgrid(3)))
      do_fft_x(:)=0 ; do_fft_y(:)=0
!      do_fft_x(1)=1 ; do_fft_y(1)=1
      do ig=1,ngtot
         g_idx(:) = modulo(g_int(:,ig),blipgrid(:))
         do_fft_x(1 + g_idx(2) + ld_bg(2)*g_idx(3)) = 1
         do_fft_y(1 + g_idx(3)) = 1
         map_igk_to_fft (ig) = 1 + g_idx(1) + ld_bg(1)*(g_idx(2) + ld_bg(2)*g_idx(3))
         if(blipreal<0)then ! gamma_only
            g_idx(:) = modulo(-g_int(:,ig),blipgrid(:))
            do_fft_x(1 + g_idx(2) + ld_bg(2)*g_idx(3)) = 1
            do_fft_y(1 + g_idx(3)) = 1
            map_minus_igk_to_fft (ig) = 1 + g_idx(1) + ld_bg(1)*(g_idx(2) + ld_bg(2)*g_idx(3))
         elseif(blipreal>0)then
            if(all(g_int(:,ig)==0))then
               map_neg_igk(ig)=ig
            elseif(unique_igk(ig))then
               do ig2=ig,ngtot
                  if(all(g_int(:,ig)+g_int(:,ig2)==0))then
                     unique_igk(ig2)=.false.
                     map_neg_igk(ig)=ig2
                     map_neg_igk(ig2)=ig
                     exit
                  endif
               enddo
            endif
         endif
      enddo
      if(blipreal>0)then !.not.gamma_only
         if(any(map_neg_igk(:)==0))then
            do ig=1,ngtot
               write(0,*)ig,g_int(:,ig),map_neg_igk(ig),unique_igk(ig)
            enddo
            CALL errore( 'pw2blip_init','G points do not pair up correctly',0)
         endif
!          if(any(unique_igk(map_neg_igk(2:)).eqv.unique_igk(2:)))then
!             do ig=1,ngtot
!                write(0,*)ig,g_int(:,ig),map_neg_igk(ig),unique_igk(ig)
!             enddo
!             CALL errore( 'pw2blip_init','G points do not pair up correctly',1)
!          endif
         if(any(map_neg_igk(map_neg_igk(:))/=(/(ig,ig=1,ngtot)/)))then
            do ig=1,ngtot
               write(0,*)ig,g_int(:,ig),map_neg_igk(ig),unique_igk(ig)
            enddo
            CALL errore( 'pw2blip_init','G points do not pair up correctly',2)
         endif
      endif

! Set up blipgrid
      allocate(psic(bg_vol)) ! local FFT grid for transform

! Calculating gamma.
      allocate(gamma(ngtot))
      gamma(:) = 1.d0
      da(1:3)=2.d0*pi/dble( blipgrid(:) )
      if(gamma_approx==1)then
         do ig=1,ngtot
            do d=1,3
               if(g_int(d,ig)/=0)then
                  k=da(d)*dble(g_int(d,ig)) ; cosk=cos(k) ; k2=k*k ; k4=k2*k2
                  gamma(ig)=gamma(ig)*k4/(6.d0*((cosk-2.d0)*cosk+1.d0))
               else
                  gamma(ig)=gamma(ig)*2.d0/3.d0
               endif
            enddo
         enddo ! ig
      elseif(gamma_approx==2)then
         do ig=1,ngtot
            gamma(ig)=1.d0/(&
               & (1.d0+0.5d0*cos(da(1)*g_vec(1,ig))) &
               &*(1.d0+0.5d0*cos(da(2)*g_vec(2,ig))) &
               &*(1.d0+0.5d0*cos(da(3)*g_vec(3,ig))) &
               &)
         enddo ! ig
      else
         write(6,*)'Bug: bad gamma_approx.' ; stop
      endif ! gamma_approx

   END SUBROUTINE pw2blip_init

   SUBROUTINE pw2blip_cleanup
      deallocate(psic,gamma,g_int)
      deallocate(map_igk_to_fft,do_fft_x,do_fft_y)
      if(blipreal<0)deallocate(map_minus_igk_to_fft) ! gammaonly
   END SUBROUTINE pw2blip_cleanup

! get_phase: find complex phase factor to rotate this orbital to the real plane
! PROBLEM: for degenerate orbitals, phase rotations are not sufficient to
! obtain orthogonal real blip wave function. A unitary transform of the whole
! eigenspace would be necessary.
! In practice, simply choosing real or imaginary part of each orbital has never
! caused serious problems.

! original version, as implemented in the blip.f90 converter tool
! destroys the normalization but seems to work reasonably well in most cases
   COMPLEX(DP) FUNCTION get_phase(psi,resqr,imsqr)
      COMPLEX(DP), INTENT(in) :: psi(ngtot)
      REAL(DP),INTENT(out) :: resqr,imsqr

      resqr = sum(abs(psi(:)+conjg(psi(map_neg_igk(:))))**2)
      imsqr = sum(abs(psi(:)-conjg(psi(map_neg_igk(:))))**2)

      if(resqr>imsqr)then
         get_phase = (1.d0,0.d0)
      else
         get_phase = (0.d0,-1.d0)
      endif
   END FUNCTION

! new version by Norbert Nemec:
! works perfectly for non-degenerate cases (finds the exact phase rotation)
! it is not known how it compares to the original version in the degenerate
! case
!    COMPLEX(DP) FUNCTION get_phase(psi)
!       COMPLEX(DP), INTENT(in) :: psi(ngtot)
!       COMPLEX(DP) :: phase
!
!       write(6,*)"complex wfn at G=0: ",psi(1)
!       phase = sqrt(sum(psi(:)*psi(map_neg_igk(:))))
!       write(6,*)"wfn phase factor ",phase
!       phase = abs(phase)/phase
!       write(6,*)"phase corrected wfn at G=0: ",psi(1)*phase
!       write(6,*)"wfn1 maximum imaginary deviation after phase correction",&
!          &maxval(abs(aimag((psi(:)+psi(map_neg_igk(:)))*phase)))
!       get_phase = phase
!    END FUNCTION


   SUBROUTINE pw2blip_transform(psi)
      USE fft_scalar, ONLY: cfft3ds

      COMPLEX(DP), INTENT(in) :: psi(ngtot)

      if(blipreal<0)then ! gamma_only
         blipreal = -1

         phase1 = (1.d0,0.d0)

         psic (:) = (0.d0, 0.d0)
         psic (map_igk_to_fft (1:ngtot)) = psi(1:ngtot)*gamma(1:ngtot)
         psic (map_minus_igk_to_fft (1:ngtot)) = conjg(psi(1:ngtot))*gamma(1:ngtot)

      elseif(blipreal>0)then ! real wfn
         blipreal = 1

         phase1 = get_phase(psi,norm_real(1),norm_imag(1))

         psic (:) = (0.d0, 0.d0)
         psic (map_igk_to_fft (:)) = (0.5d0,0.d0)*(phase1*psi(:)+conjg(phase1*psi(map_neg_igk(:))))*gamma(:)

      else ! complex wfn
         phase1 = (1.d0,0.d0)

         psic (:) = (0.d0, 0.d0)
         psic (map_igk_to_fft (1:ngtot)) = psi(1:ngtot)*gamma(1:ngtot)
      endif
      phase2 = (0.d0,0.d0)

      ! perform the transformation
      call cfft3ds (psic,blipgrid(1),blipgrid(2),blipgrid(3),&
       &ld_bg(1),ld_bg(2),ld_bg(3),+1,do_fft_x(:),do_fft_y(:))
   END SUBROUTINE

   SUBROUTINE pw2blip_transform2(psi1,psi2)
      USE fft_scalar, ONLY: cfft3ds

      COMPLEX(DP), INTENT(in) :: psi1(ngtot),psi2(ngtot)

      if(blipreal<0)then ! gamma_only
         blipreal = -2

         phase1 = (1.d0,0.d0)
         phase2 = (1.d0,0.d0)

         psic (:) = (0.d0, 0.d0)
         psic (map_igk_to_fft (1:ngtot)) = (psi1(1:ngtot)+(0.d0,1.d0)*psi2(1:ngtot))*gamma(1:ngtot)
         psic (map_minus_igk_to_fft (1:ngtot)) = conjg((psi1(1:ngtot)-(0.d0,1.d0)*psi2(1:ngtot)))*gamma(1:ngtot)
      elseif(blipreal>0)then ! real wfn
         blipreal = 2

         phase1 = get_phase(psi1,norm_real(1),norm_imag(1))
         phase2 = get_phase(psi2,norm_real(2),norm_imag(2))

         psic (:) = (0.d0, 0.d0)
         psic (map_igk_to_fft (:)) = (&
            &(0.5d0,0.d0)*(phase1*psi1(:)+conjg(phase1*psi1(map_neg_igk(:))))&
            & + (0.d0,0.5d0)*(phase2*psi2(:)+conjg(phase2*psi2(map_neg_igk(:))))&
            &)*gamma(:)
      else !
         call errore("pw2blip_transform2","BUG: can only perform one complex FFT at a time",3)
      endif


      ! perform the transformation
      call cfft3ds (psic,blipgrid(1),blipgrid(2),blipgrid(3),&
       &ld_bg(1),ld_bg(2),ld_bg(3),+1,do_fft_x(:),do_fft_y(:))
   END SUBROUTINE

   SUBROUTINE pw2blip_get(node)
      INTEGER,INTENT(in) :: node
      IF(ionode_id /= node)then
         CALL mp_get(psic,psic,me_pool,ionode_id,node,2498,intra_pool_comm)
         CALL mp_get(blipreal,blipreal,me_pool,ionode_id,node,2314,intra_pool_comm)
         CALL mp_get(norm_real(:),norm_real(:),me_pool,ionode_id,node,4532,intra_pool_comm)
         CALL mp_get(norm_imag(:),norm_imag(:),me_pool,ionode_id,node,1235,intra_pool_comm)
      endif
   END SUBROUTINE pw2blip_get

   SUBROUTINE pw2blip_stat(node,i)
      INTEGER,INTENT(in) :: node,i

      if(blipreal>0)then ! one real wfn
         if(ionode)write(6,*)"ratio (resqr:imsqr) "&
            &,norm_real(i)/(norm_real(i)+norm_imag(i)),norm_imag(i)/(norm_real(i)+norm_imag(i))
      endif
   END SUBROUTINE

   COMPLEX(dp) FUNCTION cavc(i1,i2,i3)
      INTEGER,INTENT(in) :: i1,i2,i3
      cavc = psic(1+i1+ld_bg(1)*(i2+ld_bg(2)*i3))
   END FUNCTION cavc

   REAL(dp) FUNCTION avc1(i1,i2,i3)
      INTEGER,INTENT(in) :: i1,i2,i3
      avc1 = real(psic(1+i1+ld_bg(1)*(i2+ld_bg(2)*i3)))
   END FUNCTION avc1

   REAL(dp) FUNCTION avc2(i1,i2,i3)
      INTEGER,INTENT(in) :: i1,i2,i3
      avc2 = aimag(psic(1+i1+ld_bg(1)*(i2+ld_bg(2)*i3)))
   END FUNCTION avc2


   SUBROUTINE blipeval(r,rpsi,grad,lap)
!----------------------------------------------------------------------------!
! This subroutine evaluates the value of a function, its gradient and its    !
! Laplacian at a vector point r, using the overlapping of blip functions.    !
! The blip grid is defined on a cubic cell, so r should always be given in   !
! units of the crystal lattice vectors.                                      !
!----------------------------------------------------------------------------!
      IMPLICIT NONE
      DOUBLE PRECISION,INTENT(in) :: r(3)
      COMPLEX(dp),INTENT(out) :: rpsi,grad(3),lap

      REAL(dp) t(3)
      INTEGER i(3),idx(3,4),jx,jy,jz
      REAL(dp) x(3),tx(3,4),dtx(3,4),d2tx(3,4)
      COMPLEX(dp) sderiv(6),C

      rpsi=(0.d0,0.d0) ; grad(:)=(0.d0,0.d0) ; sderiv(:)=(0.d0,0.d0)

      t(:) = r(:)*rnr(:)
      i(:) = modulo(floor(t(:)),nr(:))

      idx(:,1) = modulo(i(:)-1,nr(:))
      idx(:,2) = i(:)
      idx(:,3) = modulo(i(:)+1,nr(:))
      idx(:,4) = modulo(i(:)+2,nr(:))

      x(:)=t(:)-dble(idx(:,2)-1)
      tx(:,1)=2.d0+x(:)*(-3.d0+x(:)*(1.5d0-0.25d0*x(:))) ! == (8+x*(-12+x*(6-x)))/4 == (2-x)(4-2x+x2)/4
      dtx(:,1)=(-3.d0+x(:)*(3.d0-0.75d0*x(:)))*rnr(:)    ! == (-12+x*(12-3*x))r/4 == (2-x)(x-2)3r/4
      d2tx(:,1)=(3.d0-1.5d0*x(:))*rnr2(:)                ! == (2-x)3r2/2
      x(:)=t(:)-dble(idx(:,2))
      tx(:,2)=1.d0+x(:)*x(:)*(-1.5d0+0.75d0*x(:))        ! == (4-3x2(2-x))/4
      dtx(:,2)=x(:)*(-3.d0+2.25d0*x(:))*rnr(:)           ! == -x(12-9x)r/4
      d2tx(:,2)=(-3.d0+4.5d0*x(:))*rnr2(:)               ! == -(6-9x)r2/2
      x(:)=t(:)-dble(idx(:,2)+1)
      tx(:,3)=1.d0+x(:)*x(:)*(-1.5d0-0.75d0*x(:))        ! == (4-3x2(2+x))/4
      dtx(:,3)=x(:)*(-3.d0-2.25d0*x(:))*rnr(:)           ! == -x(12+9x)r/4
      d2tx(:,3)=(-3.d0-4.5d0*x(:))*rnr2(:)               ! == -(6+9x)r2/2
      x(:)=t(:)-dble(idx(:,2)+2)
      tx(:,4)=2.d0+x(:)*(3.d0+x(:)*(1.5d0+0.25d0*x(:)))  ! == (8+x*(12+x*(6+x)))/4 == (2+x)(4+2x+x2)/4
      dtx(:,4)=(3.d0+x(:)*(3.d0+0.75d0*x(:)))*rnr(:)     ! == (12+x*(12+3*x))r/4 == (2+x)(x+2)3r/4
      d2tx(:,4)=(3.d0+1.5d0*x(:))*rnr2(:)                ! == (2+x)3r2/2

      do jx=1,4
         do jy=1,4
            do jz=1,4
               C = cavc(idx(1,jx),idx(2,jy),idx(3,jz))
               rpsi = rpsi + C * tx(1,jx)*tx(2,jy)*tx(3,jz)
               grad(1) = grad(1) + C * dtx(1,jx)*tx(2,jy)*tx(3,jz)
               grad(2) = grad(2) + C * tx(1,jx)*dtx(2,jy)*tx(3,jz)
               grad(3) = grad(3) + C * tx(1,jx)*tx(2,jy)*dtx(3,jz)
               sderiv(1) = sderiv(1) + C * d2tx(1,jx)*tx(2,jy)*tx(3,jz)
               sderiv(2) = sderiv(2) + C * tx(1,jx)*d2tx(2,jy)*tx(3,jz)
               sderiv(3) = sderiv(3) + C * tx(1,jx)*tx(2,jy)*d2tx(3,jz)
               sderiv(4) = sderiv(4) + C * dtx(1,jx)*dtx(2,jy)*tx(3,jz)
               sderiv(5) = sderiv(5) + C * tx(1,jx)*dtx(2,jy)*dtx(3,jz)
               sderiv(6) = sderiv(6) + C * dtx(1,jx)*tx(2,jy)*dtx(3,jz)
            enddo
         enddo
      enddo

! Transformation of gradient to the Cartesian grid
      grad(1:3)=matmul(bg/alat,grad(1:3))

! The Laplacian: summing all contributions with appropriate transformation
      lap= sum(sderiv(:)*lvp(:))/alat**2

   END SUBROUTINE blipeval


   SUBROUTINE inve(v,inv)
!-----------------------!
! Inverts 3x3 matrices. !
!-----------------------!
   IMPLICIT NONE
      REAL(dp),INTENT(in) :: v(3,3)
      REAL(dp),INTENT(out) :: inv(3,3)
      REAL(dp) d
      d=v(1,1)*(v(2,2)*v(3,3)-v(2,3)*v(3,2))+ &
         &v(2,1)*(v(3,2)*v(1,3)-v(1,2)*v(3,3))+ &
         &v(3,1)*(v(1,2)*v(2,3)-v(1,3)*v(2,2))
      if(d==0.d0)then
         write(6,*)'Trying to invert a singular determinant.'
         stop
      endif
      d=1.d0/d
      inv(1,1)=(v(2,2)*v(3,3)-v(2,3)*v(3,2))*d
      inv(1,2)=(v(3,2)*v(1,3)-v(1,2)*v(3,3))*d
      inv(1,3)=(v(1,2)*v(2,3)-v(1,3)*v(2,2))*d
      inv(2,1)=(v(3,1)*v(2,3)-v(2,1)*v(3,3))*d
      inv(2,2)=(v(1,1)*v(3,3)-v(3,1)*v(1,3))*d
      inv(2,3)=(v(2,1)*v(1,3)-v(1,1)*v(2,3))*d
      inv(3,1)=(v(2,1)*v(3,2)-v(2,2)*v(3,1))*d
      inv(3,2)=(v(3,1)*v(1,2)-v(1,1)*v(3,2))*d
      inv(3,3)=(v(1,1)*v(2,2)-v(1,2)*v(2,1))*d
   END SUBROUTINE inve

END MODULE
