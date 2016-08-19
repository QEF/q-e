MODULE pw2blip
   USE kinds, ONLY: DP

   USE io_global, ONLY: ionode, ionode_id
   USE mp_pools, ONLY: me_pool,nproc_pool,intra_pool_comm
   USE mp, ONLY: mp_get
   USE control_flags, ONLY: gamma_only
   USE constants, ONLY: tpi
   USE cell_base, ONLY: at,alat
   USE fft_support, ONLY: allowed, good_fft_dimension

   PRIVATE
   PUBLIC pw2blip_init,pw2blip_cleanup,pw2blip_transform,pw2blip_transform2,&
    &blipgrid,cavc,avc1,avc2,pw2blip_get,blipeval,blip3dk,g_int

   INTEGER,PUBLIC :: blipreal = 0
   ! blipreal == 0 -- complex wfn1
   ! blipreal == 1 -- one real wfn (gamma_only)
   ! blipreal == 2 -- two real wfn (gamma_only)

   INTEGER :: ngtot
   COMPLEX(dp),ALLOCATABLE :: psic(:),cavc_flat(:)
   INTEGER :: blipgrid(3),ld_bg(3),bg_vol

   REAL(dp),ALLOCATABLE :: gamma(:)

   INTEGER,PARAMETER :: gamma_approx = 1
   REAL(dp),PARAMETER :: pi = 3.14159265358979324d0

   INTEGER,ALLOCATABLE :: map_igk_to_fft(:)
   INTEGER,ALLOCATABLE :: map_minus_igk_to_fft(:) ! gamma_only
   INTEGER,ALLOCATABLE :: do_fft_z(:),do_fft_y(:)

   INTEGER :: nr(3)
   INTEGER,ALLOCATABLE :: g_int(:,:)
   REAL(dp) :: rnr(3),rnr2(3),bg(3,3),lvp(6)

CONTAINS

   SUBROUTINE pw2blip_init(ngtot_in,g_vec,multiplicity)
      INTEGER,INTENT(in) :: ngtot_in
      REAL(dp),INTENT(in) :: g_vec(3,ngtot_in)
      REAL(dp),INTENT(in) :: multiplicity
      REAL(dp) :: da(3),k,k2,k4,cosk
      INTEGER :: ig,ig2,d,g_idx(3)
      INTEGER,PARAMETER :: nmax = 5000

      ngtot = ngtot_in

      ALLOCATE(g_int(3,ngtot))
      DO ig=1,ngtot
         g_int(1,ig) = nint (sum(g_vec(:,ig) * at (:,1)))
         g_int(2,ig) = nint (sum(g_vec(:,ig) * at (:,2)))
         g_int(3,ig) = nint (sum(g_vec(:,ig) * at (:,3)))
      ENDDO

      IF(any(g_int(:,1)/=0))THEN
         CALL errore('pw2blip_init','first G vector is not zero',0)
      ENDIF

      ! choose size of blip grid in real space
      DO d=1,3
         blipgrid(d) = 2*ceiling(dble(maxval(abs(g_int(d,:))))*multiplicity)+2
         DO WHILE(.not.allowed(blipgrid(d)))
            blipgrid(d) = blipgrid(d) + 1
         ENDDO
         IF (blipgrid(d)>nmax) &
            CALL errore ('pw2blip_init', 'blipgrid is unreasonably large', blipgrid(d))
      ENDDO

      nr(:) = blipgrid(:)
      rnr(:) = dble(nr(:))
      rnr2(:) = rnr(:)*rnr(:)

      CALL inve(at,bg)
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
      ALLOCATE(map_igk_to_fft(ngtot))
!      map_igk_to_fft(1) = 1
      IF(gamma_only)THEN
         ALLOCATE(map_minus_igk_to_fft(ngtot))
!         map_minus_igk_to_fft(1) = 1
      ENDIF
      ALLOCATE(do_fft_z(ld_bg(1)*ld_bg(2)),do_fft_y(blipgrid(1)))
      do_fft_z(:)=0 ; do_fft_y(:)=0
!      do_fft_x(1)=1 ; do_fft_y(1)=1
      DO ig=1,ngtot
         g_idx(:) = modulo(g_int(:,ig),blipgrid(:))
         do_fft_z(1 + g_idx(1) + ld_bg(1)*g_idx(2)) = 1
         do_fft_y(1 + g_idx(1)) = 1
         map_igk_to_fft (ig) = 1 + g_idx(1) + ld_bg(1)*(g_idx(2) + ld_bg(2)*g_idx(3))
         IF(gamma_only)THEN ! gamma_only
            g_idx(:) = modulo(-g_int(:,ig),blipgrid(:))
            do_fft_z(1 + g_idx(1) + ld_bg(1)*g_idx(2)) = 1
            do_fft_y(1 + g_idx(1)) = 1
            map_minus_igk_to_fft (ig) = 1 + g_idx(1) + ld_bg(1)*(g_idx(2) + ld_bg(2)*g_idx(3))
         ENDIF
      ENDDO

! Set up blipgrid
      ALLOCATE(psic(bg_vol)) ! local FFT grid for transform

! Calculating gamma.
      ALLOCATE(gamma(ngtot))
      gamma(:) = 1.d0
      da(1:3)=2.d0*pi/dble( blipgrid(:) )
      IF(gamma_approx==1)THEN
         DO ig=1,ngtot
            DO d=1,3
               IF(g_int(d,ig)/=0)THEN
                  k=da(d)*dble(g_int(d,ig)) ; cosk=cos(k) ; k2=k*k ; k4=k2*k2
                  gamma(ig)=gamma(ig)*k4/(6.d0*((cosk-2.d0)*cosk+1.d0))
               ELSE
                  gamma(ig)=gamma(ig)*2.d0/3.d0
               ENDIF
            ENDDO
         ENDDO ! ig
      ELSEIF(gamma_approx==2)THEN
         DO ig=1,ngtot
            gamma(ig)=1.d0/(&
               & (1.d0+0.5d0*cos(da(1)*g_vec(1,ig))) &
               &*(1.d0+0.5d0*cos(da(2)*g_vec(2,ig))) &
               &*(1.d0+0.5d0*cos(da(3)*g_vec(3,ig))) &
               &)
         ENDDO ! ig
      ELSE
         WRITE(6,*)'Bug: bad gamma_approx.' ; STOP
      ENDIF ! gamma_approx

   END SUBROUTINE pw2blip_init

   SUBROUTINE pw2blip_cleanup
      DEALLOCATE(psic,gamma,g_int)
      DEALLOCATE(map_igk_to_fft,do_fft_z,do_fft_y)
      IF(gamma_only)DEALLOCATE(map_minus_igk_to_fft)
   END SUBROUTINE pw2blip_cleanup

   SUBROUTINE pw2blip_transform(psi)
      USE fft_scalar, ONLY: cfft3ds

      COMPLEX(DP), INTENT(in) :: psi(ngtot)

      psic (:) = (0.d0, 0.d0)
      psic (map_igk_to_fft (1:ngtot)) = psi(1:ngtot)*gamma(1:ngtot)

      IF(gamma_only)THEN
         psic (map_minus_igk_to_fft (1:ngtot)) = conjg(psi(1:ngtot))*gamma(1:ngtot)
         blipreal = 1
      ENDIF

      ! perform the transformation
      CALL cfft3ds (psic,blipgrid(1),blipgrid(2),blipgrid(3),&
       &ld_bg(1),ld_bg(2),ld_bg(3),1,+1,do_fft_z(:),do_fft_y(:))
   END SUBROUTINE

   SUBROUTINE pw2blip_transform2(psi1,psi2)
      USE fft_scalar, ONLY: cfft3ds

      COMPLEX(DP), INTENT(in) :: psi1(ngtot),psi2(ngtot)

      IF(.not.gamma_only)THEN
         CALL errore("pw2blip_transform2","BUG: can only perform one complex FFT at a time",3)
      ENDIF

      blipreal = 2

      psic (:) = (0.d0, 0.d0)
      psic (map_igk_to_fft (1:ngtot)) = (psi1(1:ngtot)+(0.d0,1.d0)*psi2(1:ngtot))*gamma(1:ngtot)
      psic (map_minus_igk_to_fft (1:ngtot)) = conjg((psi1(1:ngtot)-(0.d0,1.d0)*psi2(1:ngtot)))*gamma(1:ngtot)

      ! perform the transformation
      CALL cfft3ds (psic,blipgrid(1),blipgrid(2),blipgrid(3),&
       &ld_bg(1),ld_bg(2),ld_bg(3),1,+1,do_fft_z(:),do_fft_y(:))
   END SUBROUTINE

   SUBROUTINE pw2blip_get(node)
      INTEGER,INTENT(in) :: node
      IF(ionode_id /= node)THEN
         CALL mp_get(psic,psic,me_pool,ionode_id,node,2498,intra_pool_comm)
         CALL mp_get(blipreal,blipreal,me_pool,ionode_id,node,2314,intra_pool_comm)
      ENDIF
   END SUBROUTINE pw2blip_get

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

      DO jx=1,4
         DO jy=1,4
            DO jz=1,4
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
            ENDDO
         ENDDO
      ENDDO

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
      IF(d==0.d0)THEN
         WRITE(6,*)'Trying to invert a singular determinant.'
         STOP
      ENDIF
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
