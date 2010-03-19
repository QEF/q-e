MODULE pw2blip
   USE kinds, ONLY: DP

   USE io_global, ONLY: ionode, ionode_id
   USE mp_global, ONLY: me_pool,nproc_pool,intra_pool_comm
   USE mp, ONLY: mp_get

   PRIVATE
   PUBLIC pw2blip_init,pw2blip_cleanup,pw2blip_transform
   PUBLIC blipgrid,cavc,pw2blip_get

   INTEGER :: ngtot
   COMPLEX(dp),ALLOCATABLE :: psic(:),cavc_flat(:)
   INTEGER :: blipgrid(3),ld_bg(3),bg_vol
   REAL(dp),ALLOCATABLE :: gamma(:)

   INTEGER,PARAMETER :: gamma_approx = 1
   REAL(dp),PARAMETER :: pi = 3.14159265358979324d0

   INTEGER,ALLOCATABLE :: map_igk_to_fft(:)
   INTEGER,ALLOCATABLE :: do_fft_x(:),do_fft_y(:)

CONTAINS

   SUBROUTINE pw2blip_init(ngtot_in,g_vec,multiplicity)
      USE cell_base, ONLY: at
      USE fft_scalar, ONLY: allowed, good_fft_dimension

      INTEGER,INTENT(in) :: ngtot_in
      REAL(dp),INTENT(in) :: g_vec(3,ngtot_in)
      REAL(dp),INTENT(in) :: multiplicity
      REAL(dp) :: da(3),k,k2,k4,cosk
      INTEGER :: ig,d,g_int(3,ngtot_in),g_idx(3),dummy(0)
      INTEGER,PARAMETER :: nmax = 5000

      ngtot = ngtot_in

      do ig=1,ngtot
         g_int(1,ig) = nint (sum(g_vec(:,ig) * at (:,1)))
         g_int(2,ig) = nint (sum(g_vec(:,ig) * at (:,2)))
         g_int(3,ig) = nint (sum(g_vec(:,ig) * at (:,3)))
      enddo

      ! choose size of blip grid in real space
      do d=1,3
         blipgrid(d) = 2*ceiling(dble(maxval(abs(g_int(d,:))))*multiplicity)+2
         do while(.not.allowed(blipgrid(d)))
            blipgrid(d) = blipgrid(d) + 1
         enddo
         if (blipgrid(d)>nmax) &
            call errore ('set_fft_dim', 'blipgrid is unreasonably large', blipgrid(d))
      enddo

      ! set up leading dimensions of fft data array
      ld_bg(1) = good_fft_dimension(blipgrid(1))
      ld_bg(2) = blipgrid(2)
      ld_bg(3) = blipgrid(3)
      bg_vol = ld_bg(1)*ld_bg(2)*ld_bg(3)

! Set up indices to fft grid: map_igk_to_fft
      allocate(map_igk_to_fft(ngtot))
      allocate(do_fft_x(blipgrid(3)*ld_bg(2)),do_fft_y(blipgrid(3)))
      do_fft_x(:)=0
      do_fft_y(:)=0
      do ig=1,ngtot
         g_idx(:) = modulo(g_int(:,ig),blipgrid(:))
         do_fft_x(1 + g_idx(2) + ld_bg(2)*g_idx(3)) = 1
         do_fft_y(1 + g_idx(3)) = 1
         map_igk_to_fft (ig) = 1 + g_idx(1) + ld_bg(1)*(g_idx(2) + ld_bg(2)*g_idx(3))
      enddo

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
      deallocate(psic,gamma)
      deallocate(map_igk_to_fft,do_fft_x,do_fft_y)
   END SUBROUTINE pw2blip_cleanup

   SUBROUTINE pw2blip_transform(psi)
      USE fft_scalar, ONLY: cfft3ds

      COMPLEX(DP), INTENT(in) :: psi(ngtot)

      psic (:) = (0.d0, 0.d0)
      psic (map_igk_to_fft (1:ngtot)) = psi(1:ngtot)*gamma(1:ngtot)

      ! perform the transformation
      call cfft3ds (psic,blipgrid(1),blipgrid(2),blipgrid(3),&
       &ld_bg(1),ld_bg(2),ld_bg(3),+1,do_fft_x(:),do_fft_y(:))
   END SUBROUTINE

   SUBROUTINE pw2blip_get(node)
      INTEGER,INTENT(in) :: node
      IF(ionode_id /= node)&
         &CALL mp_get(psic,psic,me_pool,ionode_id,node,1234,intra_pool_comm)
   END SUBROUTINE pw2blip_get

   COMPLEX(dp) FUNCTION cavc(i1,i2,i3)
      INTEGER,INTENT(in) :: i1,i2,i3
      cavc = psic(i1+ld_bg(2)*(i2-1+ld_bg(3)*(i3-1)))
   END FUNCTION cavc

END MODULE
