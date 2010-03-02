MODULE pw2blip
   USE kinds, ONLY: DP

   ! choose dense fft grid
   USE fft_base, ONLY: dfft => dfftp
   USE gvect, ONLY: map_igk_to_fft => nl,&
                    NGvec => ngm

   ! choose smooth fft grid
!   USE fft_base, ONLY: dfft => dffts
!   USE gsmooth, ONLY: map_igk_to_fft => nls,&
!                      NGvec => ngms

   ! choose coarse fft grid
!   USE gcoarse, ONLY: dfft => dfftc
!   USE gcoarse, ONLY: map_igk_to_fft => nlc

   USE io_global, ONLY: ionode, ionode_id

   PRIVATE
   PUBLIC pw2blip_init,pw2blip_cleanup,pw2blip_transform,pw2blip_gather
   PUBLIC blipgrid,cavc

   COMPLEX(dp),ALLOCATABLE :: psic(:),cavc_flat(:)
   INTEGER :: blipgrid(3)
   REAL(dp),ALLOCATABLE :: gamma(:)

   INTEGER,PARAMETER :: gamma_approx = 1
   REAL(dp),PARAMETER :: pi = 3.14159265358979324d0

CONTAINS

   SUBROUTINE pw2blip_init
      USE gvect, ONLY: ig1,ig2,ig3

      REAL(dp) :: da(3),k,k2,k4,cosk
      INTEGER :: ig,d,g_int(3)

      allocate(psic(dfft%nnr)) ! local FFT grid for transform
      if(ionode)allocate(cavc_flat(dfft%nr1x*dfft%nr2x*dfft%nr3x)) ! global FFT grid for gather
      blipgrid = (/ dfft%nr1,dfft%nr2,dfft%nr3 /)

! Calculating gamma.
      allocate(gamma(NGvec))
      gamma(:) = 1.d0
      da(1:3)=2.d0*pi/dble( (/ dfft%nr1,dfft%nr2,dfft%nr3 /) )
      if(gamma_approx==1)then
         do ig=1,NGvec
            g_int(:) = (/ig1(ig),ig2(ig),ig3(ig)/)
            do d=1,3
               if(g_int(d)/=0)then
                  k=da(d)*dble(g_int(d)) ; cosk=cos(k) ; k2=k*k ; k4=k2*k2
                  gamma(ig)=gamma(ig)*k4/(6.d0*((cosk-2.d0)*cosk+1.d0))
               else
                  gamma(ig)=gamma(ig)*2.d0/3.d0
               endif
            enddo
         enddo ! ig
      elseif(gamma_approx==2)then
         do ig=1,NGvec
            gamma(ig)=1.d0/((1.d0+0.5d0*cos(da(1)*g(1,ig))) &
               &*(1.d0+0.5d0*cos(da(2)*g(2,ig)))*(1.d0+0.5d0*cos(da(3)*g(3,ig))))
         enddo ! ig
      else
         write(6,*)'Bug: bad gamma_approx.' ; stop
      endif ! gamma_approx

   END SUBROUTINE pw2blip_init

   SUBROUTINE pw2blip_cleanup
      deallocate(psic,gamma)
      if(ionode)deallocate(cavc_flat)
   END SUBROUTINE pw2blip_cleanup

   SUBROUTINE pw2blip_transform(psi)
      USE wvfct, ONLY: npwx ,& ! allocation size of psi
                       npw  ,& ! number of plane waves at this k point
                       igk     ! map from plane waves in psi to k-independent indices
      USE fft_parallel, ONLY: tg_cft3s

      COMPLEX(DP), INTENT(in) :: psi(npwx)
      INTEGER :: ipw,ig

      ! reciprocal space is distributed in x-y sticks
      psic (:) = (0.d0, 0.d0)
      psic (map_igk_to_fft (igk(1:npw))) = psi(1:npw)*gamma(igk(1:npw))

      ! perform the transformation
      call tg_cft3s (psic, dfft, +2)
      ! real space is distributed in z planes
   END SUBROUTINE

   SUBROUTINE pw2blip_gather
      USE mp_global, ONLY: nproc_pool, intra_pool_comm
      USE mp, ONLY: mp_gather

      INTEGER :: proc, info
      INTEGER :: displs(0:nproc_pool-1), recvcount(0:nproc_pool-1)

      do proc = 0, ( nproc_pool - 1 )
         recvcount(proc) = dfft%nnp * dfft%npp(proc+1)
         if( proc == 0 )then
            displs(proc) = 0
         else
            displs(proc) = displs(proc-1) + recvcount(proc-1)
         endif
      enddo
      info = 0
      call mp_gather( psic, cavc_flat, &
       & recvcount, displs, ionode_id, intra_pool_comm )
   END SUBROUTINE

   COMPLEX(dp) FUNCTION cavc(i1,i2,i3)
      INTEGER,INTENT(in) :: i1,i2,i3
      cavc = cavc_flat(i1+dfft%nr1x*(i2-1+dfft%nr3x*(i3-1)))
   END FUNCTION cavc

END MODULE
