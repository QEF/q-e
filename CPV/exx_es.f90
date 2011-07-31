
      SUBROUTINE exx_es(nfi, c, cv)
!===============================================================
! modified from exact_exchange.f90 written by Zhaofeng and Xifan.
! Lingzhu Kong
!===============================================================

      USE kinds,                   ONLY  : DP
      USE mp,                      ONLY  : mp_barrier 
      USE mp_global,               ONLY  : nproc_image, me_image, root_image, intra_image_comm
      USE parallel_include
      USE io_global,               ONLY  : stdout
      USE cell_base,               ONLY  : at, alat, omega
      use cell_base,               only  : r_to_s, s_to_r, ainv, h
      USE electrons_base,          ONLY  : nbsp, nbspx, f, nspin, ispin
      USE gvecw,                   ONLY  : ngw
      USE wannier_module,          ONLY  : wfc
      USE cp_main_variables,       ONLY  : my_nbspx,pairv, vpsig=>exx_potential, n_exx, lmax, clm, vwc
      USE cp_main_variables,       ONLY  : odtothd_in_sp,  thdtood_in_sp, thdtood, np_in_sp
      USE constants,               ONLY  : fpi
      USE printout_base,           ONLY  : printout_base_open, printout_base_unit, printout_base_close
      USE wannier_base,            ONLY  : neigh, dis_cutoff, vnbsp 
      USE control_flags,           ONLY  : lwfpbe0nscf
      USE fft_base,                ONLY  : dffts
      USE mp_wave,                  ONLY : redistwfr

      IMPLICIT NONE
      COMPLEX(DP)   c(ngw, nbspx), cv(ngw, vnbsp)
      
      INTEGER  :: sdispls(nproc_image), sendcount(nproc_image)
      INTEGER  :: rdispls(nproc_image), recvcount(nproc_image)
#ifdef __MPI
      INTEGER  :: istatus(MPI_STATUS_SIZE)
#endif
      INTEGER     ir, i, j,nfi, ierr, nnrtot
      INTEGER     njj(nbsp),  nj_max
      INTEGER     overlap3(neigh,nbsp)
      REAl(DP)    wc(3, nbsp), middle(3,neigh), alength(3)
      REAl(DP)    ha, hb, hc, sa1, hcub
      
      REAL(DP),    ALLOCATABLE ::   vpsil(:,:)
      REAL(DP),    ALLOCATABLE ::   rho(:),rho_in_sp(:),v(:)
      REAL(DP),    ALLOCATABLE ::   psi(:,:)
      REAL(DP),    ALLOCATABLE ::   psi_v(:,:), psi_pair(:,:)
      INTEGER,     ALLOCATABLE ::   my_nbsp(:), my_vnbsp(:), my_nxyz(:)
      INTEGER,     ALLOCATABLE ::   index_my_nbsp (:, :), rk_of_vobtl (:), lindex_of_vobtl(:)

      INTEGER   iobtl, gindex_of_iobtl, irank, rk_of_obtl_trcv, rk_of_obtl_tbs
      INTEGER   obtl_tbs, lindex_obtl_tbs, obtl_trcv, lindex_obtl_trcv, me
      REAL(DP)  totalenergy, totalenergyg, tot_energy(nbsp), tmp(3), tmp2(3)

      INTEGER, allocatable  :: irecvreq(:,:)

      INTEGER   tran(3), proc, tmp_iobtl
!=============================================================================================
! General variables used in this subroutine
      nnrtot = dffts%nr1 * dffts%nr2 * dffts%nr3
      alength(1) = sqrt( at(1,1)**2 + at(2,1)**2 + at(3,1)**2 ) * alat
      alength(2) = sqrt( at(1,2)**2 + at(2,2)**2 + at(3,2)**2 ) * alat
      alength(3) = sqrt( at(1,3)**2 + at(2,3)**2 + at(3,3)**2 ) * alat
      ha = alength(1) / dffts%nr1  !nr1s in the parallel case
      hb = alength(2) / dffts%nr2  !nr2s in the parallel case
      hc = alength(3) / dffts%nr3  !nr3s in the parallel case
      hcub = omega / DBLE(nnrtot) !nnrtot in parallel
      sa1 = 1.0d0/omega

      print *, 'entering exx_es', n_exx, nfi
      if(n_exx == 0)then
         call exx_setup_nscf( nnrtot, lmax, clm, fpi, wc, vwc, nbsp, vnbsp)
      end if

!-------------------------------------------------------------------------

      if (n_exx /= 0) then
         wc(:, :) = wfc(:, :)
         do ir = 1, nbsp
            tmp = wc(:,ir)
            call r_to_s(tmp,  tmp2, ainv)
!           call pbcs(tmp2, tmp,1)
            do i = 1, 3
               tmp(i) = tmp2(i) - int(tmp2(i))
               if(tmp(i) < 0)then
                  tmp(i) = tmp(i) + 1
               endif
            enddo

            call s_to_r(tmp,  tmp2, h)
            wc(:,ir) = tmp2(:)
         end do
      endif

! initialize the output as zero
      vpsig(:, :) = 0.0d0
      totalenergy = 0.0d0

      ALLOCATE( my_nxyz ( nproc_image ) )
      ALLOCATE( my_nbsp ( nproc_image ) )
      ALLOCATE( my_vnbsp( nproc_image ) )

      my_nbsp(:) = nbsp / nproc_image

      IF( MOD(nbsp, nproc_image) /= 0)THEN
         DO i = 1, nproc_image
            IF( (i-1) < MOD(nbsp, nproc_image) ) my_nbsp(i) = my_nbsp(i)+1
         END DO
      ENDIF

      my_vnbsp(:) = vnbsp / nproc_image
      DO i = 1, nproc_image
         IF( (i-1) < MOD(vnbsp, nproc_image) )my_vnbsp(i)=my_vnbsp(i)+1
      END DO

      my_nxyz(:) = dffts%nr1x*dffts%nr2x*dffts%npp

!     print *, 'me_nbsp  = ',  my_nbsp
!     print *, 'me_vnbsp = ', my_vnbsp
!     print *, 'my_nxyz  = ',  my_nxyz

      ALLOCATE( index_my_nbsp (my_nbspx, nproc_image) )
      ALLOCATE( rk_of_vobtl ( vnbsp ) )
      ALLOCATE( lindex_of_vobtl( vnbsp ) )

      index_my_nbsp(:, :) = nbsp + 1
      do irank = 1, nproc_image
         do iobtl = 1, my_nbsp(irank)
            gindex_of_iobtl = iobtl
            do proc = 1, irank - 1, 1
               gindex_of_iobtl = gindex_of_iobtl + my_nbsp(proc)
            enddo
            if( gindex_of_iobtl > nbsp)exit
            index_my_nbsp(iobtl, irank) = gindex_of_iobtl
         enddo
      enddo
        
      do iobtl = 1, vnbsp 
         rk_of_vobtl(iobtl) = 0
         tmp_iobtl = iobtl
         do proc = 1, nproc_image
            tmp_iobtl = tmp_iobtl - my_vnbsp(proc)
            if(tmp_iobtl <= 0)THEN
              rk_of_vobtl(iobtl) = proc - 1
              exit
            endif
         enddo
      enddo

      do iobtl = 1, vnbsp
         lindex_of_vobtl(iobtl) = iobtl
         do proc = 1, nproc_image
            if(lindex_of_vobtl(iobtl) <= my_vnbsp(proc))exit
            lindex_of_vobtl(iobtl) = lindex_of_vobtl(iobtl) - my_vnbsp(proc)
         enddo
      enddo

      n_exx = n_exx + 1
      me = me_image + 1
!=========================================================================

      allocate ( psi_v(nnrtot, my_vnbsp(me) ) )
      allocate ( psi(  nnrtot, my_nbsp(me ) ) )
      allocate ( vpsil(nnrtot, my_nbsp(me ) ) )
      allocate ( v(nnrtot) )
      allocate ( rho(nnrtot), rho_in_sp(np_in_sp) )

      call start_clock('r_orbital')
      call exx_psi(cv, psi_v, nnrtot, my_vnbsp, my_nxyz, vnbsp) 
      call exx_psi(c,  psi,   nnrtot, my_nbsp , my_nxyz,  nbsp) 
      call stop_clock('r_orbital')
!                                
!=========================================================================
!                              PAIR POTENTIAL
!=========================================================================

      call exx_index_pair_nv(wc, overlap3, njj, nj_max)
      print *, 'nj_max =', nj_max
      allocate( irecvreq(nj_max,0:nproc_image-1) )
      allocate( psi_pair(nnrtot, nj_max ) , stat=ierr )
      if(ierr /= 0)write(*,*)"allocation error for psi_pair"

      vpsil(:,:) = 0.d0
      do iobtl = 1, my_nbspx
 
         print *, 'iobtl =', iobtl
         psi_pair(:,:)=0.d0

         call mp_barrier( intra_image_comm )
#ifdef __MPI
         call start_clock('send_psi')
         do j = 1, nj_max
            do irank = 1, nproc_image

               gindex_of_iobtl = index_my_nbsp(iobtl, irank)
               if( gindex_of_iobtl <= nbsp)then

               rk_of_obtl_trcv = irank - 1
               obtl_tbs        = overlap3(j, gindex_of_iobtl)

               if(obtl_tbs .ne. 0)then
                  rk_of_obtl_tbs  = rk_of_vobtl(obtl_tbs)
                  lindex_obtl_tbs = lindex_of_vobtl(obtl_tbs)

                  if( (me_image .eq. rk_of_obtl_trcv) .and. (me_image .eq. rk_of_obtl_tbs ))then
                     psi_pair(:,j) = psi_v(:, lindex_obtl_tbs)      ! local copy

                  elseif( me_image .eq. rk_of_obtl_tbs )then
                     CALL MPI_SEND( psi_v(1, lindex_obtl_tbs), nnrtot, MPI_DOUBLE_PRECISION, & 
 &                                  rk_of_obtl_trcv, j*irank, intra_image_comm,ierr )
                
                  elseif( me_image .eq. rk_of_obtl_trcv )then
!                    CALL MPI_RECV( psi_pair(1,j),           nnrtot, MPI_DOUBLE_PRECISION, &
!&                                  rk_of_obtl_tbs,  j*irank, intra_image_comm, istatus,ierr)
                     CALL MPI_IRECV( psi_pair(1,j),           nnrtot, MPI_DOUBLE_PRECISION, &
                                    rk_of_obtl_tbs,  j*irank, intra_image_comm, irecvreq(j,me_image),ierr)
                  endif
               endif
               endif
            enddo  !irank
         enddo  ! j

         do j = 1, nj_max
            do irank = 1, nproc_image
               gindex_of_iobtl = index_my_nbsp(iobtl, irank)

               if( gindex_of_iobtl <= nbsp)then
               rk_of_obtl_trcv = irank - 1
               obtl_tbs        = overlap3(j, gindex_of_iobtl)

               if(obtl_tbs .ne. 0)then
                  rk_of_obtl_tbs  = rk_of_vobtl(obtl_tbs)
                  lindex_obtl_tbs = lindex_of_vobtl(obtl_tbs)

                  if( (me_image .eq. rk_of_obtl_trcv) .and. (me_image .ne. rk_of_obtl_tbs) )then
                     CALL MPI_WAIT(irecvreq(j,me_image), istatus, ierr)
                  endif
               endif
               endif
            enddo
         enddo
         call stop_clock('send_psi')
#endif
!=======================================================================

         middle(:,:)=0.d0
         gindex_of_iobtl = index_my_nbsp(iobtl, me) 
         if( gindex_of_iobtl <= nbsp)then

         call start_clock('getpairv')
         do j = 1, njj( gindex_of_iobtl )

            IF(overlap3(j,gindex_of_iobtl) /= 0)THEN

               call getmiddlewc(wc(1,gindex_of_iobtl),vwc(1,overlap3(j,gindex_of_iobtl)), middle(1,j) )

               v(:) = 0.0d0
               call getsftv( dffts%nr1, dffts%nr2, dffts%nr3, ha, hb, hc, &
                             middle(1, j), tran)

               call start_clock('getrho_ind')
               call getrho(  nnrtot, psi(1, iobtl), psi_pair(1, j), rho, rho_in_sp, tran,sa1)

               call stop_clock('getrho_ind')

               call start_clock('getexxv')
               call getvofr( nnrtot, hcub, n_exx, rho_in_sp, v, pairv(1,1,j,iobtl), pairv(1,2,j,iobtl), tran )
               call stop_clock('getexxv')

               do ir = 1, nnrtot
                  vpsil(ir,iobtl) = vpsil(ir,iobtl) - 0.25D0*v(ir) * psi_pair(ir,j)
               end do

            END IF
         END DO !for j
         endif
         print *,'done iobtl =', iobtl
         call stop_clock('getpairv')
      end do ! iobtl

      call start_clock('vl2vg')
      call redistwfr( vpsig, vpsil, my_nxyz, my_nbsp, intra_image_comm, -1 )
      call stop_clock('vl2vg')
      print *, 'leaving exx_es'
!==============================================================================
      deallocate(rho, rho_in_sp, v, vpsil, psi, psi_v, irecvreq, psi_pair )
      deallocate(my_nbsp, my_vnbsp, my_nxyz, index_my_nbsp,rk_of_vobtl , lindex_of_vobtl )

      return
   
      END SUBROUTINE exx_es

      SUBROUTINE getsftv(nr1s, nr2s, nr3s, ha, hb, hc, wc, tran)
      USE kinds,     ONLY  : DP
      use cell_base, ONLY  : r_to_s, ainv, at, alat
      IMPLICIT   none

      INTEGER  nr1s, nr2s, nr3s, tran(3)
      REAL(DP) wc(3), ha, hb, hc, wclat(3), alength(3)

      INTEGER  i, bcm(3), wcm(3)

      alength(1) = sqrt( at(1,1)**2 + at(2,1)**2 + at(3,1)**2 ) * alat
      alength(2) = sqrt( at(1,2)**2 + at(2,2)**2 + at(3,2)**2 ) * alat
      alength(3) = sqrt( at(1,3)**2 + at(2,3)**2 + at(3,3)**2 ) * alat

      ! convert to lattice coordinates
      call r_to_s(wc, wclat, ainv)

      bcm(1) = INT( nr1s/2)
      wcm(1) = INT( wclat(1)*alength(1)/ha  ) + 1

      bcm(2) = INT( nr2s/2)
      wcm(2) = INT( wclat(2)*alength(2)/hb  ) + 1

      bcm(3) = INT(  nr3s/2 )
      wcm(3) = INT( wclat(3)*alength(3)/hc  ) + 1

      DO i = 1, 3
         tran(i) = bcm(i) - wcm(i)
      ENDDO

      RETURN
      END

!====================================================================================
      SUBROUTINE getrho( ntot, psi, psi2, rho, rho_in_sp, tran, sa1)

      USE kinds, ONLY  : DP
      USE fft_base,                ONLY  : dffts
      USE cp_main_variables,       ONLY  :  odtothd_in_sp,  thdtood_in_sp, thdtood
      USE cp_main_variables,       ONLY  :  np_in_sp 
      USE constants,               ONLY  : fpi

      IMPLICIT   none

      INTEGER  ntot, tran(3)
      REAL(DP) psi(ntot), psi2(ntot), rho(ntot), rho_in_sp(np_in_sp),sa1

      INTEGER  ir, ip, i, j, k, ii, jj, kk

      DO ip = 1, ntot
         rho(ip) = psi(ip) * psi2(ip) * sa1
      ENDDO

      DO ir = 1, np_in_sp
         i  = odtothd_in_sp(1, ir)
         j  = odtothd_in_sp(2, ir)
         k  = odtothd_in_sp(3, ir)

         ii = i - tran(1)
         jj = j - tran(2)
         kk = k - tran(3)

! get ii between (1, n)
         if( ii .gt. dffts%nr1) ii = ii - int(ii/dffts%nr1)*dffts%nr1
         if( jj .gt. dffts%nr2) jj = jj - int(jj/dffts%nr2)*dffts%nr2
         if( kk .gt. dffts%nr3) kk = kk - int(kk/dffts%nr3)*dffts%nr3

         if( ii .lt. 1) ii = ii + dffts%nr1 - int( ii/dffts%nr1 )*dffts%nr1
         if( jj .lt. 1) jj = jj + dffts%nr2 - int( jj/dffts%nr2 )*dffts%nr2
         if( kk .lt. 1) kk = kk + dffts%nr3 - int( kk/dffts%nr3 )*dffts%nr3

         ip = thdtood(ii, jj, kk)
         rho_in_sp( ir ) = rho(ip)
      ENDDO

      RETURN
      END


      !=========================================================================
      subroutine getmiddlewc(wc1, wc2, mid)
      USE kinds,      ONLY  : DP
      USE cell_base,  ONLY  : r_to_s, s_to_r, h, ainv

      IMPLICIT none

      real(DP)  wc1(3), wc2(3), mid(3)
      real(DP)  diff(3), diffs(3), mids(3)

      integer i

      do i = 1, 3
         mid(i)  = wc1(i) + wc2(i)
         diff(i) = wc1(i) - wc2(i)
      enddo

      call r_to_s(diff, diffs, ainv)
      call r_to_s(mid , mids , ainv)

      do i = 1, 3
         mids (i) = 0.5d0* (mids(i) - ABS(ANINT(diffs(i))))
      enddo

      call s_to_r(mids, mid, h)

      return
      end
