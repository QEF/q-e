
      SUBROUTINE exx_es(nfi, c, cv)
      !===============================================================
      ! modified from exact_exchange.f90 written by Zhaofeng and Xifan.
      ! Lingzhu Kong
      !===============================================================
      ! Note:  From this code exx_potential is returned after multiplying mixing parameter exxalfa.
      !        Later the exx_potential is added with GGA potential in forces.f90.
      !        In the future, full exx_potential should be returned and the mixing parameter exxalfa
      !        should be multiplied in forces.f90.
      !===============================================================

      USE kinds,                   ONLY  : DP
      USE mp,                      ONLY  : mp_barrier 
      USE mp_global,               ONLY  : nproc_image, me_image, root_image, intra_image_comm
      USE parallel_include
      USE io_global,               ONLY  : stdout
      USE cell_base,               ONLY  : h, omega
      USE electrons_base,          ONLY  : nbsp, nbspx, f, nspin, ispin
      USE gvecw,                   ONLY  : ngw
      USE wannier_module,          ONLY  : wfc
      USE exx_module,              ONLY  : pairv, pair_dist, vpsig=>exx_potential,n_exx, lmax, clm, vwc
      USE exx_module,              ONLY  : odtothd_in_sp,  thdtood_in_sp, thdtood, np_in_sp => np_in_sp_s !HK: to be fixed
      USE exx_module,              ONLY  : my_nbspx, my_nbsp, my_nxyz, rk_of_obtl, lindex_of_obtl, index_my_nbsp
      USE exx_module,              ONLY  : exx_setup_nscf, getnpinsp
      USE exx_module,              ONLY  : exxalfa
      USE constants,               ONLY  : fpi
      USE printout_base,           ONLY  : printout_base_open, printout_base_unit, printout_base_close
      USE wannier_base,            ONLY  : neigh, dis_cutoff, vnbsp
      USE control_flags,           ONLY  : lwfpbe0nscf
      USE fft_base,                ONLY  : dffts,dfftp
      USE mp_wave,                 ONLY  : redistwfr

      IMPLICIT NONE
      COMPLEX(DP)   c(ngw, nbspx), cv(ngw, vnbsp)
      
      INTEGER  :: sdispls(nproc_image), sendcount(nproc_image)
      INTEGER  :: rdispls(nproc_image), recvcount(nproc_image)
#if defined(__MPI)
      ! HK: sequential support
      INTEGER  :: istatus(MPI_STATUS_SIZE)
#endif

      INTEGER     ir, i, j,nfi, ierr, nnrtot,nr1s,nr2s,nr3s
      INTEGER     njj(nbsp),  nj_max
      INTEGER     overlap3(neigh,nbsp)
      REAl(DP)    wc(3, nbsp), middle(3,neigh)
      REAl(DP)    a(3),ha, hb, hc, sa1
      REAl(DP)    hcub, centerx, centery, centerz
      ! HK: pair distance
      REAL(DP)    d_pair(neigh/2)
      
      REAL(DP),    ALLOCATABLE ::   vpsil(:,:)
      REAL(DP),    ALLOCATABLE ::   rho(:),rho_in_sp(:),v(:)
      REAL(DP),    ALLOCATABLE ::   psi(:,:)
      REAL(DP),    ALLOCATABLE ::   psi_v(:,:), psi_pair(:,:)
      INTEGER,     ALLOCATABLE ::   my_vnbsp(:)
      INTEGER,     ALLOCATABLE ::   rk_of_vobtl (:), lindex_of_vobtl(:)

      INTEGER   iobtl, gindex_of_iobtl, irank, rk_of_obtl_trcv, rk_of_obtl_tbs
      INTEGER   obtl_tbs, lindex_obtl_tbs, obtl_trcv, lindex_obtl_trcv, me
      REAL(DP)  totalenergy, totalenergyg, tot_energy(nbsp)

      INTEGER, allocatable  :: irecvreq(:,:)

      INTEGER   tran(3), proc, tmp_iobtl
!=============================================================================================
! General variables used in this subroutine
      nr1s=dfftp%nr1; nr2s=dfftp%nr2; nr3s=dfftp%nr3 
     
      a(1)=DSQRT(h(1,1)*h(1,1)+h(2,1)*h(2,1)+h(3,1)*h(3,1))   ! lattice 1 
      a(2)=DSQRT(h(1,2)*h(1,2)+h(2,2)*h(2,2)+h(3,2)*h(3,2))   ! lattice 2 
      a(3)=DSQRT(h(1,3)*h(1,3)+h(2,3)*h(2,3)+h(3,3)*h(3,3))   ! lattice 3 
     
      ha = a(1) / nr1s  
      hb = a(2) / nr2s  
      hc = a(3) / nr3s  
     
      nnrtot = nr1s * nr2s * nr3s
      hcub = omega / DBLE(nnrtot) !nnrtot in parallel
     
      centerx = 0.5_DP * a(1)
      centery = 0.5_DP * a(2)
      centerz = 0.5_DP * a(3)
     
      sa1 = 1.0_DP/omega

      print *, 'entering exx_es', n_exx, nfi
      if(n_exx == 0)then
         call exx_setup_nscf( nnrtot, lmax, clm, fpi, wc, vwc, nbsp, vnbsp)
      end if

!-------------------------------------------------------------------------

      if (n_exx /= 0) then
         wc(:, :) = wfc(:, :)
         do ir = 1, nbsp
            if (wc(1, ir) < 0) then
               wc(1,ir) = wc(1,ir) + a(1)
            end if
            if (wc(2, ir) < 0) then
               wc(2,ir) = wc(2,ir) + a(2)
            end if
            if (wc(3, ir) < 0) then
               wc(3,ir) = wc(3,ir) + a(3)
            end if
         end do
      endif

! initialize the output as zero
      vpsig(:, :) = 0.0d0
      totalenergy = 0.0d0

      ALLOCATE( my_vnbsp( nproc_image ) )

      my_vnbsp(:) = vnbsp / nproc_image
      DO i = 1, nproc_image
         IF( (i-1) < MOD(vnbsp, nproc_image) )my_vnbsp(i)=my_vnbsp(i)+1
      END DO

     !print *, 'me_vnbsp = ', my_vnbsp
     !print *, 'my_nxyz  = ',  my_nxyz

      ALLOCATE( rk_of_vobtl ( vnbsp ) )
      ALLOCATE( lindex_of_vobtl( vnbsp ) )

!      print *, 'index_my_nbsp = ', index_my_nbsp

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

!      print *, 'rk_of_vobtl = ', rk_of_vobtl

      do iobtl = 1, vnbsp
         lindex_of_vobtl(iobtl) = iobtl
         do proc = 1, nproc_image
            if(lindex_of_vobtl(iobtl) <= my_vnbsp(proc))exit
            lindex_of_vobtl(iobtl) = lindex_of_vobtl(iobtl) - my_vnbsp(proc)
         enddo
      enddo

!      print *, 'lindex_of_vobtl = ', lindex_of_vobtl

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

#if defined(__MPI)
         ! HK: sequential support
         call mp_barrier( intra_image_comm )
#endif

         call start_clock('send_psi')
         do j = 1, nj_max
            do irank = 1, nproc_image

               gindex_of_iobtl = index_my_nbsp(iobtl, irank)
               if( gindex_of_iobtl > nbsp)exit

               rk_of_obtl_trcv = irank - 1
               obtl_tbs        = overlap3(j, gindex_of_iobtl)

               if(obtl_tbs .ne. 0)then
                  rk_of_obtl_tbs  = rk_of_vobtl(obtl_tbs)
                  lindex_obtl_tbs = lindex_of_vobtl(obtl_tbs)

                  if( (me_image .eq. rk_of_obtl_trcv) .and. (me_image .eq. rk_of_obtl_tbs ))then
                     psi_pair(:,j) = psi_v(:, lindex_obtl_tbs)      ! local copy

                  elseif( me_image .eq. rk_of_obtl_tbs )then
#if defined(__MPI)
                    ! HK: sequential support
                     CALL MPI_SEND( psi_v(1, lindex_obtl_tbs), nnrtot, MPI_DOUBLE_PRECISION, & 
 &                                  rk_of_obtl_trcv, j*irank, intra_image_comm,ierr )
#endif
                
                  elseif( me_image .eq. rk_of_obtl_trcv )then
#if defined(__MPI)
                    ! HK: sequential support
!                    CALL MPI_RECV( psi_pair(1,j),           nnrtot, MPI_DOUBLE_PRECISION, &
!&                                  rk_of_obtl_tbs,  j*irank, intra_image_comm, istatus,ierr)
                     CALL MPI_IRECV( psi_pair(1,j),           nnrtot, MPI_DOUBLE_PRECISION, &
                                    rk_of_obtl_tbs,  j*irank, intra_image_comm, irecvreq(j,me_image),ierr)
#endif
                  endif
               endif
            enddo  !irank
         enddo  ! j

         do j = 1, nj_max
            do irank = 1, nproc_image
               gindex_of_iobtl = index_my_nbsp(iobtl, irank)

               if( gindex_of_iobtl > nbsp)exit
               rk_of_obtl_trcv = irank - 1
               obtl_tbs        = overlap3(j, gindex_of_iobtl)

               if(obtl_tbs .ne. 0)then
                  rk_of_obtl_tbs  = rk_of_vobtl(obtl_tbs)
                  lindex_obtl_tbs = lindex_of_vobtl(obtl_tbs)

                  if( (me_image .eq. rk_of_obtl_trcv) .and. (me_image .ne. rk_of_obtl_tbs) )then
#if defined(__MPI)
                    ! HK: sequential support
                     CALL MPI_WAIT(irecvreq(j,me_image), istatus, ierr)
#endif
                  endif
               endif
            enddo
         enddo

         call stop_clock('send_psi')
!=======================================================================

         middle(:,:)=0.d0
         gindex_of_iobtl = index_my_nbsp(iobtl, me) 
         if( gindex_of_iobtl > nbsp)exit

         call start_clock('getpairv')
         do j = 1, njj( gindex_of_iobtl )

            IF(overlap3(j,gindex_of_iobtl) /= 0)THEN

               call getmiddlewc(wc(1,gindex_of_iobtl),vwc(1,overlap3(j,gindex_of_iobtl)), &
&                               centerx, centery, centerz, a(1), a(2), a(3), middle(1,j) )
               ! HK: get pair distance
               call get_pair_dist(wc(1,gindex_of_iobtl),vwc(1,overlap3(j,gindex_of_iobtl)),d_pair(j))

               v(:) = 0.0d0
               call getsftv( nr1s, nr2s, nr3s, ha, hb, hc, middle(1, j), tran)

               call start_clock('getrho_ind')
               call getrho(  nnrtot, psi(1, iobtl), psi_pair(1, j), rho, rho_in_sp, tran,sa1)

               call stop_clock('getrho_ind')

               call start_clock('getexxv')
               ! HK: modidfed quad_hybrid extrapolation
               call getvofr( nnrtot, hcub, n_exx, rho_in_sp, v, pairv(1,1,j,iobtl), pairv(1,2,j,iobtl),&
                pairv(1,3,j,iobtl), tran, .FALSE., d_pair(j), pair_dist(1,j,iobtl), pair_dist(2,j,iobtl),&
                pair_dist(3,j,iobtl))
               call stop_clock('getexxv')

               do ir = 1, nnrtot
                  vpsil(ir,iobtl) = vpsil(ir,iobtl) - exxalfa*v(ir) * psi_pair(ir,j)
               end do

            END IF
         END DO !for j
         print *,'done iobtl =', iobtl
         call stop_clock('getpairv')
      end do ! iobtl

      call start_clock('vl2vg')
      call redistwfr( vpsig, vpsil, my_nxyz, my_nbsp, intra_image_comm, -1 )
      call stop_clock('vl2vg')
      print *, 'leaving exx_es'
!==============================================================================
      deallocate(rho, rho_in_sp, v, vpsil, psi, psi_v, irecvreq, psi_pair )
      deallocate(my_vnbsp, rk_of_vobtl , lindex_of_vobtl )

      return
   
      END SUBROUTINE exx_es

!====================================================================================
      SUBROUTINE getrho( ntot, psi, psi2, rho, rho_in_sp, tran, sa1)

      USE kinds, ONLY  : DP
      USE exx_module,              ONLY  :  odtothd_in_sp,  thdtood_in_sp, thdtood
      USE exx_module,              ONLY  :  np_in_sp => np_in_sp_s !HK: to be fixed
      USE constants,               ONLY  : fpi
      USE fft_base,                ONLY  : dfftp

      IMPLICIT   none

      INTEGER  ntot, tran(3), nr1s,nr2s,nr3s
      REAL(DP) psi(ntot), psi2(ntot), rho(ntot), rho_in_sp(np_in_sp),sa1

      INTEGER  ir, ip, i, j, k, ii, jj, kk

      nr1s=dfftp%nr1; nr2s=dfftp%nr2; nr3s=dfftp%nr3 

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
         if( ii .gt. nr1s)ii = ii - nr1s
         if( jj .gt. nr2s)jj = jj - nr2s
         if( kk .gt. nr3s)kk = kk - nr3s

         if( ii .lt. 1)ii = ii + nr1s
         if( jj .lt. 1)jj = jj + nr2s
         if( kk .lt. 1)kk = kk + nr3s

         ip = thdtood(ii, jj, kk)
         rho_in_sp( ir ) = rho(ip)
      ENDDO

      RETURN
      END


