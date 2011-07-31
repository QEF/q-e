
      SUBROUTINE exx_gs(nfi, c)

!===============================================================
! modified from exact_exchange.f90 written by Zhaofeng and Xifan.
! In the original code, each processor holds one wave function. 
! Here I modify it so each processor holds nbsp/nproc_image
! number of wave functions. 
! Lingzhu Kong
!===============================================================

      USE kinds,                   ONLY  : DP
      USE fft_base,                ONLY  : dffts
      USE mp,                      ONLY  : mp_barrier, mp_sum
      USE mp_global,               ONLY  : nproc_image, me_image, root_image, intra_image_comm
      USE parallel_include
      USE io_global,               ONLY  : stdout
      USE cell_base,               ONLY  : at, omega, alat
      USE cell_base,               ONLY  : r_to_s, s_to_r, ainv, h, pbcs
      USE electrons_base,          ONLY  : nbsp, nbspx, nspin
      USE gvecw,                   ONLY  : ngw
      USE input_parameters,        ONLY  : nbnd
      USE wannier_module,          ONLY  : wfc
      USE cp_main_variables,       ONLY  : my_nbspx,selfv, pairv
      USE cp_main_variables,       ONLY  : vpsig=>exx_potential
      USE cp_main_variables,       ONLY  : n_exx
      USE cp_main_variables,       ONLY  : lmax, clm
      USE energies,                ONLY  : exx
      USE constants,               ONLY  : fpi
      USE printout_base,           ONLY  : printout_base_open, printout_base_unit, printout_base_close
      USE wannier_base,            ONLY  : neigh, dis_cutoff
      USE cp_main_variables,       ONLY  : odtothd_in_sp,  thdtood_in_sp
      USE cp_main_variables,       ONLY  : thdtood, np_in_sp
      USE mp_wave,                 ONLY  : redistwfr

      IMPLICIT NONE
      COMPLEX(DP)   c(ngw, nbspx)
#ifdef __MPI
      INTEGER  :: istatus(MPI_STATUS_SIZE)
#endif
      INTEGER     ir, i, j,nfi, ierr, nnrtot
      INTEGER     njj(nbsp),  nj_max
      INTEGER     overlap3(neigh/2,nbsp)

      REAl(DP)    sa1, wannierc(3, nbsp), middle(3,neigh/2)
      REAl(DP)    ha, hb, hc, alength(3) ,hcub
      
      INTEGER,     ALLOCATABLE ::   my_nbsp(:), my_nxyz(:)
      INTEGER,     ALLOCATABLE ::   index_my_nbsp (:, :), rk_of_obtl (:), lindex_of_obtl(:)

      REAL(DP),    pointer     ::   psi_pair(:,:)=>null(), vpsiforj(:,:)=>null()
      REAL(DP),    ALLOCATABLE ::   vpsil(:,:),vpsiforj_trcv(:)
      REAL(DP),    ALLOCATABLE ::   rho(:),rho_in_sp(:),v(:)
      REAL(DP),    ALLOCATABLE ::   psi(:,:)

      INTEGER   iobtl, gindex_of_iobtl, irank, rk_of_obtl_trcv, rk_of_obtl_tbs
      INTEGER   obtl_tbs, lindex_obtl_tbs, obtl_trcv, lindex_obtl_trcv 
      REAL(DP)  totalenergy, totalenergyg, tot_energy(nbsp)
      REAL(DP)  selfe, paire(neigh/2), tmp(3), tmp2(3)

      INTEGER, allocatable  :: isendreq(:),irecvreq(:,:)

      INTEGER   tran(3), proc, tmp_iobtl, me

!=============================================================================================
! General variables used in this subroutine
      nnrtot = dffts%nr1 * dffts%nr2 * dffts%nr3
      alength(1) = sqrt( at(1,1)**2 + at(2,1)**2 + at(3,1)**2 ) * alat 
      alength(2) = sqrt( at(1,2)**2 + at(2,2)**2 + at(3,2)**2 ) * alat
      alength(3) = sqrt( at(1,3)**2 + at(2,3)**2 + at(3,3)**2 ) * alat

      ha = alength(1) / dffts%nr1  !nr1s in the parallel case
      hb = alength(2) / dffts%nr2  !nr1s in the parallel case
      hc = alength(3) / dffts%nr3  !nr1s in the parallel case
      hcub = omega / DBLE(nnrtot) !nnrtot in parallel
      sa1 = 1.0d0/omega

      print *, 'entering exx_gs', n_exx, nfi, np_in_sp
      if(n_exx == 0)then
         call exx_setup( nnrtot, lmax, clm, fpi)
      end if

!-------------------------------------------------------------------------
! Get the Wannier center

      if (n_exx == 0) then
          do ir = 1, nbsp
             read(407, *) wannierc(1,ir), wannierc(2,ir), wannierc(3,ir)
          end do
      else
         wannierc(:, :) = wfc(:, :)
      end if

      do ir = 1, nbsp
         tmp = wannierc(:,ir)
         call r_to_s(tmp, tmp2, ainv)
!        call pbcs(tmp2, tmp,1)
         do i = 1, 3
            tmp(i) = tmp2(i) - int(tmp2(i))
            if(tmp(i) < 0)then
               tmp(i) = tmp(i) + 1
            endif
         enddo
         call s_to_r(tmp,  tmp2, h)
         wannierc(:,ir) = tmp2(:)
      end do

! initialize the output as zero
      vpsig(:, :) = 0.0d0
      totalenergy = 0.0d0

      allocate( my_nxyz ( nproc_image ) )
      allocate( my_nbsp ( nproc_image ) )

      my_nbsp(:) = nbsp/nproc_image

      IF( MOD(nbsp, nproc_image) /= 0)THEN
         DO i = 1, nproc_image
            IF( (i-1) < MOD(nbsp, nproc_image) ) my_nbsp(i) = my_nbsp(i)+1
         END DO
      ENDIF

      my_nxyz(:) = dffts%nr1x*dffts%nr2x*dffts%npp
      print *, my_nbsp
      print *, my_nxyz
     
      allocate( index_my_nbsp (my_nbspx, nproc_image) )
      allocate( rk_of_obtl ( nbsp ) )
      allocate( lindex_of_obtl( nbsp ) )

      index_my_nbsp(:, :) = nbsp + 1
      do irank = 1, nproc_image
         do iobtl = 1, my_nbsp(irank)
            gindex_of_iobtl = iobtl
            do proc = 1, irank - 1, 1
               gindex_of_iobtl = gindex_of_iobtl + my_nbsp(proc)
            enddo
            if( gindex_of_iobtl <= nbsp)then
               index_my_nbsp(iobtl, irank) = gindex_of_iobtl
            endif

         enddo
      enddo

      print *, 'index_my_nbsp = ', index_my_nbsp

      do iobtl = 1, nbsp
         rk_of_obtl(iobtl) = 0
         tmp_iobtl = iobtl
         do proc = 1, nproc_image
            tmp_iobtl = tmp_iobtl - my_nbsp(proc)
            if(tmp_iobtl <= 0)THEN
              rk_of_obtl(iobtl) = proc - 1
              print *, 'lrk_of_iobtl=', proc-1, rk_of_obtl(iobtl) 
              exit
            endif
         enddo
      enddo

      print *, 'rk_of_obtl = ', rk_of_obtl

      do iobtl = 1, nbsp
         lindex_of_obtl(iobtl) = iobtl
         do proc = 1, nproc_image
            if(lindex_of_obtl(iobtl) <= my_nbsp(proc))exit
            lindex_of_obtl(iobtl) = lindex_of_obtl(iobtl) - my_nbsp(proc)
         enddo
      enddo

      print *, 'lindex_of_obtl = ', lindex_of_obtl

      me = me_image + 1
      n_exx = n_exx + 1

!=========================================================================
!obtain orbitals on each local processor, stored in psi

      allocate ( psi(  nnrtot, my_nbsp(me ) ) )
      allocate ( vpsil(nnrtot, my_nbsp(me ) ) )
      allocate ( v(nnrtot) )
      allocate ( rho(nnrtot) )
      allocate ( rho_in_sp(np_in_sp) )

      call start_clock('r_orbital')
      call exx_psi(c, psi, nnrtot, my_nbsp, my_nxyz, nbsp) 
      call stop_clock('r_orbital')
!=========================================================================
!                                
!              SELF POTENTIAL FOR EACH ORBITAL STARTS HERE
!                                
!=========================================================================

      call start_clock('self_v')
      do iobtl = 1, my_nbsp(me)

         gindex_of_iobtl =  index_my_nbsp(iobtl, me)

         call getsftv ( dffts%nr1,dffts%nr2,dffts%nr3, ha, hb, hc, &
                        wannierc(1, gindex_of_iobtl), tran)

         call getrho( nnrtot, psi(1, iobtl), psi(1, iobtl), rho, rho_in_sp, tran, sa1)

         call getvofr( nnrtot, hcub, n_exx, rho_in_sp, v, selfv(1,1,iobtl), selfv(1,2,iobtl), tran )

         do ir = 1, nnrtot
            vpsil(ir,iobtl) = -0.25d0* v(ir)*psi(ir,iobtl) ! PBE0
         end do

         call vvprod(nnrtot, rho, v, selfe)
         selfe = selfe * 0.5d0 * hcub
         totalenergy = totalenergy + selfe
         write(6,*)iobtl, 'self energy' , selfe
      enddo ! do iobtl 

      call stop_clock('self_v')

!===============================================================================
!
!                              PAIR POTENTIAL
!===============================================================================

      call exx_index_pair(wannierc, overlap3, njj, nj_max)
 
!========================================================================
!                      THE MOST OUTER LOOP STARTS:
!========================================================================

      allocate( irecvreq(nj_max,0:nproc_image-1) )
      allocate( psi_pair(nnrtot, nj_max ) )       ! largest memory allocation
      vpsiforj => psi_pair

      do iobtl = 1, my_nbspx
 
         psi_pair(:,:)=0.d0

         call mp_barrier( intra_image_comm )
#ifdef __MPI
!========================================================================
         call start_clock('send_psi')
         do j = 1, nj_max
            do irank = 1, nproc_image

               gindex_of_iobtl =  index_my_nbsp(iobtl, irank)
               if( gindex_of_iobtl <= nbsp)then
               rk_of_obtl_trcv = irank - 1
               obtl_tbs        = overlap3(j, gindex_of_iobtl)

               if(obtl_tbs .ne. 0)then
                  rk_of_obtl_tbs  = rk_of_obtl(obtl_tbs)
                  lindex_obtl_tbs = lindex_of_obtl(obtl_tbs)

                  if( (me_image .eq. rk_of_obtl_trcv) .and. (me_image .eq. rk_of_obtl_tbs ))then
                     psi_pair(:,j) = psi(:, lindex_obtl_tbs)      ! local copy

                  elseif( me_image .eq. rk_of_obtl_tbs )then
                     CALL MPI_SEND( psi(1, lindex_obtl_tbs), nnrtot, MPI_DOUBLE_PRECISION, & 
                                    rk_of_obtl_trcv, j*irank, intra_image_comm,ierr )
                
                  elseif( me_image .eq. rk_of_obtl_trcv )then
                     CALL MPI_IRECV( psi_pair(1,j),           nnrtot, MPI_DOUBLE_PRECISION, &
                                    rk_of_obtl_tbs,  j*irank, intra_image_comm, irecvreq(j,me_image),ierr)
                  endif
               endif
               endif
            enddo  !irank
         enddo  ! j
!=======================================================================

         do j = 1, nj_max
            do irank = 1, nproc_image
               gindex_of_iobtl =  index_my_nbsp(iobtl, irank)
               if( gindex_of_iobtl <= nbsp)then
               rk_of_obtl_trcv = irank - 1
               obtl_tbs        = overlap3(j, gindex_of_iobtl)
               if(obtl_tbs .ne. 0)then
                  rk_of_obtl_tbs  = rk_of_obtl(obtl_tbs)
                  if( (me_image .eq. rk_of_obtl_trcv) .and. (me_image .ne. rk_of_obtl_tbs) )then
                     CALL MPI_WAIT(irecvreq(j,me_image), istatus, ierr)
                  endif
               endif
               endif
            enddo
         enddo

         call stop_clock('send_psi')
#endif
! after this loop ( do irank ), all the processor got all the overlapping orbitals 
! for the i_obtl orbital and ready to calculate pair potential 
!=======================================================================

         middle(:,:)=0.d0
         gindex_of_iobtl = index_my_nbsp(iobtl, me)
         print *,'gindex_of_iobtl = ', gindex_of_iobtl 
         if( gindex_of_iobtl <= nbsp)then

!second loop starts: calculate overlapping potential with the j_th orbitals

         call start_clock('getpairv')
         do j = 1, njj( gindex_of_iobtl )

            IF(overlap3(j,gindex_of_iobtl) /= 0)THEN

               call getmiddlewc(wannierc(1,gindex_of_iobtl),wannierc(1,overlap3(j,gindex_of_iobtl)), middle(1,j) )

               v(:) = 0.0d0
               call getsftv( dffts%nr1, dffts%nr2, dffts%nr3, ha, hb, hc, &
                             middle(1, j), tran)
               call getrho(  nnrtot, psi(1, iobtl), psi_pair(1, j), rho, rho_in_sp, tran,sa1)
               call getvofr( nnrtot, hcub, n_exx, rho_in_sp, v, pairv(1,1,j,iobtl), pairv(1,2,j,iobtl), tran )
               do ir = 1, nnrtot
                  vpsil(ir,iobtl) = vpsil(ir,iobtl) - 0.25d0 * v(ir) * psi_pair(ir,j)
                  vpsiforj(ir,j)  =                 - 0.25d0 * v(ir) * psi(ir,iobtl) 
               end do

               call vvprod(nnrtot, rho, v, paire(j))
               paire(j) = paire(j) * 0.5d0 * hcub
               totalenergy = totalenergy + 2.d0*paire(j)

            END IF
         enddo !for j

         print *, 'pair energy  follows '
         write(*,'(5f15.8)')(paire(j),j=1,njj( gindex_of_iobtl ))

         endif !gindex_of_iobtl <= nbsp

         call stop_clock('getpairv')
!
!===============================================================================
! After this loop, each processor finished the pair potentials for the 
! iobtl orbital, and shall talk to send/recv vpsiforj
!===============================================================================

         call mp_barrier(intra_image_comm)

         allocate(vpsiforj_trcv(nnrtot))
#ifdef __MPI
         call start_clock('sendv')
         do j = 1, nj_max
            do irank = 1, nproc_image
 
               gindex_of_iobtl = index_my_nbsp(iobtl, irank)
               if( gindex_of_iobtl <= nbsp)then

               rk_of_obtl_tbs  = irank - 1
               obtl_trcv = overlap3(j, gindex_of_iobtl)
       
               if ( obtl_trcv .ne. 0)then
                  rk_of_obtl_trcv  = rk_of_obtl(obtl_trcv)
                  lindex_obtl_trcv = lindex_of_obtl(obtl_trcv)

                  if( (me_image .eq. rk_of_obtl_trcv) .and. (me_image .eq. rk_of_obtl_tbs ))then
                     do ir = 1, nnrtot
                        vpsil(ir,lindex_obtl_trcv) = vpsil(ir,lindex_obtl_trcv) + vpsiforj(ir,j)
                     end do
          
                  elseif( me_image .eq. rk_of_obtl_tbs ) then
                     CALL MPI_SEND( vpsiforj(1,j), nnrtot, MPI_DOUBLE_PRECISION, &
                                    rk_of_obtl_trcv, j*irank, intra_image_comm,ierr) 
                  elseif( me_image .eq. rk_of_obtl_trcv)then
                     CALL MPI_RECV( vpsiforj_trcv, nnrtot, MPI_DOUBLE_PRECISION, &
                                 rk_of_obtl_tbs,  j*irank, intra_image_comm, istatus,ierr)
       
                     do ir = 1, nnrtot
                        vpsil(ir,lindex_obtl_trcv) = vpsil(ir,lindex_obtl_trcv) + vpsiforj_trcv(ir)
                     end do
                  endif
               endif
               endif  ! gindex_of_iobtl <= nbsp
       
            end do ! irank  
         end do  ! j

         call stop_clock('sendv')
#endif
         deallocate( vpsiforj_trcv )

      end do ! iobtl
      deallocate(psi_pair )

!============================================================================
               !THE MOST OUTER LOOP FOR PAIR POTENTIAL ENDS HERE
!============================================================================

#ifdef __MPI
!============================================================================
!collect and redistribure in format of vpsig(dffts%nnrx,nbsp)
      CALL MPI_ALLREDUCE(totalenergy, totalenergyg, 1, MPI_DOUBLE_PRECISION, &
&                        MPI_SUM, intra_image_comm, ierr)
#endif
      exx = totalenergyg
      if(nspin .eq. 1)exx = exx + totalenergyg
      write(stdout, '(a, f12.6)')'    EXX Energy' , exx

      call start_clock('vl2vg')
      call redistwfr( vpsig, vpsil, my_nxyz, my_nbsp, intra_image_comm, -1 )
      call stop_clock('vl2vg')

      print *, 'leaving exx_gs'
!==============================================================================
      deallocate( rho )
      deallocate( rho_in_sp )
      deallocate( v )
      deallocate( vpsil )
      deallocate( psi )
      deallocate( irecvreq )
      deallocate( my_nbsp )
      deallocate( my_nxyz )
      deallocate( index_my_nbsp )
      deallocate( rk_of_obtl )
      deallocate( lindex_of_obtl )
      return
   
      END SUBROUTINE exx_gs

      SUBROUTINE vvprod(n, v1, v2, prod)
      USE kinds, ONLY  : DP
      IMPLICIT   none

      INTEGER  n
      REAL(DP) prod, v1(n), v2(n)

      INTEGER  i

      prod = 0.0d0
      do i = 1, n
         prod = prod + v1(i) * v2(i)
      end do

      return
      end

