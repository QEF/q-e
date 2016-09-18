!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

!routines for rapid evaluation of Fock operators

MODULE exchange_custom

  USE kinds, ONLY: DP
  USE fft_custom_gwl

  IMPLICIT NONE

  TYPE exchange_cus
!data structure for general Fock operator
     REAL(kind=DP), DIMENSION(:,:,:), POINTER :: wfc!valence wfcs in real space
     REAL(kind=DP) :: dual!dual for defining the grid on real space
     REAL(kind=DP) :: cutoff!G space cutoff in Ryberg
     REAL(kind=DP) :: truncation_radius!for Coulomb interaction in Bohr
     REAL(kind=DP), DIMENSION(:), POINTER  :: fac!factors on G space of Coulomb interaction
     INTEGER :: nbndv(2)!number of valence states
     TYPE(fft_cus) :: fft_g2r!from wfcs to real space
     TYPE(fft_cus) :: fft_r2g!from real space to restricted G space
     TYPE(fft_cus) :: fft_small!for periodic calculations
!the following for small 
     REAL(kind=DP), DIMENSION(:), POINTER  :: fac_s
     INTEGER, DIMENSION(:,:,:), POINTER :: r2s_xy
     INTEGER, DIMENSION(:,:,:), POINTER :: s2r_xy
     INTEGER :: n(3),m(3)!I have the relation n*edge=diameter*m
     LOGICAL :: l_localized!if true consider valence wfcs as localized 
     REAL(kind=DP) :: thrs_loc!threshold for localized valence wfcs
     INTEGER, DIMENSION (:,:), POINTER :: n_loc!number of points above thrs
     INTEGER, DIMENSION (:,:,:), POINTER :: tab_loc!table for points above threshold
     INTEGER :: nspin!number of spin channels
     REAL(kind=DP), DIMENSION(:,:,:,:,:),POINTER :: wfc_red!valence wfcs in real space  ready for the reduced grid
  END TYPE exchange_cus

  SAVE

!to be use in pw.x

  LOGICAL :: l_exchange_fast=.false.
  REAL(kind=DP) :: exchange_fast_dual=2.d0
  REAL(kind=DP) :: exchange_fast_cutoff=40.d0
  REAL(kind=DP) :: exchange_fast_radius=10.d0
  INTEGER :: exchange_fast_nbndv(2)
  LOGICAL :: l_exchange_turbo=.false.
  INTEGER :: exchange_m(3)
  INTEGER :: exchange_n(3)
  LOGICAL :: l_exchange_localized=.false.
  REAL(kind=DP) :: exchange_thrs_loc
  CONTAINS

      
   SUBROUTINE periodic_fock_cus(ispin,psi,xpsi,exx_cus)
!apply Fock operator to a wavefunction                 
!experimental version work just with factor 1/2

      USE io_global, ONLY : stdout, ionode,ionode_id
      USE mp_global, ONLY : me_pool,intra_pool_comm
      USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2,bg
      USE constants, ONLY : e2, pi, tpi, fpi, RYTOEV
      USE wavefunctions_module, ONLY : psic
      USE mp, ONLY : mp_sum
      USE mp_world, ONLY : world_comm, nproc
      USE wvfct,    ONLY : npwx, npw, wg
      USE gvect
      USE mp_wave, ONLY : mergewf,splitwf

      implicit none
      
      INTEGER, INTENT(in) :: ispin! spin channel 
      COMPLEX(kind=DP), INTENT(in) :: psi(npw)
      COMPLEX(kind=DP), INTENT(inout) :: xpsi(npw)

      TYPE(exchange_cus), INTENT(in) :: exx_cus

      INTEGER, ALLOCATABLE :: r2s_xy(:) !real to small XY index
    
      
      INTEGER :: i,j,k,n,ii,jj,kk,ig,iv,jv
      INTEGER :: ix,iy,iz, ix_s,iy_s,iz_s,iz_eff
      INTEGER :: rz_start,rz_end,iqq,iqq_s,rz_start_s,rz_end_s
      REAL(kind=DP), ALLOCATABLE :: fac(:),prodr(:),prods(:,:),planes(:),vexc(:)
      REAL(kind=DP) :: qq_fact, sca
      INTEGER :: iorig, idest
      INTEGER,ALLOCATABLE :: z2proc_s(:),z2proc(:)
      INTEGER :: req, ierr
#if defined(__MPI)
      INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
      COMPLEX(kind=DP), ALLOCATABLE :: prods_g(:,:)
      COMPLEX(kind=DP), ALLOCATABLE :: psi_t(:),evc_g(:),vexc_g(:)
      REAL(kind=DP), ALLOCATABLE :: psi_r(:),psi_r_red(:,:,:),plane(:)

      INTEGER :: i_mod, j_mod,k_mod,iplane
      INTEGER :: jd,jdmax,ir,nr3small,nr3small_max,nplane
      REAL(kind=DP), ALLOCATABLE :: prod_r_red(:,:,:)
      REAL(kind=DP), ALLOCATABLE :: vexc_red(:,:,:)
      INTEGER :: ip_todo,token,ip_delta,tag,ip
      REAL(kind=DP), ALLOCATABLE :: b_plane(:,:)
      INTEGER, ALLOCATABLE :: b_iplane(:),b_z(:)
      INTEGER, ALLOCATABLE :: proc_list(:)
      INTEGER :: offset

      !write(stdout,*) 'periodic_fock'
      !FLUSH(stdout)
    
      CALL start_clock('periodic_fock')
   
      !setup correspondence grids
  
      rz_start=0
      rz_end =0
      do ii=1,me_pool + 1
         rz_start=rz_end+1
         rz_end=rz_end+exx_cus%fft_g2r%dfftt%npp(ii)
      end do

      rz_start_s=0
      rz_end_s=0
      do ii=1,me_pool + 1
         rz_start_s=rz_end_s+1
         rz_end_s=rz_end_s+exx_cus%fft_small%dfftt%npp(ii)
      end do
      
      nr3small=rz_end_s-rz_start_s+1

      allocate(z2proc_s(exx_cus%fft_small%nr3t))
      allocate(z2proc(exx_cus%fft_g2r%nr3t))
      allocate(vexc(exx_cus%fft_g2r%nrxxt))
      allocate( evc_g( exx_cus%fft_g2r%ngmt_g ) )
      allocate(vexc_g(exx_cus%fft_g2r%npwt))
      j=0
      k=0
      do ii=1,nproc
         j=k+1
         k=k+exx_cus%fft_small%dfftt%npp(ii)
         z2proc_s(j:k)=ii-1
      end do

      j=0
      k=0
      do ii=1,nproc
         j=k+1
         k=k+exx_cus%fft_g2r%dfftt%npp(ii)
         z2proc(j:k)=ii-1
      end do



   
 
      allocate(fac(exx_cus%fft_small%ngmt))

      !setup fac
      do ig=1,exx_cus%fft_small%ngmt

         qq_fact = exx_cus%fft_small%gt(1,ig)**2.d0 + exx_cus%fft_small%gt(2,ig)**2.d0 + exx_cus%fft_small%gt(3,ig)**2.d0

         if (qq_fact > 1.d-8) then
            fac(ig)=(e2*fpi/(exx_cus%fft_small%tpiba2_t*qq_fact))*(1.d0-dcos(dsqrt(qq_fact)*&
                 &(exx_cus%truncation_radius*dble(exx_cus%n(1))/dble(exx_cus%m(1)))*exx_cus%fft_small%tpiba_t))
         else
            fac(ig)=e2*fpi*((exx_cus%truncation_radius*dble(exx_cus%n(1))/dble(exx_cus%m(1)))**2.d0/2.d0)

         endif

      end do
      fac(:)=fac(:)/omega/(dble(exx_cus%n(1)*exx_cus%n(2)*exx_cus%n(3)))



      allocate(prodr(exx_cus%fft_g2r%nrxxt))
      allocate(prods(exx_cus%fft_small%nrxxt,2))
      allocate(planes(exx_cus%fft_small%nrx1t*exx_cus%fft_small%nrx2t))
      allocate(prods_g(exx_cus%fft_small%ngmt,2))
      allocate(plane(exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t))


      allocate(psi_t(exx_cus%fft_g2r%npwt))
      allocate(psi_r(exx_cus%fft_g2r%nrxxt))
      allocate(psi_r_red(exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t,nr3small,exx_cus%m(3)))
      allocate(prod_r_red(exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t,nr3small,exx_cus%m(3)))
      allocate(vexc_red(exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t,nr3small,exx_cus%m(3)))


      !loop on KS states
     
      vexc(:)=0.d0
      vexc_red=0.d0

      CALL start_clock('pf_mergesplit')
      if(exx_cus%fft_g2r%dual_t==4.d0) then
         psi_t(1:exx_cus%fft_g2r%npwt)=psi(1:exx_cus%fft_g2r%npwt)
      else
         call mergewf(psi(:),evc_g,npw,ig_l2g,me_pool,nproc,ionode_id,intra_pool_comm)
         call splitwf(psi_t(:),evc_g,exx_cus%fft_g2r%npwt,exx_cus%fft_g2r%ig_l2gt,&
              &me_pool,nproc,ionode_id,intra_pool_comm)
      endif
      CALL stop_clock('pf_mergesplit')
   
      psic(:)=(0.d0,0.d0)
      psic(exx_cus%fft_g2r%nlt(1:exx_cus%fft_g2r%npwt))  = psi_t(1:exx_cus%fft_g2r%npwt)
      psic(exx_cus%fft_g2r%nltm(1:exx_cus%fft_g2r%npwt)) = CONJG( psi_t(1:exx_cus%fft_g2r%npwt) )

      CALL start_clock('pf_fftext')
      CALL cft3t( exx_cus%fft_g2r, psic, exx_cus%fft_g2r%nr1t, exx_cus%fft_g2r%nr2t, exx_cus%fft_g2r%nr3t, &
           &exx_cus%fft_g2r%nrx1t, exx_cus%fft_g2r%nrx2t, exx_cus%fft_g2r%nrx3t, 2 )
      psi_r(1:exx_cus%fft_g2r%nrxxt)= DBLE(psic(1:exx_cus%fft_g2r%nrxxt))
      CALL stop_clock('pf_fftext')






!put the psi wavefunction already on the small z distribution among processor
!!!!!!!!!!!
!first the internal case
      do k=1,exx_cus%n(3)
         do iz=rz_start,rz_end
            iz_s=mod(iz-1+(k-1)*exx_cus%fft_g2r%nr3t,exx_cus%fft_small%nr3t)+1
            iplane=(iz-1+(k-1)*exx_cus%fft_g2r%nr3t)/exx_cus%fft_small%nr3t+1
            idest=z2proc_s(iz_s)
               !put plane on small plane
            if(me_pool==idest) then
               do iqq=1,exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t
                  psi_r_red(iqq,iz_s-rz_start_s+1,iplane)=&
                       &psi_r((iz-rz_start)*(exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t)+iqq)
               enddo
            endif

         enddo
      enddo

      nr3small_max=nr3small
#if defined(__MPI)
      CALL MPI_ALLREDUCE( nr3small, nr3small_max,1,MPI_INTEGER, MPI_MAX,intra_pool_comm, req,IERR )           
#endif
      
      allocate(b_plane(exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t,nr3small_max))
      allocate(b_iplane(nr3small_max),b_z(nr3small_max))
      allocate(proc_list(nproc))
!loop on task delta
      do ip_delta=1,nproc-1

         if(mod(ip_delta,2)==0) then
              !if(mod(me_pool+1,2)==0) then!even
              !   if(mod((me_pool+1)/2,2)==0) then
              !      token=0
              !   else
              !      token=1
              !   endif
 
            !  else
            !     
            !     if(mod((me_pool+2)/2,2)==0) then
            !        token=0
            !     else
            !        token=1
            !     endif
            !
            !  endif
            proc_list=0
            do ip=1,nproc
               if(proc_list(ip)==0) then
                  if(proc_list(mod(ip+ip_delta-1,nproc)+1)==0) then
                     proc_list(ip)=-1
                     proc_list(mod(ip+ip_delta-1,nproc)+1)=1
                  else
                  
                  endif
               endif
            enddo

            if(proc_list(me_pool+1) ==-1) then
               token=0
            else
               token=1
            endif
           else
              if(mod(me_pool+1,2)==0) then
                 token=0
              else
                 token=1
              endif
             
           endif
           do ip_todo=1,2
              if(mod(ip_todo+token,2)==0) then
!if I am a sender
  !loop on my data to see if and what I have to send
                 nplane=0
                 do k=1,exx_cus%n(3)
                    do iz=rz_start,rz_end
                       iz_s=mod(iz-1+(k-1)*exx_cus%fft_g2r%nr3t,exx_cus%fft_small%nr3t)+1
                       iplane=(iz-1+(k-1)*exx_cus%fft_g2r%nr3t)/exx_cus%fft_small%nr3t+1
                       idest=z2proc_s(iz_s)
                       if(idest==mod(me_pool+ip_delta,nproc)) then
                          nplane=nplane+1
                          do iqq=1,exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t            
                             b_plane(iqq,nplane)=&
                                  &psi_r((iz-rz_start)*(exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t)+iqq)                      
                          enddo
                          b_z(nplane)=iz_s
                          b_iplane(nplane)=iplane

                       endif
                    enddo
                 enddo
!send nplane
#if defined(__MPI)
                 idest=mod(me_pool+ip_delta,nproc)
                 CALL MPI_ISEND( nplane,1, MPI_INTEGER, idest, 0, intra_pool_comm, req,IERR ) 
                 CALL MPI_WAIT(req,istatus,ierr)

                 if(nplane>0) then
                    CALL MPI_ISEND( b_plane,exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t*nplane, MPI_DOUBLE_PRECISION, &   
                         &idest, 1, intra_pool_comm, req,IERR )
                    CALL MPI_WAIT(req,istatus,ierr)
                    CALL MPI_ISEND( b_z,nplane, MPI_INTEGER,idest, 2, intra_pool_comm, req,IERR )
                    CALL MPI_WAIT(req,istatus,ierr)
                    CALL MPI_ISEND( b_iplane,nplane, MPI_INTEGER,idest, 3, intra_pool_comm, req,IERR )
                    CALL MPI_WAIT(req,istatus,ierr)
                endif
#endif
          
              else
!if I am receiver 
  !see if and what I have to receive
#if defined(__MPI)
                 iorig=me_pool-ip_delta
                 if(iorig<0) iorig=iorig+nproc
                 CALL MPI_RECV( nplane,1, MPI_INTEGER, iorig, 0, intra_pool_comm, istatus,IERR )
                
                 if(nplane>0) then
                    CALL MPI_RECV( b_plane,exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t*nplane, MPI_DOUBLE_PRECISION, &
                         &iorig, 1, intra_pool_comm, istatus,IERR )
                    CALL MPI_RECV( b_z,nplane, MPI_INTEGER,iorig, 2, intra_pool_comm, istatus,IERR )
                    CALL MPI_RECV( b_iplane,nplane, MPI_INTEGER,iorig, 3, intra_pool_comm, istatus,IERR )
                    


                    do ii=1,nplane
                       do iqq=1,exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t
                          psi_r_red(iqq,b_z(ii)-rz_start_s+1,b_iplane(ii))=b_plane(iqq,ii)
                       enddo
                    enddo
                 endif
#endif
              endif
           enddo
        enddo
        deallocate(b_plane,b_iplane,b_z,proc_list)

!!!!!!!!!
!      do k=1,exx_cus%n(3)
!         do iz=1,exx_cus%fft_g2r%nr3t
!            if(iz >= rz_start .and. iz <= rz_end) then
!               !if Z is mine determine owner and send it                                                                      
!               iz_s=mod(iz-1+(k-1)*exx_cus%fft_g2r%nr3t,exx_cus%fft_small%nr3t)+1
!               iplane=(iz-1+(k-1)*exx_cus%fft_g2r%nr3t)/exx_cus%fft_small%nr3t+1
!               idest=z2proc_s(iz_s)
!               !put plane on small plane                                                                                         
!               if(me_pool==idest) then
!                  do iqq=1,exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t
!                     psi_r_red(iqq,iz_s-rz_start_s+1,iplane)=&
!                          &psi_r((iz-rz_start)*(exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t)+iqq)
!                  enddo
!               else
!                  do iqq=1,exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t
!                     plane(iqq)=psi_r((iz-rz_start)*(exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t)+iqq)
!                  enddo
!                  CALL MPI_ISEND( plane,exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t, MPI_DOUBLE_PRECISION, &
!                       &idest, iz, intra_pool_comm, req,IERR )
!                  CALL MPI_WAIT(req,istatus,ierr)
!               
!               endif
!
!            else
!               iz_s=mod(iz-1+(k-1)*exx_cus%fft_g2r%nr3t,exx_cus%fft_small%nr3t)+1
!               iplane=(iz-1+(k-1)*exx_cus%fft_g2r%nr3t)/exx_cus%fft_small%nr3t+1
!            !if Z o  small cell is mine receive it                                                   
!               if(z2proc_s(iz_s)==me_pool) then
!                  iorig=z2proc(iz)
!                  CALL MPI_RECV( plane, exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t, MPI_DOUBLE_PRECISION, &
!                             &iorig, iz, intra_pool_comm, istatus, IERR )
!                  do iqq=1,exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t
!                     psi_r_red(iqq,iz_s-rz_start_s+1,iplane)=plane(iqq)
!                  enddo
!                  
!                  
!
!               endif
!            endif
!         enddo
!      enddo



      !loop on KS valence states
      
      CALL start_clock('pf_inner')
      do k=1,exx_cus%n(3)
         do j=1,exx_cus%n(2)
            do i=1,exx_cus%n(1)
               do jv=1,exx_cus%nbndv(ispin),2!loop on bands
      !do product
                  if(jv<exx_cus%nbndv(ispin)) then
                     jdmax=1
                  else
                     jdmax=0
                  endif
                  prods=0.d0
                  do jd=0,jdmax
!!!!!!!!!!!!!!!
                     CALL start_clock('pf_product')
                     if(.not.exx_cus%l_localized) then
                        do ii=1,exx_cus%m(3)
                           do iz=1,nr3small
                              do iqq=1,exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t
                                 prod_r_red(iqq,iz,ii)=psi_r_red(iqq,iz,ii)*exx_cus%wfc_red(iqq,iz,ii,jv+jd,ispin)
                              enddo
                           enddo
                        enddo
                     else
                        write(stdout,*) 'to be implmented yet'
                        FLUSH(stdout)
                        stop
                     endif
                     CALL stop_clock('pf_product')


                     CALL start_clock('pf_plane')
                     do iz_s=1,nr3small
                        planes(:)=0.d0
                        do ii=1,exx_cus%m(3)
                           !do iy=1,exx_cus%fft_g2r%nr2t
                            !  do ix=1,exx_cus%fft_g2r%nr1t
                               !  iqq=(iy-1)*exx_cus%fft_g2r%nrx1t+ix!ATTENZIONE TO BE VERIFY
                           do iqq=1,exx_cus%fft_g2r%nr1t*exx_cus%fft_g2r%nr2t
                              planes(exx_cus%r2s_xy(iqq,i,j))=planes(exx_cus%r2s_xy(iqq,i,j))+&
                             &prod_r_red(iqq,iz_s,ii)
                           enddo
                             ! enddo
                           !enddo
                        enddo
                        offset=(iz_s-1)*(exx_cus%fft_small%nrx1t*exx_cus%fft_small%nrx2t)
                        do iqq_s=1,exx_cus%fft_small%nrx1t*exx_cus%fft_small%nrx2t
                           prods(offset+iqq_s,jd+1)=&
                   &   planes(iqq_s)
                        enddo


                     enddo
                     CALL stop_clock('pf_plane')

 

                  enddo!on jd
              
      ! do fft
                  CALL start_clock('pf_cache')
                  psic(1:exx_cus%fft_small%nrxxt)=cmplx(prods(1:exx_cus%fft_small%nrxxt,1),prods(1:exx_cus%fft_small%nrxxt,2))
                  CALL stop_clock('pf_cache')
                  CALL start_clock('pf_fftinner')
                  CALL cft3t( exx_cus%fft_small, psic, exx_cus%fft_small%nr1t, exx_cus%fft_small%nr2t, exx_cus%fft_small%nr3t, &
              &exx_cus%fft_small%nrx1t, exx_cus%fft_small%nrx2t, exx_cus%fft_small%nrx3t, -1 )
                  CALL stop_clock('pf_fftinner')
                  if(jdmax==0) then
                     CALL start_clock('pf_cache')
                     prods_g(1:exx_cus%fft_small%ngmt,1) = psic(exx_cus%fft_small%nlt(1:exx_cus%fft_small%ngmt))
                     CALL stop_clock('pf_cache')
      !apply fac
                     CALL start_clock('pf_product')
                     prods_g(1:exx_cus%fft_small%ngmt,1)=fac(1:exx_cus%fft_small%ngmt)*prods_g(1:exx_cus%fft_small%ngmt,1)
                     CALL stop_clock('pf_product')
                  else
                     CALL start_clock('pf_cache')
                     prods_g(1:exx_cus%fft_small%ngmt,1)=0.5d0*(psic(exx_cus%fft_small%nlt(1:exx_cus%fft_small%ngmt))+&
                          &conjg( psic(exx_cus%fft_small%nltm(1:exx_cus%fft_small%ngmt))))
                     prods_g(1:exx_cus%fft_small%ngmt,2)=(0.d0,-0.5d0)*(psic(exx_cus%fft_small%nlt(1:exx_cus%fft_small%ngmt))-&
                          &conjg( psic(exx_cus%fft_small%nltm(1:exx_cus%fft_small%ngmt))))
                     CALL stop_clock('pf_cache')
                     CALL start_clock('pf_product')
                     prods_g(1:exx_cus%fft_small%ngmt,1)=fac(1:exx_cus%fft_small%ngmt)*prods_g(1:exx_cus%fft_small%ngmt,1)
                     prods_g(1:exx_cus%fft_small%ngmt,2)=fac(1:exx_cus%fft_small%ngmt)*prods_g(1:exx_cus%fft_small%ngmt,2)
                     CALL stop_clock('pf_product')
                  endif
            !put back on large cell
                  psic=0.d0
                  CALL start_clock('pf_cache')
                  if(jdmax==0) then
                     psic(exx_cus%fft_small%nlt(1:exx_cus%fft_small%ngmt))  = prods_g(1:exx_cus%fft_small%ngmt,1)
                     psic(exx_cus%fft_small%nltm(1:exx_cus%fft_small%ngmt)) = CONJG( prods_g(1:exx_cus%fft_small%ngmt,1))
                  else
                       psic(exx_cus%fft_small%nlt(1:exx_cus%fft_small%ngmt))= prods_g(1:exx_cus%fft_small%ngmt,1)+&
                            &(0.d0,1.d0)*prods_g(1:exx_cus%fft_small%ngmt,2)
                       psic(exx_cus%fft_small%nltm(1:exx_cus%fft_small%ngmt)) = CONJG( prods_g(1:exx_cus%fft_small%ngmt,1) )+&
                         &(0.d0,1.d0)*CONJG( prods_g(1:exx_cus%fft_small%ngmt,2) )

                  endif
                  CALL stop_clock('pf_cache')
                  CALL start_clock('pf_fftinner')
                  CALL cft3t( exx_cus%fft_small, psic, exx_cus%fft_small%nr1t, exx_cus%fft_small%nr2t, exx_cus%fft_small%nr3t, &
                       &exx_cus%fft_small%nrx1t, exx_cus%fft_small%nrx2t, exx_cus%fft_small%nrx3t, +1 )
                  CALL stop_clock('pf_fftinner')
                  CALL start_clock('pf_cache')
                  prods(1:exx_cus%fft_small%nrxxt,1)=dble(psic(1:exx_cus%fft_small%nrxxt))
                  if(jdmax==1) then
                     prods(1:exx_cus%fft_small%nrxxt,2)=dimag(psic(1:exx_cus%fft_small%nrxxt))
                  endif
                  CALL stop_clock('pf_cache')
!!!!!!!!!!

                  do jd=0,jdmax
                     CALL start_clock('pf_plane2')
                     do iz_s=1,nr3small
                        plane(:)=0.d0
                        do jj=1,exx_cus%m(2)
                           do kk=1,exx_cus%m(1)
                        !      do iy_s=1,exx_cus%fft_small%nr2t
                        !         do ix_s=1,exx_cus%fft_small%nr1t
                        !            iy=mod(iy_s+exx_cus%fft_small%nr2t*(jj-1)-1,exx_cus%fft_g2r%nr2t)+1
                        !            ix=mod(ix_s+exx_cus%fft_small%nr1t*(kk-1)-1,exx_cus%fft_g2r%nr1t)+1
                        !            iqq=(iy-1)*exx_cus%fft_g2r%nrx1t+ix
                        !            iqq_s=(iy_s-1)*exx_cus%fft_small%nrx1t+ix_s
                        !            plane(iqq)=prods((iz_s-1)*(exx_cus%fft_small%nrx1t*exx_cus%fft_small%nrx2t)+iqq_s,jd+1)
                        !   
                        !         enddo
                        !      enddo
                              do iqq_s=1,exx_cus%fft_small%nrx1t*exx_cus%fft_small%nrx2t!problem if nrx1t /= nr1t
                                 plane(exx_cus%s2r_xy(iqq_s,kk,jj))=&
                                      &prods((iz_s-1)*(exx_cus%fft_small%nrx1t*exx_cus%fft_small%nrx2t)+iqq_s,jd+1)
                              enddo
                           enddo
                        enddo
                        do ii=1,exx_cus%m(3)
                           do iqq=1,exx_cus%fft_g2r%nr1t*exx_cus%fft_g2r%nr2t
                              prod_r_red(iqq,iz_s,ii)=plane(iqq)
                           enddo
                        enddo
                     enddo
                     CALL stop_clock('pf_plane2')
                     CALL start_clock('pf_product')
                     if(.not.exx_cus%l_localized) then
                        do ii=1,exx_cus%m(3)
                           do iz_s=1,nr3small
                              do iqq=1,exx_cus%fft_g2r%nr1t*exx_cus%fft_g2r%nr2t
                                 vexc_red(iqq,iz_s,ii)=vexc_red(iqq,iz_s,ii)+prod_r_red(iqq,iz_s,ii)*&
                &exx_cus%wfc_red(iqq,iz_s,ii,jv+jd,ispin)*wg(jv+jd,ispin)*dble(exx_cus%nspin)/2.d0
                              enddo
                           enddo
                        enddo
                     else
                        write(stdout,*) 'Not implemented yet'
                        stop
                     endif
                     CALL stop_clock('pf_product')
                  enddo
                  

               enddo!on nbndv
            enddo
      !end loop
         enddo
      enddo
      CALL stop_clock('pf_inner')
!!!!!!!!!!!!!   
!!!!!!!!!!!
!first the internal case

      do iz_s=rz_start_s,rz_end_s
         do ii=1,exx_cus%m(3)
            iz=mod(iz_s+exx_cus%fft_small%nr3t*(ii-1)-1,exx_cus%fft_g2r%nrx3t)+1
            idest=z2proc(iz)
            if(me_pool==idest) then
               do iqq=1,exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t
                  vexc((iz-rz_start)*exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t+iqq)=vexc_red(iqq,iz_s-rz_start_s+1,ii)
               enddo
            endif
         enddo
      enddo
      allocate(proc_list(nproc))
      
      allocate(b_plane(exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t,nr3small_max*exx_cus%m(3)))
      allocate(b_z(nr3small_max*exx_cus%m(3)))
!loop on task delta
      do ip_delta=1,nproc-1
         if(mod(ip_delta,2)==0) then
            !if(mod(me_pool+1,2)==0) then!even
            !   if(mod((me_pool+1)/2,2)==0) then
            !      token=0
            !   else
            !      token=1
            !   endif
            !   
            !else
            !   
            !   if(mod((me_pool+2)/2,2)==0) then
            !      token=0
            !   else
            !      token=1
            !   endif
            !   
            !endif
            proc_list=0
            do ip=1,nproc
               if(proc_list(ip)==0) then
                  if(proc_list(mod(ip+ip_delta-1,nproc)+1)==0) then
                     proc_list(ip)=-1
                     proc_list(mod(ip+ip_delta-1,nproc)+1)=1
                  else
                    
                  endif
               endif
            enddo

            if(proc_list(me_pool+1) ==-1) then
               token=0
            else
               token=1
            endif

            
         else
            if(mod(me_pool+1,2)==0) then
               token=0
            else
               token=1
            endif
           
         endif
         do ip_todo=1,2
            if(mod(ip_todo+token,2)==0) then
!if I am a sender
  !loop on my data to see if and what I have to send
               nplane=0
               do iz_s=rz_start_s,rz_end_s
                  do ii=1,exx_cus%m(3)
                     iz=mod(iz_s+exx_cus%fft_small%nr3t*(ii-1)-1,exx_cus%fft_g2r%nrx3t)+1
                     idest=z2proc(iz)
                     if(idest==mod(me_pool+ip_delta,nproc)) then
                        nplane=nplane+1
                        do iqq=1,exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t            
                           b_plane(iqq,nplane)=vexc_red(iqq,iz_s-rz_start_s+1,ii)
                        enddo
                        b_z(nplane)=iz
                     endif
                  enddo
               enddo
!send nplane
#if defined(__MPI)
               idest=mod(me_pool+ip_delta,nproc)
               CALL MPI_ISEND( nplane,1, MPI_INTEGER, idest, 0, intra_pool_comm, req,IERR ) 
               CALL MPI_WAIT(req,istatus,ierr)
              
               if(nplane>0) then
                  CALL MPI_ISEND( b_plane,exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t*nplane, MPI_DOUBLE_PRECISION, &   
                       &idest, 1, intra_pool_comm, req,IERR )
                  CALL MPI_WAIT(req,istatus,ierr)
                  CALL MPI_ISEND( b_z,nplane, MPI_INTEGER,idest, 2, intra_pool_comm, req,IERR )
                  CALL MPI_WAIT(req,istatus,ierr)
                  
               endif
#endif

            else
!if I am receiver 
  !see if and what I have to receive
#if defined(__MPI)
               iorig=me_pool-ip_delta
               if(iorig<0) iorig=iorig+nproc
               CALL MPI_RECV( nplane,1, MPI_INTEGER, iorig, 0, intra_pool_comm, istatus,IERR )
              
               if(nplane>0) then
                  CALL MPI_RECV( b_plane,exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t*nplane, MPI_DOUBLE_PRECISION, &
                       &iorig, 1, intra_pool_comm, istatus,IERR )
                  CALL MPI_RECV( b_z,nplane, MPI_INTEGER,iorig, 2, intra_pool_comm, istatus,IERR )
                
                  do ii=1,nplane
                     do iqq=1,exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t
                        vexc((b_z(ii)-rz_start)*exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t+iqq)=&
                             &b_plane(iqq,ii)
                     enddo
                  enddo
               endif
#endif
            endif
         enddo
      enddo
      deallocate(b_plane,b_z,proc_list)

      






!!!!!!!!!!!!!
!now find vexc from vexc_red
!      do iz_s=1,exx_cus%fft_small%nr3t
!         !send and receive z and z+alat                                                            
!         do ii=1,exx_cus%m(3)
!            iz=mod(iz_s+exx_cus%fft_small%nr3t*(ii-1)-1,exx_cus%fft_g2r%nrx3t)+1
!           
!           
!            !do periodic replica
!            if(iz_s >= rz_start_s .and. iz_s <= rz_end_s) then
!               !if Z is mine determine owner and send               
!               idest=z2proc(iz)
!               !put plane on small plane
!               if(me_pool==idest) then
!                  do iqq=1,exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t
!                     vexc((iz-rz_start)*exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t+iqq)=vexc_red(iqq,iz_s-rz_start_s+1,ii)
!                  enddo
!               else
!                  do iqq=1,exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t
!                     plane(iqq)=vexc_red(iqq,iz_s-rz_start_s+1,ii)
!                  enddo
!                  CALL MPI_ISEND( plane, exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t, MPI_DOUBLE_PRECISION, &
!                       &idest, iz, intra_pool_comm, req,IERR )
!                  CALL MPI_WAIT(req,istatus,ierr)
!                  
!               endif
!            else
!                  !if Z on  large cell is mine receive it 
!               if(z2proc(iz)==me_pool) then
!                  iorig=z2proc_s(iz_s)
!                  CALL MPI_RECV( plane, exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t, MPI_DOUBLE_PRECISION, &
!                       &iorig, iz, intra_pool_comm, istatus, IERR )
!                  do iqq=1,exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t
!                     vexc((iz-rz_start)*exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t+iqq)=plane(iqq)
!                  enddo
!               endif
!            endif
!            
!         enddo
!      enddo

!!!!!!!!!
     
!do fft and back to standard ordering 
!do scalar producs
      psic(1:exx_cus%fft_g2r%nrxxt)=dcmplx(vexc(1:exx_cus%fft_g2r%nrxxt),0.d0)
      CALL start_clock('pf_fftext')
      CALL cft3t( exx_cus%fft_g2r, psic, exx_cus%fft_g2r%nr1t, exx_cus%fft_g2r%nr2t, exx_cus%fft_g2r%nr3t, &
           &exx_cus%fft_g2r%nrx1t, exx_cus%fft_g2r%nrx2t, exx_cus%fft_g2r%nrx3t, -2 )
      CALL stop_clock('pf_fftext')
      vexc_g(1:exx_cus%fft_g2r%npwt) = psic(exx_cus%fft_g2r%nlt(1:exx_cus%fft_g2r%npwt))



!put in the order or wfcs                   
      CALL start_clock('pf_mergesplit')
      if(exx_cus%fft_g2r%dual_t==4.d0) then
         xpsi(1:exx_cus%fft_g2r%npwt)=vexc_g(1:exx_cus%fft_g2r%npwt)
      else
         call mergewf(vexc_g,evc_g,exx_cus%fft_g2r%npwt,exx_cus%fft_g2r%ig_l2gt,me_pool,nproc,ionode_id,intra_pool_comm)
         call splitwf(xpsi,evc_g,npw,ig_l2g,me_pool,nproc,ionode_id,intra_pool_comm)
      endif
      CALL stop_clock('pf_mergesplit')


     
      deallocate(fac)
      deallocate(prodr,prods)
      deallocate(planes)
      deallocate(z2proc_s)
      deallocate(z2proc)
      deallocate(prods_g)
      deallocate(vexc)
      deallocate(psi_t,psi_r,evc_g)
      deallocate(psi_r_red,plane,prod_r_red)
      deallocate(vexc_red)
      CALL stop_clock('periodic_fock')
      return
    END SUBROUTINE periodic_fock_cus



  


    SUBROUTINE periodic_dft_exchange(nbnds,psi,exx_cus)
!experimental version work just with factor 1/2

      USE io_global, ONLY : stdout, ionode,ionode_id
      USE mp_global, ONLY : me_pool,intra_pool_comm
      USE cell_base, ONLY : at, alat, tpiba, omega, tpiba2,bg
      USE constants, ONLY : e2, pi, tpi, fpi, RYTOEV
      USE wavefunctions_module, ONLY : psic
      USE mp,        ONLY : mp_sum
      USE mp_world,  ONLY : world_comm, nproc
      USE wvfct,     ONLY : npwx, npw
      USE gvect
      USE mp_wave,   ONLY : mergewf,splitwf

      implicit none

      INTEGER, INTENT(in) :: nbnds!total number of states
      COMPLEX(kind=DP), INTENT(in) :: psi(npwx,nbnds)!wavefunctions   
      TYPE(exchange_cus), INTENT(in) :: exx_cus

      TYPE(fft_cus) :: fft_small
      INTEGER, ALLOCATABLE :: r2s_xy(:) !real to small XY index
      INTEGER, ALLOCATABLE :: r2s_z(:)!real to small Z index (large)
      
      INTEGER :: i,j,k,n,ii,jj,kk,ig,iv,jv
      INTEGER :: ix,iy,iz, ix_s,iy_s,iz_s,iz_eff
      INTEGER :: rz_start,rz_end,iqq,iqq_s,rz_start_s,rz_end_s
      REAL(kind=DP), ALLOCATABLE :: fac(:),prodr(:),prods(:),planes(:),vexc(:)
      REAL(kind=DP) :: qq_fact, sca
      INTEGER :: iorig, idest
      INTEGER,ALLOCATABLE :: z2proc_s(:),z2proc(:)
      INTEGER :: req, ierr
#if defined(__MPI)
      INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
      COMPLEX(kind=DP), ALLOCATABLE :: prods_g(:)
      COMPLEX(kind=DP), ALLOCATABLE :: psi_t(:),evc_g(:)
      REAL(kind=DP), ALLOCATABLE :: psi_r(:)


      !setup small grid
      fft_small%at_t(1:3,1:3)=exx_cus%fft_g2r%at_t(1:3,1:3)
      fft_small%bg_t(1:3,1:3)=exx_cus%fft_g2r%bg_t(1:3,1:3)
      fft_small%alat_t=exx_cus%fft_g2r%alat_t/2.d0
      fft_small%omega_t=exx_cus%fft_g2r%omega_t/8.d0
      fft_small%tpiba_t=exx_cus%fft_g2r%tpiba_t*2.d0
      fft_small%tpiba2_t=exx_cus%fft_g2r%tpiba2_t*4.d0
     
      fft_small%ecutt=exx_cus%fft_g2r%ecutt
      fft_small%dual_t=exx_cus%fft_g2r%dual_t

      call initialize_fft_custom_cell(fft_small)

      write(stdout,*) 'Dimensions of real cell'
      write(stdout,*) exx_cus%fft_g2r%nr1t,exx_cus%fft_g2r%nr2t,exx_cus%fft_g2r%nr3t
      write(stdout,*) exx_cus%fft_g2r%nrx1t,exx_cus%fft_g2r%nrx2t,exx_cus%fft_g2r%nrx3t

      write(stdout,*) 'Dimensions of small cell'

      write(stdout,*) fft_small%nr1t,fft_small%nr2t,fft_small%nr3t
      write(stdout,*) fft_small%nrx1t,fft_small%nrx2t,fft_small%nrx3t  
    

      
      allocate(r2s_xy(exx_cus%fft_g2r%nrxxt))
      
      allocate(r2s_z(exx_cus%fft_g2r%nrxxt))
      !setup correspondence grids
  
      rz_start=0
      rz_end =0
      do ii=1,me_pool + 1
         rz_start=rz_end+1
         rz_end=rz_end+exx_cus%fft_g2r%dfftt%npp(ii)
      end do

      rz_start_s=0
      rz_end_s=0
      do ii=1,me_pool + 1
         rz_start_s=rz_end_s+1
         rz_end_s=rz_end_s+fft_small%dfftt%npp(ii)
      end do

      allocate(z2proc_s(fft_small%nr3t))
      allocate(z2proc(exx_cus%fft_g2r%nr3t))
      allocate(vexc(exx_cus%fft_g2r%nrxxt))
      allocate( evc_g( exx_cus%fft_g2r%ngmt_g ) )
      j=0
      k=0
      do ii=1,nproc
         j=k+1
         k=k+fft_small%dfftt%npp(ii)
         z2proc_s(j:k)=ii-1
      end do

      j=0
      k=0
      do ii=1,nproc
         j=k+1
         k=k+exx_cus%fft_g2r%dfftt%npp(ii)
         z2proc(j:k)=ii-1
      end do



      r2s_xy(:)=0
      do iz=1,exx_cus%fft_g2r%dfftt%npp(me_pool+1)
         do iy=1,exx_cus%fft_g2r%nr2t
            do ix=1,exx_cus%fft_g2r%nr1t
               iqq=(iz-1)*(exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t)+(iy-1)*exx_cus%fft_g2r%nrx1t+ix
               iy_s=mod(iy-1,fft_small%nr2t)+1
               ix_s=mod(ix-1,fft_small%nr1t)+1
               iqq_s=(iy_s-1)*fft_small%nrx1t+ix_s
               r2s_xy(iqq)=iqq_s!XY correspondance only
               iz_eff=iz+rz_start-1
               iz_s=mod(iz_eff-1,fft_small%nr3t)+1
            enddo
         enddo
      enddo
 
      allocate(fac(fft_small%ngmt))

      !setup fac
      do ig=1,fft_small%ngmt

         qq_fact = fft_small%gt(1,ig)**2.d0 + fft_small%gt(2,ig)**2.d0 + fft_small%gt(3,ig)**2.d0

         if (qq_fact > 1.d-8) then
            fac(ig)=(e2*fpi/(fft_small%tpiba2_t*qq_fact))*(1.d0-dcos(dsqrt(qq_fact)*&
                 &(exx_cus%truncation_radius/2.d0)*fft_small%tpiba_t))
         else
            fac(ig)=e2*fpi*((exx_cus%truncation_radius/2.d0)**2.d0/2.d0)

         endif

      end do
      fac(:)=fac(:)/omega


      allocate(prodr(exx_cus%fft_g2r%nrxxt))
      allocate(prods(fft_small%nrxxt))
      allocate(planes(fft_small%nrx1t*fft_small%nrx2t))
      allocate(prods_g(fft_small%ngmt))


      allocate(psi_t(exx_cus%fft_g2r%npwt))
      allocate(psi_r(exx_cus%fft_g2r%nrxxt))

      !loop on KS states
      do iv=1,nbnds
         vexc(:)=0.d0


         if(exx_cus%fft_g2r%dual_t==4.d0) then
            psi_t(1:exx_cus%fft_g2r%npwt)=psi(1:exx_cus%fft_g2r%npwt,iv)
         else
            call mergewf(psi(:,iv),evc_g,npw,ig_l2g,me_pool,nproc,ionode_id,intra_pool_comm)
            call splitwf(psi_t(:),evc_g,exx_cus%fft_g2r%npwt,exx_cus%fft_g2r%ig_l2gt,&
                 &me_pool,nproc,ionode_id,intra_pool_comm)
         endif

         psic(:)=(0.d0,0.d0)
         psic(exx_cus%fft_g2r%nlt(1:exx_cus%fft_g2r%npwt))  = psi_t(1:exx_cus%fft_g2r%npwt)
         psic(exx_cus%fft_g2r%nltm(1:exx_cus%fft_g2r%npwt)) = CONJG( psi_t(1:exx_cus%fft_g2r%npwt) )

         CALL cft3t( exx_cus%fft_g2r, psic, exx_cus%fft_g2r%nr1t, exx_cus%fft_g2r%nr2t, exx_cus%fft_g2r%nr3t, &
              &exx_cus%fft_g2r%nrx1t, exx_cus%fft_g2r%nrx2t, exx_cus%fft_g2r%nrx3t, 2 )
         psi_r(1:exx_cus%fft_g2r%nrxxt)= DBLE(psic(1:exx_cus%fft_g2r%nrxxt))


      !loop on KS valence states
      
         do jv=1,exx_cus%nbndv(1)
      !do product
            prodr(1:exx_cus%fft_g2r%nrxxt)=psi_r(1:exx_cus%fft_g2r%nrxxt)*exx_cus%wfc(1:exx_cus%fft_g2r%nrxxt,jv,1)
      !put on small cell
            !loop on real Z direction
            prods(:)=0.d0
            do iz=1,exx_cus%fft_g2r%nr3t
               if(iz >= rz_start .and. iz <= rz_end) then
                  planes(:)=0.d0
                  do iy=1,exx_cus%fft_g2r%nr2t
                     do ix=1,exx_cus%fft_g2r%nr1t
                        iqq=(iz-rz_start)*(exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t)+(iy-1)*exx_cus%fft_g2r%nrx1t+ix
                        planes(r2s_xy(iqq))=planes(r2s_xy(iqq))+prodr(iqq)
                     enddo
                  enddo
            !if Z is mine determine owner and send it
                  iz_s=mod(iz-1,fft_small%nr3t)+1
                  idest=z2proc_s(iz_s)
            !put plane on small plane
                  if(me_pool==idest) then
                     do iqq_s=1,fft_small%nrx1t*fft_small%nrx2t
                        prods((iz_s-rz_start_s)*(fft_small%nrx1t*fft_small%nrx2t)+iqq_s)=&
                            & prods((iz_s-rz_start_s)*(fft_small%nrx1t*fft_small%nrx2t)+iqq_s)+planes(iqq_s)
                     enddo
                  else
#if defined(__MPI)
                     CALL MPI_ISEND( planes, fft_small%nrx1t*fft_small%nrx2t, MPI_DOUBLE_PRECISION, &
                          &idest, iz, intra_pool_comm, req,IERR )
                     CALL MPI_WAIT(req,istatus,ierr)
#endif
                  endif
                  
               else
                  iz_s=mod(iz-1,fft_small%nr3t)+1
            !if Z o  small cell is mine receive it
                  if(z2proc_s(iz_s)==me_pool) then
                     iorig=z2proc(iz)
#if defined(__MPI)
                     CALL MPI_RECV( planes, fft_small%nrx1t*fft_small%nrx2t, MPI_DOUBLE_PRECISION, &
                          &iorig, iz, intra_pool_comm, istatus, IERR )
#endif
                     do iqq_s=1,fft_small%nrx1t*fft_small%nrx2t
                        prods((iz_s-rz_start_s)*(fft_small%nrx1t*fft_small%nrx2t)+iqq_s)=&
                             &prods((iz_s-rz_start_s)*(fft_small%nrx1t*fft_small%nrx2t)+iqq_s)+planes(iqq_s)
                     enddo

                  endif
               endif
            enddo
      ! do fft
            psic(1:fft_small%nrxxt)=cmplx(prods(1:fft_small%nrxxt),0.d0)
            CALL cft3t( fft_small, psic, fft_small%nr1t, fft_small%nr2t, fft_small%nr3t, &
                 &fft_small%nrx1t, fft_small%nrx2t, fft_small%nrx3t, -1 )
            prods_g(1:fft_small%ngmt) = psic(fft_small%nlt(1:fft_small%ngmt))

      !apply fac
            prods_g(1:fft_small%ngmt)=fac(1:fft_small%ngmt)*prods_g(1:fft_small%ngmt)
            !put back on large cell
            psic=0.d0
            psic(fft_small%nlt(1:fft_small%ngmt))  = prods_g(1:fft_small%ngmt)
            psic(fft_small%nltm(1:fft_small%ngmt)) = CONJG( prods_g(1:fft_small%ngmt))
            CALL cft3t( fft_small, psic, fft_small%nr1t, fft_small%nr2t, fft_small%nr3t, &
                 &fft_small%nrx1t, fft_small%nrx2t, fft_small%nrx3t, +1 )
            prods(1:fft_small%nrxxt)=dble(psic(1:fft_small%nrxxt))
      !loop on small z grid
            do iz_s=1,fft_small%nr3t
            !send and receive z and z+alat
               do ii=1,2
                  iz=iz_s+fft_small%nr3t*(ii-1)

            !do periodic replica
                  if(iz_s >= rz_start_s .and. iz_s <= rz_end_s) then
                   
            !if Z is mine determine owner and send
                     idest=z2proc(iz)
            !put plane on small plane
                     if(me_pool==idest) then
                        do iqq_s=1,fft_small%nrx1t*fft_small%nrx2t
                           planes(iqq_s)=prods((iz_s-rz_start_s)*(fft_small%nrx1t*fft_small%nrx2t)+iqq_s)
                        enddo
!do replicas
                        do jj=1,2
                           do kk=1,2
                              do iy_s=1,fft_small%nr2t
                                 do ix_s=1,fft_small%nr1t
                                    iy=iy_s+fft_small%nr2t*(jj-1)
                                    ix=ix_s+fft_small%nr1t*(kk-1)
                                    iqq=(iz-rz_start)*(exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t)+(iy-1)*exx_cus%fft_g2r%nrx1t+ix
                                    prodr(iqq)=planes(r2s_xy(iqq))
                                 enddo
                              enddo
                            enddo
                        enddo
                     else
                        do iqq_s=1,fft_small%nrx1t*fft_small%nrx2t
                           planes(iqq_s)=prods((iz_s-rz_start_s)*(fft_small%nrx1t*fft_small%nrx2t)+iqq_s)
                        enddo
#if defined(__MPI)
                        CALL MPI_ISEND( planes, fft_small%nrx1t*fft_small%nrx2t, MPI_DOUBLE_PRECISION, &
                             &idest, iz, intra_pool_comm, req,IERR )
                        CALL MPI_WAIT(req,istatus,ierr)
#endif
                     endif

                  else
                     !if Z on  large cell is mine receive it 
                     if(z2proc(iz)==me_pool) then
                        iorig=z2proc_s(iz_s)
#if defined(__MPI)
                        CALL MPI_RECV( planes, fft_small%nrx1t*fft_small%nrx2t, MPI_DOUBLE_PRECISION, &
                             &iorig, iz, intra_pool_comm, istatus, IERR )
#endif
                        do jj=1,2
                           do kk=1,2
                              do iy_s=1,fft_small%nr2t
                                 do ix_s=1,fft_small%nr1t
                                    iy=iy_s+fft_small%nr2t*(jj-1)
                                    ix=ix_s+fft_small%nr1t*(kk-1)
                                    iqq=(iz-rz_start)*(exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t)+(iy-1)*exx_cus%fft_g2r%nrx1t+ix
                                    prodr(iqq)=planes(r2s_xy(iqq))
                                 enddo
                              enddo
                           enddo
                        enddo

                     endif

                  endif
               enddo!ii

            enddo
            !do product
            vexc(1:exx_cus%fft_g2r%nrxxt)=vexc(1:exx_cus%fft_g2r%nrxxt)+prodr(1:exx_cus%fft_g2r%nrxxt)*&
                 &exx_cus%wfc(1:exx_cus%fft_g2r%nrxxt,jv,1)
      !sum up result terms
      
      !end loop
         enddo
!do scalar producs
         sca=0.d0
         do i=1,exx_cus%fft_g2r%nrxxt
            sca=sca+psi_r(i)*vexc(i)
         enddo
         call mp_sum(sca,world_comm)
         sca=sca/dble(exx_cus%fft_g2r%nr1t*exx_cus%fft_g2r%nr2t*exx_cus%fft_g2r%nr3t)
         write(stdout,*) 'PERIODIC EXCHANGE', iv, sca
         FLUSH(stdout)
      !end loop
      enddo
      deallocate(r2s_xy,r2s_z)
      deallocate(fac)
      deallocate(prodr,prods)
      deallocate(planes)
      deallocate(z2proc_s)
      deallocate(z2proc)
      deallocate(prods_g)
      deallocate(vexc)
      deallocate(psi_t,psi_r,evc_g)
      return
    END SUBROUTINE periodic_dft_exchange


    SUBROUTINE setup_exx_cus(nspin,num_nbndv_max,num_nbndv,ks_wfcs, exx_cus, dual, cutoff, truncation_radius)
!ATTENZIONE now only for cubic cell to be extended   


      USE io_global,  ONLY : stdout, ionode, ionode_id
      USE wvfct,    ONLY : npwx, npw, nbnd
      USE gvecw, ONLY : ecutwfc
      USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2,bg
      USE constants, ONLY : e2, pi, tpi, fpi, RYTOEV
      USE wavefunctions_module, ONLY : psic
      USE mp_global, ONLY : intra_pool_comm, me_pool
      USE gvect
      USE mp_wave, ONLY : mergewf,splitwf
      USE mp, ONLY : mp_barrier, mp_sum
      USE mp_world, ONLY : world_comm, mpime, nproc

      IMPLICIT NONE

      INTEGER, INTENT(in) :: nspin!spin multiplicity
      INTEGER, INTENT(in) ::num_nbndv_max!max number of valence states for both spins
      INTEGER, INTENT(in) :: num_nbndv(2)!number of valence states
      COMPLEX(kind=DP), INTENT(in) :: ks_wfcs(npwx, num_nbndv_max,2)!KS valence wfcs
      TYPE(exchange_cus), INTENT(out) :: exx_cus!the structure to be created
      REAL(kind=DP), INTENT(in) :: dual!for defining the real space grid
      REAL(kind=DP), INTENT(in) :: cutoff !for the defining the G space grid
      REAL(kind=DP), INTENT(in) :: truncation_radius!Bohr

      REAL(kind=DP) :: qq_fact
      INTEGER :: ig,ii,i,j,iv,ir,n_max,is,k,jj,kk
      COMPLEX(kind=DP), ALLOCATABLE :: state_fc_t(:,:),evc_g(:)
      INTEGER :: rz_start,rz_end,ix,iy,iz,iqq,iqq_s,ix_s,iy_s
      INTEGER :: rz_start_s,rz_end_s,nr3small
      INTEGER :: iorig, idest,iz_s,jv
      INTEGER,ALLOCATABLE :: z2proc_s(:),z2proc(:)
      INTEGER iplane
      REAL(kind=DP), ALLOCATABLE :: plane(:)
      INTEGER :: req, ierr
#if defined(__MPI)
      INTEGER :: istatus(MPI_STATUS_SIZE)
#endif

      CALL start_clock('setup_exx')
!setup parameters

      exx_cus%nbndv(1:2)=num_nbndv(1:2)
      exx_cus%dual=dual
      exx_cus%cutoff=cutoff
      exx_cus%truncation_radius=truncation_radius

      exx_cus%l_localized=l_exchange_localized
      exx_cus%thrs_loc=exchange_thrs_loc
      exx_cus%nspin=nspin

!define grids

      exx_cus%fft_g2r%ecutt=ecutwfc
      exx_cus%fft_g2r%dual_t=dual
      call mp_barrier( world_comm )
      write(stdout,*) 'Before initialize_fft_custom',exx_cus%fft_g2r%ecutt,exx_cus%fft_g2r%dual_t
      FLUSH(stdout)
      call initialize_fft_custom(exx_cus%fft_g2r)
      write(stdout,*) "GRID G to R", exx_cus%fft_g2r%nr1t, exx_cus%fft_g2r%nr2t, exx_cus%fft_g2r%nr3t
      write(stdout,*) "GRID G to R",exx_cus%fft_g2r%npwt
      FLUSH(stdout)

      exx_cus%fft_r2g%ecutt=cutoff
      exx_cus%fft_r2g%dual_t=ecutwfc*dual/cutoff
      call initialize_fft_custom(exx_cus%fft_r2g)
      write(stdout,*) "GRID R to G", exx_cus%fft_r2g%nr1t, exx_cus%fft_r2g%nr2t, exx_cus%fft_r2g%nr3t
      write(stdout,*) "GRID R to G",exx_cus%fft_r2g%npwt
      FLUSH(stdout)
      
      if(l_exchange_turbo) then
  !setup small grid
         exx_cus%m(1:3)=exchange_m(1:3)
         exx_cus%n(1:3)=exchange_n(1:3)
         do i=1,3
            exx_cus%fft_small%at_t(1:3,i)=exx_cus%fft_g2r%at_t(1:3,i)
            exx_cus%fft_small%bg_t(1:3,i)=exx_cus%fft_g2r%bg_t(1:3,i)
         enddo
         exx_cus%fft_small%alat_t=exx_cus%fft_g2r%alat_t*dble(exchange_n(1))/dble(exchange_m(1))
         exx_cus%fft_small%omega_t=exx_cus%fft_g2r%omega_t*(dble(exchange_n(1))/dble(exchange_m(1)))**3.d0
         exx_cus%fft_small%tpiba_t=exx_cus%fft_g2r%tpiba_t/(dble(exchange_n(1))/dble(exchange_m(1)))
         exx_cus%fft_small%tpiba2_t=exx_cus%fft_g2r%tpiba2_t/(dble(exchange_n(1))/dble(exchange_m(1)))**2.d0

         exx_cus%fft_small%ecutt=exx_cus%fft_g2r%ecutt
         exx_cus%fft_small%dual_t=exx_cus%fft_g2r%dual_t

         call initialize_fft_custom_cell(exx_cus%fft_small)

         allocate(exx_cus%r2s_xy(exx_cus%fft_g2r%nrxxt,exx_cus%n(1),exx_cus%n(2)))
         allocate(exx_cus%fac_s(exx_cus%fft_small%ngmt))
     !setup correspondence grids                                                                                                                                     

         rz_start=0
         rz_end =0
         do ii=1,me_pool + 1
            rz_start=rz_end+1
            rz_end=rz_end+exx_cus%fft_g2r%dfftt%npp(ii)
         end do

         exx_cus%r2s_xy(:,:,:)=0
         do i=1,exx_cus%n(1)
            do j=1,exx_cus%n(2)
               do iz=1,exx_cus%fft_g2r%dfftt%npp(me_pool+1)
                  do iy=1,exx_cus%fft_g2r%nr2t
                     do ix=1,exx_cus%fft_g2r%nr1t
                        iqq=(iz-1)*(exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t)+(iy-1)*exx_cus%fft_g2r%nrx1t+ix
                        iy_s=mod(iy-1+(j-1)*exx_cus%fft_g2r%nr2t,exx_cus%fft_small%nr2t)+1
                        ix_s=mod(ix-1+(i-1)*exx_cus%fft_g2r%nr1t,exx_cus%fft_small%nr1t)+1
                        iqq_s=(iy_s-1)*exx_cus%fft_small%nrx1t+ix_s
                        exx_cus%r2s_xy(iqq,i,j)=iqq_s!XY correspondance only
                     enddo
                  enddo
               enddo
            enddo
         enddo

         allocate(exx_cus%s2r_xy(exx_cus%fft_small%nrx1t*exx_cus%fft_small%nrx2t,exx_cus%m(1),exx_cus%m(2)))
           do jj=1,exx_cus%m(2)
              do kk=1,exx_cus%m(1)
                 do iy_s=1,exx_cus%fft_small%nr2t
                    do ix_s=1,exx_cus%fft_small%nr1t
                       iy=mod(iy_s+exx_cus%fft_small%nr2t*(jj-1)-1,exx_cus%fft_g2r%nr2t)+1
                       ix=mod(ix_s+exx_cus%fft_small%nr1t*(kk-1)-1,exx_cus%fft_g2r%nr1t)+1
                       iqq=(iy-1)*exx_cus%fft_g2r%nrx1t+ix
                       iqq_s=(iy_s-1)*exx_cus%fft_small%nrx1t+ix_s
                       exx_cus%s2r_xy(iqq_s,kk,jj)=iqq
                    enddo
                 enddo
              enddo
           enddo
      endif
!setup fac

      allocate(exx_cus%fac(exx_cus%fft_r2g%ngmt))
      
      do ig=1,exx_cus%fft_r2g%ngmt

         qq_fact = exx_cus%fft_r2g%gt(1,ig)**2.d0 + exx_cus%fft_r2g%gt(2,ig)**2.d0 + exx_cus%fft_r2g%gt(3,ig)**2.d0

         if (qq_fact > 1.d-8) then
            exx_cus%fac(ig)=(e2*fpi/(tpiba2*qq_fact))*(1.d0-dcos(dsqrt(qq_fact)*exx_cus%truncation_radius*tpiba))
         else
            exx_cus%fac(ig)=e2*fpi*(exx_cus%truncation_radius**2.d0/2.d0)

         endif

      end do
      exx_cus%fac(:)=exx_cus%fac(:)/omega

!put wfcs in real space
      allocate(exx_cus%wfc(exx_cus%fft_g2r%nrxxt,num_nbndv_max,nspin))

      allocate(state_fc_t(exx_cus%fft_g2r%npwt,num_nbndv_max))
      allocate( evc_g( exx_cus%fft_g2r%ngmt_g ) )

      do is=1,nspin
         if(exx_cus%fft_g2r%dual_t==4.d0) then
            state_fc_t(1:exx_cus%fft_g2r%npwt,1:exx_cus%nbndv(is))=ks_wfcs(1:exx_cus%fft_g2r%npwt,1:exx_cus%nbndv(is),is)
         else
            do ii=1,exx_cus%nbndv(is)
               call mergewf(ks_wfcs(:,ii,is),evc_g,npw,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
               call splitwf(state_fc_t(:,ii),evc_g,exx_cus%fft_g2r%npwt,exx_cus%fft_g2r%ig_l2gt,&
                    &mpime,nproc,ionode_id,intra_pool_comm)
            enddo
         endif





         do ii=1,exx_cus%nbndv(is),2
            psic(:)=(0.d0,0.d0)
            if(ii==exx_cus%nbndv(is)) then
               psic(exx_cus%fft_g2r%nlt(1:exx_cus%fft_g2r%npwt))  = state_fc_t(1:exx_cus%fft_g2r%npwt,ii)
               psic(exx_cus%fft_g2r%nltm(1:exx_cus%fft_g2r%npwt)) = CONJG( state_fc_t(1:exx_cus%fft_g2r%npwt,ii) )
            else
               psic(exx_cus%fft_g2r%nlt(1:exx_cus%fft_g2r%npwt))=state_fc_t(1:exx_cus%fft_g2r%npwt,ii)+&
                    &(0.d0,1.d0)*state_fc_t(1:exx_cus%fft_g2r%npwt,ii+1)
               psic(exx_cus%fft_g2r%nltm(1:exx_cus%fft_g2r%npwt)) = CONJG( state_fc_t(1:exx_cus%fft_g2r%npwt,ii) )+&
                    &(0.d0,1.d0)*CONJG( state_fc_t(1:exx_cus%fft_g2r%npwt,ii+1) )
            endif
            CALL cft3t( exx_cus%fft_g2r, psic, exx_cus%fft_g2r%nr1t, exx_cus%fft_g2r%nr2t, exx_cus%fft_g2r%nr3t, &
                 &exx_cus%fft_g2r%nrx1t, exx_cus%fft_g2r%nrx2t, exx_cus%fft_g2r%nrx3t, 2 )
            exx_cus%wfc(1:exx_cus%fft_g2r%nrxxt,ii,is)= DBLE(psic(1:exx_cus%fft_g2r%nrxxt))
            if(ii/=exx_cus%nbndv(1)) exx_cus%wfc(1:exx_cus%fft_g2r%nrxxt,ii+1,is)=&
                 &DIMAG(psic(1:exx_cus%fft_g2r%nrxxt))
         enddo
      enddo
!now put on the reduced grid
      if(l_exchange_turbo) then
!find maximum
         rz_start=0
         rz_end =0
         do ii=1,me_pool + 1
            rz_start=rz_end+1
            rz_end=rz_end+exx_cus%fft_g2r%dfftt%npp(ii)
         end do
         
         rz_start_s=0
         rz_end_s=0
         do ii=1,me_pool + 1
            rz_start_s=rz_end_s+1
            rz_end_s=rz_end_s+exx_cus%fft_small%dfftt%npp(ii)
         end do
         nr3small=rz_end_s-rz_start_s+1
         allocate(exx_cus%wfc_red(exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t,nr3small,exx_cus%m(3),num_nbndv_max,exx_cus%nspin))
       
         allocate(z2proc_s(exx_cus%fft_small%nr3t))
         allocate(z2proc(exx_cus%fft_g2r%nr3t))
      
         j=0
         k=0
         do ii=1,nproc
            j=k+1
            k=k+exx_cus%fft_small%dfftt%npp(ii)
            z2proc_s(j:k)=ii-1
         end do
         
         j=0
         k=0
         do ii=1,nproc
            j=k+1
            k=k+exx_cus%fft_g2r%dfftt%npp(ii)
            z2proc(j:k)=ii-1
         end do

      
         allocate(plane(exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t))

         do k=1,exx_cus%n(3)
            do iz=1,exx_cus%fft_g2r%nr3t
               if(iz >= rz_start .and. iz <= rz_end) then
               !if Z is mine determine owner and send it                    
                  iz_s=mod(iz-1+(k-1)*exx_cus%fft_g2r%nr3t,exx_cus%fft_small%nr3t)+1
                  iplane=(iz-1+(k-1)*exx_cus%fft_g2r%nr3t)/exx_cus%fft_small%nr3t+1
                  idest=z2proc_s(iz_s)
               !put plane on small plane             
               if(me_pool==idest) then
                  do is=1,exx_cus%nspin
                     do jv=1,exx_cus%nbndv(is)
                        do iqq=1,exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t
                           exx_cus%wfc_red(iqq,iz_s-rz_start_s+1,iplane,jv,is)=&
                                &exx_cus%wfc((iz-rz_start)*(exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t)+iqq,jv,is) 
                        enddo
                     enddo
                  enddo
               else
                   do is=1,exx_cus%nspin
                      do jv=1,exx_cus%nbndv(is)
                         do iqq=1,exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t
                            plane(iqq)=exx_cus%wfc((iz-rz_start)*(exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t)+iqq,jv,is)
                         enddo
#if defined(__MPI)
                         CALL MPI_ISEND( plane,exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t, MPI_DOUBLE_PRECISION, &
                              &idest, iz, intra_pool_comm, req,IERR )
                         CALL MPI_WAIT(req,istatus,ierr)
#endif
                      enddo
                   enddo
               endif

            else
               iz_s=mod(iz-1+(k-1)*exx_cus%fft_g2r%nr3t,exx_cus%fft_small%nr3t)+1
               iplane=(iz-1+(k-1)*exx_cus%fft_g2r%nr3t)/exx_cus%fft_small%nr3t+1
            !if Z o  small cell is mine receive it                          
               if(z2proc_s(iz_s)==me_pool) then
                  iorig=z2proc(iz)
                  do is=1,exx_cus%nspin
                     do jv=1,exx_cus%nbndv(is)
#if defined(__MPI)
                        CALL MPI_RECV( plane, exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t, MPI_DOUBLE_PRECISION, &
                             &iorig, iz, intra_pool_comm, istatus, IERR )
#endif
                        do iqq=1,exx_cus%fft_g2r%nrx1t*exx_cus%fft_g2r%nrx2t
                           exx_cus%wfc_red(iqq,iz_s-rz_start_s+1,iplane,jv,is)=plane(iqq)
                        enddo
                        
                     enddo
                  enddo
                  
               endif
            endif
         enddo
      enddo
      deallocate(z2proc,z2proc_s,plane)
   endif
       deallocate(state_fc_t,evc_g)

       if(exx_cus%l_localized) then
          allocate(exx_cus%n_loc(num_nbndv_max,nspin))
          allocate(exx_cus%tab_loc(exx_cus%fft_g2r%nrxxt,num_nbndv_max,nspin))!memory could be reduce here
          do is=1,nspin
             do iv=1,exx_cus%nbndv(is)
                exx_cus%n_loc(iv,is)=0
                do ir=1,exx_cus%fft_g2r%nrxxt
                   if((exx_cus%wfc(ir,iv,is)**2.d0>exx_cus%thrs_loc) )then
                      exx_cus%n_loc(iv,is)=exx_cus%n_loc(iv,is)+1
                      exx_cus%tab_loc(exx_cus%n_loc(iv,is),iv,is)=ir
                   endif
                enddo
                n_max=exx_cus%n_loc(iv,is)
                call mp_sum(n_max,world_comm)
                write(stdout,*) 'Using localized wfcs for exchange:',is,iv,n_max,&
                     &exx_cus%fft_g2r%nr1t*exx_cus%fft_g2r%nr2t*exx_cus%fft_g2r%nr3t
             enddo
          enddo
       endif

       CALL stop_clock('setup_exx')
      return
    END SUBROUTINE setup_exx_cus

    SUBROUTINE free_memory_exx_cus(exx_cus)

      IMPLICIT NONE

      TYPE(exchange_cus) :: exx_cus

      deallocate(exx_cus%wfc)
      deallocate(exx_cus%fac)
      if(l_exchange_turbo) then
         deallocate(exx_cus%fac_s)
         deallocate(exx_cus%r2s_xy)
         deallocate(exx_cus%wfc_red)
         deallocate(exx_cus%s2r_xy)
      endif
      if(exx_cus%l_localized) then
         deallocate(exx_cus%n_loc)
         deallocate(exx_cus%tab_loc)
      endif

      return

    END SUBROUTINE free_memory_exx_cus


    SUBROUTINE fock_cus(psi,xpsi,exx_cus)
!apply Fock operator to a wavefunction

      USE io_global, ONLY : stdout, ionode, ionode_id
      USE wvfct,     ONLY : npwx, npw, nbnd
      USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2,bg
      USE constants, ONLY : e2, pi, tpi, fpi, RYTOEV
      USE wavefunctions_module, ONLY : psic
      USE gvect
      USE mp_pools,  ONLY : intra_pool_comm
      USE mp_world,  ONLY : mpime, nproc
      USE mp_wave,   ONLY : mergewf,splitwf


      IMPLICIT NONE

      COMPLEX(kind=DP), INTENT(in) :: psi(npw)
      COMPLEX(kind=DP), INTENT(inout) :: xpsi(npw)
      TYPE(exchange_cus), INTENT(in) :: exx_cus


      REAL(kind=DP), ALLOCATABLE :: prods_r(:,:),psi_r(:),prod_tot(:)
      COMPLEX(kind=DP), ALLOCATABLE :: prods_g(:,:),evc_g(:), psi_t(:), prod_tot_g(:)
      INTEGER :: iv,ig
      INTEGER, ALLOCATABLE :: igkt(:)

      CALL start_clock('fock')

      allocate(prods_r(exx_cus%fft_g2r%nrxxt,exx_cus%nbndv(1)))
      allocate(prod_tot(exx_cus%fft_g2r%nrxxt))
      allocate(prods_g(exx_cus%fft_r2g%npwt,exx_cus%nbndv(1)))
      allocate(prod_tot_g(exx_cus%fft_g2r%npwt))
      allocate(psi_t(exx_cus%fft_g2r%npwt))
      allocate(psi_r(exx_cus%fft_g2r%nrxxt))

!put psi on the g2r G grid
      allocate( evc_g( exx_cus%fft_g2r%ngmt_g ) )


      allocate( igkt( exx_cus%fft_g2r%npwt ) )
      do ig=1,exx_cus%fft_g2r%npwt
         igkt(ig)=ig
      enddo
      

      if(exx_cus%fft_g2r%dual_t==4.d0) then
         psi_t(1:exx_cus%fft_g2r%npwt)=psi(1:exx_cus%fft_g2r%npwt)
      else
         call mergewf(psi,evc_g,npw,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
         call splitwf(psi_t(:),evc_g,exx_cus%fft_g2r%npwt,exx_cus%fft_g2r%ig_l2gt,&
                 &mpime,nproc,ionode_id,intra_pool_comm)
      endif



!trasform psi to R space
      
      psic(:)=(0.d0,0.d0)
     
      psic(exx_cus%fft_g2r%nlt(1:exx_cus%fft_g2r%npwt))  = psi_t(1:exx_cus%fft_g2r%npwt)
      psic(exx_cus%fft_g2r%nltm(1:exx_cus%fft_g2r%npwt)) = CONJG( psi_t(1:exx_cus%fft_g2r%npwt) )
     
      CALL cft3t( exx_cus%fft_g2r, psic, exx_cus%fft_g2r%nr1t, exx_cus%fft_g2r%nr2t, exx_cus%fft_g2r%nr3t, &
           &exx_cus%fft_g2r%nrx1t, exx_cus%fft_g2r%nrx2t, exx_cus%fft_g2r%nrx3t, 2 )
      psi_r(1:exx_cus%fft_g2r%nrxxt)= DBLE(psic(1:exx_cus%fft_g2r%nrxxt))



!products with \Psi_v
      do iv=1,exx_cus%nbndv(1)
         prods_r(1:exx_cus%fft_g2r%nrxxt,iv)=psi_r(1:exx_cus%fft_g2r%nrxxt)*&
              &exx_cus%wfc(1:exx_cus%fft_g2r%nrxxt,iv,1)
      enddo

!to G r2G grid
      do iv=1,exx_cus%nbndv(1),2
         if(iv==exx_cus%nbndv(1)) then
            psic(1:exx_cus%fft_r2g%nrxxt)=dcmplx(prods_r(1:exx_cus%fft_r2g%nrxxt,iv),0.d0)
         else
            psic(1:exx_cus%fft_r2g%nrxxt)=dcmplx(prods_r(1:exx_cus%fft_r2g%nrxxt,iv),prods_r(1:exx_cus%fft_r2g%nrxxt,iv+1))
         endif
     
         CALL cft3t( exx_cus%fft_r2g, psic, exx_cus%fft_r2g%nr1t, exx_cus%fft_r2g%nr2t, exx_cus%fft_r2g%nr3t, &
              &exx_cus%fft_r2g%nrx1t, exx_cus%fft_r2g%nrx2t, exx_cus%fft_r2g%nrx3t, -2 )
         if(iv==exx_cus%nbndv(1)) then
            prods_g(1:exx_cus%fft_r2g%npwt, iv) = psic(exx_cus%fft_r2g%nlt(igkt(1:exx_cus%fft_r2g%npwt)))
         else
            prods_g(1:exx_cus%fft_r2g%npwt, iv)= 0.5d0*(psic(exx_cus%fft_r2g%nlt(igkt(1:exx_cus%fft_r2g%npwt)))+&
                 &conjg( psic(exx_cus%fft_r2g%nltm(igkt(1:exx_cus%fft_r2g%npwt)))))
            prods_g(1:exx_cus%fft_r2g%npwt, iv+1)= (0.d0,-0.5d0)*(psic(exx_cus%fft_r2g%nlt(igkt(1:exx_cus%fft_r2g%npwt))) - &
                      &conjg(psic(exx_cus%fft_r2g%nltm(igkt(1:exx_cus%fft_r2g%npwt)))))
         endif
      enddo



!multiply with fac
      do iv=1,exx_cus%nbndv(1)
         prods_g(1:exx_cus%fft_r2g%npwt,iv)=prods_g(1:exx_cus%fft_r2g%npwt,iv)*exx_cus%fac(1:exx_cus%fft_r2g%npwt)
      enddo
!to R r2G grid

      do iv=1,exx_cus%nbndv(1),2
          psic(:)=(0.d0,0.d0)
          if(iv==exx_cus%nbndv(1)) then
             psic(exx_cus%fft_r2g%nlt(1:exx_cus%fft_r2g%npwt))  = prods_g(1:exx_cus%fft_r2g%npwt,iv)
             psic(exx_cus%fft_r2g%nltm(1:exx_cus%fft_r2g%npwt)) = CONJG( prods_g(1:exx_cus%fft_r2g%npwt,iv) )
          else
             psic(exx_cus%fft_r2g%nlt(1:exx_cus%fft_r2g%npwt))=prods_g(1:exx_cus%fft_r2g%npwt,iv)+&
                  &(0.d0,1.d0)*prods_g(1:exx_cus%fft_r2g%npwt,iv+1)
             psic(exx_cus%fft_r2g%nltm(1:exx_cus%fft_r2g%npwt)) = CONJG( prods_g(1:exx_cus%fft_r2g%npwt,iv) )+&
                  &(0.d0,1.d0)*CONJG( prods_g(1:exx_cus%fft_r2g%npwt,iv+1) )
          endif
          CALL cft3t( exx_cus%fft_r2g, psic, exx_cus%fft_r2g%nr1t, exx_cus%fft_r2g%nr2t, exx_cus%fft_r2g%nr3t, &
               &exx_cus%fft_r2g%nrx1t, exx_cus%fft_r2g%nrx2t, exx_cus%fft_r2g%nrx3t, 2 )
          prods_r(1:exx_cus%fft_r2g%nrxxt,iv)= DBLE(psic(1:exx_cus%fft_r2g%nrxxt))
          if(iv/=exx_cus%nbndv(1)) prods_r(1:exx_cus%fft_r2g%nrxxt,iv+1)=&
               &DIMAG(psic(1:exx_cus%fft_r2g%nrxxt))
       enddo



!products with \Psi_v

       do iv=1,exx_cus%nbndv(1)
          prods_r(1:exx_cus%fft_g2r%nrxxt,iv)= prods_r(1:exx_cus%fft_g2r%nrxxt,iv)*exx_cus%wfc(1:exx_cus%fft_g2r%nrxxt,iv,1)
       enddo

!sum up

       prod_tot=0.d0
       do iv=1,exx_cus%nbndv(1)
          prod_tot(1:exx_cus%fft_g2r%nrxxt)=prod_tot(1:exx_cus%fft_g2r%nrxxt)+prods_r(1:exx_cus%fft_g2r%nrxxt,iv)
       enddo

!transform to G space g2r grid


       psic(1:exx_cus%fft_g2r%nrxxt)=dcmplx(prod_tot(1:exx_cus%fft_g2r%nrxxt),0.d0)

       CALL cft3t( exx_cus%fft_g2r, psic, exx_cus%fft_g2r%nr1t, exx_cus%fft_g2r%nr2t, exx_cus%fft_g2r%nr3t, &
              &exx_cus%fft_g2r%nrx1t, exx_cus%fft_g2r%nrx2t, exx_cus%fft_g2r%nrx3t, -2 )

       prod_tot_g(1:exx_cus%fft_g2r%npwt) = psic(exx_cus%fft_g2r%nlt(igkt(1:exx_cus%fft_g2r%npwt)))
       


!put in the order or wfcs
       if(exx_cus%fft_g2r%dual_t==4.d0) then
          xpsi(1:exx_cus%fft_g2r%npwt)=prod_tot_g(1:exx_cus%fft_g2r%npwt)
       else
          call mergewf(prod_tot_g,evc_g,exx_cus%fft_g2r%npwt,exx_cus%fft_g2r%ig_l2gt,mpime,nproc,ionode_id,intra_pool_comm)
          call splitwf(xpsi,evc_g,npw,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
       endif


      deallocate(prods_r,prods_g,prod_tot)
      deallocate(evc_g,psi_t,psi_r,prod_tot_g)
      deallocate(igkt)

      CALL stop_clock('fock')
      return

    END SUBROUTINE fock_cus

    SUBROUTINE fast_vexx(lda, n, m, psi, hpsi,exx_cus,exxalpha,ispin)


      IMPLICIT NONE

      INTEGER          :: lda, n, m, nqi, myrank, mysize
      COMPLEX(DP) :: psi(lda,m) 
      COMPLEX(DP) :: hpsi(lda,m)
      TYPE(exchange_cus) :: exx_cus
      REAL(kind=DP) :: exxalpha
      INTEGER, INTENT(in) :: ispin

      COMPLEX(kind=DP), ALLOCATABLE :: xpsi(:)
      INTEGER :: ii

      allocate(xpsi(lda))
      do ii=1,m
         if(.not.l_exchange_turbo) then
            call fock_cus(psi(1,ii),xpsi,exx_cus)
         else
            call periodic_fock_cus(ispin,psi(1,ii),xpsi,exx_cus)
         endif
         hpsi(1:n,ii)=hpsi(1:n,ii)-exxalpha*xpsi(1:n)
      enddo
      deallocate(xpsi)
      return

    END SUBROUTINE fast_vexx

    FUNCTION exchange_energy_fast(exx_cus,exxalpha)

      USE wvfct,    ONLY : npwx, npw, nbnd
      USE wavefunctions_module, ONLY : evc
      USE mp, ONLY : mp_sum
      USE mp_world, ONLY : world_comm
      USE gvect, ONLY : gstart
      USE io_files, ONLY : prefix, tmp_dir, nwordwfc,iunwfc

      IMPLICIT NONE

      REAL(kind=DP) :: exchange_energy_fast
      TYPE(exchange_cus) :: exx_cus
      REAL(kind=DP) :: exxalpha

      INTEGER :: ii,ig,is
      COMPLEX(kind=DP), ALLOCATABLE :: psi(:,:),xpsi(:)

      exchange_energy_fast=0.d0
      allocate(xpsi(npwx),psi(npwx,nbnd))
      
      do is=1,exx_cus%nspin
         if(exx_cus%nspin==1) then
            psi(1:npw,1:exx_cus%nbndv(is))=evc(1:npw,1:exx_cus%nbndv(is))
         else
            CALL davcio(psi,2*nwordwfc,iunwfc,is,-1)
         endif
         do ii=1,exx_cus%nbndv(is)
            if(.not.l_exchange_turbo) then
               call fock_cus(psi(:,ii),xpsi,exx_cus)
            else
               call periodic_fock_cus(is,psi(:,ii),xpsi,exx_cus)
            endif
            do ig=1,npw
               exchange_energy_fast=exchange_energy_fast+2.d0*dble(psi(ig,ii)*conjg(xpsi(ig)))
            enddo
            if(gstart==2) exchange_energy_fast=exchange_energy_fast-dble(psi(1,ii)*conjg(xpsi(1)))
         enddo
      enddo
      deallocate(xpsi,psi)
      
      call mp_sum(exchange_energy_fast,world_comm)
      if(exx_cus%nspin==1) then
         exchange_energy_fast=-exchange_energy_fast*exxalpha*2.d0!the 2 is for spin ATTENZIONE
      else
         exchange_energy_fast=-exchange_energy_fast*exxalpha
      endif
      return
    END FUNCTION exchange_energy_fast

    subroutine dft_exchange_fast(ispin,nbnd_s,psi,exx_cus)
      USE io_global,            ONLY : stdout, ionode, ionode_id
      USE wvfct,    ONLY : npwx, npw, nbnd, wg
      USE gvect
      USE mp, ONLY : mp_sum
      USE mp_world, ONLY : world_comm
      
      implicit none
      
      INTEGER, INTENT(in) :: ispin!spin channel
      INTEGER, INTENT(in) :: nbnd_s!number of states
      COMPLEX(kind=DP), INTENT(in) :: psi(npwx,nbnd_s)!wavefunctions
      TYPE(exchange_cus), INTENT(in) :: exx_cus!descriptor of exchange

      INTEGER :: ii,ig
      COMPLEX(kind=DP), ALLOCATABLE :: xpsi(:,:)
      REAL(kind=DP) :: sca

      allocate(xpsi(npwx,nbnd_s))
!loop on states

      do ii=1,nbnd_s
!apply X operator
         if(.not.l_exchange_turbo) then
            call fock_cus(psi(:,ii),xpsi(:,ii),exx_cus)
         else
            call periodic_fock_cus(ispin,psi(:,ii),xpsi(:,ii),exx_cus)
         endif
      enddo

!calculate overlap
      do ii=1,nbnd_s
         sca=0.d0
         do ig=1,npw
            sca=sca+2.d0*dble(psi(ig,ii)*conjg(xpsi(ig,ii)))
         enddo
         if(gstart==2) sca=sca-dble(psi(1,ii)*conjg(xpsi(1,ii)))
         call mp_sum(sca,world_comm)
         write(stdout,*) 'EXCHANGE FAST',ii, sca
      enddo
      FLUSH(stdout)
      deallocate(xpsi)
      return
    end subroutine dft_exchange_fast
END MODULE exchange_custom
