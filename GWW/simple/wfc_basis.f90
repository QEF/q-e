!this subroutine creates a optimal basis for the periodic functions |u_{nk}>

subroutine wfc_basis
  USE kinds, ONLY : DP
  USE io_files,             ONLY : prefix, iunwfc, nwordwfc,tmp_dir
  USE wavefunctions, ONLY : evc,psic
  USE wvfct, ONLY : nbnd,npw,npwx,et 
  USE mp, ONLY : mp_sum,mp_barrier
  USE klist, ONLY : nks,ngk,xk, igk_k
  USE becmod,        ONLY : bec_type, becp, calbec,allocate_bec_type, deallocate_bec_type
  USE uspp,     ONLY : nkb, vkb, becsum, nhtol, nhtoj, indv, okvan
  USE uspp_param, ONLY : upf, nh
  USE noncollin_module, ONLY: npol, noncolin
  USE mp_world, ONLY : world_comm
  USE spin_orb, ONLY: lspinorb
  USE ions_base,  ONLY : nat, nsp, ityp
  USE io_global, ONLY : stdout, ionode
  USE input_simple
  USE wannier_gw, ONLY : num_nbndv
  USE fft_base,         ONLY : dffts,dfftp
  USE fft_interfaces,ONLY : fwfft, invfft
  USE gvect, ONLY : ngm, gstart,gg, g
  

  implicit none


  COMPLEX(kind=DP), ALLOCATABLE :: omat(:,:),omatj(:),omatij(:,:)
  INTEGER ::  ipol
  TYPE(BEC_TYPE), ALLOCATABLE :: bec0(:),bec1(:),bec2(:)
  INTEGER :: ijkb0,nt,na,jh,ih,i,j,ikb,jkb
  INTEGER :: ik,ig,nfound
  INTEGER :: npw1,ntot,isca
  COMPLEX(kind=DP), ALLOCATABLE :: wfc0(:,:),wfc1(:,:),wfc3(:,:)
  COMPLEX(kind=DP) :: csca
  COMPLEX(kind=DP), EXTERNAL :: ZDOTC

  COMPLEX(kind=DP), ALLOCATABLE :: wfc_t(:,:)
  INTEGER, EXTERNAL :: find_free_unit
  INTEGER :: iun
  COMPLEX(kind=DP), ALLOCATABLE :: rwfc(:,:), cwfc(:),rprod(:),prod_g(:,:)
  real(kind=DP), allocatable :: grids(:),gridd(:)
  INTEGER :: ii,jj
  COMPLEX(kind=DP), ALLOCATABLE :: valbec(:),zmat(:,:),valbecj(:,:,:)
  LOGICAL :: debug=.false.
  INTEGER :: ikk, jkk, kk

  COMPLEX(kind=DP), ALLOCATABLE :: wfc_e_tmp1(:),wfc_e_tmp2(:)
  TYPE(BEC_TYPE), ALLOCATABLE :: bec_e_tmp1(:),bec_e_tmp2(:)

  INTEGER :: iq(3)
  !
  call start_clock('optimal_basis')
  
!determine npw_max
  if(nks>1) then
     !rewind (unit = iunigk)
     !READ( iunigk ) igk
     npw = ngk (1)
     npw_max=maxval(igk_k(1:npw,1))
     do ik=2,nks
        !READ( iunigk ) igk
        npw = ngk (ik)
        isca=maxval(igk_k(1:npw,ik))
        if(isca>npw_max) npw_max=isca
     end do
  else
     npw = ngk (1)
     npw_max=maxval(igk_k(1:npw,1))
     
  endif
  write(stdout,*) 'NPWX NPW_MAX', npwx, npw_max

  

 ! if (nks>1) rewind (unit = iunigk)
   npw = ngk (1)
   !IF ( nks > 1 ) READ( iunigk ) igk
   
   call davcio (evc, 2*nwordwfc, iunwfc, 1, -1)
      

!allocate basis 

   ntot=num_val+num_cond
   ntot_e=num_val+num_cond
   allocate(wfc_t(npol*npwx,ntot_e))
   

!put k=1 wfcs on basis
   wfc_t=(0.d0,0.d0)
   do ipol=0,npol-1
      wfc_t(1+ipol*npwx:npw+ipol*npwx,1:ntot_e)=evc(1+ipol*npwx:npw+ipol*npwx,num_nbndv(1)-num_val+1:num_nbndv(1)+num_cond)
   enddo
   
!check orthonormality
!   allocate(omat(ntot_e,ntot_e))
!   call ZGEMM('C','N',ntot_e,ntot_e,npol*npwx,(1.d0,0.d0),wfc_t,npol*npwx,wfc_t,npol*npwx,(0.d0,0.d0),omat,ntot_e)
!   call mp_sum(omat,world_comm)


   
!   if(ionode) then
!      do i=1,ntot_e
!         !do j=1,ntot_e
!            write(stdout,*) 'K POINT 1, ORTHONORMALITY CHECK: ', i, i, omat(i,i)
         !enddo
!      end do
!   endif
!   deallocate(omat)

   allocate (wfc0(npw_max*npol,ntot_e),wfc1(npw_max*npol,ntot_e))
   wfc0=(0.d0,0.d0)
   do ipol=0,npol-1
      do ig=1,npw
         wfc0(igk_k(ig,1)+ipol*npw_max,1:ntot_e)=evc(ig+ipol*npwx,num_nbndv(1)-num_val+1:num_nbndv(1)+num_cond)
      enddo
   enddo


!loop on k
   do ik=2,nks!DEBUG should start from 2
      write(stdout,*) 'K POINT: ', ik
      call start_clock('wfc_loop')
      
      npw1 = ngk (ik)
      npw = ngk (ik)
      
      call davcio (evc, 2*nwordwfc, iunwfc, ik, -1)
      

      wfc_t=(0.d0,0.d0)
      do ipol=0,npol-1
         wfc_t(1+ipol*npwx:npw+ipol*npwx,1:ntot)=evc(1+ipol*npwx:npw+ipol*npwx,num_nbndv(1)-num_val+1:num_nbndv(1)+num_cond)
      enddo

      wfc1=(0.d0,0.d0)
      do ipol=0,npol-1
         do ig=1,npw1
            wfc1(igk_k(ig,ik)+ipol*npw_max,1:ntot)=evc(ig+ipol*npwx,num_nbndv(1)-num_val+1:num_nbndv(1)+num_cond)
         enddo
      enddo


     
!we need the bec factors only at k-point ik

    
     
       
!!allocate  overlap matrix
       allocate(omat(ntot_e,ntot))

!!at this point wfc_e are expressed by k-states up to ik-1

!!calculate overlap
       call start_clock('zgemm')
       call ZGEMM('C','N',ntot_e,ntot,npol*npw_max,(1.d0,0.d0),wfc0,npol*npw_max,wfc1,npol*npw_max,(0.d0,0.d0),omat,ntot_e)
       call stop_clock('zgemm')
       call mp_sum(omat,world_comm)
       
       
       call mp_barrier(world_comm)
       if(ionode.and.debug) then
          do i=1,ntot_e
             do j=1,ntot
                write(stdout,*) 'K POINT I J, ORTHONORMALITY CHECK: ', i,j, omat(i,j)
             enddo
          enddo
       endif
       
!!project out
       call start_clock('zgemm')
       call ZGEMM('N','N',npw_max*npol,ntot,ntot_e,(-1.d0,0.d0),wfc0,npw_max*npol,omat,ntot_e,(1.d0,0.d0),wfc1,npw_max*npol)
       call stop_clock('zgemm')
!updates bec1
!loop on k-points up to ik

       
 

!DEBUG part check if they are really orthogonal
       !!caluclate overlap
       if(debug) then
          call start_clock('zgemm')
          call ZGEMM('C','N',ntot_e,ntot,npol*npw_max,(1.d0,0.d0),wfc0,npol*npw_max,wfc1,npol*npw_max,(0.d0,0.d0),omat,ntot_e)
          call stop_clock('zgemm')
          call mp_sum(omat,world_comm)
         
          
          
          if(ionode.and.debug) then
             do i=1,ntot_e
                do j=1,ntot
                   write(stdout,*) 'K POINT I J, ORTHONORMALITY CHECK2: ', i,j, omat(i,j)
                enddo
             enddo
          endif
       endif!debug
       
!!orhtonormalize them

!!add to basis and updates arrays and counters
       deallocate(omat)


       call start_clock('wfc_optimal')
       call optimal_gram_schmidt_z(ntot,wfc1,1,s_bands,nfound)
       
       call stop_clock('wfc_optimal')

!update wfc
!
       allocate(wfc3(npw_max*npol,ntot_e+nfound))
       wfc3=(0.d0,0.d0)
       wfc3(1:npw_max*npol,1:ntot_e)=wfc0(1:npw_max*npol,1:ntot_e)
       wfc3(1:npw_max*npol,ntot_e+1:ntot_e+nfound)=wfc1(1:npw_max*npol,1:nfound)
       deallocate(wfc0)
       ntot_e=ntot_e+nfound
       allocate(wfc0(npw_max*npol,ntot_e))
       wfc0=(0.d0,0.0)
       wfc0(1:npw_max*npol,1:ntot_e)=wfc3(1:npw_max*npol,1:ntot_e)
       deallocate(wfc3)
!do the same copying for the bec factors
     
!check orthonormality of basis

       if(debug) then
          allocate(omat(ntot_e,ntot_e))
          call ZGEMM('C','N',ntot_e,ntot_e,npol*npw_max,(1.d0,0.d0),wfc0,npol*npw_max,wfc0,npol*npw_max,(0.d0,0.d0),omat,ntot_e)
          call mp_sum(omat,world_comm)
     
          do ii=1,ntot_e
             do jj=1,ntot_e
                write(stdout,*) 'ORTHONORMALITY :', ik, omat(ii,jj)
             enddo
          enddo
          deallocate(omat)
       endif

!recalculate bec's
       write(stdout,*) 'DIMENSION OF BASIS', ntot_e
       call stop_clock('wfc_loop')
       call print_clock('wfc_loop')
       call print_clock('wfc_optimal')
       call print_clock('zgemm')
    end do

!copy results to common variable (sigh..)
    allocate(wfc_e(npw_max*npol,ntot_e))
    wfc_e(1:npw_max*npol,1:ntot_e)=wfc0(1:npw_max*npol,1:ntot_e)
 
!write data on file
    if(ionode) then
       iun=find_free_unit()
       open( unit= iun, file=trim(tmp_dir)//trim(prefix)//'.wfc_basis', &
            &status='unknown',form='unformatted')
    
!number of k points
       write(iun) nks
!number of valence states to be used
       write(iun) num_val
!number of conduction states to be used
       write(iun) num_cond
       !ntot_e
       write(iun) ntot_e
    endif
!loop on k points
    !if(nks>1)  rewind (unit = iunigk)

    allocate(omat(ntot_e,ntot))
    do ik=1,nks
       write(stdout,*) ' IK', ik
!!k-point xyz
       if(ionode) write(iun)  xk(1:3,ik)
!!weight
       npw = ngk (ik)

       call davcio (evc, 2*nwordwfc, iunwfc, ik, -1)
            
!!energies
       if(ionode) write(iun) et(num_nbndv(1)-num_val+1:num_nbndv(1)+num_cond,ik) 
!!matrix <wfc_e|wfc_ik>
       wfc1=(0.d0,0.d0)
       do ipol=0,npol-1
          do ig=1,npw
             wfc1(igk_k(ig,ik)+ipol*npw_max,1:ntot)=evc(ig+ipol*npwx,num_nbndv(1)-num_val+1:num_nbndv(1)+num_cond)
          enddo
       enddo
       call start_clock('zgemm')
       call ZGEMM('C','N',ntot_e,ntot,npol*npw_max,(1.d0,0.d0),wfc_e,npol*npw_max,wfc1,npol*npw_max,(0.d0,0.d0),omat,ntot_e)
       call stop_clock('zgemm')
       call mp_sum(omat,world_comm)
     
!check for completness
       if(debug) then
          allocate(zmat(ntot,ntot))
          call start_clock('zgemm')
          call ZGEMM('C','N',ntot,ntot,ntot_e,(1.d0,0.d0),omat,ntot_e,omat,ntot_e,(0.d0,0.d0),zmat,ntot)
          call stop_clock('zgemm')
          do ii=1,ntot
             do jj=1,ntot
                write(stdout,*) 'BASIS CHECK ik ii jj', ik,ii,zmat(ii,jj)
             enddo
          enddo
          deallocate(zmat)
       endif
!write omat on disk
       if(ionode) then
          do i=1,ntot
             write(iun) omat(1:ntot_e,i)
          enddo
       endif
!if required  calculates product terms for being checked
       if(debug) then
          ALLOCATE(rwfc(dfftp%nnr*npol,ntot))
          ALLOCATE(cwfc(1:dffts%nnr*npol))
          allocate(grids(dffts%nnr),gridd(dfftp%nnr))
          ALLOCATE(rprod(dfftp%nnr*npol))
          allocate(prod_g(npw_max*npol,ntot))
          allocate(zmat(ntot,ntot))
          do ii=1,ntot
             psic(1:dffts%nnr)=0.d0
             psic(dffts%nl(1:npw_max))=wfc1(1:npw_max,ii)
             CALL invfft ('Wave', psic, dffts)
             cwfc(1:dffts%nnr)=psic(1:dffts%nnr)

             if(npol>1) then
                psic(1:dffts%nnr)=0.d0
                psic(dffts%nl(1:npw_max))=wfc1(npw_max+1:npw_max+npw_max,ii)
                CALL invfft ('Wave', psic, dffts)
                cwfc(dffts%nnr+1:2*dffts%nnr)=psic(1:dffts%nnr)
             endif


             do ipol=0,npol-1
                grids(1:dffts%nnr)=dble(cwfc(1+ipol*dffts%nnr:dffts%nnr+ipol*dffts%nnr))
                call interpolate (gridd, grids, 1)
                rwfc(1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr,ii)=gridd(1:dfftp%nnr)
                grids(1:dffts%nnr)=dimag(cwfc(1+ipol*dffts%nnr:dffts%nnr+ipol*dffts%nnr))
                call interpolate (gridd, grids, 1)
                rwfc(1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr,ii)=rwfc(1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr,ii)&
                     &+(0.d0,1.d0)*gridd(1:dfftp%nnr)
             enddo

          enddo

          deallocate(cwfc,grids,gridd)
       
       


          do ii=1,ntot
             do jj=1,ntot
!form product in real space
               
                rprod(1:dfftp%nnr*npol)=conjg(rwfc(1:dfftp%nnr*npol,ii))*&
                     &rwfc(1:dfftp%nnr*npol,jj)

                do ipol=0,npol-1
                   psic(1:dfftp%nnr)=rprod( 1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr)
                   CALL fwfft ('Rho', psic, dfftp)
                   do ig=1,npw_max
                      prod_g(ig+npw_max*ipol,jj)=psic(dfftp%nl(ig))
                   enddo
                  
                enddo
             enddo
!calculate squared modulus of product
             call start_clock('zgemm')
             call ZGEMM('C','N',ntot,ntot,npw_max*npol,(1.d0,0.d0),prod_g,npw_max*npol,prod_g,npw_max*npol,&
                  &(0.d0,0.d0),zmat,ntot)
             call stop_clock('zgemm')
             call mp_sum(zmat, world_comm)
             if(debug) then
                do jj=1,ntot
                   write(stdout,*) 'DEBUG PRODUCTS ik, i ,j', ik,ii,jj,zmat(jj,jj)
                enddo
             endif
          enddo

!do the same starting from e vectors
          allocate(wfc_e_tmp1(npw_max*npol))
          allocate(wfc_e_tmp2(npw_max*npol))
          allocate(bec_e_tmp1(nks),bec_e_tmp2(nks))
          do ikk=1,nks
             call allocate_bec_type (nkb,1,bec_e_tmp1(ikk))
             call allocate_bec_type (nkb,1,bec_e_tmp2(ikk))
          enddo


          do ii=1,ntot
             wfc_e_tmp1=0.d0
             do kk=1,ntot_e
                wfc_e_tmp1(1:npw_max*npol)=wfc_e_tmp1(1:npw_max*npol)+wfc_e(1:npw_max*npol,kk)*omat(kk,ii)
             enddo
               
                
            
             

             do jj=1,ntot
                wfc_e_tmp2=0.d0
                do kk=1,ntot_e
                   wfc_e_tmp2(1:npw_max*npol)= wfc_e_tmp2(1:npw_max*npol)+wfc_e(1:npw_max*npol,kk)*omat(kk,jj)
                enddo
             
                
                ALLOCATE(cwfc(1:dffts%nnr*npol))
                allocate(grids(dffts%nnr),gridd(dfftp%nnr))
         
         
                psic(1:dffts%nnr)=0.d0
                psic(dffts%nl(1:npw_max))=wfc_e_tmp1(1:npw_max)
                CALL invfft ('Wave', psic, dffts)
                cwfc(1:dffts%nnr)=psic(1:dffts%nnr)

                if(npol>1) then
                   psic(1:dffts%nnr)=0.d0
                   psic(dffts%nl(1:npw_max))=wfc_e_tmp1(npw_max+1:npw_max+npw_max)
                   CALL invfft ('Wave', psic, dffts)
                   cwfc(dffts%nnr+1:2*dffts%nnr)=psic(1:dffts%nnr)
                endif


                do ipol=0,npol-1
                   grids(1:dffts%nnr)=dble(cwfc(1+ipol*dffts%nnr:dffts%nnr+ipol*dffts%nnr))
                   call interpolate (gridd, grids, 1)
                   rwfc(1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr,ii)=gridd(1:dfftp%nnr)
                   grids(1:dffts%nnr)=dimag(cwfc(1+ipol*dffts%nnr:dffts%nnr+ipol*dffts%nnr))
                   call interpolate (gridd, grids, 1)
                   rwfc(1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr,ii)=rwfc(1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr,ii)&
                        &+(0.d0,1.d0)*gridd(1:dfftp%nnr)
                enddo


                psic(1:dffts%nnr)=0.d0
                psic(dffts%nl(1:npw_max))=wfc_e_tmp2(1:npw_max)
                CALL invfft ('Wave', psic, dffts)
                cwfc(1:dffts%nnr)=psic(1:dffts%nnr)

                if(npol>1) then
                   psic(1:dffts%nnr)=0.d0
                   psic(dffts%nl(1:npw_max))=wfc_e_tmp2(npw_max+1:npw_max+npw_max)
                      CALL invfft ('Wave', psic, dffts)
                   cwfc(dffts%nnr+1:2*dffts%nnr)=psic(1:dffts%nnr)
                endif


                do ipol=0,npol-1
                   grids(1:dffts%nnr)=dble(cwfc(1+ipol*dffts%nnr:dffts%nnr+ipol*dffts%nnr))
                   call interpolate (gridd, grids, 1)
                   rwfc(1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr,jj)=gridd(1:dfftp%nnr)
                   grids(1:dffts%nnr)=dimag(cwfc(1+ipol*dffts%nnr:dffts%nnr+ipol*dffts%nnr))
                   call interpolate (gridd, grids, 1)
                   rwfc(1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr,jj)=rwfc(1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr,jj)&
                        &+(0.d0,1.d0)*gridd(1:dfftp%nnr)
                enddo


                rprod(1:dfftp%nnr*npol)=conjg(rwfc(1:dfftp%nnr*npol,ii))*&
                     &rwfc(1:dfftp%nnr*npol,jj)

                do ipol=0,npol-1
                   psic(1:dfftp%nnr)=rprod( 1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr)
                   CALL fwfft ('Rho', psic, dfftp)
                   do ig=1,npw_max
                      prod_g(ig+npw_max*ipol,1)=psic(dfftp%nl(ig))
                   enddo
                                   
                   
                enddo

                csca=ZDOTC(npw_max*npol,prod_g,1,prod_g,1)
                call mp_sum(csca,world_comm)

                if(debug) write(stdout,*) 'DEBUG PRODUCTS2 ik, i ,j', ik,ii,jj,csca




                deallocate(cwfc,grids,gridd)
             enddo
          enddo


          deallocate(wfc_e_tmp1,wfc_e_tmp2)
          deallocate(bec_e_tmp1,bec_e_tmp2)

          deallocate(rwfc,rprod,prod_g,zmat)
    
       endif!debug
    enddo!ik
    write(stdout,*) ' '
    write(stdout,*) 'TOTAL NUMBER OF OPTIMAL BASIS VECTORS :', ntot_e
    write(stdout,*) ' '
!chech for completness of basis
    if(debug) then
       allocate(zmat(ntot_e,ntot_e))
       call start_clock('zgemm')
       call ZGEMM('C','N',ntot_e,ntot_e,npol*npw_max,(1.d0,0.d0),wfc_e,npol*npw_max,wfc_e,npol*npw_max,(0.d0,0.d0),zmat,ntot_e)
       call stop_clock('zgemm')
       call mp_sum(zmat,world_comm)
     
       do ii=1,ntot_e
          do jj=1,ntot_e
             write(stdout,*) 'CHECK OPTIMAL BASIS:', ii,jj, zmat(ii,jj)
          enddo
       enddo
       deallocate(zmat)
    endif!debug
   
    deallocate(omat)
    if(ionode) close(iun)

    

   deallocate(wfc0,wfc1)

  call stop_clock('optimal_basis')

  return


end subroutine wfc_basis




SUBROUTINE optimal_gram_schmidt_z(num_in,wfcs,ithres,thres,num_out)
!this subroutine performs a gram_schmidt orthonormalization and retains
!vectors which are above the give threshold

  USE kinds,                ONLY : DP
  USE mp_world, ONLY : world_comm, mpime, nproc
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE io_global,            ONLY : stdout, ionode,ionode_id
  USE noncollin_module, ONLY: npol, noncolin
  USE input_simple, ONLY : npw_max,vkb_max
  USE becmod,        ONLY : bec_type,calbec,allocate_bec_type, deallocate_bec_type
  USE uspp,     ONLY : nkb, vkb, becsum, nhtol, nhtoj, indv, okvan
  USE uspp_param, ONLY : upf, nh
  USE noncollin_module, ONLY: npol, noncolin
  USE spin_orb, ONLY: lspinorb
  USE ions_base,  ONLY : nat, nsp, ityp



 implicit none

  INTEGER, INTENT(in) :: num_in!number of initial vectors
  COMPLEX(kind=DP), INTENT(inout) :: wfcs(npw_max*npol,num_in)!in input non-orthonormal in output optimal basis
  INTEGER, INTENT(in) :: ithres!kind of threshold
  REAL(kind=DP), INTENT(in) :: thres!thrshold for the optimal basis
  INTEGER, INTENT(out) :: num_out!final number of orthonormal basis functions



  INTEGER :: i,j
  COMPLEX(kind=DP), ALLOCATABLE :: prod(:)
  COMPLEX(kind=DP) :: csca
  COMPLEX(kind=DP), EXTERNAL :: zdotc
  REAL(kind=DP) :: sca

  TYPE(BEC_TYPE) :: bec0
  TYPE(BEC_TYPE), ALLOCATABLE :: bec1(:)
  INTEGER :: ijkb0,nt,na,jh,ih,ikb,jkb,ipol



  allocate(prod(num_in))
  num_out=0


  do i=1,num_in
     if(num_out >0) then
        call zgemv('C',npw_max*npol,num_out,(1.d0,0.d0), wfcs,npw_max*npol,wfcs(1,i),1,(0.d0,0.d0),prod,1)
        call mp_sum(prod(1:num_out),world_comm)
        call start_clock('zgemm')
        call zgemm('N','N',npw_max*npol,1,num_out,(-1.d0,0.d0),wfcs,npw_max*npol,prod,num_in,(1.d0,0.d0),wfcs(1,i),npw_max*npol)
        call stop_clock('zgemm')
     endif
     csca = zdotc(npw_max*npol,wfcs(1,i),1,wfcs(1,i),1)
     call mp_sum(csca,world_comm)



     if(dble(csca) >= thres) then
        num_out=num_out+1
        sca=dsqrt(dble(csca))
        wfcs(1:npw_max*npol,num_out)=wfcs(1:npw_max*npol,i)/sca
     endif
  enddo


  deallocate(prod)
  return
END SUBROUTINE optimal_gram_schmidt_z
SUBROUTINE debug_wfc(npw)
   integer :: npw
   call errore('debug_wfc','not implemented',abs(npw))
END SUBROUTINE debug_wfc



