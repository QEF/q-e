!this routine forms a common basis for products of wave functions
!it write it to disk
!it forms product with basis for wave-functions and write it to disk

subroutine product_basis
  USE kinds, ONLY : DP
  USE io_files,             ONLY : prefix, iunwfc, nwordwfc, tmp_dir
  USE wavefunctions_module, ONLY : psic
  USE wvfct, ONLY : npw,npwx,et
  USE mp, ONLY : mp_sum
  USE klist, ONLY : nks,ngk,xk
  USE becmod,        ONLY : bec_type, becp, calbec,allocate_bec_type, deallocate_bec_type
  USE uspp,     ONLY : nkb, vkb, becsum, nhtol, nhtoj, indv, okvan
  USE uspp_param, ONLY : upf, nh
  USE noncollin_module, ONLY: npol, noncolin
  USE mp_world, ONLY : world_comm
  USE spin_orb, ONLY: lspinorb
  USE ions_base,  ONLY : nat, nsp, ityp
  USE io_global, ONLY : stdout, ionode
  USE input_simple
  USE gvect, ONLY : ngm, gstart,gg, g
  USE fft_base,         ONLY : dffts,dfftp
  USE fft_interfaces,ONLY : fwfft, invfft
  



  implicit none

  INTEGER :: ii,jj,ipol,ig,ir,num,num_out,iun
  COMPLEX(kind=DP), ALLOCATABLE :: rwfc(:,:), cwfc(:),rprod(:),prod_g(:,:)
  real(kind=DP), allocatable :: grids(:),gridd(:)
  COMPLEX(kind=DP), ALLOCATABLE :: valbec(:),gprod(:),sprod(:),valbecm(:,:,:,:,:,:)
  LOGICAL, PARAMETER :: ldebug=.false.
  COMPLEX(kind=DP):: csca
  COMPLEX(kind=DP), ALLOCATABLE :: omat(:,:),prod_tmp(:,:),zmat(:,:),omatij(:,:),omatm(:,:)
  INTEGER, EXTERNAL :: find_free_unit
  TYPE(BEC_TYPE) :: bec1
  INTEGER :: ikk,jkk,jjj
  INTEGER :: iq(3)
  INTEGER, PARAMETER :: block_states=1
  INTEGER :: nblock,nks_prod
  write(stdout,*) 'Routine product_basis'


!put wfc in real space without augmentation part

  ALLOCATE(rwfc(dfftp%nnr*npol,ntot_e))
  ALLOCATE(cwfc(1:dffts%nnr*npol))
  allocate(grids(dffts%nnr),gridd(dfftp%nnr))
  ALLOCATE(rprod(dfftp%nnr*npol))


  nks_prod=(2*nkpoints(1)+1)*(2*nkpoints(2)+1)*(2*nkpoints(3)+1)
  rwfc=(0.d0,0.d0)
  do ii=1,ntot_e
     psic(1:dffts%nnr)=0.d0
     psic(dffts%nl(1:npw_max))=wfc_e(1:npw_max,ii)
     CALL invfft ('Wave', psic, dffts)
     cwfc(1:dffts%nnr)=psic(1:dffts%nnr)

     if(npol>1) then
        psic(1:dffts%nnr)=0.d0
        psic(dffts%nl(1:npw_max))=wfc_e(npw_max+1:npw_max+npw_max,ii)
        CALL invfft ('Wave', psic, dffts)
        cwfc(dffts%nnr+1:2*dffts%nnr)=psic(1:dffts%nnr)
     endif


     do ipol=0,npol-1
!        grids(1:dffts%nnr)=dble(cwfc(1+ipol*dffts%nnr:dffts%nnr+ipol*dffts%nnr))
        !call interpolate (gridd, grids, 1)
        !rwfc(1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr,ii)=gridd(1:dfftp%nnr) ! DEBUG ATTENZIONE DA COMPLETARE     
       ! grids(1:dffts%nnr)=dimag(cwfc(1+ipol*dffts%nnr:dffts%nnr+ipol*dffts%nnr))
       !  call interpolate (gridd, grids, 1)
       !  rwfc(1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr,ii)=rwfc(1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr,ii)&
       !       &+(0.d0,1.d0)*gridd(1:dfftp%nnr)

        rwfc(1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr,ii)=cwfc(1+ipol*dffts%nnr:dffts%nnr+ipol*dffts%nnr)
   enddo

   enddo

   deallocate(cwfc,grids,gridd)

 


!calculate bec's for US case
!IMPORTANT BEC_E MUST BE ALREADY ALLOCATED AND PROPERLY INITIALISED



!loop on wfc states
   nprod_e=0
   allocate(prod_g(npw_max*npol,ntot_e))
    
   do ii=1,ntot_e
      call start_clock('Product_ciclo')
      write(stdout,*) 'Forming products with ', ii
      do jjj=1,ntot_e,block_states
!!form products with all the remaining i>=j
         do jj=jjj,min(jjj+block_states-1,ntot_e)
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
         


      enddo
  
!!project out basis already built
!!gram schmidt with threshold

!check for consistency
      if(ldebug) then
         ALLOCATE(cwfc(1:dffts%nnr*npol))
         allocate(sprod(dffts%nnr*npol))
          !do jj=1,ntot_e-ii+1
             do jj=1,ntot_e
             psic(1:dffts%nnr)=0.d0
             psic(dffts%nl(1:npw_max))=prod_g(1:npw_max,jj)
             CALL invfft ('Wave', psic, dffts)
             cwfc(1:dffts%nnr)=psic(1:dffts%nnr)

             if(npol>1) then
                psic(1:dffts%nnr)=0.d0
                psic(dffts%nl(1:npw_max))=prod_g(npw_max+1:npw_max+npw_max,jj)
                CALL invfft ('Wave', psic, dffts)
                cwfc(dffts%nnr+1:2*dffts%nnr)=psic(1:dffts%nnr)
             endif
             sprod(1:dffts%nnr*npol)=cwfc(1:dffts%nnr*npol)

             csca=(0.d0,0.d0)
             do ir=1,dffts%nnr*npol
                csca=csca+sprod(ir)
             enddo
             call mp_sum(csca, world_comm)
             csca=csca/(dffts%nr1*dffts%nr2*dffts%nr3)
             write(stdout,*) 'VERIFICA :' , ii, jj, csca
          enddo

          deallocate(sprod,cwfc)
       endif
       
       write(stdout,*) 'Projecting out'
!!project out basis already built
       num=ntot_e
       if(nprod_e>0) then
          allocate(omat(nprod_e,num))
          call start_clock('Product_zgemm')
          call ZGEMM('C','N',nprod_e,num,npw_max*npol,(1.0d0,0.d0),prod_e,npw_max*npol,prod_g,npw_max*npol,(0.d0,0.d0),omat,nprod_e)
          call stop_clock('Product_zgemm')
          call mp_sum(omat,world_comm)
          call start_clock('Product_zgemm')
          call ZGEMM('N','N',npw_max*npol,num,nprod_e,(-1.d0,0.d0),prod_e,npw_max*npol,omat,nprod_e,(1.d0,0.d0),prod_g,npw_max*npol)
          call stop_clock('Product_zgemm')
          deallocate(omat)
       endif
       write(stdout,*) 'Calling Gram Schmidt'
!gram schmidt with threshold
       call start_clock('Product_gram')
       call optimal_gram_schmidt_nc(num,prod_g,s_product,num_out)
       call stop_clock('Product_gram')
       if(num_out>0 ) then
          if(nprod_e>0) then 
             allocate(prod_tmp(npw_max*npol,nprod_e))
             prod_tmp(1:npw_max*npol,1:nprod_e)=prod_e(1:npw_max*npol,1:nprod_e)
             deallocate(prod_e)
             allocate(prod_e(npw_max*npol,nprod_e+num_out))
             prod_e(1:npw_max*npol,1:nprod_e)=prod_tmp(1:npw_max*npol,1:nprod_e)
             prod_e(1:npw_max*npol,nprod_e+1:nprod_e+num_out)=prod_g(1:npw_max*npol,1:num_out)
             nprod_e=nprod_e+num_out
             deallocate(prod_tmp)
          else
             allocate(prod_e(npw_max*npol,nprod_e+num_out))
             prod_e(1:npw_max*npol,1:num_out)=prod_g(1:npw_max*npol,1:num_out)
             nprod_e=num_out
          endif
       endif
       write(stdout, *) 'PRODUCT BASIS: TOT NEW ',nprod_e,num_out
       call stop_clock('Product_ciclo')
       call print_clock('Product_ciclo')
       call print_clock('Product_zgemm')
       call print_clock('Product_uszgemm')
       call print_clock('Product_gram')
       call print_clock('values_mat')
    enddo! on ii
!check for orthonormality
    if(ldebug) then
  
       allocate(zmat(nprod_e,nprod_e))
       call ZGEMM('C','N',nprod_e,nprod_e,npol*npw_max,(1.d0,0.d0),prod_e,npol*npw_max,prod_e,npol*npw_max,&
            &(0.d0,0.d0),zmat,nprod_e)
       call mp_sum(zmat,world_comm)
       do ii=1,nprod_e
          do jj=1,nprod_e
             write(stdout,*) 'PRODUCT BESIS ORHTONORMALITY:', ii,jj,zmat(ii,jj)
          enddo
       enddo
  
       deallocate(zmat)
    endif
!form products of product basis with couple of wfc basis
!consider cubic matrices for simplicity
!and write them on disk
   if(ionode) then
      iun=find_free_unit()
      open( unit= iun, file=trim(tmp_dir)//trim(prefix)//'.product_basis', &
            &status='unknown',form='unformatted')
!product number
      write(iun) nprod_e
!wfc number
      write(iun) ntot_e
   endif


    allocate(omat(nprod_e,ntot_e))
    do ii=1,ntot_e
       do jjj=1,ntot_e,block_states
!!form products with all the remaining i>=j 
          do jj=jjj,min(jjj+block_states-1,ntot_e)
          
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
          
       enddo

       if(ldebug) then
          do jj=1,ntot_e
             csca=0.d0
             do ig=1,npw_max*npol
                csca=csca+conjg(prod_g(ig,jj))*prod_g(ig,jj)
             enddo
             call mp_sum(csca, world_comm)
             write(stdout,*) 'DEBUG E* E PRODUCT', ii, jj, csca
          enddo
       endif
       

       call ZGEMM('C','N',nprod_e,ntot_e,npw_max*npol,(1.d0,0.d0),prod_e,npw_max*npol,prod_g,npw_max*npol,(0.d0,0.d0),omat,nprod_e)
       call mp_sum(omat,world_comm)
       if(ionode) then
          do jj=1,ntot_e
             write(iun) omat(1:nprod_e,jj)
          enddo
       endif
    enddo

    
    if(ionode) close(iun)
    deallocate(omat)
  

   deallocate(rprod,prod_g)
   deallocate(rwfc)

 end subroutine product_basis



!!To be updated in a norm-conserving way

SUBROUTINE optimal_gram_schmidt_nc(num_in,wfcs,thres,num_out)
!this subroutine performs a gram_schmidt orthonormalization and retains
!vectors which are above the give threshold

  USE kinds,                ONLY : DP
  USE mp_world, ONLY : world_comm, mpime, nproc
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE io_global,            ONLY : stdout, ionode,ionode_id
  USE noncollin_module, ONLY: npol, noncolin
  USE input_simple, ONLY : npw_max,vkb_max
  USE noncollin_module, ONLY: npol, noncolin
  USE spin_orb, ONLY: lspinorb
  USE ions_base,  ONLY : nat, nsp, ityp



 implicit none

  INTEGER, INTENT(in) :: num_in!number of initial vectors
  COMPLEX(kind=DP), INTENT(inout) :: wfcs(npw_max*npol,num_in)!in input non-orthonormal in output optimal basis
  REAL(kind=DP), INTENT(in) :: thres!threshold for the optimal basis          
  INTEGER, INTENT(out) :: num_out!final number of orthonormal basis functions


  INTEGER :: i,j
  COMPLEX(kind=DP), ALLOCATABLE :: prod(:)
  COMPLEX(kind=DP) :: csca
  COMPLEX(kind=DP), EXTERNAL :: zdotc
  REAL(kind=DP) :: sca


  num_out=0


  allocate(prod(num_in))
  do i=1,num_in
     if(num_out >0) then

        call zgemv('C',npw_max*npol,num_out,(1.d0,0.d0), wfcs,npw_max*npol,wfcs(1,i),1,(0.d0,0.d0),prod,1)
        call mp_sum(prod(1:num_out),world_comm)

        call zgemm('N','N',npw_max*npol,1,num_out,(-1.d0,0.d0),wfcs,npw_max*npol,prod,num_in,(1.d0,0.d0),wfcs(1,i),npw_max*npol)
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
END SUBROUTINE optimal_gram_schmidt_nc



