  !subroutine for the matrix elements <e_i|p|e_j>
subroutine epe

  use kinds, only : DP
  use input_simple
  use mp, only : mp_sum,mp_barrier
  use mp_world, only : world_comm
  use io_global, only : ionode,stdout
  use io_files, only : prefix, tmp_dir
  USE klist, ONLY : nks,ngk,xk
  USE becmod,        ONLY : bec_type, becp, calbec,allocate_bec_type, deallocate_bec_type
  USE uspp,     ONLY : nkb, vkb, becsum, nhtol, nhtoj, indv, okvan
  USE uspp_param, ONLY : upf, nh
  USE noncollin_module, ONLY: npol, noncolin
  USE wavefunctions_module, ONLY : psic
  USE fft_base,         ONLY : dffts,dfftp
  USE fft_interfaces,ONLY : fwfft, invfft
  USE gvect, ONLY : ngm, gstart,gg, g
  USE cell_base, ONLY : tpiba


  implicit none
  !wfc_e(npw_max,ntot_e)
  INTEGER, EXTERNAL :: find_free_unit
  integer :: i,j,a,iun
  complex(kind=DP), allocatable :: omat(:,:,:),ge(:,:,:)
  COMPLEX(kind=DP), ALLOCATABLE :: rwfc(:), cwfc(:),rprod(:),prod_g(:)
  real(kind=DP), allocatable :: grids(:),gridd(:)
  INTEGER :: ipol, ii, ig, ikk, jkk, ir
  COMPLEX(kind=DP), ALLOCATABLE ::  phase(:),wfc_us_e(:,:)
  REAL(kind=DP) :: sca
  INTEGER :: iq(3),nks_prod
 

  write(stdout,*)'Subroutine epe'

  allocate(omat(ntot_e,ntot_e,3))
  allocate(ge(npw_max*npol,ntot_e,3))

  !Neglect non-local terms
  !applico G ad e

  if(.not.okvan) then


     ge=0.d0
     do a=1,3!3 direzioni di G
        do ipol=1,npol
           do ig=1,npw_max
              do i=1,ntot_e
                 ge(ig+(ipol-1)*npw_max,i,a) = g(a,ig)*wfc_e(ig+(ipol-1)*npw_max,i)*2.d0*tpiba
              end do
           end do
        enddo
     end do
     do a=1,3
        call zgemm('C','N',ntot_e,ntot_e,npw_max*npol,(1.0d0,0.0d0),wfc_e,npw_max*npol,ge(:,:,a),&
                      &  npw_max*npol,(0.0d0,0.0d0),omat(:,:,a),ntot_e)
     end do
     call mp_sum(omat,world_comm)
     
  else
 !we use the following approximation also for SOC
 !first we calculate |\psi_i(r)|^2=\rho_i(r)
 !then \psi_i'(r)=sqrt(rho(r))*eps(i\phi(r)) with the same phase of the smooth part of \psi_i
 !COULD BE USED EVEN FOR NC 

     allocate(rwfc(dfftp%nnr*npol))
     allocate(cwfc(1:dffts%nnr*npol))
     allocate(grids(dffts%nnr),gridd(dfftp%nnr))
     allocate(rprod(dfftp%nnr*npol))
     allocate(prod_g(npw_max*npol))
     allocate(phase(dfftp%nnr*npol))
     allocate(wfc_us_e(npw_max*npol,ntot_e))

     !call set_ntotbec
    
     
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
           grids(1:dffts%nnr)=dble(cwfc(1+ipol*dffts%nnr:dffts%nnr+ipol*dffts%nnr))
           !call interpolate (gridd, grids, 1) DEBUG ATTENZIONE DA COMPLETARE
           rwfc(1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr)=gridd(1:dfftp%nnr)
           grids(1:dffts%nnr)=dimag(cwfc(1+ipol*dffts%nnr:dffts%nnr+ipol*dffts%nnr))
           !call interpolate (gridd, grids, 1) DEBUG ATTENZIONE DA COMPLETARE     
           rwfc(1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr)=rwfc(1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr)&
                &+(0.d0,1.d0)*gridd(1:dfftp%nnr)
        enddo
        
        rprod(1:dfftp%nnr*npol)=conjg(rwfc(1:dfftp%nnr*npol))*&
             &rwfc(1:dfftp%nnr*npol)

        do ipol=0,npol-1
           psic(1:dfftp%nnr)=rprod( 1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr)
           CALL fwfft ('Rho', psic, dfftp)
           do ig=1,npw_max
              prod_g(ig+npw_max*ipol)=psic(dfftp%nl(ig))
           enddo
        enddo
       
     
!determine phase
        do ipol=0,npol-1
           do ir=1,dfftp%nnr
              sca=sqrt(dble(conjg(rwfc(ir+ipol*dfftp%nnr))*rwfc(ir+ipol*dfftp%nnr)))
              if(sca>=1d-10) then
                 phase(ir+ipol*dfftp%nnr)=rwfc(ir+ipol*dfftp%nnr)/sca
              else
                 phase(ir+ipol*dfftp%nnr)=(1.d0,0.d0)
              endif
           enddo
        enddo


 !fft back 


        do ipol=0,npol-1
           psic(1:dfftp%nnr)=0.d0
           psic(dfftp%nl(1:npw_max))=prod_g(npw_max+1:npw_max+npw_max)
           CALL invfft ('Wave', psic, dfftp)
           rprod(dfftp%nnr*ipol+1:dfftp%nnr+dfftp%nnr*ipol)=psic(1:dfftp%nnr)
           
        enddo
!add phase
        do ipol=0,npol-1
           do ir=1,dfftp%nnr
              sca=sqrt(abs(dble(rprod(ir+ipol*dfftp%nnr))))
              rprod(ir+ipol*dfftp%nnr)=phase(ir+ipol*dfftp%nnr)*sca
           enddo
        enddo
!fft forth
        do ipol=0,npol-1
           psic(1:dfftp%nnr)=rprod( 1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr)
           CALL fwfft ('Rho', psic, dfftp)
           do ig=1,npw_max
              wfc_us_e(ig+npw_max*ipol,ii)=psic(dfftp%nl(ig))
           enddo
        enddo


     enddo


     do a=1,3!3 direzioni di G
        do i=1,ntot_e
           do ipol=0,npol-1
              ge(1+ipol*npw_max:npw_max+ipol*npw_max,i,a)=g(a,1:npw_max)*wfc_us_e(1+ipol*npw_max:npw_max+ipol*npw_max,i)&
                   &*2.d0*tpiba
           enddo
        enddo
     enddo
       
     do a=1,3
        call zgemm('C','N',ntot_e,ntot_e,npw_max*npol,(1.0d0,0.0d0),wfc_us_e,npw_max*npol,ge(1,1,a),&
              & npw_max*npol,(0.0d0,0.0d0), omat(1,1,a),ntot_e)
     end do
     call mp_sum(omat,world_comm)



     deallocate(cwfc,grids,gridd)
     deallocate(rwfc,prod_g,rprod)

    
     deallocate(phase,wfc_us_e)


  endif

  if(ionode) then
     iun=find_free_unit()
     open( unit= iun, file=trim(tmp_dir)//trim(prefix)//'.epe', &
             &status='unknown',form='unformatted')
     write(iun) ntot_e
     do a=1,3
        do i=1,ntot_e
           write(iun) omat(1:ntot_e,i,a)
        end do
     end do
     
     close(iun)
  end if


  deallocate(ge)
  deallocate(omat)
  return
end subroutine epe
