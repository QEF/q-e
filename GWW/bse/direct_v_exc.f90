subroutine direct_v_exc(a_in,fc,a_out)
! computes the v part of the direct term of the exc Hamiltonian 
USE fft_custom_gwl
use bse_basic_structures
use exciton
USE wavefunctions_module, ONLY :  psic
USE gvect,                 ONLY : ig_l2g
USE io_global, ONLY : stdout, ionode, ionode_id
USE mp_world, ONLY : mpime, nproc
USE mp_pools, ONLY: intra_pool_comm
USE mp_wave, ONLY : mergewf,splitwf
USE io_global, ONLY : stdout,ionode
USE lsda_mod, ONLY :nspin
USE gvect, ONLY : gstart
USE mp, ONLY : mp_sum
USE mp_world,             ONLY : world_comm



implicit none
type(exc), intent(in) :: a_in
type(exc):: a_out
type(exc_r):: a_in_rt
type(exc_r):: a_tmp_rt
type(ii_mat) :: iimat
type(vww_prod) :: vww
type(fft_cus) :: fc


COMPLEX(kind=DP), allocatable :: vwwg_t(:,:)
COMPLEX(kind=DP), ALLOCATABLE :: evc_g(:)


!real(kind=dp), allocatable :: phivwwr(:,:)
COMPLEX(kind=DP) :: csca

integer ii, iv,ispin, iimax

logical debug

call start_clock('direct_v_exc')
debug=.false.

!if(debug) then 
!   if(ionode) write(stdout,*) 'Direct_v_exc #1'
!endif


! allocate tmp matrix
call initialize_exc_r(a_tmp_rt)
a_tmp_rt%nrxxt=fc%nrxxt 
a_tmp_rt%numb_v=a_in%numb_v
a_tmp_rt%label=12
allocate(a_tmp_rt%ar(a_tmp_rt%nrxxt,a_tmp_rt%numb_v))


! FFT a_in to real space (dual grid)
call initialize_exc_r(a_in_rt)
call fft_a_exc(a_in,fc,a_in_rt)


! read iimat, that tells us for every w_iv which other valence band ivp is
! overlapping, and read the corresponding v*w_iv*w_ivp(G) products

call initialize_imat(iimat)

do ispin=1,nspin
! note that for spin-polarized case, probably this do-loop will have to include
! also the rest of the subroutine, or something like that
   call read_iimat(iimat,ispin)   
enddo


if(debug) then 
   if(ionode) write(stdout,*) 'Direct_v_exc #5'
   if(ionode) write(stdout,*) 'a_in%numb_v=',a_in%numb_v
   if(ionode) write(stdout,*) 'a_in%npw=',a_in%npw
   if(ionode) write(stdout,*) 'iimat%np_max=',iimat%np_max
endif

call initialize_vww_prod(vww)
call read_vww_prod(1,a_in%numb_v,a_in%npw,iimat%np_max,iimat,vww)

!if(debug) then 
!   if(ionode) write(stdout,*) 'Direct_v_exc #6'
!endif

! for every element iv of the excitonic wavefunction vector, here we FFT all
! the available v*w_iv*w_ivp(G) products, multiply by a_in_rt%ar(:,ivp)
! sum over ivp, and FFT back

allocate(vwwg_t(fc%npwt,iimat%np_max))
!allocate(phivwwr(fc%nrxxt,a_in%numb_v))
allocate(evc_g(fc%ngmt_g ))



a_tmp_rt%ar(1:a_tmp_rt%nrxxt,1:a_tmp_rt%numb_v) =0.d0
do iv=1,a_in%numb_v
    
   vwwg_t(1:fc%npwt,1:iimat%np_max)=dcmplx(0.d0,0.d0)
   iimax=0
   do ii=1,iimat%np_max
      if (iimat%iimat(ii,iv)>0) then
         iimax=iimax+1
      else
         exit
      endif
   enddo
   if(iimax>0) then
      call reorderwfp_col(iimax,vww%npw,fc%npwt,vww%vww(1,1,iv),vwwg_t, vww%npw,fc%npwt, &
           & ig_l2g,fc%ig_l2gt,fc%ngmt_g,mpime, nproc,intra_pool_comm )
   endif
!!!!!!!!!!!!!!!!
!   do ii=1,iimat%np_max
!      if (iimat%iimat(ii,iv)==0) exit
!      if(fc%dual_t==4.d0) then
!         vwwg_t(1:fc%npwt,ii)= vww%vww(1:fc%npwt,ii,iv)
!      else
!        call reorderwfp_col(1,vww%npw,fc%npwt,vww%vww(1,ii,iv),vwwg_t(1,ii), vww%npw,fc%npwt, &
!           & ig_l2g,fc%ig_l2gt,fc%ngmt_g,mpime, nproc,intra_pool_comm )
!
!      endif
!
!      
!   enddo

   
   do ii=1,iimat%np_max,2
      if(debug) then
         if(ionode) write(stdout,*) 'ii,iv,iimat', ii,iv,iimat%iimat(ii,iv)
      endif 
      if (iimat%iimat(ii,iv)==0) exit
      psic(:)=(0.d0,0.d0)
      if ((ii==iimat%np_max).or.(iimat%iimat(ii+1,iv)==0)) then
          psic(fc%nlt(1:fc%npwt))  = vwwg_t(1:fc%npwt,ii)
          psic(fc%nltm(1:fc%npwt)) = CONJG( vwwg_t(1:fc%npwt,ii) )
      else      
          psic(fc%nlt(1:fc%npwt))=vwwg_t(1:fc%npwt,ii)+(0.d0,1.d0)*vwwg_t(1:fc%npwt,ii+1)
          psic(fc%nltm(1:fc%npwt))=CONJG(vwwg_t(1:fc%npwt,ii))+(0.d0,1.d0)*CONJG(vwwg_t(1:fc%npwt,ii+1))
      endif
      CALL cft3t( fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, 2 )
      a_tmp_rt%ar(1:fc%nrxxt,iv)= &
            &a_tmp_rt%ar(1:fc%nrxxt,iv)+DBLE(psic(1:fc%nrxxt))*a_in_rt%ar(1:a_in_rt%nrxxt,iimat%iimat(ii,iv))
      if ((ii/=iimat%np_max).and.(iimat%iimat(ii+1,iv)/=0)) then
         a_tmp_rt%ar(1:fc%nrxxt,iv)=&
         &a_tmp_rt%ar(1:fc%nrxxt,iv)+DIMAG(psic(1:fc%nrxxt))*a_in_rt%ar(1:a_in_rt%nrxxt,iimat%iimat(ii+1,iv))
      endif
   enddo
enddo

!if(debug) then 
!   if(ionode) write(stdout,*) 'Direct_v_exc #7'
!endif


call free_memory_exc_a_r(a_in_rt)

call fftback_a_exc(a_tmp_rt,fc,a_out)

! free memory
call free_memory_exc_a_r(a_tmp_rt)
call free_imat(iimat)
call free_vww_prod(vww)
deallocate(vwwg_t)
!deallocate(vwwr_t)
deallocate(evc_g)


call stop_clock('direct_v_exc')
end subroutine


