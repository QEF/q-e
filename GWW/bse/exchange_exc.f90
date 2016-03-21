subroutine exchange_exc(a_in,v,v_rt,fc,a_out)
! this subroutine applies the exchange term of the Hamiltonian to the excitonic
! wavefuntion vector a 
!if l_gtrick==.true. uses vonly and does not use v_rt     

USE kinds,            ONLY : DP
USE exciton
use bse_basic_structures
use bse_wannier,      ONLY: l_truncated_coulomb, &
                             truncation_radius, &
                             dual_bse
USE gvect
USE constants,        ONLY : e2, fpi
USE cell_base,        ONLY : tpiba,omega,tpiba2
!USE mp_wave,          ONLY : mergewf,splitwf
USE fft_custom_gwl
USE gvecw,              ONLY : ecutwfc
USE io_global, ONLY : stdout, ionode, ionode_id
USE mp_world, ONLY : mpime, nproc
USE mp_pools, ONLY: intra_pool_comm
USE mp_world,             ONLY : world_comm
USE wavefunctions_module, ONLY :  psic
USE mp,          ONLY :mp_barrier

!USE io_files,         ONLY : find_free_unit


implicit none
INTEGER, EXTERNAL :: find_free_unit
type(exc) :: a_in,a_out
type(exc_r) :: a_rt
type(v_state) :: v
type(v_state_r) :: v_rt
type(fft_cus) :: fc


REAL(kind=DP), ALLOCATABLE :: fac(:)
COMPLEX(kind=DP), ALLOCATABLE :: fac_t(:)
COMPLEX(kind=DP), ALLOCATABLE :: cfac(:)
COMPLEX(kind=DP), ALLOCATABLE :: evc_g(:)
REAL(kind=DP) :: qq
integer ig,iv,iunu
COMPLEX(kind=DP), allocatable :: psiv_phiv(:)

logical :: debug

call start_clock('exchange_exc')
debug=.false.

if(debug) then 
  write(*,*) 'Starting to compute the exchange term'
  call mp_barrier(world_comm)
endif

allocate(fac(a_in%npw))
allocate(cfac(a_in%npw))

if(l_truncated_coulomb) then
   do ig=1,a_in%npw
     qq = g(1,ig)**2.d0 + g(2,ig)**2.d0 + g(3,ig)**2.d0
     if (qq > 1.d-8) then
        fac(ig)=(e2*fpi/(tpiba2*qq))*(1.d0-dcos(dsqrt(qq)*truncation_radius*tpiba))
     else
        fac(ig)=e2*fpi*(truncation_radius**2.d0/2.d0)
     endif
   enddo
   fac(:)=fac(:)/omega
else
   fac(:)=0.d0
   fac(1:a_in%npw)=vg_q(1:a_in%npw)
   !set the 0 component of v to 0 to have vbar (see S. Albrecht phD thesis or Onida,Reining,and Rubio RMP)
   if (v%gstart==2) fac(1)=0.d0
endif

cfac(1:a_in%npw)=dcmplx(fac(1:a_in%npw))


if(debug) then 
  write(*,*) 'vbar built'
  call mp_barrier(world_comm)
endif

! now distribuite the G vector in the dual grid order

allocate(fac_t(fc%npwt))
allocate(evc_g(fc%ngmt_g))
call reorderwfp_col(1,a_in%npw,fc%npwt,cfac(1),fac_t(1),a_in%npw,fc%npwt, &
     & ig_l2g,fc%ig_l2gt,fc%ngmt_g,mpime, nproc,intra_pool_comm )

!call  mergewf(cfac,evc_g,a_in%npw,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
!call  splitwf(fac_t,evc_g,fc%npwt,fc%ig_l2gt,mpime,nproc,ionode_id,intra_pool_comm)
deallocate(evc_g)


if(debug) then 
  write(*,*) 'vbar distributed'
  call mp_barrier(world_comm)
endif

!FFT the excitonic wavefunction vector to real space (dual grid)
call initialize_exc_r(a_rt)
call fft_a_exc(a_in,fc,a_rt)

if(debug) then 
  write(*,*) 'fft exc state'
  call mp_barrier(world_comm)
endif
!compute the psi_iv(r)*phi_iv(r) product in r-space (dual grid) for every iv 
!FFT back to G space (dual grid order), multiply with vbar, and sum over iv
 
!be careful! for the moment there is no spin structure in the a_in_r object so this works
!only for nspin==1

allocate(psiv_phiv(fc%npwt))

psiv_phiv(1:fc%npwt)=(0.d0,0.d0)

do iv=1,a_in%numb_v,2
   if(debug) then 
     write(*,*) 'exchange term main loop iv=',iv
     write(*,*) 'a_in%numb_v=',a_in%numb_v
     write(*,*) 'a_rt%nrxxt=',a_rt%nrxxt
!     call mp_barrier
   endif
   a_rt%ar(1:a_rt%nrxxt,iv) = a_rt%ar(1:a_rt%nrxxt,iv)*v_rt%wfnrt(1:v_rt%nrxxt,iv,1)
   if(iv/=a_in%numb_v) a_rt%ar(1:a_rt%nrxxt,iv+1) = a_rt%ar(1:a_rt%nrxxt,iv+1)*v_rt%wfnrt(1:v_rt%nrxxt,iv+1,1)
   if (iv==a_in%numb_v) then 
      psic(1:fc%nrxxt)=dcmplx(a_rt%ar(1:fc%nrxxt,iv),0.d0)
   else
      psic(1:fc%nrxxt)=dcmplx(a_rt%ar(1:fc%nrxxt,iv),a_rt%ar(1:fc%nrxxt,iv+1))
   endif 
   if(debug) then 
     write(*,*) 'before fft'
     write(*,*) 'fc%nr1t',fc%nr1t
     write(*,*) 'fc%nr2t',fc%nr2t
     write(*,*) 'fc%nr3t',fc%nr3t
     write(*,*) 'fc%nrx1t',fc%nrx1t
     write(*,*) 'fc%nrx2t',fc%nrx2t
     write(*,*) 'fc%nrx3t',fc%nrx3t
!     call mp_barrier
   endif
   CALL cft3t(fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, -2 )
   if(debug) then 
     write(*,*) 'after fft'
     call mp_barrier(world_comm)
   endif
   if (iv==a_in%numb_v) then 
      psiv_phiv(1:fc%npwt)= psiv_phiv(1:fc%npwt)+psic(fc%nlt(1:fc%npwt))*fac_t(1:fc%npwt)
   else
      psiv_phiv(1:fc%npwt)= psiv_phiv(1:fc%npwt)+&
           &0.5*(psic(fc%nlt(1:fc%npwt))+conjg( psic(fc%nltm(1:fc%npwt))))*fac_t(1:fc%npwt)+&
           &(0.d0,-0.5d0)*(psic(fc%nlt(1:fc%npwt))-conjg(psic(fc%nltm(1:fc%npwt))))*fac_t(1:fc%npwt)
   endif
   if(debug) then 
     write(*,*) 'end of loop'
     call mp_barrier(world_comm)
   endif
enddo

if(fc%gstart_t==2) psiv_phiv(1)= (0.d0,0.d0)

!FFT to real space and multiply by valence state wavefunction vector to create
!the components of the excitonic vector and then FFT back 

psic(:)=0.d0
psic(fc%nlt(1:fc%npwt))  = psiv_phiv(1:fc%npwt)
psic(fc%nltm(1:fc%npwt)) = CONJG(psiv_phiv(1:fc%npwt))

CALL cft3t( fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, 2 )

do iv=1,a_in%numb_v
    a_rt%ar(1:a_rt%nrxxt,iv)=v_rt%wfnrt(1:v_rt%nrxxt,iv,1)*dble(psic(1:fc%nrxxt))
enddo

call fftback_a_exc(a_rt,fc,a_out)

! project into the conduction band manifold

call pc_operator_exc(a_out,v,1)

deallocate(fac_t)
call free_memory_exc_a_r(a_rt)

call stop_clock('exchange_exc')

return
end subroutine

