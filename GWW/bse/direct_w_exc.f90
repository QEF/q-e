subroutine direct_w_exc(a_in,fc,a_out)
! this subroutine computes the w part of the direct term of the exc Hamiltonian
 
USE fft_custom_gwl
use bse_basic_structures
use exciton
USE gvect
use bse_wannier, ONLY: l_truncated_coulomb, &
           truncation_radius
USE constants,        ONLY : e2, fpi
USE cell_base,        ONLY : tpiba,omega,tpiba2
!USE io_files,             ONLY : find_free_unit, prefix, diropn
USE io_files,             ONLY :  prefix, diropn
USE wavefunctions_module, ONLY :  psic
USE io_global, ONLY : stdout, ionode, ionode_id
USE mp_world, ONLY : mpime, nproc
USE mp_pools, ONLY: intra_pool_comm
USE mp_wave, ONLY : mergewf,splitwf
USE polarization
USE lsda_mod, ONLY :nspin
USE io_global, ONLY : stdout,ionode
USE mp,          ONLY :mp_barrier
USE mp_world,             ONLY : world_comm





implicit none
INTEGER, EXTERNAL :: find_free_unit
type(bse_z) :: z
type(polaw) :: pw
type(exc):: a_in
type(exc):: a_out
type(exc_r):: a_in_rt
type(exc_r):: a_tmp_rt
type(fft_cus) :: fc
type(ii_mat) :: iimat



REAL(kind=DP), ALLOCATABLE :: fac(:)
COMPLEX(kind=DP), ALLOCATABLE :: p_basis(:,:)
COMPLEX(kind=DP), ALLOCATABLE :: p_basis_t(:,:)
REAL(kind=DP), ALLOCATABLE :: p_basis_r(:,:)
REAL(kind=DP), ALLOCATABLE :: zvphi(:)
REAL(kind=dp), ALLOCATABLE :: zvv(:)
COMPLEX(kind=DP), ALLOCATABLE :: evc_g(:)

INTEGER ::iungprod
INTEGER :: ig,ii,iv,ispin
REAL(kind=DP) :: qq
LOGICAL       :: exst


INTEGER :: vpmax,k
REAL(kind=DP), allocatable :: zp(:,:)
REAL(kind=DP), allocatable :: pizeta(:,:)
REAL(kind=DP), allocatable :: vphipizeta(:,:)

logical debug

call start_clock('direct_w_exc')

debug=.false.



! read  iimat
call initialize_imat(iimat)

do ispin=1,nspin
   call read_iimat(iimat,ispin) 
enddo

! read z terms
call initialize_bse_z(z)
call read_z(1,iimat,z)


! get Coulomb potential
allocate(fac(a_in%npw))
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
endif


! read polarization basis and multiply per V

iungprod = find_free_unit()
allocate(p_basis(a_in%npw,z%numw_prod))
CALL diropn( iungprod, 'wiwjwfc_red', a_in%npw*2, exst )

do ii=1,z%numw_prod
   call davcio(p_basis(:,ii),a_in%npw*2,iungprod,ii,-1)
   p_basis(1:a_in%npw,ii)=p_basis(1:a_in%npw,ii)*dcmplx(fac(1:a_in%npw))
enddo

call mp_barrier(world_comm)

close(iungprod)

! FFT to real space (dual grid)
allocate(p_basis_t(fc%npwt,z%numw_prod)) 
allocate(p_basis_r(fc%nrxxt,z%numw_prod))
allocate(evc_g(fc%ngmt_g ))

if(fc%dual_t==4.d0) then
   p_basis_t(1:fc%npwt,1:z%numw_prod)=p_basis(1:a_in%npw,1:z%numw_prod)
else
    call reorderwfp_col(z%numw_prod,a_in%npw,fc%npwt,p_basis(1,1),p_basis_t(1,1),a_in%npw,fc%npwt, &
           & ig_l2g,fc%ig_l2gt,fc%ngmt_g,mpime, nproc,intra_pool_comm )

!   do ii=1,z%numw_prod
!      call mergewf(p_basis(:,ii),evc_g,a_in%npw,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
!      call splitwf(p_basis_t(:,ii),evc_g,fc%npwt,fc%ig_l2gt,mpime,nproc,ionode_id,intra_pool_comm)
!   enddo
endif

deallocate(evc_g)
deallocate(p_basis)

call start_clock('direct_w_cft3t')
do ii=1,z%numw_prod,2
   psic(1:fc%nrxxt)=(0.d0,0.d0)
   if (ii==z%numw_prod) then
      psic(fc%nlt(1:fc%npwt))  = p_basis_t(1:fc%npwt,ii)
      psic(fc%nltm(1:fc%npwt)) = CONJG( p_basis_t(1:fc%npwt,ii) )
   else
      psic(fc%nlt(1:fc%npwt))=p_basis_t(1:fc%npwt,ii)+(0.d0,1.d0)*p_basis_t(1:fc%npwt,ii+1)
      psic(fc%nltm(1:fc%npwt))=CONJG(p_basis_t(1:fc%npwt,ii))+(0.d0,1.d0)*CONJG(p_basis_t(1:fc%npwt,ii+1))
   endif
      CALL cft3t( fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, 2 )
       p_basis_r(1:fc%nrxxt,ii)= DBLE(psic(1:fc%nrxxt))
      if(ii/=z%numw_prod) p_basis_r(1:fc%nrxxt,ii+1)= DIMAG(psic(1:fc%nrxxt))
enddo
call stop_clock('direct_w_cft3t')

deallocate(p_basis_t) 

!read P
call initialize_polaw(pw)
call read_polaw_global(0, pw)


call mp_barrier(world_comm)

!FFT the input excitonic vector to real space (dual grid)
call initialize_exc_r(a_in_rt)
call fft_a_exc(a_in,fc,a_in_rt)

call mp_barrier(world_comm)

! allocate tmp matrix
call initialize_exc_r(a_tmp_rt)
a_tmp_rt%nrxxt=fc%nrxxt 
a_tmp_rt%numb_v=a_in%numb_v
a_tmp_rt%label=12
allocate(a_tmp_rt%ar(a_tmp_rt%nrxxt,a_tmp_rt%numb_v))

!compute line by line the output excitonic vector

!!!!!!!!!!!!!!!!!dgemm subroutine!!!!!!!!!!!!!!!!!!!!!
call start_clock('direct_w_dgemv')
allocate(zp(z%numw_prod,iimat%np_max))
!!allocate(pizeta(z%numw_prod,iimat%np_max)) 
!!allocate(vphipizeta(fc%nrxxt,iimat%np_max)) 
!
a_tmp_rt%ar(1:a_tmp_rt%nrxxt,1:a_tmp_rt%numb_v) =0.d0
!
do iv=1, a_in%numb_v
   zp(1:z%numw_prod,1:iimat%np_max)=0.d0
   vpmax=0
 
   call start_clock('dgemv1')
   do ii=1, iimat%np_max
      if (iimat%iimat(ii,iv)==0) cycle
      vpmax=vpmax+1
      do k=1, z%numw_prod
         zp(k,ii)=z%z(k,ii,iv)
      enddo
   enddo

   allocate(pizeta(z%numw_prod,vpmax)) 
   allocate(vphipizeta(fc%nrxxt,vpmax)) 

   call stop_clock('dgemv1')


   call start_clock('dgemv2')
   call dgemm('N','N', z%numw_prod,vpmax, z%numw_prod,1.d0,pw%pw,z%numw_prod,zp(1,1),&
              z%numw_prod,0.d0,pizeta(1,1),z%numw_prod)
   call stop_clock('dgemv2')


   call start_clock('dgemv3')
   call dgemm('N','N', fc%nrxxt, vpmax, z%numw_prod,1.d0,p_basis_r(1,1),fc%nrxxt,&
              pizeta(1,1), z%numw_prod, 0.d0, vphipizeta(1,1),fc%nrxxt)
   call stop_clock('dgemv3')

!!     sum up
   call start_clock('dgemv4')
   do ii=1,vpmax
      a_tmp_rt%ar(1:fc%nrxxt,iv)= &
            &a_tmp_rt%ar(1:fc%nrxxt,iv)+a_in_rt%ar(1:fc%nrxxt,iimat%iimat(ii,iv))*vphipizeta(1:fc%nrxxt,ii)
   enddo
!
   call stop_clock('dgemv4')
!
!
   deallocate(pizeta)
   deallocate(vphipizeta)
call mp_barrier(world_comm)
enddo
call stop_clock('direct_w_dgemv')

deallocate(zp)

!!!!!!!!!!!!!!!!!!dgemv subroutine!!!!!!!!!!!!!!!!!!!!!
!call start_clock('direct_w_dgemv')
!allocate(zvv(z%numw_prod))
!allocate(zvphi(fc%nrxxt))
!a_tmp_rt%ar(1:a_tmp_rt%nrxxt,1:a_tmp_rt%numb_v) =0.d0
!!
!do iv=1, a_in%numb_v
!   do ii=1, iimat%np_max !ii is ivp
!      if (iimat%iimat(ii,iv)==0) exit
!      call start_clock('dgemv1')
!      call dgemv('N',z%numw_prod,z%numw_prod,1.d0,pw%pw,z%numw_prod,z%z(1,ii,iv),1,0.d0,zvv,1)
!      call stop_clock('dgemv1')
!      call start_clock('dgemv2')
!      call dgemv('N',fc%nrxxt,z%numw_prod,1.d0,p_basis_r,fc%nrxxt,zvv,1,0.d0,zvphi,1)
!      call stop_clock('dgemv2')
!!     sum up
!      call start_clock('dgemv3')
!      a_tmp_rt%ar(1:fc%nrxxt,iv)= &
!            &a_tmp_rt%ar(1:fc%nrxxt,iv)+a_in_rt%ar(1:a_in_rt%nrxxt,iimat%iimat(ii,iv))*zvphi(1:fc%nrxxt)
!      call stop_clock('dgemv3')
!   enddo
!enddo
!call stop_clock('direct_w_dgemv')

!deallocate(zvv)
!deallocate(zvphi)

call free_memory_exc_a_r(a_in_rt)
call free_bse_z(z)
call free_memory_polaw(pw)
call free_imat(iimat)

call start_clock('wdirect_fftback')
!FFT back to provide the output excitonic wave vector in G-space
call fftback_a_exc(a_tmp_rt,fc,a_out)
call stop_clock('wdirect_fftback')

call free_memory_exc_a_r(a_tmp_rt)

call stop_clock('direct_w_exc')

return
end subroutine

