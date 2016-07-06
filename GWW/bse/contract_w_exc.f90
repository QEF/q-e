MODULE contract_w

USE kinds, ONLY : DP,sgl
use bse_basic_structures, ONLY : ii_mat

SAVE

type(ii_mat) :: iimat_contract
REAL(kind=sgl), POINTER  :: vphipizeta_save(:,:)
REAL(kind=sgl), POINTER  :: vww_save(:,:)
COMPLEX(kind=sgl), POINTER  :: vphipizeta_save_g(:,:)
COMPLEX(kind=sgl), POINTER  :: vww_save_g(:,:)
INTEGER, POINTER ::  vpmax_ii(:),vpmax_ii_start(:),vpmax_ii_end(:)
INTEGER :: vpmax_tot

CONTAINS

subroutine free_memory_contrac_w
  use bse_wannier, ONLY : l_gtrick
  implicit none
  if(.not.l_gtrick) then
     deallocate(vphipizeta_save)
  else
     deallocate(vphipizeta_save_g)
  endif
  deallocate(vpmax_ii)
  deallocate(vpmax_ii_start,vpmax_ii_end)
  if(.not.l_gtrick) then
     deallocate(vww_save)
  else
     deallocate(vww_save_g)
  endif
end subroutine free_memory_contrac_w

subroutine contract_w_build(fc)
! this subroutine computes the w part of the direct term of the exc Hamiltonian

USE wvfct,                 ONLY :  npw 
USE fft_custom_gwl
use bse_basic_structures
use exciton
USE gvect
use bse_wannier, ONLY: l_truncated_coulomb, &
           truncation_radius, l_gtrick
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
logical :: debug=.false.

type(bse_z) :: z
type(polaw) :: pw

type(fft_cus) :: fc
type(ii_mat) :: iimat



REAL(kind=DP), ALLOCATABLE :: fac(:)
COMPLEX(kind=DP), ALLOCATABLE :: p_basis(:,:)
COMPLEX(kind=DP), ALLOCATABLE :: p_basis_t(:,:)
REAL(kind=DP), ALLOCATABLE :: p_basis_r(:,:)
REAL(kind=DP), ALLOCATABLE :: zvphi(:)
REAL(kind=DP), ALLOCATABLE :: zvv(:)
COMPLEX(kind=DP), ALLOCATABLE :: evc_g(:)

INTEGER ::iungprod
INTEGER :: ig,ii,iv,ispin
REAL(kind=DP) :: qq
LOGICAL       :: exst


INTEGER :: vpmax,k
REAL(kind=DP), allocatable :: zp(:,:)
REAL(kind=DP), allocatable :: pizeta(:,:)
REAL(kind=DP), allocatable :: vphipizeta(:,:)

INTEGER ::  kilobytes


call start_clock('direct_w_exc')
!CALL memstat( kilobytes )
!write(stdout,*) 'memory0', kilobytes
!FLUSH(stdout)




! read  iimat
call initialize_imat(iimat)

do ispin=1,nspin
   call read_iimat(iimat,ispin) 
enddo

! read z terms
call initialize_bse_z(z)
call read_z(1,iimat,z)

FLUSH( stdout ) 

! get Coulomb potential
allocate(fac(npw))
if(l_truncated_coulomb) then
   do ig=1,npw
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
   fac(1:npw)=vg_q(1:npw)
endif


! read polarization basis and multiply per V

iungprod = find_free_unit()
allocate(p_basis(npw,z%numw_prod))
CALL diropn( iungprod, 'wiwjwfc_red', npw*2, exst )

do ii=1,z%numw_prod
   call davcio(p_basis(:,ii),npw*2,iungprod,ii,-1)
   p_basis(1:npw,ii)=p_basis(1:npw,ii)*dcmplx(fac(1:npw))
enddo

call mp_barrier(world_comm)

close(iungprod)
!CALL memstat( kilobytes )
!write(stdout,*) 'memory1', kilobytes
!FLUSH(stdout)

! FFT to real space (dual grid)
allocate(p_basis_t(fc%npwt,z%numw_prod)) 
allocate(p_basis_r(fc%nrxxt,z%numw_prod))
allocate(evc_g(fc%ngmt_g ))


if(fc%dual_t==4.d0) then
   p_basis_t(1:fc%npwt,1:z%numw_prod)=p_basis(1:npw,1:z%numw_prod)
else
    call reorderwfp_col(z%numw_prod,npw,fc%npwt,p_basis(1,1),p_basis_t(1,1),npw,fc%npwt, &
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

!CALL memstat( kilobytes )
!write(stdout,*) 'memory2', kilobytes
!FLUSH(stdout)


! allocate tmp matrix


!compute line by line the output excitonic vector

!!!!!!!!!!!!!!!!!dgemm subroutine!!!!!!!!!!!!!!!!!!!!!
call start_clock('direct_w_dgemv')
!write(stdout,*) 'memory2',z%numw_prod,iimat%np_max
allocate(zp(z%numw_prod,iimat%np_max))
!!allocate(pizeta(z%numw_prod,iimat%np_max)) 
!!allocate(vphipizeta(fc%nrxxt,iimat%np_max)) 
!calculate vpmax_tot
vpmax_tot=0
allocate(vpmax_ii( iimat%numb_v))
allocate(vpmax_ii_start( iimat%numb_v))
allocate(vpmax_ii_end( iimat%numb_v))
do iv=1, iimat%numb_v

   vpmax=0
   vpmax_ii_start(iv)=vpmax_tot+1
   do ii=1, iimat%np_max
      if (iimat%iimat(ii,iv)==0) cycle
      vpmax=vpmax+1
   enddo
   vpmax_ii(iv)=vpmax
   vpmax_tot=vpmax_tot+vpmax
   vpmax_ii_end(iv)=vpmax_tot
enddo
write(stdout,*)  'VPHIZETA_SAVE :', fc%nrxxt,vpmax_tot
FLUSH(stdout)
if(.not.l_gtrick) then
   allocate(vphipizeta_save(fc%nrxxt,vpmax_tot))
else
   allocate(vphipizeta_save_g(fc%npwt,vpmax_tot))
endif
write(stdout,*)  'VPHIZETA_SAVE :', fc%nrxxt,vpmax_tot
!
!CALL memstat( kilobytes )
!write(stdout,*) 'memory3', kilobytes
!FLUSH(stdout)
do iv=1, iimat%numb_v
   zp(1:z%numw_prod,1:iimat%np_max)=0.d0
   vpmax=0
 
   call start_clock('dgemv1')
   do ii=1, iimat%np_max
      if (iimat%iimat(ii,iv)==0) cycle
      vpmax=vpmax+1
      do k=1, z%numw_prod
         zp(k,vpmax)=z%z(k,ii,iv)!ATTENZIONE era zp(ii)
      enddo
   enddo

  if(vpmax>0) then
     allocate(pizeta(z%numw_prod,vpmax)) 
     allocate(vphipizeta(fc%nrxxt,vpmax)) 
  else
     allocate(pizeta(z%numw_prod,1))
     allocate(vphipizeta(fc%nrxxt,1))
  endif
   call stop_clock('dgemv1')


   call start_clock('dgemv2')


   call dgemm('N','N', z%numw_prod,vpmax, z%numw_prod,1.d0,pw%pw,z%numw_prod,zp,&
              z%numw_prod,0.d0,pizeta,z%numw_prod)
   call stop_clock('dgemv2')


   call start_clock('dgemv3')
   call dgemm('N','N', fc%nrxxt, vpmax, z%numw_prod,1.d0,p_basis_r(1,1),fc%nrxxt,&
              pizeta(1,1), z%numw_prod, 0.d0, vphipizeta(1,1),fc%nrxxt)
   call stop_clock('dgemv3')
   if(.not.l_gtrick) then
      vphipizeta_save(1:fc%nrxxt,vpmax_ii_start(iv):vpmax_ii_end(iv))= real(vphipizeta(1:fc%nrxxt,1:vpmax))
   else
      do ii=1,vpmax,2
         if(ii==vpmax) then
            psic(1:fc%nrxxt)=dcmplx(vphipizeta(1:fc%nrxxt,ii),0.d0)
         else
            psic(1:fc%nrxxt)=dcmplx(vphipizeta(1:fc%nrxxt,ii),vphipizeta(1:fc%nrxxt,ii+1))
         endif
         CALL cft3t(fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, -2 )
         if(ii==vpmax) then
            vphipizeta_save_g(1:fc%npwt,vpmax_ii_start(iv)+ii-1)=cmplx(psic(fc%nlt(1:fc%npwt)))
         else
            vphipizeta_save_g(1:fc%npwt,vpmax_ii_start(iv)+ii-1)=&
                 &cmplx(0.5d0*(psic(fc%nlt(1:fc%npwt))+conjg( psic(fc%nltm(1:fc%npwt)))))
            vphipizeta_save_g(1:fc%npwt,vpmax_ii_start(iv)+ii-1+1)=&
                 &cmplx((0.d0,-0.5d0)*(psic(fc%nlt(1:fc%npwt)) - conjg(psic(fc%nltm(1:fc%npwt)))))
         endif
      enddo
   endif

   !
!
   deallocate(pizeta)
   deallocate(vphipizeta)
call mp_barrier(world_comm)
enddo
call stop_clock('direct_w_dgemv')

deallocate(zp)
deallocate(p_basis_r)

call free_bse_z(z)
call free_memory_polaw(pw)
call free_imat(iimat)


!CALL memstat( kilobytes )
!write(stdout,*) 'memory4', kilobytes
!FLUSH(stdout)
call stop_clock('direct_w_exc')

return
end subroutine contract_w_build



subroutine contract_w_apply(a_in,fc,a_out)
! this subroutine computes the w part of the direct term of the exc Hamiltonian
 
USE fft_custom_gwl
use bse_basic_structures
use exciton
USE gvect
use bse_wannier, ONLY: l_truncated_coulomb, &
           truncation_radius,l_gtrick
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

type(polaw) :: pw
type(exc):: a_in
type(exc):: a_out
type(exc_r):: a_in_rt
type(exc_r):: a_tmp_rt
type(fft_cus) :: fc





INTEGER ::iungprod
INTEGER :: ig,ii,iv,ispin
REAL(kind=DP) :: qq
LOGICAL       :: exst


INTEGER ::k


logical debug

call start_clock('direct_w_contract')

debug=.false.
! read  iimat

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
call start_clock('contract_w_dgemv')


!
a_tmp_rt%ar(1:a_tmp_rt%nrxxt,1:a_tmp_rt%numb_v) =0.d0
!
do iv=1, a_in%numb_v

   if(.not.l_gtrick) then
      do ii=1,vpmax_ii(iv)
         a_tmp_rt%ar(1:fc%nrxxt,iv)= &
              &a_tmp_rt%ar(1:fc%nrxxt,iv)+a_in_rt%ar(1:fc%nrxxt,iimat_contract%iimat(ii,iv))*&
              &dble(vphipizeta_save(1:fc%nrxxt,vpmax_ii_start(iv)-1+ii))
      enddo
   else
      do ii=1,vpmax_ii(iv),2
         psic(1:fc%nrxxt)=(0.d0,0.d0)
         if (ii==vpmax_ii(iv)) then
            psic(fc%nlt(1:fc%npwt))  = dcmplx(vphipizeta_save_g(1:fc%npwt,vpmax_ii_start(iv)-1+ii))
            psic(fc%nltm(1:fc%npwt)) = dcmplx(CONJG( vphipizeta_save_g(1:fc%npwt,vpmax_ii_start(iv)-1+ii) ))
         else
            psic(fc%nlt(1:fc%npwt))= dcmplx(vphipizeta_save_g(1:fc%npwt,vpmax_ii_start(iv)-1+ii))+&
                 &(0.0,1.0)* dcmplx(vphipizeta_save_g(1:fc%npwt,vpmax_ii_start(iv)-1+ii+1))
            psic(fc%nltm(1:fc%npwt))=DCMPLX(CONJG( vphipizeta_save_g(1:fc%npwt,vpmax_ii_start(iv)-1+ii))+&
                 &(0.0,1.0)*CONJG( vphipizeta_save_g(1:fc%npwt,vpmax_ii_start(iv)-1+ii+1)))
         endif
         CALL cft3t( fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, 2 )
         a_tmp_rt%ar(1:fc%nrxxt,iv)= &
              &a_tmp_rt%ar(1:fc%nrxxt,iv)+a_in_rt%ar(1:fc%nrxxt,iimat_contract%iimat(ii,iv))*&
              &dble(psic(1:fc%nrxxt))
         if (ii/=vpmax_ii(iv)) then
             a_tmp_rt%ar(1:fc%nrxxt,iv)= &
                  &a_tmp_rt%ar(1:fc%nrxxt,iv)+a_in_rt%ar(1:fc%nrxxt,iimat_contract%iimat(ii+1,iv))*&
                  &dimag(psic(1:fc%nrxxt))

         endif
      enddo
   endif



enddo
call stop_clock('contract_w_dgemv')




call free_memory_exc_a_r(a_in_rt)



call start_clock('wdirect_fftback')
!FFT back to provide the output excitonic wave vector in G-space
call fftback_a_exc(a_tmp_rt,fc,a_out)
call stop_clock('wdirect_fftback')

call free_memory_exc_a_r(a_tmp_rt)

FLUSH( stdout )
call stop_clock('direct_w_contract')

return
end subroutine contract_w_apply


subroutine contract_v_build(fc)
!TO BE CALLED AFTER CONTRACT_W_BUILD  
!IMPORTANT: REQUIRES IIMAT_CONTRACT
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
USE wvfct,                 ONLY :  npw
USE bse_wannier, ONLY : l_gtrick

implicit none


type(vww_prod) :: vww
type(fft_cus) :: fc


COMPLEX(kind=DP), allocatable :: vwwg_t(:,:)
COMPLEX(kind=DP), ALLOCATABLE :: evc_g(:)


!real(kind=dp), allocatable :: phivwwr(:,:)
COMPLEX(kind=DP) :: csca

integer ii, iv,ispin, iimax

logical debug

write(stdout,*) 'Routine contract_v_build'

debug=.false.

!if(debug) then 
!   if(ionode) write(stdout,*) 'Direct_v_exc #1'
!endif
write(stdout,*) 'VWW_SAVE :', fc%nrxxt,vpmax_tot
if(.not.l_gtrick) then
   allocate(vww_save(fc%nrxxt,vpmax_tot))
else
   allocate(vww_save_g(fc%npwt,vpmax_tot))
endif


call initialize_vww_prod(vww)
call read_vww_prod(1,iimat_contract%numb_v,npw,iimat_contract%np_max,iimat_contract,vww)

!if(debug) then 
!   if(ionode) write(stdout,*) 'Direct_v_exc #6'
!endif

! for every element iv of the excitonic wavefunction vector, here we FFT all
! the available v*w_iv*w_ivp(G) products, multiply by a_in_rt%ar(:,ivp)
! sum over ivp, and FFT back

allocate(vwwg_t(fc%npwt,iimat_contract%np_max))
!allocate(phivwwr(fc%nrxxt,a_in%numb_v))
allocate(evc_g(fc%ngmt_g ))


write(stdout,*) 'ATT1'
do iv=1,iimat_contract%numb_v
    
   vwwg_t(1:fc%npwt,1:iimat_contract%np_max)=dcmplx(0.d0,0.d0)
   iimax=0
   do ii=1,iimat_contract%np_max
      if (iimat_contract%iimat(ii,iv)>0) then
         iimax=iimax+1
      else
         exit
      endif
   enddo
   if(iimax>0) then
      call reorderwfp_col(iimax,vww%npw,fc%npwt,vww%vww(1,1,iv),vwwg_t, vww%npw,fc%npwt, &
           & ig_l2g,fc%ig_l2gt,fc%ngmt_g,mpime, nproc,intra_pool_comm )
   endif
   
   if(.not.l_gtrick) then
      do ii=1,vpmax_ii(iv),2
        
      
         psic(1:fc%nrxxt)=(0.d0,0.d0)
         if (ii==vpmax_ii(iv)) then
            psic(fc%nlt(1:fc%npwt))  = vwwg_t(1:fc%npwt,ii)
            psic(fc%nltm(1:fc%npwt)) = CONJG( vwwg_t(1:fc%npwt,ii) )
         else      
            psic(fc%nlt(1:fc%npwt))=vwwg_t(1:fc%npwt,ii)+(0.d0,1.d0)*vwwg_t(1:fc%npwt,ii+1)
            psic(fc%nltm(1:fc%npwt))=CONJG(vwwg_t(1:fc%npwt,ii))+(0.d0,1.d0)*CONJG(vwwg_t(1:fc%npwt,ii+1))
         endif
         CALL cft3t( fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, 2 )
         vww_save(1:fc%nrxxt,vpmax_ii_start(iv)+ii-1)=real(psic(1:fc%nrxxt))
         if (ii/=vpmax_ii(iv)) then
            vww_save(1:fc%nrxxt,vpmax_ii_start(iv)+ii)=aimag(psic(1:fc%nrxxt))
         endif

      enddo
   else
      do ii=1,vpmax_ii(iv)
         vww_save_g(1:fc%npwt,vpmax_ii_start(iv)+ii-1)=cmplx(vwwg_t(1:fc%npwt,ii))
      end do
   endif
enddo

!if(debug) then 
!   if(ionode) write(stdout,*) 'Direct_v_exc #7'
!endif

FLUSH( stdout )


call free_vww_prod(vww)
deallocate(vwwg_t)
!deallocate(vwwr_t)
deallocate(evc_g)


end subroutine


subroutine contract_v_apply(a_in,fc,a_out)
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
USE bse_wannier, ONLY : l_gtrick


implicit none
type(exc), intent(in) :: a_in
type(exc):: a_out
type(exc_r):: a_in_rt
type(exc_r):: a_tmp_rt


type(fft_cus) :: fc


COMPLEX(kind=DP), allocatable :: vwwg_t(:,:)
COMPLEX(kind=DP), ALLOCATABLE :: evc_g(:)

COMPLEX(kind=DP) :: csca

integer ii, iv,ispin, iimax

logical debug

call start_clock('direct_v_contract')
debug=.false.

call initialize_exc_r(a_tmp_rt)
a_tmp_rt%nrxxt=fc%nrxxt
a_tmp_rt%numb_v=a_in%numb_v
a_tmp_rt%label=12
allocate(a_tmp_rt%ar(a_tmp_rt%nrxxt,a_tmp_rt%numb_v))


! FFT a_in to real space (dual grid)
call initialize_exc_r(a_in_rt)
call fft_a_exc(a_in,fc,a_in_rt)





a_tmp_rt%ar(1:a_tmp_rt%nrxxt,1:a_tmp_rt%numb_v) =0.d0
do iv=1,a_in%numb_v
  

   if(.not.l_gtrick) then
      do ii=1,vpmax_ii(iv)
         a_tmp_rt%ar(1:fc%nrxxt,iv)= &
              &a_tmp_rt%ar(1:fc%nrxxt,iv)+DBLE(vww_save(1:fc%nrxxt,vpmax_ii_start(iv)+ii-1))*&
              &a_in_rt%ar(1:a_in_rt%nrxxt,iimat_contract%iimat(ii,iv))
      enddo
   else
      do ii=1,vpmax_ii(iv),2
         psic(1:fc%nrxxt)=(0.d0,0.d0)
         if (ii==vpmax_ii(iv)) then
            psic(fc%nlt(1:fc%npwt))  = dcmplx(vww_save_g(1:fc%npwt,vpmax_ii_start(iv)+ii-1))
            psic(fc%nltm(1:fc%npwt)) = dcmplx(CONJG(  vww_save_g(1:fc%npwt,vpmax_ii_start(iv)+ii-1) ))
         else
            psic(fc%nlt(1:fc%npwt))=dcmplx(vww_save_g(1:fc%npwt,vpmax_ii_start(iv)+ii-1)+&
                 &(0.0,1.0)*vww_save_g(1:fc%npwt,vpmax_ii_start(iv)+ii-1+1))
            psic(fc%nltm(1:fc%npwt))=dcmplx(CONJG(vww_save_g(1:fc%npwt,vpmax_ii_start(iv)+ii-1))&
                 &+(0.0,1.0)*CONJG(vww_save_g(1:fc%npwt,vpmax_ii_start(iv)+ii-1+1)))
         endif
         call start_clock('d_v_fft')
         CALL cft3t( fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, 2 )
         call stop_clock('d_v_fft')
         a_tmp_rt%ar(1:fc%nrxxt,iv)= &
              &a_tmp_rt%ar(1:fc%nrxxt,iv)+DBLE(psic(1:fc%nrxxt))*&
              &a_in_rt%ar(1:a_in_rt%nrxxt,iimat_contract%iimat(ii,iv))
         if (ii/=vpmax_ii(iv)) then
            a_tmp_rt%ar(1:fc%nrxxt,iv)= &
                 &a_tmp_rt%ar(1:fc%nrxxt,iv)+DIMAG(psic(1:fc%nrxxt))*&
                 &a_in_rt%ar(1:a_in_rt%nrxxt,iimat_contract%iimat(ii+1,iv))
         endif
      enddo
   endif
enddo



call free_memory_exc_a_r(a_in_rt)

call fftback_a_exc(a_tmp_rt,fc,a_out)

! free memory
call free_memory_exc_a_r(a_tmp_rt)



call stop_clock('direct_v_contract')
end subroutine







END MODULE contract_w
