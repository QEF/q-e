subroutine find_eig(vstate,vstate_r,cstate,wcstate,fc)
!this subroutine finds the lowest n_eig eigenvectors and eigenvalues 
!through the  conjugate gradient or steepest descent minimization scheme
!for each eigenvector it computes the optical amplitude and eventually the
!excitonic wavefunction

USE exciton
USE io_global, ONLY : stdout,ionode

USE bse_basic_structures
USE fft_custom_gwl
USE wvfct,    ONLY : npw,npwx,nbnd
USE bse_wannier, ONLY:num_nbndv,eps,lambda,maxit,n_eig,l_cgrad,l_plotexc,&
                      plotn_min,plotn_max,l_plotaverage,l_restart,l_verbose,&
                      n_eig_start,l_finite,l_contraction,l_gtrick,l_dielectric
USE mp,          ONLY : mp_barrier
USE mp_world,             ONLY : world_comm
USE constants,   ONLY : RYTOEV
USE contract_w
USE lsda_mod, ONLY :nspin
!USE eqv,         ONLY : dpsi, dvpsi, eprec



implicit none
!type(exc), allocatable :: bse_spectrum(:)
type(v_state) :: vstate
type(v_state_r) :: vstate_r
type(c_state) :: cstate
type(c_state) :: wcstate
type(fft_cus) :: fc

type(spectrum) :: bse_sp

type(v_state) :: psibar(3) ! formula (43) of Rev. Mod. Phys. 73, 515 
                           ! for the three polarization directions computed through DFPT 
integer :: i,nstart,ipol,ispin
complex(DP), allocatable:: dvpsi(:,:)

real(DP) :: EVTORY



EVTORY=1.d0/RYTOEV

call start_clock('find_eig')

write(stdout,*) 'Routine find_eig'
FLUSH(stdout)

call initialize_spectrum(bse_sp)

if(l_contraction) then
   write(stdout,*) 'CALL contract_w_build'
   FLUSH(stdout)
   call contract_w_build(fc)
endif
bse_sp%neig=n_eig
allocate(bse_sp%en(bse_sp%neig))
allocate(bse_sp%a(bse_sp%neig,3))

allocate(bse_spectrum(n_eig))

if(l_contraction) then
! read  iimat
   call initialize_imat(iimat_contract)
   do ispin=1,nspin
      call read_iimat(iimat_contract,ispin)
   enddo
   write(stdout,*) 'CALL contract_v_build'
   FLUSH(stdout)
   call contract_v_build(fc)
endif


do i=1,n_eig
   call initialize_exc(bse_spectrum(i))
   bse_spectrum(i)%npw=npw
   bse_spectrum(i)%numb_v=num_nbndv(1)
   allocate(bse_spectrum(i)%a(bse_spectrum(i)%npw,bse_spectrum(i)%numb_v)) 
   bse_spectrum(i)%label=i
enddo




if(l_restart==1)then
   nstart=n_eig_start
   do i=1,n_eig_start-1
      call read_exc(i, bse_spectrum(i),l_verbose)
      bse_sp%en(i)=bse_spectrum(i)%e*EVTORY
   enddo
else
   nstart=1
endif



!compute the eigenfunction and eigenvalues
if(l_restart<2) then
   do i=nstart,n_eig
      if(l_cgrad) then
         call conjgrad(i,vstate,vstate_r,cstate,wcstate,fc,bse_sp%en(i)) 
      else
         call sdescent(i,vstate,vstate_r,cstate,wcstate,fc,bse_sp%en(i))
      endif
      call write_exc(bse_spectrum(i))
   enddo
else if(l_restart==2) then
   do i=1,n_eig
      call read_exc(i, bse_spectrum(i),l_verbose)
      bse_sp%en(i)=bse_spectrum(i)%e*EVTORY
   enddo
endif


call mp_barrier(world_comm)

if(l_gtrick)  call v_wfng_to_wfnr(vstate,fc,vstate_r)

!compute the optical amplitudes 

if(l_dielectric) then
!compute the  |psibar(iv)>
   if(.not.l_finite) then
      allocate (dvpsi ( npwx , num_nbndv(1)))
      do ipol=1,3 
         call initialize_v_state(psibar(ipol))
         psibar(ipol)%nspin= vstate%nspin
         psibar(ipol)%numb_v(:)=vstate%numb_v(:)
         psibar(ipol)%npw=npw
         psibar(ipol)%gstart=vstate%gstart
         
         allocate( psibar(ipol)%wfn(psibar(ipol)%npw,psibar(ipol)%numb_v(1),psibar(ipol)%nspin))
     
         call dvpsi_e (1, ipol,dvpsi(1,1),.false.)
         do i=1,num_nbndv(1)
            psibar(ipol)%wfn(1:npw,i,1)&
                 & = dvpsi(1:npw,i)
         enddo
      enddo
      deallocate (dvpsi)
   endif
   call mp_barrier(world_comm)

   do ipol=1,3
      do i=1,n_eig
         call absorption(vstate_r,psibar(ipol)%wfn(1,1,1),fc,i,bse_sp%a(i,ipol),ipol)
!      if(ionode) write(stdout,*)'Eigv#',i,'E',bse_spectrum(i)%e, 'Amp',bse_sp%a(i)  
!      if(ionode) write(stdout,*)'Eigv#',i,'E',bse_sp%en(i), 'Amp',bse_sp%a(i)  
      enddo
   enddo
   call mp_barrier(world_comm)

!build up the spectrum
   do ipol=1,3
      call build_spectrum(bse_sp%a(1,ipol),bse_sp%en(1),ipol)
   enddo

   deallocate(bse_spectrum)
   call free_memory_spectrum(bse_sp)
   do ipol=1,3
      call free_v_state(psibar(ipol))
   enddo

endif

!plot the excitonic wfn
if(l_plotexc) then
   if(l_plotaverage) then
      call plot_excwfn(plotn_min,plotn_max,vstate_r,fc)
   else
      do i=plotn_min,plotn_max
         call plot_excwfn(i,i,vstate_r,fc)
      enddo 
   endif
endif


do i=1,n_eig
   call free_memory_exc_a(bse_spectrum(i))
enddo


if(l_contraction) then
   call free_memory_contrac_w
   call free_imat(iimat_contract)
endif

call stop_clock('find_eig')
return



end subroutine
