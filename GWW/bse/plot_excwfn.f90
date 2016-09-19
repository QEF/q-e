subroutine plot_excwfn(nstart,nend,vstate_r,fc)
! this subroutine computes and writes on file the nplot-th excitonic wavefunction
! to be read by pp.x
! note that this subroutine is working only for gamma only calculations 

USE exciton
USE bse_wannier, ONLY:num_nbndv,&
                      r_hole,l_plotaverage
use bse_basic_structures
USE pwcom
USE fft_custom_gwl
USE io_global, ONLY : stdout,ionode,ionode_id
USE io_files, ONLY : tmp_dir,prefix
USE mp_world, ONLY : mpime, nproc
USE mp, ONLY: mp_sum
USE mp_world,             ONLY : world_comm
!USE io_files, ONLY : find_free_unit
USE ions_base,      ONLY : nat, tau, atm,ityp
 
implicit none
INTEGER, EXTERNAL :: find_free_unit

integer :: nplot,nstart,nend
type(v_state_r) :: vstate_r
type(exc_r) :: a_rt
type(fft_cus) :: fc

integer ::nxh,nyh,nzh,nh
INTEGER :: nr3s_start, nr3s_end
real(kind=dp), allocatable :: psi_exc(:)
real(kind=dp), allocatable :: psi_excio(:)
real(kind=dp), allocatable :: psi_excsum(:)
real(kind=dp), allocatable :: v_rh(:)

integer :: iv,ii,iplane,iz,ounit,ix,iy,iip
logical :: debug

CHARACTER(5) :: nfile

call start_clock('plot_excwfn')
debug=.true.

!check if all variables are ok from read_file in main 
if (debug) then
   if(ionode) then
      !bg(:,i) are the reciprocal lattice vectors, b_i,
      !in tpiba=2pi/alat units: b_i(:) = bg(:,i)/tpiba
      !at(:,i) are the lattice vectors of the simulation cell, a_i,
      !in alat units: a_i(:) = at(:,i)/alat
      write(stdout,*) 'plotexcwfn bg(:,1)=',bg(1,1),bg(2,1),bg(3,1)
      write(stdout,*) 'plotexcwfn bg(:,2)=',bg(1,2),bg(2,2),bg(3,2)
      write(stdout,*) 'plotexcwfn bg(:,3)=',bg(1,3),bg(2,3),bg(3,3)
      write(stdout,*) 'plotexcwfn alat=',alat
   endif
endif

!find FFT grid point (dual grid) closer to r_hole (given in alat units)

nxh = nint ( (r_hole(1)*bg(1,1) + r_hole(2)*bg(2,1) + r_hole(3)*bg(3,1) )*fc%nr1t) + 1
nyh = nint ( (r_hole(1)*bg(1,2) + r_hole(2)*bg(2,2) + r_hole(3)*bg(3,2) )*fc%nr2t) + 1
nzh = nint ( (r_hole(1)*bg(1,3) + r_hole(2)*bg(2,3) + r_hole(3)*bg(3,3) )*fc%nr3t) + 1

allocate(v_rh(num_nbndv(1)))
v_rh(:)=0.d0


!get the valence wavefunctions at the nxh,nyh,nzh (only one processor has it!) 
#if !defined(__MPI)
nh=(nzh-1)*fc%nrx1t*fc%nrx2t+(nyh-1)*fc%nrx1t+nxh
v_rh(:)=vstate_r%wfnrt(nh,:,1)
#else
nr3s_start=0
nr3s_end =0
do ii=1,mpime+1
   nr3s_start=nr3s_end+1
   nr3s_end=nr3s_end+fc%dfftt%npp(ii)
enddo


do iplane=1,fc%dfftt%npp(mpime+1)
   iz=nr3s_start+iplane-1
   if (iz==nzh) then
      nh=(iplane-1)*fc%nrx1t*fc%nrx2t+(nyh-1)*fc%nrx1t+nxh
      v_rh(:)=vstate_r%wfnrt(nh,:,1)
   endif
enddo
call mp_sum(v_rh,world_comm)
#endif


if (debug) then
   if(ionode) write(stdout,*) 'plotexcwfn qui'
endif
!stop

!allocate and initialize the excitonic wavefunction 
allocate(psi_exc(fc%nrxxt))
psi_exc(1:fc%nrxxt)=0.d0

allocate(psi_excsum(fc%nrx1t*fc%nrx2t*fc%nrx3t))
psi_excsum(1:fc%nrx1t*fc%nrx2t*fc%nrx3t)=0.d0

do nplot=nstart,nend
if (debug) then
   if(ionode) write(stdout,*) 'plotexcwfn qui2',nplot
endif
!
!FFT the excitonic wavefunction vector to real space (dual grid)
   call initialize_exc_r(a_rt)
   call fft_a_exc(bse_spectrum(nplot),fc,a_rt)

!now compute the exitonic wavefunction
   do iv=1,num_nbndv(1)
      psi_exc(1:fc%nrxxt)=psi_exc(1:fc%nrxxt)+v_rh(iv)*&
                                           a_rt%ar(1:a_rt%nrxxt,iv)
   enddo

if (debug) then
   if(ionode) write(stdout,*) 'plotexcwfn qui3',nplot
endif
!square modulus
   psi_exc(1:fc%nrxxt)=psi_exc(1:fc%nrxxt)**2

   if(debug) then
      if(ionode) write(stdout,*) 'fc%nr1t, fc%nr2t, fc%nr3t', fc%nr1t, fc%nr2t, fc%nr3t
      if(ionode) write(stdout,*) 'fc%nrx1t, fc%nrx2t, fc%nrx3t', fc%nrx1t, fc%nrx2t, fc%nrx3t
   endif

!Now gather the excitonic wavefunction from all the processors
!and sum for the l_plotaverage case (when nstart is different from nend) 

   allocate(psi_excio(fc%nrx1t*fc%nrx2t*fc%nrx3t))

   psi_excio(1:fc%nrx1t*fc%nrx3t*fc%nrx3t)=0.d0
   do iplane=1,fc%dfftt%npp(mpime+1)
     iz=nr3s_start+iplane-1
     do iy=1,fc%nr2t
        do ix=1,fc%nr1t
           ii=(iz-1)*(fc%nrx1t*fc%nrx2t)+(iy-1)*fc%nrx1t+ix
           iip=(iplane-1)*fc%nrx1t*fc%nrx2t+(iy-1)*fc%nrx1t+ix
           psi_excio(ii)=psi_exc(iip) 
        enddo
     enddo 
   enddo
   call mp_sum(psi_excio,world_comm)

   psi_excsum(1:fc%nrx1t*fc%nrx3t*fc%nrx3t)=psi_excsum(1:fc%nrx1t*fc%nrx3t*fc%nrx3t)+&
                      psi_excio(1:fc%nrx1t*fc%nrx3t*fc%nrx3t)/(real(nend)-real(nstart)+1.d0)

   if (debug) then
      if(ionode) write(stdout,*) 'plotexcwfn qui3',nplot
   endif
   call free_memory_exc_a_r(a_rt)
   deallocate(psi_excio)
enddo ! nplot


!
! XCRYSDEN FORMAT
!
if(ionode) then
   ounit=find_free_unit()
!   open(ounit,file='exc_average.xsf',form='formatted')
   if(l_plotaverage)   open(ounit,file='exc_average.xsf',form='formatted')
   if(.not.l_plotaverage)   then
       write(nfile,'(5i1)') &
        & nstart/10000,mod(nstart,10000)/1000,mod(nstart,1000)/100,mod(nstart,100)/10,mod(nstart,10)
      open(ounit,file=trim(tmp_dir)//trim(prefix)//'.exc.xsf'//nfile,form='formatted')
   endif
   CALL xsf_struct (alat, at, nat, tau, atm, ityp, ounit)
   CALL xsf_fast_datagrid_3d &
             (psi_excsum, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, at, alat, ounit)
!   close(ounit)   
endif


deallocate(v_rh,psi_exc)

deallocate(psi_excsum)

call stop_clock('plot_excwfn')
end subroutine plot_excwfn
