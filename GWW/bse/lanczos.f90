subroutine lanczos(vstate,vstate_r,cstate,wcstate,fc)
!this subroutine computes the absorption spectrum through the lanczos procedure
USE exciton

USE bse_basic_structures
USE fft_custom_gwl
USE bse_wannier, ONLY: nit_lcz, l_contraction
USE contract_w
USE lsda_mod, ONLY :nspin
USE io_global, ONLY : stdout

implicit none
type(v_state),   intent(in) :: vstate
type(v_state_r), intent(in) :: vstate_r
type(c_state),   intent(in) :: cstate
type(c_state),   intent(in) :: wcstate
type(fft_cus),   intent(in) :: fc


real(kind=DP), allocatable :: a(:,:),b(:,:)
integer :: ispin

call start_clock('lanczos')
allocate (a(nit_lcz,3))
allocate (b(nit_lcz,3))
if(l_contraction) then
   write(stdout,*) 'CALL contract_w_build'
   FLUSH(stdout)
   call contract_w_build(fc)
   call initialize_imat(iimat_contract)
   do ispin=1,nspin
      call read_iimat(iimat_contract,ispin)
   enddo
   write(stdout,*) 'CALL contract_v_build'
   FLUSH(stdout)
   call contract_v_build(fc)

endif


!perform lanczos iterations

call lanczos_iterations(vstate,vstate_r,cstate,wcstate,fc,a(1,1),b(1,1))

!build the continuum fraction

call lanczos_cf(a(1,1),b(1,1))

if(l_contraction) then
   call free_memory_contrac_w
   call free_imat(iimat_contract)
endif

deallocate(a)
deallocate(b)

call stop_clock('lanczos')
return
end subroutine

subroutine lanczos_iterations(vstate,vstate_r,cstate,wcstate,fc,a,b)
!this subroutine computes the lanczos iteration to get the a(i) and b(i) for the
!continued fraction

USE exciton

USE bse_basic_structures
USE fft_custom_gwl
USE wvfct,       ONLY : npw,npwx,nbnd
USE bse_wannier, ONLY: num_nbndv, nit_lcz,l_restart_lcz, nlcz_restart
USE mp,          ONLY : mp_barrier,mp_bcast
USE mp_world,    ONLY : world_comm,mpime
USE io_global,   ONLY : ionode,ionode_id
use io_files,    ONLY : tmp_dir, prefix


implicit none
INTEGER, EXTERNAL :: find_free_unit
type(v_state),   intent(in) :: vstate
type(v_state_r), intent(in) :: vstate_r
type(c_state),   intent(in) :: cstate
type(c_state),   intent(in) :: wcstate
type(fft_cus),   intent(in) :: fc

integer :: i,ipol,j,is,iunab,iuni,nstart
complex(DP), allocatable:: dvpsi(:,:) !formula (43) of Rev. Mod. Phys. 73, 515 


type(exc) im1_s ! |i-1> 
type(exc) i_s ! |i>
type(exc) ip1_s ! |i+1>

real(kind=DP) :: bi,bim1
real(kind=DP), intent(inout) :: a(nit_lcz,3),b(nit_lcz,3)
CHARACTER(5) :: nproc
CHARACTER(5) :: nfile

logical :: debug

call start_clock('lanczos_iterations')
debug=.true.

call initialize_exc(im1_s)
im1_s%label=1
im1_s%npw=npw
im1_s%numb_v=num_nbndv(1)
allocate(im1_s%a(im1_s%npw,im1_s%numb_v))

call initialize_exc(i_s)
i_s%label=1
i_s%npw=npw
i_s%numb_v=num_nbndv(1)
allocate(i_s%a(i_s%npw,i_s%numb_v))

call initialize_exc(ip1_s)
ip1_s%label=1
ip1_s%npw=npw
ip1_s%numb_v=num_nbndv(1)
allocate(ip1_s%a(ip1_s%npw,ip1_s%numb_v))

allocate (dvpsi ( npwx , num_nbndv(1)))

if(l_restart_lcz) then
   if(ionode) then
      if(debug) write(*,*) 'Restarting lanczos'
      iunab = find_free_unit()
      open(unit=iunab, file=trim(tmp_dir)//trim(prefix)//'.lczrestart_ab.dat', status='unknown', form='unformatted')
      read(iunab) a(1:nlcz_restart,1), a(1:nlcz_restart,2), a(1:nlcz_restart,3)
      read(iunab) b(1:nlcz_restart,1), b(1:nlcz_restart,2), b(1:nlcz_restart,3)
      close(iunab)
   endif
   call mp_bcast(a(1:(nlcz_restart),1:3),ionode_id, world_comm)
   call mp_bcast(b(1:(nlcz_restart),1:3),ionode_id, world_comm)
endif

do ipol=1,3 

   if(.not.l_restart_lcz) then
!     compute the  |psibar(iv)> and set it as initial excitonic state |i-1> for 
!     the lanczos procedure
      if(debug) write(*,*) 'before dvpsi'
      call dvpsi_e (1, ipol,dvpsi(1,1),.false.)
      if(debug) write(*,*) 'after dvpsi'
      do i=1,num_nbndv(1)
         im1_s%a(1:npw,i)= dvpsi(1:npw,i)
      enddo

      call normalize_exc(im1_s)
      if(debug) write(*,*) 'after normalize_exc'
 
!     apply the exc Hamiltonian  
      call exc_h_a(im1_s,i_s,vstate,vstate_r,cstate,wcstate,fc) 
      if(debug) write(*,*) 'after exc_h_a'

!     a(1)= <1|H|1>
      call sproduct_exc(im1_s,i_s,a(1,ipol))

      do i=1,num_nbndv(1)
         i_s%a(1:npw,i)=i_s%a(1:npw,i)-dcmplx(a(1,ipol),0.d0)*im1_s%a(1:npw,i)
      enddo

!     b(1)=bim1=|H|1>-a(1)|1>|
      call sproduct_exc(i_s,i_s,bim1)
      b(1,ipol)=sqrt(bim1)

!     project into the conduction manifold
      do is = 1,vstate%nspin
         call pc_operator_exc(i_s,vstate,is)
      enddo
 
!     and normalize 
      call normalize_exc(i_s)


!     apply the exc Hamiltonian  
      call exc_h_a(i_s,ip1_s,vstate,vstate_r,cstate,wcstate,fc) 
   
!     a(2)= <2|H|2>
      call sproduct_exc(i_s,ip1_s,a(2,ipol))

      do i=1,num_nbndv(1)
         ip1_s%a(1:npw,i)=ip1_s%a(1:npw,i)-dcmplx(a(2,ipol),0.d0)*i_s%a(1:npw,i)-dcmplx(b(1,ipol),0.d0)*im1_s%a(1:npw,i)
      enddo

!    b(2)=bi=|H|2>-a(2)|2>-b(1)|1>|
      call sproduct_exc(ip1_s,ip1_s,bi)
      b(2,ipol)=sqrt(bi)

!     project into the conduction manifold
      do is = 1,vstate%nspin
         call pc_operator_exc(ip1_s,vstate,is)
      enddo
 
!     and normalize 
      call normalize_exc(ip1_s)
      nstart=3
    else
!   read starting excitonic vector
      nstart=nlcz_restart+1
      iuni = find_free_unit()
      write(nproc,'(5i1)') &
              & mpime/10000,mod(mpime,10000)/1000,mod(mpime,1000)/100,mod(mpime,100)/10,mod(mpime,10)
      write(nfile,'(5i1)') &
           & ipol/10000,mod(ipol,10000)/1000,mod(ipol,1000)/100,mod(ipol,100)/10,mod(ipol,10)
      open( unit=iuni, file=trim(tmp_dir)//trim(prefix)//'.lcz_is.'// nfile //'.'// nproc , status='unknown',form='unformatted')
      read(iuni) i_s%label
      read(iuni) i_s%npw
      read(iuni) i_s%numb_v
      read(iuni) i_s%e
      do j=1,i_s%numb_v
         read(iuni)  i_s%a(1:i_s%npw,j)
      enddo
      close(iuni)

      iuni = find_free_unit()
      write(nproc,'(5i1)') &
              & mpime/10000,mod(mpime,10000)/1000,mod(mpime,1000)/100,mod(mpime,100)/10,mod(mpime,10)
      write(nfile,'(5i1)') &
           & ipol/10000,mod(ipol,10000)/1000,mod(ipol,1000)/100,mod(ipol,100)/10,mod(ipol,10)
      open( unit=iuni, file=trim(tmp_dir)//trim(prefix)//'.lcz_ip1s.'// nfile //'.'// nproc , status='unknown',form='unformatted')
      read(iuni) ip1_s%label
      read(iuni) ip1_s%npw
      read(iuni) ip1_s%numb_v
      read(iuni) ip1_s%e
      do j=1,ip1_s%numb_v
         read(iuni)  ip1_s%a(1:ip1_s%npw,j)
      enddo
      close(iuni)
    endif

!   Now start lanczos iteration

    do j=nstart,nit_lcz

      if(ionode.and.(mod(j,10)==0)) write(*,*) 'lanczos iteration #', j
      do i=1,num_nbndv(1)
         im1_s%a(1:npw,i)=i_s%a(1:npw,i) ! |j-1>=|j>
      enddo

      do i=1,num_nbndv(1)
         i_s%a(1:npw,i)=ip1_s%a(1:npw,i) ! |j>=|j+1>
      enddo
     
!     apply the exc Hamiltonian  
      call exc_h_a(i_s,ip1_s,vstate,vstate_r,cstate,wcstate,fc) 
   
!     a(j)= <j|H|j>
      call sproduct_exc(i_s,ip1_s,a(j,ipol))

      do i=1,num_nbndv(1)
         ip1_s%a(1:npw,i)=ip1_s%a(1:npw,i)-dcmplx(a(j,ipol),0.d0)*i_s%a(1:npw,i)&
                          -dcmplx(b(j-1,ipol),0.d0)*im1_s%a(1:npw,i)
      enddo

!     b(j)=|H|j>-a(j)|j>-b(j-1)|j-1>|
      call sproduct_exc(ip1_s,ip1_s,bi)
      b(j,ipol)=sqrt(bi)

!     project into the conduction manifold
      do is = 1,vstate%nspin
         call pc_operator_exc(ip1_s,vstate,is)
      enddo
 
!     and normalize 
      call normalize_exc(ip1_s)
     
      call mp_barrier(world_comm)
    enddo ! end of lanczos iterations

!    write restart information on file

   iuni = find_free_unit()
   write(nproc,'(5i1)') &
           & mpime/10000,mod(mpime,10000)/1000,mod(mpime,1000)/100,mod(mpime,100)/10,mod(mpime,10)
   write(nfile,'(5i1)') &
        & ipol/10000,mod(ipol,10000)/1000,mod(ipol,1000)/100,mod(ipol,100)/10,mod(ipol,10)
   open( unit=iuni, file=trim(tmp_dir)//trim(prefix)//'.lcz_is.'// nfile //'.'// nproc , status='unknown',form='unformatted')
   write(iuni) i_s%label
   write(iuni) i_s%npw
   write(iuni) i_s%numb_v
   write(iuni) i_s%e
   do j=1,i_s%numb_v
      write(iuni)  i_s%a(1:i_s%npw,j)
   enddo
   close(iuni)

   iuni = find_free_unit()
   write(nproc,'(5i1)') &
           & mpime/10000,mod(mpime,10000)/1000,mod(mpime,1000)/100,mod(mpime,100)/10,mod(mpime,10)
   write(nfile,'(5i1)') &
        & ipol/10000,mod(ipol,10000)/1000,mod(ipol,1000)/100,mod(ipol,100)/10,mod(ipol,10)
   open( unit=iuni, file=trim(tmp_dir)//trim(prefix)//'.lcz_ip1s.'// nfile //'.'// nproc , status='unknown',form='unformatted')
   write(iuni) ip1_s%label
   write(iuni) ip1_s%npw
   write(iuni) ip1_s%numb_v
   write(iuni) ip1_s%e
   do j=1,ip1_s%numb_v
      write(iuni)  ip1_s%a(1:ip1_s%npw,j)
   enddo
   close(iuni)

enddo !ipol

!    write restart information on file
if(ionode) then
      iunab = find_free_unit()
      open(unit=iunab, file=trim(tmp_dir)//trim(prefix)//'.lczrestart_ab.dat', status='unknown', form='unformatted')
      write(iunab) a(1:nit_lcz,1), a(1:nit_lcz,2), a(1:nit_lcz,3)
      write(iunab) b(1:nit_lcz,1), b(1:nit_lcz,2), b(1:nit_lcz,3)
      close(iunab)
endif

if(debug) then
   if(ionode) then
     do ipol=1,3
       do j=1,nit_lcz
         write(*,*) 'ipol, it, a', ipol, j, a(j,ipol)
       enddo
     enddo
     do ipol=1,3
       do j=1,nit_lcz
         write(*,*) 'ipol, it, b', ipol, j, b(j,ipol)
       enddo
     enddo
    endif
endif 

!free memory
deallocate (dvpsi)
call free_memory_exc_a(im1_s)
call free_memory_exc_a(i_s)
call free_memory_exc_a(ip1_s)

call stop_clock('lanczos_iterations')
return
end subroutine


subroutine lanczos_cf(a,b)

USE bse_wannier, ONLY: nit_lcz, spectra_e_min,spectra_e_max,spectra_nstep
USE constants,   ONLY : RYTOEV, PI
USE kinds,       ONLY: DP
USE io_global, ONLY : ionode

implicit none
real(DP), intent(in) :: a(nit_lcz,3),b(nit_lcz,3)

complex(DP), allocatable      :: comega_g(:)
complex(DP), allocatable      :: den(:)
real(DP), allocatable         :: abss(:,:)    ! epsilon2
real(DP), allocatable         :: rp_abss(:,:) ! epsilon1
real(DP)                      :: eta,step,e_start

integer :: ipol,j,i
logical :: debug, im

call start_clock('lanczos_cf')
debug=.false.
eta=0.001d0

allocate(comega_g(spectra_nstep))
allocate(den(spectra_nstep))
allocate(abss(spectra_nstep,3))
allocate(rp_abss(spectra_nstep,3))

!build the energy grid (including a small imaginary part eta)
step=(spectra_e_max-spectra_e_min)/(dble(spectra_nstep-1)*RYTOEV)
e_start=spectra_e_min/RYTOEV

do i=0, spectra_nstep-1
   comega_g(i+1)=dcmplx((e_start+dble(i)*step),eta)
enddo

do ipol=1,3
!  build the continued fraction
   
   den(1:spectra_nstep)=comega_g(1:spectra_nstep)-dcmplx(a(nit_lcz,ipol),0.d0)

   if((debug).and.(ionode)) then
            write(*,*) 'ipol, den'
         do i=1,spectra_nstep
            write(*,'(I1,I6,2F8.4)') ipol, i, real(den(i)), aimag(den(i))
      enddo
   endif 

   do j=nit_lcz-1,1,-1
      den(1:spectra_nstep)= comega_g(1:spectra_nstep)&
                            -dcmplx(a(j,ipol),0.d0)&
                            -dcmplx(b(j,ipol)**2.d0,0.d0)/den(1:spectra_nstep)

      if((debug).and.(ionode)) then
            write(*,*) 'ipol, den'
         do i=1,spectra_nstep
            write(*,'(I1,I6,2F8.4)') ipol, i, real(den(i)), aimag(den(i))
         enddo
      endif 

   enddo

   abss(1:spectra_nstep,ipol)=-4*PI*aimag(dcmplx(1.d0,0.d0)/den(1:spectra_nstep))
   rp_abss(1:spectra_nstep,ipol)=1.d0-4*PI*real(dcmplx(1.d0,0.d0)/den(1:spectra_nstep))

   if((debug).and.(ionode)) then
     write(*,*) 'ABSORPTION ipol', ipol
     do i=1,spectra_nstep
        write(*,'(I6,F8.4,F12.4)') i, real(comega_g(i))*RYTOEV, abss(i,ipol)
     enddo
   endif 


enddo !ipol
im=.true.
call print_spectrum(abss,im)
im=.false.
call print_spectrum(rp_abss,im)


deallocate(comega_g)
deallocate(den)
deallocate(abss)

call stop_clock('lanczos_cf')
return
end subroutine

