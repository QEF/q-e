!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!



   PROGRAM gww_fit

   USE  input_gw,            ONLY : input_options, read_input_gww
   USE  io_global,           ONLY : stdout, ionode
   USE  self_energy_storage
   USE  expansion
   USE  energies_gww
   USE  start_end
   USE  mp_world,           ONLY : mpime
   USE  para_gww
   USE  times_gw
   USE  w_divergence
   USE constants,            ONLY : RYTOEV
   USE io_files,  ONLY :  prefix,tmp_dir


   implicit none
   INTEGER, EXTERNAL :: find_free_unit
   TYPE(input_options)   :: options
   TYPE(self_storage)    :: ss
   TYPE(self_expansion)  :: se
   TYPE(self_on_real)    :: sr
   TYPE(quasi_particles) :: qp
   TYPE(times_freqs)     :: tf
   TYPE(gv_time)         :: gt


   INTEGER :: ii, jj, iun1,iun2,iun,idumm, it
   REAL(kind=DP) :: rdumm1,rdumm2,rdumm3,ss_r,ss_i,ddumm
   CHARACTER(5) :: nfile
   COMPLEX(kind=DP) :: zz, sigmac,dsigmac
   REAL(kind=DP)  :: offset
   INTEGER :: i_homo
   INTEGER :: ifil!DEBUG
   INTEGER, ALLOCATABLE :: order(:)

!setup MPI environment


   call startup



   write(stdout,*)  'PROGRAM GWW-FIT: Version 0.2'



!initialize arrays
   call initialize_quasi_particle(qp)
   call initialize_self_storage(ss)
   call initialize_self_on_real(sr)

!  read in input structure

   call read_input_gww(options)
   FLUSH(stdout)
   call setup_para_gww(options%n, options%max_i, options%i_min, options%i_max)
! setup time/frequency grid if required
   
   

   call setup_timefreq(tf,options)
   FLUSH(stdout)
!  call test_fft(tf)

!read in LDA and HF energies from bands.dat
   qp%max_i=options%max_i
   qp%nspin=options%nspin
   allocate(qp%ene_dft_ks(options%max_i,qp%nspin))
   allocate(qp%ene_hf(options%max_i,qp%nspin))
   allocate(qp%ene_gw(options%max_i,qp%nspin))
   allocate(qp%ene_gw_pert(options%max_i,qp%nspin))
   allocate(order(options%max_i))
   open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'bands.dat', status='old',form='formatted')
   read(iun,*) i_homo
   do ii=1,options%max_i
      read(iun,*) idumm,qp%ene_dft_ks(ii,1),rdumm1,rdumm2,qp%ene_hf(ii,1)
   enddo
   close(iun)
   qp%ene_dft_ks(:,1)=qp%ene_dft_ks(:,1)/RYTOEV
   qp%ene_hf(:,1)=qp%ene_hf(:,1)/RYTOEV
!setup ss and read in date from fit graphs
   ss%ontime=.false.
   ss%whole_s=.false.
   ss%n=options%n
   ss%n_grid_fit=options%n_grid_fit
   ss%max_i=options%max_i
   ss%i_min=options%i_min
   ss%i_max=options%i_max
   ss%tau=options%tau
   ss%nspin=1

   if(options%l_order) then
      open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'order.dat', status='old',form='formatted')
      do ii=1,options%max_i
         read(iun,*) idumm,order(ii),ddumm
      enddo
      close(iun)
   else
      do ii=1,options%max_i
         order(ii)=ii
      enddo
   endif

   allocate(ss%diag_freq_fit(ss%max_i,2*ss%n_grid_fit+1,1))
   do ii=ss%i_min,ss%i_max
      write(nfile,'(5i1)') &  
           & order(ii)/10000,mod(order(ii),10000)/1000,mod(order(ii),1000)/100,mod(order(ii),100)/10,mod(order(ii),10)
      !ifil=ii-4
      !write(nfile,'(5i1)') &
      !     & ifil/10000,mod(ifil,10000)/1000,mod(ifil,1000)/100,mod(ifil,100)/10,mod(ifil,10)
      iun1 = find_free_unit()
      open( unit=iun1, file=trim(tmp_dir)//trim(prefix)//'-'//'re_on_im'// nfile, status='old',form='formatted')
      open( unit=iun2, file=trim(tmp_dir)//trim(prefix)//'-'//'im_on_im'// nfile, status='old',form='formatted')
      do jj=-ss%n_grid_fit,ss%n_grid_fit
          read(iun1,*) rdumm1,rdumm2,ss_r,rdumm3
          read(iun2,*) rdumm1,rdumm2,ss_i,rdumm3
          ss%diag_freq_fit(ii,jj+ss%n_grid_fit+1,1)=cmplx(ss_r,ss_i)
      enddo
      close(iun1)
      close(iun2)
   enddo


   if(options%n_real_axis>=0) then
      call create_self_on_real(options, sr)
   endif
   call create_self_energy_fit( tf, se, ss, options,sr,.false.)  
   call free_memory_self_storage(ss)
   call free_memory_self_on_real(sr)

   
  
   offset=-(qp%ene_dft_ks(i_homo+1,1)+qp%ene_dft_ks(i_homo,1))/2.d0
   !offset=-(6.07889+6.75583)/2.d0/13.606d0!DEBUG

   do ii=options%i_min,options%i_max
      call value_on_frequency(se,ii,qp%ene_dft_ks(ii,1)+offset,sigmac,1)
      call derivative_on_frequency(se,ii,qp%ene_dft_ks(ii,1)+offset,dsigmac,1)
      
      write(stdout,*) 'value, zeta:',ii,sigmac,dsigmac,offset
      zz=(1.d0,0.d0)-dsigmac
     
      qp%ene_gw(ii,1)=qp%ene_dft_ks(ii,1)+offset +&
         (sigmac+qp%ene_hf(ii,1)-qp%ene_dft_ks(ii,1))/zz
     
      write(stdout,*) 'GW-PERT energy', ii,real(qp%ene_gw(ii,1)-offset)*RYTOEV

      qp%ene_gw_pert(ii,1)=qp%ene_gw(ii,1)-offset
!self-consistency loop
      do it=1,10
         call value_on_frequency_complex(se,ii,qp%ene_gw(ii,1),sigmac,1)
         write(stdout,*) 'Iteration energy',it,sigmac
      
         qp%ene_gw(ii,1)=qp%ene_dft_ks(ii,1)+offset+sigmac+qp%ene_hf(ii,1)-qp%ene_dft_ks(ii,1) 
      
      enddo
      qp%ene_gw(ii,1)= qp%ene_gw(ii,1)-offset
   enddo

   !call create_quasi_particles(options,qp,se)
   call free_memory_self_expansion(se)
      
    if(ionode) then
     write(stdout,*) 'QUASI-PARTICLES ENERGIES IN Ev:'
     do ii=options%i_min,options%i_max
        write(stdout,'(''State:'',i5,''LDA  :'',f10.5,'' GW-PERT  :'',f10.5,'' GW  :'',f10.5, &
     & '' HF-pert :'',f10.5)') &
             & ii,qp%ene_dft_ks(ii,1)*RYTOEV, real(qp%ene_gw_pert(ii,1))*RYTOEV, &
             & real(qp%ene_gw(ii,1))*RYTOEV,qp%ene_hf(ii,1)*RYTOEV
     enddo
     write(stdout,*) 'IMAGINARY ENERGIES IN Ev:'
     do ii=options%i_min,options%i_max
        write(stdout,'(''State:'',i5,'' GW (Im) :'',f10.5)') ii,aimag(qp%ene_gw(ii,1))*RYTOEV
     enddo
  endif

   

  deallocate(order)


!stops MPI
   call free_memory_times_freqs(tf)
   call free_memory_para_gww
   call stop_run

   stop
 END PROGRAM gww_fit

