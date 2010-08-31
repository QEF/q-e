!P.Umari  program GWW
! Modified by G. Stenuit
!
!this is the main part of the GWW program

   PROGRAM gww_fit

   USE  input_gw,            ONLY : input_options, read_input_gww
   USE  io_global,           ONLY : stdout, ionode
   USE  io_files,            ONLY : find_free_unit
   USE  self_energy_storage
   USE  expansion
   USE  energies_gww
   USE  start_end
   USE  mp_global,           ONLY : mpime
   USE  para_gww
   USE  times_gw
   USE  w_divergence
   USE constants,            ONLY : RYTOEV


   implicit none

   TYPE(input_options)   :: options
   TYPE(self_storage)    :: ss
   TYPE(self_expansion)  :: se
   TYPE(quasi_particles) :: qp
   TYPE(times_freqs)     :: tf
   TYPE(gv_time)         :: gt


   INTEGER :: ii, jj, iun1,iun2,iun,iunu,idumm, it
   REAL(kind=DP) :: rdumm1,rdumm2,rdumm3
   CHARACTER(5) :: nfile
   COMPLEX(kind=DP) :: zz, sigmac,dsigmac
   REAL(kind=DP)  :: offset
   INTEGER :: i_homo
   REAL(kind=DP), ALLOCATABLE :: ene_dft_u(:)
   REAL(kind=DP), ALLOCATABLE :: ss_r(:),ss_i(:)


!setup MPI environment


   call startup



   write(stdout,*)  'PROGRAM GWW-FIT: Version 0.2'



!initialize arrays
   call initialize_quasi_particle(qp)


!  read in input structure

   call read_input_gww(options)
   call flush_unit(stdout)
   call setup_para_gww(options%n, options%max_i, options%i_min, options%i_max)
! setup time/frequency grid if required



   call setup_timefreq(tf,options)
   call flush_unit(stdout)
!  call test_fft(tf)

!read in LDA and HF energies from bands.dat
   qp%max_i=options%max_i
   allocate(qp%ene_dft_ks(options%max_i))
   allocate(qp%ene_hf(options%max_i))
   allocate(qp%ene_gw(options%max_i))
   allocate(qp%ene_gw_pert(options%max_i))

!!! this file bands.dat can be created by doing
! 1. grep LDA out_00001 and remove the first two lines
! 2. awk '{print NR,$4,$7,$10,$13}' of the previous output
! 3. put in the first line the numero of the HOMO level
   iun = find_free_unit()
   open( unit=iun, file='bands.dat', status='old',form='formatted')
   read(iun,*) i_homo
   do ii=1,options%max_i
      read(iun,*) idumm,qp%ene_dft_ks(ii),rdumm1,rdumm2,qp%ene_hf(ii)
   enddo
   close(iun)
   qp%ene_dft_ks(:)=qp%ene_dft_ks(:)/RYTOEV
   qp%ene_hf(:)=qp%ene_hf(:)/RYTOEV

!!!! NEW if remainder=10 => read the file .hubbard_u and modifies the hf energies
!!!! this is done here only because in the pw.x run the hubbard correction in the
!!! first-order calculation was not done.
!!! in the code, NOW, everything is done correctly and so remainder=10 should
!!! NO longer be used !
   if(options%remainder==10) then
      write(stdout,*) 'Read the file .hubbard_u'
      call flush_unit(stdout)
      allocate(ene_dft_u(options%max_i))
      iunu = find_free_unit()
      open(unit=iunu,file=trim(options%prefix)//'.hubbard_u',status='old',form='unformatted')
      read(iunu) idumm
      do ii=1,options%max_i
         read(iunu) ene_dft_u(ii)
      enddo
      do ii=1,options%max_i
         write(stdout,*) 'ene_dft_u(', ii, ')=', ene_dft_u(ii)*RYTOEV, ' eV'
         call flush_unit(stdout)
      enddo
      qp%ene_hf(:)=qp%ene_hf(:)-ene_dft_u(:)
      deallocate(ene_dft_u)
      close(iunu)
   endif


!setup ss and read in date from fit graphs
   ss%ontime=.false.
   ss%whole_s=.false.
   ss%n=options%n
   ss%n_grid_fit=options%n_grid_fit
   ss%max_i=options%max_i
   ss%i_min=options%i_min
   ss%i_max=options%i_max
   ss%tau=options%tau


   allocate(ss%diag_freq_fit(ss%max_i,2*ss%n_grid_fit+1))
   do ii=ss%i_min,ss%i_max
      !
      allocate(ss_r(-ss%n_grid_fit:ss%n_grid_fit))
      allocate(ss_i(-ss%n_grid_fit:ss%n_grid_fit))
      !
      write(nfile,'(5i1)') &
           & ii/10000,mod(ii,10000)/1000,mod(ii,1000)/100,mod(ii,100)/10,mod(ii,10)
      iun1 = find_free_unit()
      write(stdout,*) 'iun1=', iun1
      call flush_unit(stdout)
      open( unit=iun1, file='re_on_im'// nfile, status='old',form='formatted')
      do jj=-ss%n_grid_fit,ss%n_grid_fit
          read(iun1,*) rdumm1,rdumm2,ss_r(jj),rdumm3
      enddo

      iun2 = find_free_unit()
      write(stdout,*) 'iun1=', iun1, ' and iun2=', iun2
      call flush_unit(stdout)
      open( unit=iun2, file='im_on_im'// nfile, status='old',form='formatted')
      do jj=-ss%n_grid_fit,ss%n_grid_fit
          read(iun2,*) rdumm1,rdumm2,ss_i(jj),rdumm3
      enddo

      do jj=-ss%n_grid_fit,ss%n_grid_fit
          !!! write(stdout,*) 'ss_r=', ss_r(jj), ' and ss_i=', ss_i(jj)
          ss%diag_freq_fit(ii,jj+ss%n_grid_fit+1)=cmplx(ss_r(jj),ss_i(jj))
      enddo
      !call flush_unit(stdout)
      !
      close(iun1)
      close(iun2)
      !
      deallocate(ss_r,ss_i)
   enddo


   call create_self_energy_fit( tf, se, ss, options)
   call free_memory_self_storage(ss)



   offset=-(qp%ene_dft_ks(i_homo+1)+qp%ene_dft_ks(i_homo))/2.d0
   do ii=options%i_min,options%i_max
      call value_on_frequency(se,ii,qp%ene_dft_ks(ii)+offset,sigmac)
      call derivative_on_frequency(se,ii,qp%ene_dft_ks(ii)+offset,dsigmac)

      write(stdout,*) 'value, zeta:',ii,sigmac,dsigmac,offset
      zz=(1.d0,0.d0)-dsigmac

      qp%ene_gw(ii)=qp%ene_dft_ks(ii)+offset +&
         (sigmac+qp%ene_hf(ii)-qp%ene_dft_ks(ii))/zz

      write(stdout,*) 'GW-PERT energy', ii,real(qp%ene_gw(ii)-offset)*RYTOEV

      qp%ene_gw_pert(ii)=qp%ene_gw(ii)-offset
!self-consistency loop
      do it=1,10
         call value_on_frequency_complex(se,ii,qp%ene_gw(ii),sigmac)
         write(stdout,*) 'Iteration energy',it,sigmac

         qp%ene_gw(ii)=qp%ene_dft_ks(ii)+offset+sigmac+qp%ene_hf(ii)-qp%ene_dft_ks(ii)

      enddo
      qp%ene_gw(ii)= qp%ene_gw(ii)-offset
   enddo

   !call create_quasi_particles(options,qp,se)
   call free_memory_self_expansion(se)
   call printout_quasi_particles(qp)

    !if(ionode) then
    !  write(stdout,*) 'QUASI-PARTICLES ENERGIES IN Ev:'
    !  do ii=options%i_min,options%i_max
    !     write(stdout,'(''State:'',i5,''LDA  :'',f10.5,'' GW-PERT  :'',f10.5,'' GW  :'',f10.5, &
    !                   '' HF-pert :'',f10.5)') &
    !                   & ii,qp%ene_dft_ks(ii)*RYTOEV, real(qp%ene_gw_pert(ii))*RYTOEV, &
    !                   & real(qp%ene_gw(ii))*RYTOEV,qp%ene_hf(ii)*RYTOEV
    !  enddo
    !  write(stdout,*) 'IMAGINARY ENERGIES IN Ev:'
    !  do ii=options%i_min,options%i_max
    !     write(stdout,'(''State:'',i5,'' GW (Im) :'',f10.5)') ii,aimag(qp%ene_gw(ii))*RYTOEV
    !  enddo
    !endif






!stops MPI
   call free_memory_times_freqs(tf)
   call free_memory_para_gww
   call stop_run

   stop
 END PROGRAM gww_fit

