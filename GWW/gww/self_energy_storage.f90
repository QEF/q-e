!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

MODULE self_energy_storage
!this modules contains the structure and subroutines
!to store the expectation values of the self-energy 
!and to perform ffts and fits
!in parallel version the calculations on times are parallelized

  USE kinds, ONLY  : DP 
 
  TYPE self_storage
!descriptor of <Psi_i|\Sigma|\Psi_j>
     LOGICAL :: ontime!if .true. data is on imaginary time , otherwise imaginary frequency
     LOGICAL :: whole_s!if .true. also the off-diagonal elements are considered
     INTEGER :: n!number of sample on positive and on negative times (total of 2*n+1 samples)
     INTEGER :: n_grid_fit!number of sample on positive and on negative frequencies for fit (total of 2*n+1 samples)
     INTEGER :: max_i!number of states considered
     INTEGER :: i_min!minimum state to be calculated
     INTEGER :: i_max!maximum state to be calculated
     INTEGER :: nspin!spin multiplicity
     REAL(kind=DP) :: tau!max time (on imaginary axes)
     COMPLEX(kind=DP), DIMENSION(:,:,:), POINTER :: diag !values <Psi_i|\Sigma|\Psi_i>,time_j
     COMPLEX(kind=DP), DIMENSION(:,:,:,:), POINTER :: whole !values  <Psi_i|\Sigma|\Psi_j>,time_k 
     COMPLEX(kind=DP), DIMENSION(:,:,:), POINTER :: diag_freq_fit !values <Psi_i|\Sigma|\Psi_i>,on frequency for fit 
     COMPLEX(kind=DP), DIMENSION(:,:,:,:), POINTER :: whole_freq_fit !values  <Psi_i|\Sigma|\Psi_j>, on frequency for fit
     REAL(kind=DP), POINTER, DIMENSION(:,:) :: ene_remainder!for storing remainders
     INTEGER :: i_min_whole!minimum state to be calculated for off-diagonal elements
     INTEGER :: i_max_whole!maximum state to be calculated for off-diagonal elements
  END  TYPE self_storage


  TYPE self_on_real
!descriptor of <Psi_i|\Sigma|\Psi_j> on an arbitrary grid on the real self_energy axis (inside HOMO-LUMO gap)
!or in general on an arbitrary grid in complex plane
     INTEGER :: n!number of samples
     INTEGER :: max_i!number of states considered
     INTEGER :: i_min!minimum state to be calculated                     
     INTEGER :: i_max!maximum state to be calculated
     INTEGER :: nspin!spin multiplicity 
     COMPLEX(kind=DP), DIMENSION(:), POINTER :: grid!grid point
     COMPLEX(kind=DP), DIMENSION(:,:,:), POINTER :: diag!diagonal expectation values
  END TYPE self_on_real

CONTAINS

  SUBROUTINE initialize_self_on_real(sr)
    implicit none
    TYPE(self_on_real) :: sr
    nullify(sr%grid)
    nullify(sr%diag)
    return
  END SUBROUTINE initialize_self_on_real

  SUBROUTINE free_memory_self_on_real(sr)
    implicit none
    TYPE(self_on_real) :: sr
    if(associated(sr%grid)) deallocate(sr%grid)
    nullify(sr%grid)
     if(associated(sr%diag)) deallocate(sr%diag)
    nullify(sr%diag)
  END SUBROUTINE free_memory_self_on_real


  SUBROUTINE initialize_self_storage(ss)
    implicit none
    TYPE(self_storage) :: ss
  
    
    nullify(ss%diag)
    nullify(ss%whole)
    nullify(ss%ene_remainder)
    nullify(ss%diag_freq_fit)
    nullify(ss%whole_freq_fit)

    return

  END SUBROUTINE initialize_self_storage

  SUBROUTINE create_self_on_real(options, sr)
!this subroutine create the object self_on_real reading the data from disk
    USE io_global,          ONLY : stdout, ionode,ionode_id
    USE input_gw,           ONLY : input_options
    USE mp,                 ONLY : mp_bcast
    USE mp_world,           ONLY : world_comm
    USE io_files,             ONLY : prefix,tmp_dir

    implicit none
    INTEGER, EXTERNAL :: find_free_unit
    TYPE(input_options) :: options
    TYPE(self_on_real)  :: sr
    INTEGER :: ii,iw,iun
    CHARACTER(5) :: nfile
    REAL(kind=DP) :: x, y1, y2

    sr%max_i=options%max_i
    sr%i_min=options%i_min
    sr%i_max=options%i_max
    sr%n = options%n_real_axis
    sr%nspin=1

    call initialize_self_on_real(sr)
    allocate(sr%grid(options%n_real_axis))
    allocate(sr%diag(options%n_real_axis,options%max_i,1))
    sr%grid(:)=(0.d0,0.d0)
    sr%diag(:,:,:)=(0.d0,0.d0)
    do ii=options%i_min, options%i_max
       write(nfile,'(5i1)') &
            & ii/10000,mod(ii,10000)/1000,mod(ii,1000)/100,mod(ii,100)/10,mod(ii,10)
       if(ionode) then
          iun = find_free_unit()
          open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'self_on_real'// nfile, status='unknown',form='formatted')
          do iw=1,options%n_real_axis
             read(iun,*)  x, y1, y2
             if(x<-0.315) y1=y1+0.266184-0.004408
             sr%grid(iw)=dcmplx(x,0.d0)
             sr%diag(iw,ii,1)=dcmplx(y1,y2)
          enddo
          close(iun)
       endif
       call mp_bcast(sr%diag(:,ii,1),ionode_id,world_comm)
    enddo
    call mp_bcast(sr%grid,ionode_id,world_comm)
    return
  END SUBROUTINE create_self_on_real

  SUBROUTINE write_self_storage_ondisk(ss, options)
!this subroutine writes the green function on disk
!the file name is taken from the label

    USE io_global,          ONLY : stdout, ionode
    USE input_gw,           ONLY : input_options
    USE mp,                 ONLY : mp_barrier
    USE mp_world,           ONLY : world_comm
    USE io_files,             ONLY : prefix,tmp_dir
    implicit none
    INTEGER, EXTERNAL :: find_free_unit
    TYPE(self_storage) :: ss!the self_energy descriptor to be written on file
    TYPE(input_options) :: options!for debug flag

    INTEGER :: iw, jw, kw, iun, is

    if(ionode) then
       iun = find_free_unit()
     
       open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'storage', status='unknown',form='unformatted')
       
       write(iun) ss%ontime
       write(iun) ss%whole_s
       write(iun) ss%n
       write(iun) ss%max_i
       write(iun) ss%i_min
       write(iun) ss%i_max
       write(iun) ss%tau
       write(iun) ss%n_grid_fit
       write(iun) ss%i_min_whole
       write(iun) ss%i_max_whole
       write(iun) ss%nspin
       do is=1,ss%nspin
          do iw=1,2*ss%n+1
             write(iun) ss%diag(1:ss%max_i,iw,is)
          end do
          if(ss%whole_s) then
             do iw=1,2*ss%n+1
                write(iun) ss%whole(ss%i_min_whole:ss%i_max_whole,1:ss%max_i,iw,is)
             end do
          endif
          do iw=1,2*ss%n_grid_fit+1
             write(iun) ss%diag_freq_fit(1:ss%max_i,iw,is)
          end do
          if(ss%whole_s) then
             do iw=1,2*ss%n_grid_fit+1
                write(iun) ss%whole_freq_fit(ss%i_min_whole:ss%i_max_whole,1:ss%max_i,iw,is)
             end do
          endif
       enddo
       close(iun)
    endif
    call mp_barrier( world_comm )
  END SUBROUTINE write_self_storage_ondisk

  SUBROUTINE read_self_storage_ondisk(ss, options)
!this subroutine writes the green function on disk
!the file name is taken from the label

    USE io_global,          ONLY : stdout, ionode, ionode_id
    USE input_gw,           ONLY : input_options
    USE mp,                 ONLY : mp_barrier, mp_bcast
    USE mp_world,           ONLY : world_comm
     USE io_files,             ONLY : prefix,tmp_dir
    implicit none
    INTEGER, EXTERNAL :: find_free_unit
    TYPE(self_storage) :: ss!the self_energy descriptor to be read from file
    TYPE(input_options) :: options!for debug flag
    
    INTEGER :: iw, jw, kw, iun,is

    if(ionode) then
       iun = find_free_unit()
       open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'storage', status='old',form='unformatted')
    endif

  !  call free_memory_self_storage(ss)

    if(ionode) then
    
       read(iun) ss%ontime
       read(iun) ss%whole_s
       read(iun) ss%n
       read(iun) ss%max_i
       read(iun) ss%i_min
       read(iun) ss%i_max
       read(iun) ss%tau
       read(iun) ss%n_grid_fit
       read(iun) ss%i_min_whole
       read(iun) ss%i_max_whole
       read(iun) ss%nspin
       
    endif

    call mp_bcast(ss%ontime, ionode_id,world_comm)
    call mp_bcast(ss%whole_s, ionode_id,world_comm)
    call mp_bcast(ss%n, ionode_id,world_comm)
    call mp_bcast(ss%max_i, ionode_id,world_comm)
    call mp_bcast(ss%i_min, ionode_id,world_comm)
    call mp_bcast(ss%i_max, ionode_id,world_comm)
    call mp_bcast(ss%tau, ionode_id,world_comm)
    call mp_bcast(ss%n_grid_fit, ionode_id,world_comm)
    call mp_bcast(ss%i_min_whole, ionode_id,world_comm)
    call mp_bcast(ss%i_max_whole, ionode_id,world_comm)
    call mp_bcast(ss%nspin, ionode_id,world_comm)

!check for consistency
    if(ss%max_i/=options%max_i) then
      write(stdout,*) 'Routine  read_self_storage_ondisk max_i wrong'
      stop
    endif

!allocates
    if(ss%whole_s) then
      allocate(ss%whole(ss%i_min_whole:ss%i_max_whole,ss%max_i,2*ss%n+1,ss%nspin))    
   else
      nullify(ss%whole)
   endif
   allocate(ss%diag(ss%max_i,2*ss%n+1,ss%nspin))
     
 
    


    allocate(ss%ene_remainder(ss%max_i,ss%nspin))
    
    if(ss%whole_s) then
      allocate(ss%whole_freq_fit(ss%i_min_whole:ss%i_max_whole,ss%max_i,2*ss%n_grid_fit+1,ss%nspin))
    else
       nullify(ss%whole_freq_fit)
    endif
    allocate(ss%diag_freq_fit(ss%max_i,2*ss%n_grid_fit+1,ss%nspin))
     



    if(ionode) then
     
       do is=1,ss%nspin
          do iw=1,2*ss%n+1
             read(iun) ss%diag(1:ss%max_i,iw,is)
          end do
          if(ss%whole_s) then
             do iw=1,2*ss%n+1
                read(iun) ss%whole(ss%i_min_whole:ss%i_max_whole,1:ss%max_i,iw,is)
             end do
          endif
          
          
          do iw=1,2*ss%n_grid_fit+1
             read(iun) ss%diag_freq_fit(1:ss%max_i,iw,is)
          end do
          if(ss%whole_s) then
             do iw=1,2*ss%n_grid_fit+1
                read(iun) ss%whole_freq_fit(ss%i_min_whole:ss%i_max_whole,1:ss%max_i,iw,is)
             end do
          endif
       enddo
       close(iun)
    endif
   
    call mp_bcast(ss%diag, ionode_id,world_comm)
    if(ss%whole_s) then
       call mp_bcast(ss%whole, ionode_id,world_comm)
    endif
   
   
    call mp_bcast(ss%diag_freq_fit, ionode_id,world_comm)
    if(ss%whole_s) then
       call mp_bcast(ss%whole_freq_fit, ionode_id,world_comm)
    endif

    return
   END SUBROUTINE read_self_storage_ondisk


  SUBROUTINE free_memory_self_storage(ss)
!deallocate if allocated
    implicit none
    TYPE(self_storage) :: ss
    if(associated(ss%diag)) deallocate(ss%diag)
    nullify(ss%diag)
    if(associated(ss%whole)) deallocate(ss%whole)
    nullify(ss%whole)
    if(associated(ss%ene_remainder)) deallocate(ss%ene_remainder)
    nullify(ss%ene_remainder)
    if(associated(ss%diag_freq_fit)) deallocate(ss%diag_freq_fit)
    nullify(ss%diag_freq_fit)
    if(associated(ss%whole_freq_fit)) deallocate(ss%whole_freq_fit)
    nullify(ss%whole_freq_fit)

  END SUBROUTINE


  subroutine write_self_on_real(sr,ifile)
!this subroutine writes the self-energy  function on disk
!the file name is taken from the label 

    USE io_global,          ONLY : stdout, ionode
    USE mp,                 ONLY : mp_barrier
    USE mp_world,           ONLY : world_comm
     USE io_files,             ONLY : prefix,tmp_dir
    implicit none
    INTEGER, EXTERNAL :: find_free_unit
    TYPE(self_on_real),INTENT(in) :: sr!the self_energy descriptor to be written on file  
    INTEGER,INTENT(in) :: ifile!0 for integration part 1 for total part
    
    INTEGER :: iun

    if(ionode) then
       iun = find_free_unit()
       if(ifile==0) then
          open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'self_on_realA', status='unknown',form='unformatted')
       else
          open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'self_on_realB', status='unknown',form='unformatted')
       endif
       write(iun) sr%n
       write(iun) sr%max_i
       write(iun) sr%i_min
       write(iun) sr%i_max
       write(iun) sr%nspin
       write(iun) sr%grid(1:sr%n)
       write(iun) sr%diag(1:sr%n,1:sr%max_i,1:sr%nspin)
       close(iun)
    endif


    return
  end subroutine write_self_on_real

  subroutine read_self_on_real(sr,ifile)
!this subroutine reads the self-energy  function on disk
!the file name is taken from the label  

    USE io_global,          ONLY : stdout, ionode, ionode_id
    USE mp,                 ONLY : mp_bcast
    USE mp_world,           ONLY : world_comm
     USE io_files,             ONLY : prefix,tmp_dir
    implicit none
    INTEGER, EXTERNAL :: find_free_unit
    TYPE(self_on_real),INTENT(out) :: sr!the self_energy descriptor to be written on file                
    INTEGER,INTENT(in) :: ifile!0 for integration part 1 for total part                                              

    INTEGER :: iun

    if(ionode) then
       iun = find_free_unit()
       if(ifile==0) then
          open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'self_on_realA', status='old',form='unformatted')
       else
          open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'self_on_realB', status='old',form='unformatted')
       endif
       read(iun) sr%n
       read(iun) sr%max_i
       read(iun) sr%i_min
       read(iun) sr%i_max
       read(iun) sr%nspin
    endif
    call mp_bcast(sr%n, ionode_id,world_comm)
    call mp_bcast(sr%max_i,ionode_id,world_comm)
    call mp_bcast(sr%i_min,ionode_id,world_comm)
    call mp_bcast(sr%i_max,ionode_id,world_comm)
    call mp_bcast(sr%nspin,ionode_id,world_comm)
    allocate(sr%grid(sr%n))
    allocate(sr%diag(sr%n,sr%max_i,sr%nspin))
    if(ionode) then
       read(iun) sr%grid(1:sr%n)
       read(iun) sr%diag(1:sr%n,1:sr%max_i,1:sr%nspin)
       close(iun)
    endif
    call mp_bcast(sr%grid,ionode_id,world_comm)
    call mp_bcast(sr%diag,ionode_id,world_comm)

    return
  end subroutine read_self_on_real


  subroutine do_self_on_real(options,tf,ss,sr)
!this subroutine calculate the integral part of the self energy on real frequency axis
!at the end it calculates also the self_eenergy on imaginary axis
    USE io_global, ONLY : stdout
    USE input_gw,           ONLY : input_options
    USE times_gw,          ONLY : times_freqs


    implicit none
    
    TYPE(input_options) :: options
    TYPE(times_freqs)   :: tf
    TYPE(self_storage)  :: ss
    TYPE(self_on_real)  :: sr

    INTEGER :: is,ii,iw
    REAL(kind=DP) :: freq,energy

!set up sr
    sr%n=options%n_real_axis
    sr%max_i=options%max_i
    sr%i_min=options%i_min
    sr%i_max=options%i_max
    sr%nspin=options%nspin

    allocate(sr%grid(sr%n))
    allocate(sr%diag(sr%n,sr%max_i,sr%nspin))
    do iw=0,sr%n-1
       freq=(options%real_energy_max-options%real_energy_min)/dble(sr%n)*dble(iw)+options%real_energy_min
       sr%grid(iw+1)=dcmplx(freq,0.d0)
    enddo
       
!loop on frequencies
    do iw=1,sr%n
       energy=dble(sr%grid(iw))
       call do_self_lanczos_time(ss, tf ,options,.true.,energy)
!do fft
       call fft_storage_grid_fit(tf, ss)       
!extract data from ss
       do is=1,sr%nspin
          do ii=sr%i_min,sr%i_max
             sr%diag(iw,ii,is)=ss%diag_freq_fit(ii,ss%n_grid_fit+1,is)
          enddo
       enddo
       call free_memory_self_storage(ss)
    enddo
!anlytic continuation case
!energy should be at the middle of homo-lumo gap
    call do_self_lanczos_time(ss, tf ,options,.false.,energy)
    call fft_storage_grid_fit(tf, ss)
    return
  end subroutine do_self_on_real



  SUBROUTINE set_remainder(ss, qp)
!this subroutine simply copy the self-energy remainders
!from ss to qp, in order to allow restarting 
    USE energies_gww, ONLY : quasi_particles

    implicit none

    TYPE(self_storage) :: ss
    TYPE(quasi_particles) :: qp

    if(.not.associated(qp%ene_remainder))  allocate(qp%ene_remainder(ss%max_i,1))
    qp%ene_remainder(:,1)=ss%ene_remainder(:,1)
    
    return
  END SUBROUTINE set_remainder


  SUBROUTINE create_self_ontime(tf, ss,options,qp)
!this subroutine creates the structure self_storege
!on imaginary time
   USE  constants,         ONLY : eps8
   USE io_global,          ONLY : stdout, ionode 
   USE input_gw,           ONLY : input_options
   USE basic_structures,   ONLY : q_mat, wannier_u, wp_psi,v_pot,free_memory
   USE green_function,     ONLY : green,read_green,free_memory_green, initialize_green
   USE polarization,       ONLY : polaw,free_memory_polaw,read_polaw, initialize_polaw, &
                                   &invert_v_pot,distribute_v_pot, collect_v_pot
   USE compact_product
   USE mp,                 ONLY : mp_sum, mp_barrier
   USE para_gww,           ONLY : is_my_time, is_my_pola
   USE energies_gww,           ONLY : quasi_particles
   USE times_gw,           ONLY : times_freqs
   USE w_divergence
   USE mp_world,            ONLY : world_comm,nproc,mpime

   implicit none


   TYPE(times_freqs), INTENT(in) :: tf!for times grid
   TYPE(input_options), INTENT(in) :: options! for imaginary time range and number of samples
   TYPE(self_storage) :: ss!
   TYPE(quasi_particles), INTENT(in) :: qp!for the HF energies if required

   TYPE(green)     :: gg,gm!green function
   TYPE(q_mat)     :: qm!overlap of orthonormalized wannier products with wannier products
   TYPE(polaw)     :: ww!dressed interaction
   TYPE(wannier_u) :: uu!transformation matrix ks to wannier
   TYPE(contraction) :: cr!to speed up calculation
   TYPE(contraction_index) :: cri! index of contraction
   TYPE(contraction_state) :: crs!state contraction data
   TYPE(wp_psi) :: wp!for remainder calculations
   TYPE(gv_time) :: gt!for the treatment of the G=0,G=0 divergence of W
   TYPE(v_pot) :: vp,vpi,vpid



   REAL(kind=DP) :: time,dt
   INTEGER       :: iw,ii,jj
   REAL(kind=DP) :: offset
   COMPLEX(kind=DP) :: sene
   INTEGER :: l_blk, nbegin,nend
   REAL(kind=DP), ALLOCATABLE :: wtemp(:,:)


   nullify(vp%vmat)
   nullify(vpi%vmat)
   nullify(vpid%vmat)

   if(options%l_self_from_pola .or. options%l_self_beta) then
      if(options%w_divergence == 2) then
         call read_data_pw_v(vp,options%prefix,options%debug,0,.true.)
      else
         call read_data_pw_v(vp,options%prefix,options%debug,0,.false.)
      endif
      call invert_v_pot(vp,vpi)
      call free_memory(vp)
      call distribute_v_pot(vpi,vpid)
      call free_memory(vpi)
   endif

!set self_energy descriptor

!  call free_memory_self_storage(ss)
  ss%ontime=.true.
  ss%max_i=options%max_i
  ss%i_min=options%i_min
  ss%i_max=options%i_max
  ss%n=options%n
  ss%tau=options%tau
  ss%whole_s=options%whole_s
  ss%nspin=1

  if(tf%grid_fit/=0) then
     ss%n_grid_fit=tf%n_grid_fit
  else
     ss%n_grid_fit=tf%n
  endif

  if(ss%whole_s) then
   allocate(ss%whole(ss%max_i,ss%max_i,2*ss%n+1,1))
   ss%whole(:,:,:,:)=(0.d0,0.d0)
   allocate(ss%whole_freq_fit(ss%max_i,ss%max_i,2*ss%n_grid_fit+1,1))
   ss%whole_freq_fit(:,:,:,:)=(0.d0,0.d0)
   nullify(ss%diag)
   nullify(ss%diag_freq_fit)
  else
   allocate(ss%diag(ss%max_i,2*ss%n+1,1))
   ss%diag(:,:,:)=(0.d0,0.d0)
   nullify(ss%whole)
   allocate(ss%diag_freq_fit(ss%max_i,2*ss%n_grid_fit+1,1))
   ss%diag_freq_fit(:,:,:)=(0.d0,0.d0)
   nullify(ss%whole_freq_fit)
  endif

!set up self-energy remainders

  allocate(ss%ene_remainder(ss%max_i,1))
 
  ss%ene_remainder(:,1)=0.d0
  


  if(.not.options%lvcprim_file .and. .not.options%l_self_beta) then
  
!read U matrix
     call read_data_pw_u(uu,options%prefix) 
!read overlap matrix Q
     call read_data_pw_q(qm,options%prefix, options%l_self_from_pola)
     dt = ss%tau/real(ss%n)
     if(options%use_contractions) then
        if(.not.options%l_contraction_single_state) then
           write(stdout,*) 'call do_contraction'!ATTENZIONE
           call  do_contraction(qm,uu,cr, options%max_i)
           write(stdout,*) 'done do_contraction'!ATTENZIONE
           call  write_contraction(cr,options)
           write(stdout,*) 'done do_contraction'!ATTENZIONE
        else
!contraction index and states already available on disk
           call read_contraction_index(cri, options)
        endif
     endif
     !loop   
     call initialize_green(gg)
     call initialize_polaw(ww)

     l_blk= (2*ss%n+1)/nproc
     if(l_blk*nproc < (2*ss%n+1)) l_blk = l_blk+1
     nbegin=mpime*l_blk+1 -(ss%n+1)
     nend=nbegin+l_blk-1
     if(nend > ss%n) nend = ss%n

!     do iw=-ss%n,ss%n
!        if(is_my_time(iw)) then

     do iw=nbegin,nbegin+l_blk-1
        if(iw <= ss%n)  then

           write(stdout,*) 'Time :',iw!ATTENZIONE
           FLUSH(stdout)
           time=dt*real(iw)
!read dressed interaction
!we take care of the symmetry t ==> -t          
           call read_polaw(abs(iw),ww,options%debug,options%l_verbose)

!some controls
           if(.not. ww%ontime) then
              write(stdout,*) 'Routine create_self_ontime: imaginary time required'
              stop
           endif
           if(tf%l_fft_timefreq) then
              if(abs(time-ww%time) >= eps8) then
                 write(stdout,*) 'Routine create_self_ontime: imaginary time does not correspond'
                 stop
              endif
           endif
           
           if(options%l_self_from_pola) then
!if required  obtains the dressed polarization
              call collect_v_pot(vpi,vpid)
              allocate(wtemp(ww%numpw,ww%numpw))
              call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,vpi%vmat,ww%numpw,ww%pw,ww%numpw,&
                    &0.d0, wtemp,ww%numpw)
              call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,wtemp,ww%numpw,vpi%vmat,ww%numpw,&
                   &0.d0,ww%pw,ww%numpw)
              call free_memory(vpi)
              deallocate(wtemp)
           endif
           call read_green(iw,gg,options%debug,.false.)



!some controls
           if(.not. gg%ontime) then
              write(*,*) 'Routine create_self_ontime: imaginary time required'
              stop
           endif
           if(tf%l_fft_timefreq) then
              if(abs(time-gg%time) >= eps8) then
                 write(*,*) 'Routine create_self_ontime: imaginary time does not correspond'
                 stop
              endif
           endif
       !calculate elements
           if(ss%whole_s) then
              do ii=ss%i_min,ss%i_max
                 do jj=ss%i_min,ss%i_max
                    if(.not.options%use_contractions) then
                       call self_energy(ii,jj,ss%whole(ii,jj,iw+ss%n+1,1),time,qm,uu,gg,ww)
                    else
                       call self_energy_contraction(ii,jj,ss%whole(ii,jj,iw+ss%n+1,1),time,cr,gg,ww)
                    endif
                 enddo
              enddo
           else
              do ii=ss%i_min,ss%i_max
                 write(stdout,*) 'State:', ii
                 FLUSH(stdout)
                 if(.not.options%use_contractions) then
                    call self_energy(ii,ii,ss%diag(ii,iw+ss%n+1,1),time,qm,uu,gg,ww)
                 else
                    if(.not.options%l_contraction_single_state) then
                       call self_energy_contraction(ii,ii,ss%diag(ii,iw+ss%n+1,1),time,cr,gg,ww)
                    else
                       crs%state=ii
                       call read_contraction_state(cri,crs,options)
                       call self_energy_contraction_state(ii,ii,ss%diag(ii,iw+ss%n+1,1),time,cri,crs,gg,ww)
                       call free_memory_contraction_state(crs)
                    endif
                 endif
              enddo
           endif
           
!at zero time 1/2 positive G and 1/2 negative
           
           if(iw==0) then
              
              do ii=1,ss%max_i
                 ss%diag(ii,iw+ss%n+1,1)=0.5d0*ss%diag(ii,iw+ss%n+1,1)
              enddo
              call read_green(iw,gg,options%debug,.true.)
              do ii=ss%i_min,ss%i_max
                  write(stdout,*) 'State:', ii
                 FLUSH(stdout)
                 if(.not.options%use_contractions) then
                    call self_energy(ii,ii,sene,time,qm,uu,gg,ww)
                 else
                    if(.not.options%l_contraction_single_state) then
                       call self_energy_contraction(ii,ii,sene,time,cr,gg,ww)
                    else
                       crs%state=ii
                       call read_contraction_state(cri,crs,options)
                       call self_energy_contraction_state(ii,ii,sene,time,cri,crs,gg,ww)
                       call free_memory_contraction_state(crs)
                    endif
                 endif
                 ss%diag(ii,iw+ss%n+1,1)=ss%diag(ii,iw+ss%n+1,1)+0.5d0*sene
              enddo
           endif
        else
           if(options%l_self_from_pola) then
              call collect_v_pot(vpi,vpid)
              call free_memory(vpi)
           endif


        endif!on is_my_time
     enddo
  
     call free_memory(vpid)


    
    
    if(ss%whole_s) then
       call mp_sum(ss%whole(:,:,:,:),world_comm)
    else
       call mp_sum(ss%diag(:,:,:),world_comm)
    end if
    
    call free_memory(uu)
    call free_memory(qm)
    call free_memory_polaw(ww)
    call free_memory_green(gg)
    if(.not.options%l_contraction_single_state) &
         & call free_memory_contraction(cr)
  



 else
   !FROM VCPRIM FILE
    call selfenergy_ontime_file(ss,tf,options)
 endif

!if required add coulomb-like term for the treatment of the (G=0,G=0) divergence of W
 if(options%w_divergence == 2 ) then
    call initialize_gv_time(gt)
    call read_gv_time(gt)
    !consistency check
     if(options%max_i /= gt%max_i) then
        write(stdout,*) 'max_i not correct'
        stop
     endif
     call setup_gv_time(gt)
     
     do iw=1,2*gt%n+1
        ss%diag(:,iw,1)=ss%diag(:,iw,1)+gt%ex(:,iw)
     enddo
     call free_memory_gv_time(gt)
  endif

  return
 
END SUBROUTINE create_self_ontime


  SUBROUTINE write_storage(tf,ss)
!this subroutine write on standard output
!the values of write_storage
    USE io_global,      ONLY : stdout, ionode
    USE constants,      ONLY : pi
    USE mp,             ONLY : mp_barrier
    USE mp_world,       ONLY : world_comm
    USE times_gw,       ONLY : times_freqs

    implicit none

    TYPE(times_freqs), INTENT(in) :: tf!for time grid
    TYPE(self_storage), INTENT(in) :: ss
    INTEGER :: iw,ii,jj
    REAL(kind=DP) :: time,dt,totalfrequency,totalperiod,omega

    if(ionode) then
       if(ss%ontime) then
          write(stdout,*) '--------Sigma on imaginary time----------'
          dt=ss%tau/real(ss%n)
          do iw=-ss%n,ss%n
             if(tf%l_fft_timefreq) then
                time=dt*real(iw)
             else
                time=tf%times(iw)
             endif
             if(ss%whole_s) then
                do ii=1,ss%max_i
                   do jj=1,ss%max_i
                      write(stdout,*) time,ii,jj,ss%whole(ii,jj,iw+ss%n+1,1)
                   enddo
                enddo
             else
                do ii=1,ss%max_i
                   write(stdout,*) iw, time,ii, ss%diag(ii,iw+ss%n+1,1)
                enddo
             endif
          enddo
       else
          write(stdout,*) '--------Sigma on imaginary frequency----------'
          totalperiod=2.d0*ss%tau+2.d0*ss%tau/real(ss%n)
          totalfrequency=(2.d0*pi/totalperiod)
          
          do iw=-ss%n,ss%n
             if(tf%l_fft_timefreq) then
                omega=totalfrequency*real(iw)
             else
                omega=tf%freqs(iw)
             endif
             if(ss%whole_s) then
                do ii=1,ss%max_i
                   do jj=1,ss%max_i
                      write(stdout,*) omega,ii,jj,ss%whole(ii,jj,iw+ss%n+1,1)
                   enddo
                enddo
             else
                do ii=1,ss%max_i
                   write(stdout,*) omega,ii, ss%diag(ii,iw+ss%n+1,1)
                enddo
             endif
          enddo
       endif
    endif
    return

  END SUBROUTINE write_storage

  SUBROUTINE fft_storage_grid(tf,ss)
!this subroutine performs a FFT on the storage data
!inverse or direct determined by ontime
!uses grid

   USE io_global, ONLY : stdout
   USE constants, ONLY : pi
   USE times_gw,  ONLY : times_freqs


   implicit none

   TYPE(times_freqs), INTENT(in) :: tf!for time frequency grids
   TYPE(self_storage), INTENT(inout) :: ss!input data

   COMPLEX(kind=DP), DIMENSION(:), ALLOCATABLE :: ss_old, tmpc
   COMPLEX(kind=DP),  DIMENSION(:,:), ALLOCATABLE :: factors
   INTEGER :: ii,jj, is,js,kk


   INTEGER, PARAMETER :: nmesh=30

   REAL(kind=DP) :: b_p,b_m,r_p,r_m
   COMPLEX(kind=DP) :: a_p,a_m, cor_1,cor_2
   REAL(kind=DP), ALLOCATABLE :: x(:),w(:)
   COMPLEX(kind=DP), ALLOCATABLE ::  fij(:,:), fp(:),fm(:)


   allocate(ss_old(2*tf%n+1), tmpc(2*tf%n+1))
   allocate(factors(-tf%n:tf%n, -tf%n:tf%n))

!setup factors for every time posistion
   do ii=-tf%n, tf%n
      if(ss%ontime) then!time to frequency transform
         do jj=-tf%n,tf%n
             factors(jj,ii)=tf%weights_time(jj)*exp((0.d0,-1.d0)*tf%freqs(ii)*tf%times(jj))
         enddo
         factors(:,ii)=factors(:,ii)*(0.d0,-1.d0)
      else!frequency to time transform
         do jj=-tf%n,tf%n
            factors(jj,ii)=tf%weights_freq(jj)*exp((0.d0,1.d0)*tf%times(ii)*tf%freqs(jj))
         enddo
         factors(:,ii)=factors(:,ii)*(0.d0,1.d0)/(2.d0*pi)
      endif
   enddo
   if(ss%whole_s) then!full matrix
      do is=1,ss%max_i
         do js=1,ss%max_i
!copy array to be transformed
            ss_old(:)=ss%whole(is,js,:,1)
!transform
            do ii=-tf%n,tf%n
               do kk=-tf%n,tf%n
                  tmpc(kk+tf%n+1)=ss_old(kk+tf%n+1)*factors(kk,ii)
               enddo
               ss%whole(is,js,ii+tf%n+1,1)=sum(tmpc(1:2*tf%n+1))
            enddo
         enddo
      enddo
   else
      if(tf%l_fourier_fit_time .and. ss%ontime) then
         allocate(fij(-tf%n:tf%n,nmesh))
         allocate(fp(nmesh),fm(nmesh))
         allocate(x(nmesh),w(nmesh))
         x(:)=0.d0
         w(:)=0.d0
         call legzo(nmesh,x,w)

         !x(:)=x(:)*tf%tau/2.d0
         !x(:)=x(:)+tf%tau/2.d0
         !w(:)=w(:)*tf%tau/2.d0

         x(:)=x(:)*(tf%times(tf%n)-tf%tau)/2.d0
         x(:)=x(:)+(tf%times(tf%n)-tf%tau)/2.d0+tf%tau
         w(:)=w(:)*(tf%times(tf%n)-tf%tau)/2.d0

         do ii=-tf%n,tf%n
            do jj=1,nmesh
               fij(ii,jj)=exp((0.d0,-1.d0)*tf%freqs(ii)*x(jj))
            enddo
         enddo
      endif
      do is=1,ss%max_i
 !copy array to be transformed
         ss_old(:)=ss%diag(is,:,1)
!transform
         do ii=-tf%n,tf%n
            do kk=-tf%n,tf%n
               tmpc(kk+tf%n+1)=ss_old(kk+tf%n+1)*factors(kk,ii)
            enddo
            ss%diag(is,ii+tf%n+1,1)=sum(tmpc(1:2*tf%n+1))
         enddo

         if(tf%l_fourier_fit_time .and. ss%ontime) then
            r_p=dble(ss_old(2*tf%n)/ss_old(2*tf%n+1))
            write(stdout,*) 'RP',ss_old(2*tf%n),ss_old(2*tf%n-5)
            if(r_p <= 1.d0) r_p = tf%g_tau
            b_p=log(r_p)/(tf%times(tf%n)-tf%times(tf%n-1))
            a_p=ss_old(2*tf%n)/(exp(-b_p*tf%times(tf%n-1)))
            if(r_p == tf%g_tau) a_p=0.d0

            if(abs(ss_old(2)) > 1.d-10 .and. abs(ss_old(1)) > 1.d-10) then
               r_m=dble(ss_old(2)/ss_old(1))
               if(r_m <= 1.d0) r_m = tf%g_tau
               b_m=log(r_m)/(tf%times(-tf%n+1)-tf%times(-tf%n))
               a_m=ss_old(2)/(exp(b_m*tf%times(-tf%n+1)))
               if(r_m == tf%g_tau) a_m=0.d0
            else
               r_m=0.d0
               a_m=(0.d0,0.d0)
               b_m=0.d0
            endif

           
            do jj=1,nmesh
               fp(jj)=a_p*exp(-b_p*x(jj))*w(jj)
            enddo
            if(r_m /=0.d0) then
               do jj=1,nmesh
                  fm(jj)=a_m*exp(-b_m*x(jj))*w(jj)
               enddo
            endif

            do ii=-tf%n,tf%n
             !  cor_1=(0.d0,-1.d0)*(a_p/(b_p+(0.d0,1.d0)*tf%freqs(ii)))
             !  if(r_m /= 0.d0) then
             !     cor_1=cor_1+(0.d0,-1.d0)*(a_m/(b_m-(0.d0,1.d0)*tf%freqs(ii)))
             !  endif
               cor_2=0.d0
               do jj=1,nmesh
                  cor_2=cor_2-fij(ii,jj)*fp(jj)
                  if(r_m /=0.d0) then
                     cor_2=cor_2-conjg(fij(ii,jj))*fm(jj)
                  endif
               enddo
               cor_2=cor_2*(0.d0,-1.d0)
               ss%diag(is,ii+tf%n+1,1)=ss%diag(is,ii+tf%n+1,1)!-cor_2!+cor_1+cor_2
               write(stdout,*) 'COR2' , cor_2
            enddo
         endif
      enddo
      if(tf%l_fourier_fit_time .and.ss%ontime) deallocate(fij,fp,fm,x,w)
   endif
   if(ss%ontime) then
      ss%ontime=.false.
   else
      ss%ontime=.true.
   endif


   deallocate(ss_old,tmpc)
   deallocate(factors)

   return
 END SUBROUTINE fft_storage_grid

  SUBROUTINE fft_storage(ss)
!this subroutine performs a FFT on the storage data
!inverse or direct determined by ontime
   USE io_global, ONLY : stdout
   USE constants, ONLY : pi
   USE fft_scalar, ONLY : cft_1z

   implicit none

   TYPE(self_storage) :: ss!input data



   REAL(kind=DP) :: totalperiod,omega,time,totalfrequency
   INTEGER :: iw,ii,ipos
   COMPLEX(kind=DP), ALLOCATABLE :: inz(:),outz(:)
   COMPLEX(kind=DP) :: fact
   INTEGER*8 :: plan

   totalperiod=2.d0*ss%tau+2.d0*ss%tau/real(ss%n)
   totalfrequency=(2.d0*pi/totalperiod)*real(2*ss%n+2)
   allocate(inz(2*ss%n+2),outz(2*ss%n+2))
  
   if(.not.ss%whole_s) then
     if(ss%ontime) then!time to frequency transformation
       ss%ontime=.false.
!loop on states
       do ii=1,ss%max_i

         inz(:)=(0.d0,0.d0)
         do iw=-ss%n,ss%n
           ipos=iw+ss%n+2
           inz(ipos)=ss%diag(ii,iw+ss%n+1,1)
         enddo
         inz(1)=inz(2)

         call cft_1z(inz,1,2*ss%n+2,2*ss%n+2, -1,outz)
         outz(:)=outz(:)*dble(2*ss%n+2)

         do iw=0,2*ss%n+2-1
           if(iw <= (2*ss%n+1)) then
             omega=(2.d0*pi/totalperiod)*real(iw)
           else
             omega=(2.d0*pi/totalperiod)*real(iw-2*ss%n-2)
           endif
           fact=exp((0.d0,-1.d0)*omega*totalperiod/2.d0)*(0.d0,-1.d0)*(ss%tau/real(ss%n))
           outz(iw+1)=outz(iw+1)*fact
         enddo

         do iw=0,2*ss%n+1
           if(iw/=(ss%n+1)) then
             if(iw < (ss%n+1)) then
               ss%diag(ii,ss%n+iw+1,1)=outz(iw+1)
             else
               ss%diag(ii,iw-ss%n-2+1,1)=outz(iw+1)
             endif
           endif
         enddo
         write(*,*) 'ELIMINATO:', outz(ss%n+1)
       enddo

   else !frequency to time transform
       ss%ontime=.true.
!loop on states
       do ii=1,ss%max_i

         inz(:)=(0.d0,0.d0)
         do iw=-ss%n,ss%n
           ipos=iw+ss%n+2
           inz(ipos)=ss%diag(ii,iw+ss%n+1,1)
         enddo



         call cft_1z(inz,1,2*ss%n+2,2*ss%n+2, 1,outz)

         do iw=0,2*ss%n+2-1
           if(iw <= (2*ss%n+1)) then
             time=(ss%tau/real(ss%n))*real(iw)
           else
             time=(ss%tau/real(ss%n))*real(iw-2*ss%n-2)
           endif
           fact=exp((0.d0,+1.d0)*time*totalfrequency/2.d0)*(0.d0,+1.d0)/totalperiod
           outz(iw+1)=outz(iw+1)*fact
         enddo


         do iw=0,2*ss%n+1
           if(iw/=(ss%n+1)) then
             if(iw < (ss%n+1)) then
               ss%diag(ii,ss%n+iw+1,1)=outz(iw+1)
             else
               ss%diag(ii,iw-ss%n-2+1,1)=outz(iw+1)
             endif
           endif
         enddo

       enddo


   endif

   else
     write(stdout,*) 'ENTIRE SIGMA NOT IMPLEMENTED YET'
   endif
  deallocate(inz,outz)
  return
    
  END SUBROUTINE 

  SUBROUTINE test_fft(tf)
!just a fft test

  USE times_gw, ONLY : times_freqs

  implicit none

  TYPE(times_freqs), INTENT(in) :: tf

  TYPE(self_storage) :: ss
  INTEGER :: n,iw
  REAL(kind=DP) :: tau, lambda

  n=100
  tau=25.
  lambda=2.


  ss%ontime=.true.
  ss%whole_s=.false.
  ss%n=tf%n
  ss%tau=tf%tau
  ss%max_i=1
  n=tf%n
  allocate(ss%diag(1,2*n+1,1))
  nullify(ss%whole)
  
  do iw=-n,n
    ss%diag(1,iw+n+1,1)=exp(-(real(iw)*tau/real(n)/lambda)**2.)
  enddo
 
  call write_storage(tf,ss)

  call fft_storage(ss)
  
  call write_storage(tf,ss)

  call fft_storage(ss)

  call write_storage(tf,ss)

  call free_memory_self_storage(ss)
  return
  END SUBROUTINE  

  SUBROUTINE addconduction_self_ontime(ss, options)
!this subroutine adds to the self_energy of conduction states
!on negative imaginary times, the part due to terms \Psi_c'\Psic\w_P

    USE io_global,     ONLY : stdout, ionode
    USE input_gw,      ONLY : input_options
    USE basic_structures,  ONLY : v_pot,wannier_u_prim, v_pot_prim,free_memory, ortho_polaw
    USE green_function,    ONLY : green, read_green, free_memory_green, initialize_green
    USE polarization,      ONLY : polaw, free_memory_polaw, read_polaw, invert_v_pot, invert_ortho_polaw,&
         & orthonormalize_inverse, initialize_polaw, orthonormalize_vpot, distribute_ortho_polaw, collect_ortho_polaw,&
         &  distribute_v_pot, collect_v_pot
    USE mp,                ONLY : mp_sum
    USE para_gww,          ONLY : is_my_pola
    USE mp_world,          ONLY : world_comm,nproc,mpime
  
    implicit none

    TYPE(input_options) :: options
    TYPE(self_storage)  :: ss

    TYPE(v_pot) :: vp,vpi,vpid
    TYPE(ortho_polaw) :: op,opi, opd, opid
    TYPE(polaw) :: ww!dressed interaction
    TYPE(wannier_u_prim) :: wup
    TYPE(v_pot_prim) :: vpp
    TYPE(green) :: gg

    INTEGER iw,jw,kw,it,ii
    REAL(kind=DP), ALLOCATABLE :: wtemp(:,:)
    REAL(kind=DP), ALLOCATABLE :: cp(:,:,:) !arrys for contraction c',c, numpw
    REAL(kind=DP), ALLOCATABLE :: qg(:,:)
    COMPLEX(kind=DP), ALLOCATABLE :: sene(:,:)
    REAL(kind=DP), ALLOCATABLE :: gf_t(:,:)
    REAL(kind=DP), ALLOCATABLE :: pwcp_t(:,:)
    REAL(kind=DP), EXTERNAL :: ddot

    INTEGER :: l_blk, nbegin,nend
    INTEGER  :: i_first

    nullify(vp%vmat)
    nullify(vpi%vmat)
    nullify(op%on_mat)
    nullify(opi%on_mat)
    nullify(wup%umat)
    nullify(vpp%ij)
    nullify(vpp%vmat)
    nullify(opd%on_mat)
    nullify(opid%on_mat)
    nullify(vpid%vmat)

    call initialize_green(gg)
    call initialize_polaw(ww)
    write(stdout,*) 'addconduction_self_ontime OLD 1'!ATTENZIONE
    FLUSH(stdout)


!read coulombian potential and calculate inverse

    if(ss%whole_s) then
       write(stdout,*) 'Whole s not implemented YET'
       stop
    endif




    call read_data_pw_u_prim(wup,options%prefix)

    i_first=max(ss%i_min,wup%nums_occ+1)

    write(stdout,*) 'addconduction_self_ontime1_2'!ATTENZIONE
    FLUSH(stdout)

    if(options%w_divergence==2) then
       call read_data_pw_v_pot_prim(vpp, options%prefix,.true.)
    else
       call read_data_pw_v_pot_prim(vpp, options%prefix,.false.)
    endif


    write(stdout,*) 'addconduction_self_ontime1_3'!ATTENZIONE
     FLUSH(stdout)

   allocate(sene(-ss%n:0,ss%i_max-wup%nums_occ))
   sene(:,:)=(0.d0,0.d0)


!set up contraction array \sum_j U^{C'}_ij Vjkl

  !  allocate(cp(vpp%numpw, wup%nums-wup%nums_occ,options%max_i-wup%nums_occ))
   allocate(cp(vpp%numpw, wup%nums-wup%nums_occ,i_first:ss%i_max))
   cp(:,:,:)=0.d0
    do iw=1,vpp%numpw_prim
       do ii=i_first,ss%i_max
          do kw=1,vpp%numpw
             cp(kw,vpp%ij(2,iw)-wup%nums_occ,ii)=cp(kw,vpp%ij(2,iw)-wup%nums_occ,ii)+&
                  &dble(wup%umat(ii-wup%nums_occ,vpp%ij(1,iw)))*vpp%vmat(iw,kw)
          enddo
       enddo
    enddo

    call free_memory(vpp)
    call free_memory(wup)!in this way only the data is deallocated

     write(stdout,*) 'addconduction_self_ontime1_4'!ATTENZIONE
     FLUSH(stdout)

    if(options%w_divergence == 2) then
       call read_data_pw_v(vp,options%prefix,options%debug,0,.true.)
    else
       call read_data_pw_v(vp,options%prefix,options%debug,0,.false.)
    endif
    if(options%lnonorthogonal) then
       call read_data_pw_ortho_polaw(op,options%prefix)
       call orthonormalize_vpot(op,vp)
    endif
    call invert_v_pot(vp,vpi)
    call free_memory(vp)
    write(stdout,*) 'addconduction_self_ontime1_45'
    call distribute_v_pot(vpi,vpid)
    call free_memory(vpi)
    if(options%lnonorthogonal) then
       call invert_ortho_polaw(op,opi)
       write(stdout,*) 'addconduction_self_ontime1_5 op',op%numpw!ATTENZIONE
       call distribute_ortho_polaw(op,opd)
       call free_memory(op)
       write(stdout,*) 'addconduction_self_ontime1_6 opd',opd%numpw!ATTENZIONE

       call distribute_ortho_polaw(opi,opid)
       call free_memory(opi)
    endif


    l_blk= (ss%n+1)/nproc
     if(l_blk*nproc < (ss%n+1)) l_blk = l_blk+1
     nbegin=mpime*l_blk+1 -(ss%n+1)
     nend=nbegin+l_blk-1
     if(nend > 0) nend = 0




    write(stdout,*) 'addconduction_self_ontime5',nbegin,l_blk!ATTENZIONE
    FLUSH(stdout)
!loop on negative imaginary times
    do it=nbegin,nbegin+l_blk-1
       if(it <= 0)  then
          write(stdout,*) 'addconduction_self_ontime time', it!ATTENZIONE
          FLUSH(stdout)
!we take care of the symmetru t ==> -t

          call read_polaw(abs(it),ww,options%debug,options%l_verbose)
          write(stdout,*) 'addconduction_self_ontime6 ww', ww%numpw!ATTENZIONE
           if(options%lnonorthogonal) then 
              call collect_ortho_polaw(opi,opid)
              write(stdout,*) 'dimensions', opi%numpw, opid%numpw
              call orthonormalize_inverse(opi,ww)
              call free_memory(opi)
           endif
           write(stdout,*) 'addconduction_self_ontime7'!ATTENZIONE
          allocate(wtemp(ww%numpw,ww%numpw))

          call collect_v_pot(vpi,vpid)

          call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,vpi%vmat,ww%numpw,ww%pw,ww%numpw,&
               &0.d0, wtemp,ww%numpw)



          call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,wtemp,ww%numpw,vpi%vmat,ww%numpw,&
              &0.d0,ww%pw,ww%numpw)

          call free_memory(vpi)

          deallocate(wtemp)
          
          if(options%lnonorthogonal) then 
             call collect_ortho_polaw(op,opd)
             call orthonormalize_inverse(op,ww)
             call free_memory(op)
          endif
           write(stdout,*) 'addconduction_self_ontime8'!ATTENZIONE
           FLUSH(stdout)
          call read_green(it,gg,options%debug,.true.)

          allocate(gf_t(wup%nums-wup%nums_occ,wup%nums-wup%nums_occ))
          do iw=1,(wup%nums-wup%nums_occ)
             do jw=1,(wup%nums-wup%nums_occ)
                gf_t(jw,iw)=gg%gf_p(jw+wup%nums_occ, iw+wup%nums_occ,1)
             enddo
          enddo

          do ii=i_first,ss%i_max
             write(stdout,*) 'II' , ii 
             FLUSH(stdout)
             allocate(qg(ww%numpw,wup%nums-wup%nums_occ))
             call dgemm('N','N',ww%numpw,wup%nums-wup%nums_occ,wup%nums-wup%nums_occ,1.d0,&
                & cp(:,:,ii),ww%numpw,gf_t,wup%nums-wup%nums_occ,0.d0,qg,ww%numpw)


             allocate(pwcp_t(ww%numpw,wup%nums-wup%nums_occ))
             call dgemm('N','N',ww%numpw,wup%nums-wup%nums_occ,ww%numpw,1.d0,&
              &ww%pw,ww%numpw,cp(:,:,ii),ww%numpw,0.d0,pwcp_t,ww%numpw)

             do iw=1,(wup%nums-wup%nums_occ)
                sene(it,ii-wup%nums_occ)=sene(it,ii-wup%nums_occ)+&
                     &ddot(ww%numpw,qg(:,iw),1,pwcp_t(:,iw),1)*gg%factor*ww%factor
             enddo

             deallocate(pwcp_t)
             deallocate(qg)

             sene(it,ii-wup%nums_occ)=sene(it,ii-wup%nums_occ)*(0.d0,1.d0)
             if(it==0) sene(it,ii-wup%nums_occ)=sene(it,ii-wup%nums_occ)*0.5d0
             write(stdout,*) 'Conduction contribution', it,ii, sene(it,ii-wup%nums_occ)
          enddo
          deallocate(gf_t)
       else
          if(options%lnonorthogonal) then 
             call collect_ortho_polaw(opi,opid)
             call free_memory(opi)
          endif
          call collect_v_pot(vpi,vpid)
          call free_memory(vpi)
          if(options%lnonorthogonal) then 
             call collect_ortho_polaw(op,opd)
             call free_memory(op)          
          endif
       endif
    enddo
    call mp_sum(sene(-ss%n:0,:),world_comm)
    do ii=1,ss%i_max-wup%nums_occ
       do it=-ss%n,0
          ss%diag(ii+wup%nums_occ,it+ss%n+1,1)=ss%diag(ii+wup%nums_occ, it+ss%n+1,1)+sene(it,ii)
       enddo
    enddo
!!!!!!!!!!!
    call free_memory(vpid)
    if(options%lnonorthogonal) then 
       call free_memory(opd)
       call free_memory(opi)
       call free_memory(opid)
    endif
    call free_memory_polaw(ww)
    call free_memory_green(gg)
    deallocate(cp)
    deallocate(sene)
    return
  END SUBROUTINE addconduction_self_ontime

  SUBROUTINE fft_storage_grid_fit(tf,ss)
!this subroutine performs a FFT from time to frequency on the storage data
!from W,P grid to fit grid
!in case also for diagonal elements

   USE io_global, ONLY : stdout
   USE constants, ONLY : pi
   USE times_gw,  ONLY : times_freqs


   implicit none

   TYPE(times_freqs), INTENT(in) :: tf!for time frequency grids
   TYPE(self_storage), INTENT(inout) :: ss!input data

   COMPLEX(kind=DP), DIMENSION(:), ALLOCATABLE :: tmpc
   COMPLEX(kind=DP),  DIMENSION(:,:), ALLOCATABLE :: factors
   INTEGER :: ii,jj, is,js,kk


   INTEGER, PARAMETER :: nmesh=30
   

   REAL(kind=DP) :: b_p,b_m,r_p,r_m
   COMPLEX(kind=DP) :: a_p,a_m, cor_1,cor_2
   REAL(kind=DP), ALLOCATABLE :: x(:),w(:)
   COMPLEX(kind=DP), ALLOCATABLE ::  fij(:,:), fp(:),fm(:)
   INTEGER :: ispin




   allocate(tmpc(2*tf%n+1))
   allocate(factors(-tf%n:tf%n, -tf%n_grid_fit:tf%n_grid_fit))

   


!setup factors for every time position
   do ii=-tf%n_grid_fit, tf%n_grid_fit
      do jj=-tf%n,tf%n
         factors(jj,ii)=tf%weights_time(jj)*exp((0.d0,-1.d0)*tf%freqs_fit(ii)*tf%times(jj))
      enddo
      factors(:,ii)=factors(:,ii)*(0.d0,-1.d0)
   enddo
   if(tf%l_fourier_fit_time .and. ss%ontime) then
      allocate(fij(-tf%n_grid_fit:tf%n_grid_fit,nmesh))
      allocate(fp(nmesh),fm(nmesh))
      allocate(x(nmesh),w(nmesh))
      x(:)=0.d0
      w(:)=0.d0
      call legzo(nmesh,x,w)
      
      x(:)=x(:)*(tf%times(tf%n)-tf%tau)/2.d0
      x(:)=x(:)+(tf%times(tf%n)-tf%tau)/2.d0+tf%tau
      w(:)=w(:)*(tf%times(tf%n)-tf%tau)/2.d0

      do ii=-tf%n_grid_fit,tf%n_grid_fit
         do jj=1,nmesh
            fij(ii,jj)=exp((0.d0,-1.d0)*tf%freqs_fit(ii)*x(jj))
         enddo
      enddo
   endif

   do ispin=1,ss%nspin
      do is=1,ss%max_i

!transform
         do ii=-tf%n_grid_fit,tf%n_grid_fit
            do kk=-tf%n,tf%n
               tmpc(kk+tf%n+1)=ss%diag(is,kk+tf%n+1,ispin)*factors(kk,ii)
            enddo
            ss%diag_freq_fit(is,ii+tf%n_grid_fit+1,ispin)=sum(tmpc(1:2*tf%n+1))
         enddo
      
         if(ss%whole_s) then
            do js=ss%i_min_whole,ss%i_max_whole
               do ii=-tf%n_grid_fit,tf%n_grid_fit
                  do kk=-tf%n,tf%n
                     tmpc(kk+tf%n+1)=ss%whole(js,is,kk+tf%n+1,ispin)*factors(kk,ii)
                  enddo
                  ss%whole_freq_fit(js,is,ii+tf%n_grid_fit+1,ispin)=sum(tmpc(1:2*tf%n+1))
               enddo
            enddo
         endif
         

         if(tf%l_fourier_fit_time .and. ss%ontime) then
            r_p=dble(ss%diag(is,2*tf%n,1)/ss%diag(is,2*tf%n+1,1))
            if(r_p <= 1.d0) r_p = tf%g_tau
            b_p=log(r_p)/(tf%times(tf%n)-tf%times(tf%n-1))
            a_p=ss%diag(is,2*tf%n,1)/(exp(-b_p*tf%times(tf%n-1)))
            if(r_p == tf%g_tau) a_p=0.d0
         
            if(abs(ss%diag(is,2,1)) > 1.d-10 .and. abs(ss%diag(is,1,1)) > 1.d-10) then
               r_m=dble(ss%diag(is,2,1)/ss%diag(is,1,1))
               if(r_m <= 1.d0) r_m = tf%g_tau
               b_m=log(r_m)/(tf%times(-tf%n+1)-tf%times(-tf%n))
               a_m=ss%diag(is,2,1)/(exp(b_m*tf%times(-tf%n+1)))
               if(r_m == tf%g_tau) a_m=0.d0
            else
               r_m=0.d0
               a_m=(0.d0,0.d0)
               b_m=0.d0
            endif
            

            do jj=1,nmesh
               fp(jj)=a_p*exp(-b_p*x(jj))*w(jj)
            enddo
            if(r_m /=0.d0) then
               do jj=1,nmesh
                  fm(jj)=a_m*exp(-b_m*x(jj))*w(jj)
               enddo
            endif
            
            do ii=-tf%n_grid_fit,tf%n_grid_fit
               cor_2=0.d0
               do jj=1,nmesh
                  cor_2=cor_2-fij(ii,jj)*fp(jj)
                  if(r_m /=0.d0) then
                     cor_2=cor_2-conjg(fij(ii,jj))*fm(jj)
                  endif
               enddo
               cor_2=cor_2*(0.d0,-1.d0)
               ss%diag_freq_fit(is,ii+tf%n_grid_fit+1,1)=ss%diag_freq_fit(is,ii+tf%n_grid_fit+1,1)!-cor_2!+cor_1+cor_2
            enddo
         endif
      enddo
   enddo
   if(tf%l_fourier_fit_time .and.ss%ontime) deallocate(fij,fp,fm,x,w)
   
   
   deallocate(tmpc)
   deallocate(factors)
   

   return
 END SUBROUTINE fft_storage_grid_fit

  SUBROUTINE addconduction_self_ontime_file(ss, tf ,options)
!this subroutine adds to the self_energy of conduction states
!on negative imaginary times, the part due to terms \Psi_c'\Psic\w_P
!using terms from file

    USE io_global,     ONLY : stdout, ionode
    USE input_gw,      ONLY : input_options
    USE basic_structures,  ONLY : v_pot,wannier_u,free_memory, ortho_polaw,initialize_memory,cprim_prod
    USE green_function,    ONLY : green, read_green, free_memory_green, initialize_green
    USE polarization,      ONLY : polaw, free_memory_polaw, read_polaw, invert_v_pot, invert_ortho_polaw,&
         & orthonormalize_inverse, initialize_polaw, orthonormalize_vpot, distribute_ortho_polaw, collect_ortho_polaw,&
         &  distribute_v_pot, collect_v_pot
    USE mp,                ONLY : mp_sum
    USE para_gww,          ONLY : is_my_pola
    USE mp_world,            ONLY : world_comm,nproc,mpime
    USE times_gw,  ONLY : times_freqs


   implicit none

   TYPE(times_freqs), INTENT(in) :: tf!for time frequency grids
   TYPE(input_options) :: options
   TYPE(self_storage)  :: ss

   TYPE(v_pot) :: vp,vpi,vpid
   TYPE(ortho_polaw) :: op,opi, opd, opid
   TYPE(polaw) :: ww!dressed interaction
   TYPE(wannier_u) :: wu!structure to be read and initialized
   TYPE(cprim_prod) :: cpp!the producs c' c' v  wp

   INTEGER iw,jw,kw,it,ii,jj
   INTEGER :: l_blk, nbegin,nend
   REAL(kind=DP) :: offset
   COMPLEX(kind=DP), ALLOCATABLE :: sene(:,:)
   REAL(kind=DP), ALLOCATABLE:: wtemp(:,:), vtemp(:)
   REAL(kind=DP), EXTERNAL :: ddot
   LOGICAL :: ok_read
   
   nullify(vp%vmat)
   nullify(vpi%vmat)
   nullify(op%on_mat)
   nullify(opi%on_mat)
   nullify(opd%on_mat)
   nullify(opid%on_mat)
   nullify(vpid%vmat)
   nullify(wu%umat)

   call initialize_memory(cpp)

   call read_data_pw_u(wu, options%prefix)
   deallocate(wu%umat)

   if(.not.options%l_hf_energies) then
      if(wu%nums > wu%nums_occ(1)) then
         offset=-(wu%ene(wu%nums_occ(1)+1,1)+wu%ene(wu%nums_occ(1),1))/2.d0
      else
         offset=-wu%ene(wu%nums_occ(1),1)
      endif
   else
      write(stdout,*) 'HF energies to be implemented YET'
      stop
      !if(wu%nums > wu%nums_occ(1)) then
      !   offset=-(ene_hf(wu%nums_occ(1)+1)+ene_hf(wu%nums_occ(1)))/2.d0
      !else
      !   offset=-ene_hf(wu%nums_occ(1))
      !endif
   endif



   call initialize_polaw(ww)
   write(stdout,*) 'addconduction_self_ontime_file1'!ATTENZIONE

!read coulombian potential and calculate inverse

   if(ss%whole_s) then
      write(stdout,*) 'Whole s not implemented YET'
      stop
   endif


   if(options%w_divergence == 2) then
      call read_data_pw_v(vp,options%prefix,options%debug,0,.true.)
   else
      call read_data_pw_v(vp,options%prefix,options%debug,0,.false.)
   endif
   call read_data_pw_ortho_polaw(op,options%prefix)
   call orthonormalize_vpot(op,vp)
   call invert_v_pot(vp,vpi)
   call free_memory(vp)
   write(stdout,*) 'addconduction_self_ontime1_45'
   call distribute_v_pot(vpi,vpid)
   call free_memory(vpi)
   call invert_ortho_polaw(op,opi)

   write(stdout,*) 'addconduction_self_ontime1_5 op',op%numpw!ATTENZIONE
   call distribute_ortho_polaw(op,opd)
   call free_memory(op)
   write(stdout,*) 'addconduction_self_ontime1_6 opd',opd%numpw!ATTENZIONE
   call distribute_ortho_polaw(opi,opid)
   call free_memory(opi)

   l_blk= (ss%n+1)/nproc
   if(l_blk*nproc < (ss%n+1)) l_blk = l_blk+1
   nbegin=mpime*l_blk+1 -(ss%n+1)
   nend=nbegin+l_blk-1
   if(nend > 0) nend = 0
   
   
   write(stdout,*) 'addconduction_self_ontime5',nbegin,l_blk!ATTENZIONE
   FLUSH(stdout)

   allocate(sene(-ss%n:0,options%i_min:options%i_max))
   sene(:,:)=(0.d0,0.d0)



!loop on negative imaginary times
   do it=nbegin,nbegin+l_blk-1
      if(it <= 0)  then
         write(stdout,*) 'addconduction_self_ontime time', it!ATTENZIONE
!we take care of the symmetru t ==> -t

         call read_polaw(abs(it),ww,options%debug,options%l_verbose)
         write(stdout,*) 'addconduction_self_ontime6'!ATTENZIONE
         FLUSH(stdout)
         call collect_ortho_polaw(opi,opid)
         write(stdout,*) 'addconduction_self_ontime6.1'!ATTENZIONE
         call orthonormalize_inverse(opi,ww)
         write(stdout,*) 'addconduction_self_ontime6.2'!ATTENZIONE
         call free_memory(opi)
         write(stdout,*) 'addconduction_self_ontime7'!ATTENZIONE
         FLUSH(stdout)
         allocate(wtemp(ww%numpw,ww%numpw))

         call collect_v_pot(vpi,vpid)
            
         call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,vpi%vmat,ww%numpw,ww%pw,ww%numpw,&
                 &0.d0, wtemp,ww%numpw)



         call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,wtemp,ww%numpw,vpi%vmat,ww%numpw,&
              &0.d0,ww%pw,ww%numpw)

         call free_memory(vpi)

         deallocate(wtemp)

         call collect_ortho_polaw(op,opd)
         call orthonormalize_inverse(op,ww)
         call free_memory(op)

!!now ww contains \tilde{ww}    
         write(stdout,*) 'addconduction_self_ontime8'!ATTENZIONE
         FLUSH(stdout)
        
!read in cprim_prod
!first multiplication
!second multiplication
!copy on sene
            !loop on c' states
         do ii=max(options%i_min,wu%nums_occ(1)+1),options%i_max
            cpp%cprim=ii
            call read_data_pw_cprim_prod(cpp, options%prefix,.false.,ok_read,.false.,.false.)

!loop on c
            allocate(vtemp(cpp%numpw))
            do jj=1,cpp%nums_cond
            
!multiply W_ijS_jc =T_ic
               call dgemv('N',ww%numpw,ww%numpw,1.d0,ww%pw,ww%numpw,cpp%cpmat(:,jj),1,0.d0,vtemp,1)
!multiply S_icTi_c
               sene(it,ii)=sene(it,ii)+ddot(cpp%numpw,vtemp,1,cpp%cpmat(:,jj),1)*&
                    & exp((wu%ene(jj+wu%nums_occ(1),1)+offset)*tf%times(it))*ww%factor*(0.d0,-1.d0)
         
            enddo
            sene(it,ii)=sene(it,ii)*(0.d0,1.d0)
            if(it==0) sene(it,ii)=sene(it,ii)*0.5d0
            write(stdout,*) 'Conduction contribution', it,ii, sene(it,ii)
            FLUSH(stdout)
            deallocate(vtemp)
         enddo
      else
         call collect_ortho_polaw(opi,opid)
         call free_memory(opi)
         call collect_v_pot(vpi,vpid)
         call free_memory(vpi)
         call collect_ortho_polaw(op,opd)
         call free_memory(op)
         do ii=max(options%i_min,wu%nums_occ(1)+1),options%i_max
            cpp%cprim=ii
            call read_data_pw_cprim_prod(cpp, options%prefix,.false.,ok_read,.false.,.false.)
         enddo
      endif
   enddo
   call mp_sum(sene(-ss%n:0,:),world_comm)
   do ii=max(options%i_min,wu%nums_occ(1)+1),options%i_max
       do it=-ss%n,0
          ss%diag(ii,it+ss%n+1,1)=ss%diag(ii, it+ss%n+1,1)+sene(it,ii)
       enddo
    enddo
   
 !copy sene results on ss with opportune factors
!!!!!!!!!!!
   call free_memory(vpid)
   call free_memory(opd)
   call free_memory(opi)
   call free_memory(opid)
   call free_memory_polaw(ww)
   call free_memory( cpp)
   deallocate(sene)
   return
 END SUBROUTINE addconduction_self_ontime_file



  SUBROUTINE selfenergy_ontime_file(ss, tf ,options)
!this subroutine calculates the self_energy of selected states
!using terms from file or from strategy BETA

    USE io_global,     ONLY : stdout, ionode
    USE input_gw,      ONLY : input_options
    USE basic_structures,  ONLY : v_pot,wannier_u,free_memory, ortho_polaw,initialize_memory,cprim_prod,q_mat,&
                                        & wannier_u_prim,v_pot_prim
    USE green_function,    ONLY : green, read_green, free_memory_green, initialize_green
    USE polarization,      ONLY : polaw, free_memory_polaw, read_polaw, invert_v_pot, invert_ortho_polaw,&
         & orthonormalize_inverse, initialize_polaw, orthonormalize_vpot, distribute_ortho_polaw, collect_ortho_polaw,&
         &  distribute_v_pot, collect_v_pot
    USE mp,                ONLY : mp_sum, mp_barrier
    USE para_gww,          ONLY : is_my_pola
    USE mp_world,            ONLY : world_comm,nproc,mpime
    USE times_gw,  ONLY : times_freqs


   implicit none

   TYPE(times_freqs), INTENT(in) :: tf!for time frequency grids
   TYPE(input_options) :: options
   TYPE(self_storage)  :: ss

   TYPE(v_pot) :: vp,vpi,vpid
   TYPE(ortho_polaw) :: op,opi, opd, opid
   TYPE(polaw) :: ww!dressed interaction
   TYPE(wannier_u) :: wu!structure to be read and initialized
   TYPE(cprim_prod) :: cpp,cppd!the producs c' c' v  wp
   TYPE(q_mat) :: qm, qmd!for strategy beta
   TYPE(wannier_u_prim) :: wup!for strategy beta
   TYPE(v_pot_prim) :: vpp,vppd!for strategy beta


   INTEGER iw,jw,kw,it,ii,jj
   INTEGER :: l_blk, nbegin,nend
   REAL(kind=DP) :: offset
   COMPLEX(kind=DP), ALLOCATABLE :: sene(:,:)
   REAL(kind=DP), ALLOCATABLE:: wtemp(:,:), vtemp(:,:)
   REAL(kind=DP), EXTERNAL :: ddot
   LOGICAL :: ok_read

   nullify(vp%vmat)
   nullify(vpi%vmat)
   nullify(op%on_mat)
   nullify(opi%on_mat)
   nullify(opd%on_mat)
   nullify(opid%on_mat)
   nullify(vpid%vmat)
   nullify(wu%umat)
   nullify(wup%umat)
   nullify(vpp%vmat)

   call initialize_memory(cpp)
   call initialize_memory(cppd)


   if(options%l_self_beta) ok_read=.true.
   call read_data_pw_u(wu, options%prefix)
   if(.not.options%l_self_beta) deallocate(wu%umat)

   if(.not.options%l_hf_energies) then
      if(wu%nums > wu%nums_occ(1)) then
         offset=-(wu%ene(wu%nums_occ(1)+1,1)+wu%ene(wu%nums_occ(1),1))/2.d0
      else
         offset=-wu%ene(wu%nums_occ(1),1)
      endif
   else
      write(stdout,*) 'HF energies to be implemented YET'
      stop
      !if(wu%nums > wu%nums_occ(1)) then
      !   offset=-(ene_hf(wu%nums_occ(1)+1)+ene_hf(wu%nums_occ(1)))/2.d0
      !else
      !   offset=-ene_hf(wu%nums_occ(1))
      !endif
   endif



   call initialize_polaw(ww)
   write(stdout,*) 'addconduction_self_ontime1'!ATTENZIONE
   FLUSH(stdout)

!read coulombian potential and calculate inverse

   if(ss%whole_s) then
      write(stdout,*) 'Whole s not implemented YET'
      stop
   endif


   if(options%w_divergence == 2) then
      call read_data_pw_v(vp,options%prefix,options%debug,0,.true.)
   else
      call read_data_pw_v(vp,options%prefix,options%debug,0,.false.)
   endif
   if(options%lnonorthogonal) then
      call read_data_pw_ortho_polaw(op,options%prefix)
      call orthonormalize_vpot(op,vp)
   endif
   call invert_v_pot(vp,vpi)
   call free_memory(vp)
   write(stdout,*) 'addconduction_self_ontime1_45'
   FLUSH(stdout)
   call distribute_v_pot(vpi,vpid)
   call free_memory(vpi)
   if(options%lnonorthogonal) then
      call invert_ortho_polaw(op,opi)
   endif

   write(stdout,*) 'addconduction_self_ontime1_5 op',op%numpw!ATTENZIONE
   FLUSH(stdout)
   if(options%lnonorthogonal) then
      call distribute_ortho_polaw(op,opd)
      call free_memory(op)
  
      write(stdout,*) 'addconduction_self_ontime1_6 opd',opd%numpw!ATTENZIONE
      FLUSH(stdout)
      call distribute_ortho_polaw(opi,opid)
      call free_memory(opi)
   endif
   l_blk= (2*ss%n+1)/nproc
   if(l_blk*nproc < (2*ss%n+1)) l_blk = l_blk+1
   nbegin=mpime*l_blk+1 -(ss%n+1)
   nend=nbegin+l_blk-1
   if(nend > ss%n) nend = ss%n


   write(stdout,*) 'addconduction_self_ontime5',nbegin,l_blk!ATTENZIONE
   FLUSH(stdout)
   allocate(sene(-ss%n:ss%n,options%i_min:options%i_max))
   sene(:,:)=(0.d0,0.d0)

!if required read and distribute q_mat
   if(options%l_self_beta) then
      call read_data_pw_q(qm,options%prefix,.true.)
      call distribute_qmat(qm,qmd)
      call free_memory(qm)
      if(options%i_max > wu%nums_occ(1)) then
         if(options%w_divergence == 2) then
            call read_data_pw_v_pot_prim(vpp,options%prefix, .true.)
         else
            call read_data_pw_v_pot_prim(vpp,options%prefix, .false.)
         endif
         call distribute_v_pot_prim(vpp,vppd)
         call free_memory(vpp)
         call read_data_pw_u_prim(wup,options%prefix)
      endif
   endif
      


!loop on negative imaginary times
   do it=nbegin,nbegin+l_blk-1
      if(it <= ss%n)  then
         write(stdout,*) 'addconduction_self_ontime time', it!ATTENZIONE
         FLUSH(stdout)
!we take care of the symmetru t ==> -t

         call read_polaw(abs(it),ww,options%debug,options%l_verbose)
         write(stdout,*) 'addconduction_self_ontime6'!ATTENZIONE
         FLUSH(stdout)
         if(options%lnonorthogonal) then
            call collect_ortho_polaw(opi,opid)
            write(stdout,*) 'addconduction_self_ontime6.1'!ATTENZIONE
            call orthonormalize_inverse(opi,ww)
            write(stdout,*) 'addconduction_self_ontime6.2'!ATTENZIONE
            call free_memory(opi)
         endif
         write(stdout,*) 'addconduction_self_ontime7'!ATTENZIONE
         FLUSH(stdout)
         allocate(wtemp(ww%numpw,ww%numpw))
         
         call collect_v_pot(vpi,vpid)

         call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,vpi%vmat,ww%numpw,ww%pw,ww%numpw,&
                 &0.d0, wtemp,ww%numpw)



         call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,wtemp,ww%numpw,vpi%vmat,ww%numpw,&
              &0.d0,ww%pw,ww%numpw)

         call free_memory(vpi)

         deallocate(wtemp)

         if(options%lnonorthogonal) then
            call collect_ortho_polaw(op,opd)
            call orthonormalize_inverse(op,ww)
            call free_memory(op)
         endif

!!now ww contains \tilde{ww}
         write(stdout,*) 'addconduction_self_ontime8'!ATTENZIONE
         FLUSH(stdout)

!read in cprim_prod
!first multiplication
!second multiplication
!copy on sene
            !loop on c' states
         do ii=options%i_min,options%i_max
            cpp%cprim=ii
            call mp_barrier( world_comm )
            if(.not.options%l_self_beta) then
               call read_data_pw_cprim_prod(cpp, options%prefix,.true.,ok_read,.false.,.false.)
            else
!read qmat
               call  create_vcprim(cppd, ii ,wu, qmd)
               if(ii>wu%nums_occ(1)) then
!if required adds the conduction term
                  call add_vcprim_conduction(cppd, wu, wup, vppd)
               end if
               call collect_cprim_prod(cpp,cppd)
               call free_memory(cppd)
            endif
            if(ok_read) then

!loop on c
               allocate(vtemp(cpp%numpw,max(cpp%nums_occ,cpp%nums-cpp%nums_occ)))
               if(it <= 0) then
                  call dgemm('N','N',ww%numpw,cpp%nums-cpp%nums_occ,ww%numpw,1.d0,ww%pw,ww%numpw,&
                       & cpp%cpmat(:,cpp%nums_occ+1:cpp%nums),cpp%lda,0.d0,vtemp,ww%numpw)
                  
                  do jj=cpp%nums_occ+1,cpp%nums

!multiply W_ijS_jc =T_ic
!                     call dgemv('N',ww%numpw,ww%numpw,1.d0,ww%pw,ww%numpw,cpp%cpmat(:,jj),1,0.d0,vtemp,1)
                    
!multiply S_icTi_c
                     sene(it,ii)=sene(it,ii)+ddot(cpp%numpw,vtemp(:,jj-cpp%nums_occ),1,cpp%cpmat(:,jj),1)*&
                          & exp((wu%ene(jj,1)+offset)*tf%times(it))*ww%factor*(0.d0,-1.d0)

                  enddo
                  sene(it,ii)=sene(it,ii)*(0.d0,1.d0)
                  write(stdout,*) 'Conduction contribution', it,ii, sene(it,ii)
                  FLUSH(stdout)
               endif
               if(it >= 0) then
                   call dgemm('N','N',ww%numpw,cpp%nums_occ,ww%numpw,1.d0,ww%pw,ww%numpw,&
                       &cpp%cpmat(:,1:cpp%nums_occ),cpp%lda,0.d0,vtemp,ww%numpw)

                  do jj=1,cpp%nums_occ

!multiply W_ijS_jc =T_ic
!                     call dgemv('N',ww%numpw,ww%numpw,1.d0,ww%pw,ww%numpw,cpp%cpmat(:,jj),1,0.d0,vtemp,1)
!multiply S_icTi_c
                     sene(it,ii)=sene(it,ii)+ddot(cpp%numpw,vtemp(:,jj),1,cpp%cpmat(:,jj),1)*&
                          & exp((wu%ene(jj,1)+offset)*tf%times(it))*ww%factor*(0.d0,+1.d0)

                  enddo
                  sene(it,ii)=sene(it,ii)*(0.d0,1.d0)
                  write(stdout,*) 'Conduction contribution', it,ii, sene(it,ii)
                  FLUSH(stdout)
               endif
               if(it==0) sene(it,ii)=sene(it,ii)*0.5d0
               deallocate(vtemp)
            endif
         enddo
      else
         if(options%lnonorthogonal) then
            call collect_ortho_polaw(opi,opid)
            call free_memory(opi)
         endif
         call collect_v_pot(vpi,vpid)
         call free_memory(vpi)
         if(options%lnonorthogonal) then
            call collect_ortho_polaw(op,opd)
            call free_memory(op)
         endif
         do ii=options%i_min,options%i_max
            cpp%cprim=ii
            call mp_barrier( world_comm )
            if(.not.options%l_self_beta) then
               call read_data_pw_cprim_prod(cpp, options%prefix,.true.,ok_read,.false.,.false.)
            else
!read qmat
               call  create_vcprim(cppd, ii ,wu, qmd)
               if(ii>wu%nums_occ(1)) then
!if required adds the conduction term
                  call add_vcprim_conduction(cppd, wu, wup, vppd)
               end if
               call collect_cprim_prod(cpp,cppd)
               call free_memory(cppd)
            endif
         enddo
      endif
   enddo
   call mp_sum(sene(-ss%n:ss%n,:),world_comm)
   do ii=options%i_min,options%i_max
       do it=-ss%n,ss%n
          ss%diag(ii,it+ss%n+1,1)=ss%diag(ii, it+ss%n+1,1)+sene(it,ii)
       enddo
    enddo

 !copy sene results on ss with opportune factors
!!!!!!!!!!!
    if(options%l_self_beta) call free_memory(qmd)
    if(options%l_self_beta .and. options%i_max > wu%nums_occ(1) ) then
       call free_memory(vppd)
       call free_memory(wup)
    endif
   call free_memory(vpid)
   if(options%lnonorthogonal) then
      call free_memory(opd)
      call free_memory(opi)
      call free_memory(opid)
   endif
   call free_memory_polaw(ww)
   call free_memory( cpp)
   call free_memory(cppd)
   call free_memory( wu)
   deallocate(sene)
  return
 END SUBROUTINE selfenergy_ontime_file



  SUBROUTINE selfenergy_ontime_upper(ss, tf ,options)
!this subroutine calculates the self_energy of selected states


    USE io_global,     ONLY : stdout, ionode
    USE input_gw,      ONLY : input_options
    USE basic_structures,  ONLY : v_pot,wannier_u,free_memory, ortho_polaw,initialize_memory,cprim_prod,q_mat,&
                                        & wannier_u_prim,v_pot_prim, upper_states
    USE green_function,    ONLY : green, read_green, free_memory_green, initialize_green
    USE polarization,      ONLY : polaw, free_memory_polaw, read_polaw, invert_v_pot, invert_ortho_polaw,&
         & orthonormalize_inverse, initialize_polaw, orthonormalize_vpot, distribute_ortho_polaw, collect_ortho_polaw,&
         &  distribute_v_pot, collect_v_pot
    USE mp,                ONLY : mp_sum
    USE para_gww,          ONLY : is_my_pola
    USE mp_world,            ONLY : world_comm,nproc,mpime
    USE times_gw,  ONLY : times_freqs


   implicit none

   TYPE(times_freqs), INTENT(in) :: tf!for time frequency grids
   TYPE(input_options) :: options
   TYPE(self_storage)  :: ss

   TYPE(v_pot) :: vp,vpi,vpid
   TYPE(ortho_polaw) :: op,opi, opd, opid
   TYPE(polaw) :: ww!dressed interaction
   TYPE(wannier_u) :: wu!structure to be read and initialized
   TYPE(cprim_prod) :: cpp,cppd!the producs c' c' v  wp
   TYPE(upper_states) :: us


   INTEGER iw,jw,kw,it,ii,jj
   INTEGER :: l_blk, nbegin,nend
   REAL(kind=DP) :: offset
   COMPLEX(kind=DP), ALLOCATABLE :: sene(:,:)
   REAL(kind=DP), ALLOCATABLE:: wtemp(:,:), vtemp(:,:)
   REAL(kind=DP), EXTERNAL :: ddot
   LOGICAL :: ok_read

   nullify(vp%vmat)
   nullify(vpi%vmat)
   nullify(op%on_mat)
   nullify(opi%on_mat)
   nullify(opd%on_mat)
   nullify(opid%on_mat)
   nullify(vpid%vmat)
   nullify(wu%umat)



   call initialize_memory(cpp)
   call initialize_memory(cppd)


   call read_data_pw_u(wu, options%prefix)
   deallocate(wu%umat)
   call initialize_memory(us)
   call read_data_pw_upper_states(us,options%prefix)


   if(.not.options%l_hf_energies) then
      if(wu%nums > wu%nums_occ(1)) then
         offset=-(wu%ene(wu%nums_occ(1)+1,1)+wu%ene(wu%nums_occ(1),1))/2.d0
      else
         offset=-wu%ene(wu%nums_occ(1),1)
      endif
   else
      write(stdout,*) 'HF energies to be implemented YET'
      stop
      !if(wu%nums > wu%nums_occ(1)) then
      !   offset=-(ene_hf(wu%nums_occ(1)+1)+ene_hf(wu%nums_occ(1)))/2.d0
      !else
      !   offset=-ene_hf(wu%nums_occ(1))
      !endif
   endif



   call initialize_polaw(ww)
   write(stdout,*) 'addconduction_self_upper1'!ATTENZIONE
   FLUSH(stdout)

!read coulombian potential and calculate inverse

   if(ss%whole_s) then
      write(stdout,*) 'Whole s not implemented YET'
      stop
   endif


   if(options%w_divergence == 2) then
      call read_data_pw_v(vp,options%prefix,options%debug,0,.true.)
   else
      call read_data_pw_v(vp,options%prefix,options%debug,0,.false.)
   endif
   if(options%lnonorthogonal) then
      call read_data_pw_ortho_polaw(op,options%prefix)
      call orthonormalize_vpot(op,vp)
   endif
   call invert_v_pot(vp,vpi)
   call free_memory(vp)
   write(stdout,*) 'addconduction_self_upper1_45'
   FLUSH(stdout)
   call distribute_v_pot(vpi,vpid)
   call free_memory(vpi)
   if(options%lnonorthogonal) then
      call invert_ortho_polaw(op,opi)
   endif

   write(stdout,*) 'addconduction_self_upper1_5 op',op%numpw!ATTENZIONE
   FLUSH(stdout)
   if(options%lnonorthogonal) then
      call distribute_ortho_polaw(op,opd)
      call free_memory(op)

      write(stdout,*) 'addconduction_self_upper_6 opd',opd%numpw!ATTENZIONE
      FLUSH(stdout)
      call distribute_ortho_polaw(opi,opid)
      call free_memory(opi)
   endif

   l_blk= (ss%n+1)/nproc
   if(l_blk*nproc < (ss%n+1)) l_blk = l_blk+1
   nbegin=mpime*l_blk + 1 - (ss%n+1)
   nend=nbegin+l_blk-1
   if(nend > 0) nend = 0


   write(stdout,*) 'addconduction_self_upper5',nbegin,l_blk!ATTENZIONE
   FLUSH(stdout)
   allocate(sene(-ss%n:ss%n,options%i_min:options%i_max))
   sene(:,:)=(0.d0,0.d0)




!loop on negative imaginary times
   do it=nbegin,nbegin+l_blk-1
      if(it <= ss%n)  then
         write(stdout,*) 'addconduction_self_ontime time', it!ATTENZIONE
         FLUSH(stdout)
!we take care of the symmetru t ==> -t

         call read_polaw(abs(it),ww,options%debug,options%l_verbose)
         write(stdout,*) 'addconduction_self_upper6'!ATTENZIONE
         FLUSH(stdout)
         if(options%lnonorthogonal) then
            call collect_ortho_polaw(opi,opid)
            write(stdout,*) 'addconduction_self_ontime6.1'!ATTENZIONE
            call orthonormalize_inverse(opi,ww)
            write(stdout,*) 'addconduction_self_ontime6.2'!ATTENZIONE
            call free_memory(opi)
         endif
         write(stdout,*) 'addconduction_self_upper7'!ATTENZIONE
         FLUSH(stdout)
         allocate(wtemp(ww%numpw,ww%numpw))

         call collect_v_pot(vpi,vpid)

         call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,vpi%vmat,ww%numpw,ww%pw,ww%numpw,&
                 &0.d0, wtemp,ww%numpw)



         call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,wtemp,ww%numpw,vpi%vmat,ww%numpw,&
              &0.d0,ww%pw,ww%numpw)

         call free_memory(vpi)

         deallocate(wtemp)

         if(options%lnonorthogonal) then
            call collect_ortho_polaw(op,opd)
            call orthonormalize_inverse(op,ww)
            call free_memory(op)
         endif

!!now ww contains \tilde{ww}
         write(stdout,*) 'addconduction_self_ontime8'!ATTENZIONE
         FLUSH(stdout)

!read in cprim_prod
!first multiplication
!second multiplication
!copy on sene
            !loop on c' states
         do ii=options%i_min,options%i_max
            cpp%cprim=ii
            call read_data_pw_cprim_prod(cpp, options%prefix,.true.,ok_read,.false.,.true.)
            
            if(ok_read) then

!loop on c
               allocate(vtemp(cpp%numpw,us%nums_reduced))
               if(it <= 0) then
                  call dgemm('N','N',ww%numpw,us%nums_reduced,ww%numpw,1.d0,ww%pw,ww%numpw,&
                       & cpp%cpmat,cpp%lda,0.d0,vtemp,ww%numpw)

                  do jj=1,us%nums_reduced

!multiply W_ijS_jc =T_ic
!                     call dgemv('N',ww%numpw,ww%numpw,1.d0,ww%pw,ww%numpw,cpp%cpmat(:,jj),1,0.d0,vtemp,1)

!multiply S_icTi_c
                     sene(it,ii)=sene(it,ii)+ddot(cpp%numpw,vtemp(:,jj),1,cpp%cpmat(:,jj),1)*&
                          & exp((us%ene(jj)+offset)*tf%times(it))*ww%factor*(0.d0,-1.d0)

                  enddo
                  sene(it,ii)=sene(it,ii)*(0.d0,1.d0)
                  write(stdout,*) 'Conduction contribution', it,ii, sene(it,ii)
                  FLUSH(stdout)
               endif
               if(it >= 0) then
               endif
               if(it==0) sene(it,ii)=sene(it,ii)*0.5d0
               deallocate(vtemp)
            endif
         enddo
      else
         if(options%lnonorthogonal) then
            call collect_ortho_polaw(opi,opid)
            call free_memory(opi)
         endif
         call collect_v_pot(vpi,vpid)
         call free_memory(vpi)
         if(options%lnonorthogonal) then
            call collect_ortho_polaw(op,opd)
            call free_memory(op)
         endif
         do ii=options%i_min,options%i_max
            cpp%cprim=ii
            call read_data_pw_cprim_prod(cpp, options%prefix,.true.,ok_read,.false.,.true.)
         enddo
      endif
   enddo
   call mp_sum(sene(-ss%n:ss%n,:),world_comm)
   do ii=options%i_min,options%i_max
       do it=-ss%n,ss%n
          ss%diag(ii,it+ss%n+1,1)=ss%diag(ii, it+ss%n+1,1)+sene(it,ii)
       enddo
    enddo

 !copy sene results on ss with opportune factors
!!!!!!!!!!!
   
    
   call free_memory(vpid)
   if(options%lnonorthogonal) then
      call free_memory(opd)
      call free_memory(opi)
      call free_memory(opid)
   endif
   call free_memory_polaw(ww)
   call free_memory( cpp)
   call free_memory(cppd)
   call free_memory( wu)
   call free_memory(us)
   deallocate(sene)
  return
 END SUBROUTINE selfenergy_ontime_upper


 SUBROUTINE self_on_real_print(sr)
!this subroutine writes on charcter file the
!self energy on real axis

   USE io_global, ONLY : ionode
   USE io_files,             ONLY : prefix,tmp_dir
  
   implicit none
   INTEGER, EXTERNAL :: find_free_unit
   TYPE(self_on_real) :: sr

   INTEGER :: is ,ii,ie,iun
   CHARACTER(5) :: nfile

   if(ionode) then
   !loop on spin
      do is=1,sr%nspin
!loop on states
         do ii=sr%i_min,sr%i_max
            write(nfile,'(5i1)') &
                 & ii/10000,mod(ii,10000)/1000,mod(ii,1000)/100,mod(ii,100)/10,mod(ii,10)

         !openfile
            iun = find_free_unit()
            if(is==1) then
               open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'self_on_real'// nfile, status='unknown',form='formatted')
            else
               open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'self_on_real'// nfile, status='unknown',form='formatted')
            endif
            do ie=1,sr%n
               write(iun,*) dble(sr%grid(ie)),dble(sr%diag(ie,ii,is)),dimag(sr%diag(ie,ii,is))
            enddo
            close(iun)
         enddo
      enddo
   endif
   return

 END SUBROUTINE self_on_real_print

 SUBROUTINE self_on_real_value(sr,ii,ispin,energy,value,ierr)
!this subroutine gives the self_energy at point energy (complex) using
!linear extrapolation of real part only

   implicit none

   TYPE(self_on_real) :: sr
   INTEGER, INTENT(in) :: ii!KS  state
   INTEGER, INTENT(in) :: ispin!spin channel
   COMPLEX(kind=DP), INTENT(in) :: energy
   COMPLEX(kind=DP), INTENT(out) :: value
   INTEGER, INTENT(out) :: ierr

   INTEGER :: ie
   REAL(kind=DP) :: delta,lun

   ierr=0
   if(dble(sr%grid(1))>dble(energy) .or. dble(sr%grid(sr%n))<dble(energy)) then
      ierr=1
      return
   endif

   do ie=1,sr%n-1
      if(dble(sr%grid(ie)) < dble(energy) .and. dble(energy) <= dble(sr%grid(ie+1))) then
         lun=dble(energy)-dble(sr%grid(ie))
         delta=dble(sr%grid(ie+1))-dble(sr%grid(ie))
         value=sr%diag(ie,ii,ispin)+(-sr%diag(ie,ii,ispin)+sr%diag(ie+1,ii,ispin))*lun/delta
         exit
      endif
   enddo
 


   return

 END SUBROUTINE self_on_real_value


 END MODULE
