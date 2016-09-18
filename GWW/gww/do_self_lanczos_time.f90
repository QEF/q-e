!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!



subroutine do_reducible_pola(tf ,options)
!this subroutine calculates and writes on disk the reducible polarizability from the screen interaction

  USE kinds,             ONLY : DP
  USE io_global,         ONLY : stdout, ionode, ionode_id
  USE input_gw,          ONLY : input_options
  USE basic_structures,  ONLY : v_pot,wannier_u,free_memory, initialize_memory,lanczos_chain, vt_mat_lanczos,tt_mat_lanczos,&
                                     & semicore
  USE green_function,    ONLY : green, read_green, free_memory_green, initialize_green
  USE polarization,      ONLY : polaw, free_memory_polaw, read_polaw, write_polaw,invert_v_pot, initialize_polaw, &
                                  & read_polaw_global
  USE mp,                ONLY : mp_sum, mp_bcast
  USE mp_world,          ONLY : nproc,mpime,world_comm
  USE times_gw,          ONLY : times_freqs
  USE self_energy_storage, ONLY : self_storage,write_self_storage_ondisk,free_memory_self_storage
  USE lanczos
  USE constants,          ONLY : tpi,pi
  USE start_end ! debug                                                                                                                        
  USE parallel_include

  implicit none

  TYPE(times_freqs), INTENT(in) :: tf!for time frequency grids                                                                                 
  TYPE(input_options) :: options
  TYPE(self_storage)  :: ss

  TYPE(v_pot) :: vp,vpi
  TYPE(polaw) :: ww!dressed polarization   
  INTEGER :: l_blk, nbegin,nend, nsize
  REAL(kind=DP), ALLOCATABLE:: wtemp(:,:)
  INTEGER :: iw
  REAL(kind=DP) :: v_head

  nullify(vp%vmat)
  nullify(vpi%vmat)
  call initialize_polaw(ww)



  write(stdout,*) 'Trasform W to Pgreek'
  FLUSH(stdout)
  
  if(options%w_divergence == 2) then
     call read_data_pw_v(vp,options%prefix,options%debug,0,.true.)
  else
     call read_data_pw_v(vp,options%prefix,options%debug,0,.false.)
  endif
  v_head=vp%vmat(vp%numpw,vp%numpw)
  
  call invert_v_pot(vp,vpi)
  call free_memory(vp)
  
  l_blk= (tf%n+1)/nproc
  if(l_blk*nproc < (tf%n+1)) l_blk = l_blk+1
  nbegin=mpime*l_blk
  nend=nbegin+l_blk-1


!loop on imaginary frequency i\omega                                                                                                                                                                                                                       
  do iw=nbegin,nbegin+l_blk-1
     if(iw <= tf%n) then
        
        call read_polaw(iw,ww,options%debug,options%l_verbose)
        
        allocate(wtemp(ww%numpw,ww%numpw))
        call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,vpi%vmat,ww%numpw,ww%pw,ww%numpw,&
             &0.d0, wtemp,ww%numpw)
        call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,wtemp,ww%numpw,vpi%vmat,ww%numpw,&
             &0.d0,ww%pw,ww%numpw)
        deallocate(wtemp)
        
        call write_polaw(ww,options%debug)
     endif
  enddo


  call free_memory(vpi)
  call free_memory_polaw(ww)
  write(stdout,*) 'Done'
  FLUSH(stdout)



  return
end subroutine do_reducible_pola

subroutine do_self_lanczos_time(ss, tf ,options,l_real_axis,energy)
!this subroutine calculte the self-energy on time using fourier trasfrom using the lanczos scheme

  USE kinds,             ONLY : DP
  USE io_global,         ONLY : stdout, ionode, ionode_id
  USE input_gw,          ONLY : input_options
  USE basic_structures,  ONLY : v_pot,wannier_u,free_memory, initialize_memory,lanczos_chain, vt_mat_lanczos,tt_mat_lanczos,&
                                     & semicore
  USE green_function,    ONLY : green, read_green, free_memory_green, initialize_green
  USE polarization,      ONLY : polaw, free_memory_polaw, read_polaw, write_polaw,invert_v_pot, initialize_polaw, &
                                  & read_polaw_global
  USE mp,                ONLY : mp_sum, mp_bcast
  USE mp_world,          ONLY : nproc,mpime,world_comm
  USE times_gw,          ONLY : times_freqs
  USE self_energy_storage, ONLY : self_storage,write_self_storage_ondisk,free_memory_self_storage
  USE lanczos
  USE constants,          ONLY : tpi,pi
  USE start_end ! debug
  USE parallel_include
  USE io_files,  ONLY : prefix, tmp_dir

  implicit none

  TYPE(times_freqs), INTENT(in) :: tf!for time frequency grids
  TYPE(input_options) :: options
  TYPE(self_storage)  :: ss
  LOGICAL, INTENT(in) :: l_real_axis
!if true calculates on real frequency axis at given energy 
  REAL(kind=DP), INTENT(in) :: energy!energy on real axis at which calculating the self-energy (or part of it)
                                     !only if l_real_axis == true


  TYPE(v_pot) :: vp,vpi
  TYPE(polaw) :: ww!dressed polarization
  REAL(kind=DP) :: inv_epsi,v_head
  INTEGER :: l_blk, nbegin,nend, nsize, l_blk_freq, nbegin_freq,nend_freq
  REAL(kind=DP), ALLOCATABLE:: wtemp(:,:)
  INTEGER :: iw
  TYPE(wannier_u) :: uu
  REAL(kind=DP) :: offset
  TYPE(lanczos_chain) :: lc
  REAL(kind=DP), ALLOCATABLE :: re_e_mat(:,:,:),im_e_mat(:,:,:)
  COMPLEX(kind=DP), ALLOCATABLE :: e_mat_tmp(:,:,:)
  COMPLEX(kind=DP) :: af(1)
  REAL(kind=DP), ALLOCATABLE :: pw_mat(:,:,:),pw_mat_t(:,:,:)
  INTEGER :: numpw,numt,numl
  INTEGER :: it, ii,jj
  REAL(kind=DP) :: time,factor
  COMPLEX(kind=DP) :: cfactor
  REAL(kind=DP), ALLOCATABLE :: pw_tmp(:,:),pw_dumm(:,:)
  INTEGER :: iproc_time,ierr

  TYPE(vt_mat_lanczos) :: vtl,vtl_j
  TYPE(tt_mat_lanczos) :: ttl,ttl_j

  REAL(kind=DP), ALLOCATABLE :: re_g_mat(:,:,:), im_g_mat(:,:,:),re_g_mat_t(:,:,:), im_g_mat_t(:,:,:)
  REAL(kind=DP), ALLOCATABLE :: tmp_mat(:,:),tmp_mat1(:,:),tmp_mat2(:,:)
  REAL(kind=DP), ALLOCATABLE :: g_tmp(:,:), g_dumm(:,:), re_h_mat(:,:),im_h_mat(:,:)
  REAL(kind=DP), EXTERNAL :: DDOT

  LOGICAL :: l_single=.true.!if true e_mat is saved in single precision 

  REAL(kind=4), ALLOCATABLE :: re_e_mat_single(:,:,:),im_e_mat_single(:,:,:)
  REAL(kind=DP), ALLOCATABLE :: e_mat_double(:,:) 
  INTEGER :: l_blk_t, nbegin_t,nend_t, nsize_t,in
  INTEGER, PARAMETER :: ndivt=1!10
  INTEGER :: l_blk_g, nbegin_g,nend_g, nsize_g!paremter for optional dedicated frequency grid for G
  INTEGER :: j_min, j_max, is

  TYPE(semicore) :: sc
  REAL(kind=DP), ALLOCATABLE :: tmp_vec_sc(:)
  INTEGER :: iv_sc

  LOGICAL, PARAMETER :: l_distribute_sm=.true.!if true the S matrices are distributed among mpi tasks instead of being read from disk
                                              !it requires the parameter l_single==false
  INTEGER :: l_blk_sm, nbegin_sm,nend_sm, nsize_sm,iproc
  REAL(kind=DP), ALLOCATABLE :: st_save(:,:,:), sl_save(:,:,:)
  COMPLEX(kind=DP), ALLOCATABLE :: exp_table(:,:)
  REAL(kind=DP), ALLOCATABLE :: re_e_mat_t(:,:,:),im_e_mat_t(:,:,:)!in time for fourier trasform
  REAL(kind=DP), ALLOCATABLE :: re_e_mat_part(:,:,:),im_e_mat_part(:,:,:)!in time for storing partial calculations
  REAL(kind=DP), ALLOCATABLE :: pw_part_t(:,:,:)!in time for storing partial calculations 
  INTEGER :: n_cycles, i_cycles,i_min_cycles,i_max_cycles
  INTEGER :: n_list(2),iun,iun2
  INTEGER, ALLOCATABLE :: i_list(:,:)
  INTEGER, EXTERNAL :: find_free_unit

  if(options%whole_s) then
     l_single=.false.
  endif
  write(stdout,*) 'Routine do_self_lanczos_time'
  FLUSH(stdout)

  if(options%l_big_system) then
     n_cycles=options%i_max-options%i_min+1
  else
     n_cycles=1
  endif

  if(options%l_list) then
     if(ionode) then
        iun =  find_free_unit()
        open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'list_1.dat', status='old')
        read(iun,*) n_list(1)
        if(uu%nspin==2) then
           iun2 =  find_free_unit()
           open( unit=iun2, file=trim(tmp_dir)//trim(prefix)//'-'//'list_2.dat', status='old')
           read(iun,*) n_list(2)
        else
           n_list(2)=0
        endif
     endif
     call mp_bcast(n_list,ionode_id,world_comm)
     allocate(i_list(max(n_list(1),n_list(2)),2))
     i_list=0
     if(ionode) then
        do ii=1,n_list(1)
           read(iun,*) i_list(ii,1)
        enddo
        close(iun)
        if(uu%nspin==2) then
           do ii=1,n_list(2)
              read(iun2,*) i_list(ii,2)
           enddo
           close(iun2)
        endif
     endif
     call mp_bcast(i_list,ionode_id,world_comm)
     n_cycles=n_list(1)
  endif



!keeps in memory G, P(i\tau)
  

  nullify(vp%vmat)
  nullify(vpi%vmat)
  call initialize_polaw(ww)

  
  call initialize_memory(sc)
   

!calculate offset
!read in DFT energies

  call read_data_pw_u(uu,options%prefix)

  ss%ontime=.true.
  ss%max_i=options%max_i
  ss%i_min=options%i_min
  ss%i_max=options%i_max
  ss%n=tf%n
  ss%tau=options%tau
  ss%whole_s=options%whole_s
  ss%n_grid_fit=tf%n_grid_fit
  if(ss%whole_s) then
     ss%i_min_whole=options%i_min_whole
     ss%i_max_whole=options%i_max_whole
  endif

  ss%nspin=uu%nspin
  if(ss%whole_s) then
     allocate(ss%whole(ss%i_min_whole:ss%i_max_whole,ss%max_i,2*ss%n+1,ss%nspin))
     ss%whole(:,:,:,:)=(0.d0,0.d0)
     allocate(ss%whole_freq_fit(ss%i_min_whole:ss%i_max_whole,ss%max_i,2*ss%n_grid_fit+1,ss%nspin))
     ss%whole_freq_fit(:,:,:,:)=(0.d0,0.d0)
     allocate(ss%diag(ss%max_i,2*ss%n+1,ss%nspin))
     ss%diag(:,:,:)=(0.d0,0.d0)
     allocate(ss%diag_freq_fit(ss%max_i,2*ss%n_grid_fit+1,ss%nspin))
     ss%diag_freq_fit(:,:,:)=(0.d0,0.d0)
  else
     allocate(ss%diag(ss%max_i,2*ss%n+1,ss%nspin))
     ss%diag(:,:,:)=(0.d0,0.d0)
     nullify(ss%whole)
     allocate(ss%diag_freq_fit(ss%max_i,2*ss%n_grid_fit+1,ss%nspin))
     ss%diag_freq_fit(:,:,:)=(0.d0,0.d0)
     nullify(ss%whole_freq_fit)
  endif
!for compatibility                                                                                                                      
  allocate(ss%ene_remainder(ss%max_i,1))
  ss%ene_remainder(:,1)=0.d0

!loop on spin                                                                                                                             
  do is=1,ss%nspin


!NOT_TO_BE_INCLUDED_START
     if(options%l_semicore) call read_data_pw_semicore(sc, options%prefix, is)
!NOT_TO_BE_INCLUDED_END
     if(.not.l_real_axis) then
        if(uu%nums > uu%nums_occ(is)) then
           offset=-(uu%ene(uu%nums_occ(is)+1,is)+uu%ene(uu%nums_occ(is),is))/2.d0! CUSSI XE GIUSTO DEBUG 
           !offset=-(uu%ene(uu%nums_occ(2)+1,2)+uu%ene(uu%nums_occ(1),1))/2.d0
        else
           offset=-uu%ene(uu%nums_occ(is),is)
        endif
     else
        offset=-energy
     endif


!!!!!!!!
     l_blk= (tf%n+1)/nproc
     if(l_blk*nproc < (tf%n+1)) l_blk = l_blk+1
     nbegin=mpime*l_blk
     nend=nbegin+l_blk-1
     if(nend > tf%n) nend=tf%n
     nsize=nend-nbegin+1
     

 !read polarizability matrices     
  
     if( options%l_verbose) write(stdout,*) 'Read Pgreek'
     FLUSH(stdout)

     do iw=nbegin,nend
        call read_polaw(iw,ww,options%debug,options%l_verbose)
        if(iw==nbegin) allocate(pw_mat(ww%numpw,ww%numpw,l_blk))
        pw_mat(:,:,iw-nbegin+1)=ww%pw(:,:)
        call free_memory_polaw(ww)
     enddo
     numpw=ww%numpw
     call mp_bcast(numpw, ionode_id,world_comm)
     if(nbegin > tf%n) allocate(pw_mat(numpw,numpw,l_blk))

!Fourier trasform reducible polarizability matrices to imaginary time


     write(stdout,*) 'Fourier trasform Pgreek'
     FLUSH(stdout)


     allocate(pw_tmp(numpw,numpw))
     allocate(pw_dumm(numpw,numpw))
     allocate(pw_mat_t(numpw,numpw,l_blk))
!loop on time
     do it=0,tf%n
!each procs sums up its matrices in frequency with opportune factor
        time=tf%times(it)
        pw_tmp(:,:)=0.d0
        do iw=nbegin,nend
           factor=2.d0*dble(tf%weights_freq(iw)*exp((0.d0,1.d0)*tf%times(it)*tf%freqs_eff(iw)))/(2.d0*pi)
           pw_tmp(:,:)=pw_tmp(:,:)+pw_mat(:,:,iw-nbegin+1)*factor
        enddo
#if defined(__MPI)
!the distribution of times on procs is the same of that for frequecies
        iproc_time=it/l_blk
!all processors sums to iproc_time
        if(iproc_time==mpime) then
           call MPI_REDUCE(pw_tmp,pw_mat_t(1,1,it-nbegin+1),numpw*numpw,MPI_DOUBLE_PRECISION,MPI_SUM,iproc_time,world_comm,ierr)
        else
           call MPI_REDUCE(pw_tmp,pw_dumm,numpw*numpw,MPI_DOUBLE_PRECISION,MPI_SUM,iproc_time,world_comm,ierr)
        endif

#else
        pw_mat_t(:,:,it+1)=pw_tmp

#endif


!mp_sum to processer owing that time
!this processor put on opportune array
     enddo
     deallocate(pw_tmp)
     deallocate(pw_dumm)
     deallocate(pw_mat)

     l_blk_g= (tf%n_g+1)/nproc
     if(l_blk_g*nproc < (tf%n_g+1)) l_blk_g = l_blk_g+1
     nbegin_g=mpime*l_blk_g
     nend_g=nbegin_g+l_blk_g-1
     if(nend_g > tf%n_g) nend_g=tf%n_g
     nsize_g=nend_g-nbegin_g+1
!allocate and compute table for Fourier trasform
     allocate(exp_table(tf%n+1,l_blk_g))

     do it=0,tf%n
        do iw=nbegin_g,nend_g
           exp_table(it+1,iw-nbegin_g+1)=exp((0.d0,1.d0)*tf%times(it)*tf%freqs_g_eff(iw))
        enddo
     enddo




     do i_cycles=1,n_cycles
!calculates G

     call initialize_memory(lc)
     
     if(.not.options%l_list) then
        call read_data_pw_lanczos_chain(lc, i_cycles, options%prefix, .false.,is)
     else
        call read_data_pw_lanczos_chain(lc, i_list(i_cycles,is), options%prefix, .false.,is)
     endif
     write(stdout,*) 'Lanczos dimensions', lc%numt,lc%num_steps,is
     FLUSH(stdout)


  
 

     if(.not.l_single) then
        allocate(re_e_mat(lc%numt,lc%numt,l_blk_g))
        allocate(im_e_mat(lc%numt,lc%numt,l_blk_g))
    
        allocate(e_mat_tmp(lc%numt,lc%numt,1))
        do iw=nbegin_g,nbegin_g+l_blk_g-1
           if(iw <= tf%n_g) then
              af(1)=dcmplx(offset,-tf%freqs_g(iw))
              call solve_lanczos_complex(1,af,e_mat_tmp,lc)
              re_e_mat(:,:,iw-nbegin_g+1)=dble(e_mat_tmp(:,:,1))
              im_e_mat(:,:,iw-nbegin_g+1)=dimag(e_mat_tmp(:,:,1))
           else
              call solve_lanczos_fake_complex(lc)
           endif
     
        end do
        deallocate(e_mat_tmp)
!for entire self-energy store 
!!!!!!!!!!!!!!
        if( ss%whole_s .and.l_distribute_sm ) then
           if( options%l_verbose) write(stdout,*) 'before allocate',lc%numt,l_blk
           FLUSH(stdout)
           allocate(re_e_mat_t(lc%numt,lc%numt,l_blk), im_e_mat_t(lc%numt,lc%numt,l_blk))
           if( options%l_verbose) write(stdout,*) 'after '
           allocate(g_tmp(lc%numt,lc%numt))
           allocate(g_dumm(lc%numt,lc%numt))
           if( options%l_verbose) write(stdout,*) 'ATT1'
           FLUSH(stdout)

           numt=lc%numt
!loop on time                                                                       
           do it=0,tf%n
!each procs sums up its matrices in frequency with opportune factor
              time=tf%times(it)
              g_tmp(:,:)=0.d0
              do iw=nbegin_g,nend_g
                 factor=2.d0*dble(tf%weights_freq_g(iw)*exp_table(it+1,iw-nbegin_g+1))/(2.d0*pi)
                 g_tmp(1:numt,1:numt)=g_tmp(1:numt,1:numt)+re_e_mat(1:numt,1:numt,iw-nbegin_g+1)*factor
              enddo
#if defined(__MPI)
!the distribution of times on procs is the same of that for frequecies
              iproc_time=it/l_blk
!all processors sums to iproc_time
              if(iproc_time==mpime) then
                 call MPI_REDUCE(g_tmp,re_e_mat_t(1,1,it-nbegin+1),numt*numt,&
                      &MPI_DOUBLE_PRECISION,MPI_SUM,iproc_time,world_comm,ierr)
              else
                 call MPI_REDUCE(g_tmp,g_dumm,numt*numt,MPI_DOUBLE_PRECISION,&
                      &MPI_SUM,iproc_time,world_comm,ierr)
              endif

#else
              re_e_mat_t(1:numt,1:numt,it+1)=g_tmp(1:numt,1:numt)

#endif

           enddo
!loop on time
           do it=0,tf%n
!each procs sums up its matrices in frequency with opportune factor
              time=tf%times(it)
              g_tmp(:,:)=0.d0
              do iw=nbegin_g,nend_g
                 factor=-2.d0*dimag(tf%weights_freq_g(iw)*exp_table(it+1,iw-nbegin_g+1))/(2.d0*pi)
                 g_tmp(1:numt,1:numt)=g_tmp(1:numt,1:numt)+im_e_mat(1:numt,1:numt,iw-nbegin_g+1)*factor
              enddo
#if defined(__MPI)
!the distribution of times on procs is the same of that for frequecies
              iproc_time=it/l_blk
!all processors sums to iproc_time
              if(iproc_time==mpime) then
                 call MPI_REDUCE(g_tmp,im_e_mat_t(1,1,it-nbegin+1),numt*numt,&
                      &MPI_DOUBLE_PRECISION,MPI_SUM,iproc_time,world_comm,ierr)
              else
                 call MPI_REDUCE(g_tmp,g_dumm,numt*numt,MPI_DOUBLE_PRECISION,&
                      &MPI_SUM,iproc_time,world_comm,ierr)
              endif

#else
              im_e_mat_t(1:numt,1:numt,it+1)=g_tmp(1:numt,1:numt)

#endif

           enddo

           deallocate(g_tmp,g_dumm)
        endif
        if( options%l_verbose) write(stdout,*) 'ATT2' !DEBUG
        FLUSH(stdout)
  
     else
       if(i_cycles==1) then
          allocate(re_e_mat_single(lc%numt,lc%numt,l_blk_g))
          allocate(im_e_mat_single(lc%numt,lc%numt,l_blk_g))
       endif
        do iw=nbegin_g,nbegin_g+l_blk_g-1
           if(iw <= tf%n_g) then
              af(1)=dcmplx(offset,-tf%freqs_g(iw))
              call solve_lanczos_single(af(1),re_e_mat_single(1,1,iw-nbegin_g+1),im_e_mat_single(1,1,iw-nbegin_g+1),lc)
           else
              call solve_lanczos_fake_single(lc)
           endif

        end do
        l_blk_t= (lc%numt)/ndivt
        if(l_blk_t*ndivt < (lc%numt)) l_blk_t = l_blk_t+1


        if(i_cycles==1) allocate(e_mat_double(lc%numt,l_blk_t))

     endif
     if( options%l_verbose) write(stdout,*) 'Done'
     FLUSH(stdout)

     call initialize_memory(vtl)
     call initialize_memory(ttl)
     call initialize_memory(vtl_j)
     call initialize_memory(ttl_j)

!if required read all S matrices and distribute among mpi tasks

     if( ss%whole_s .and.l_distribute_sm ) then
        l_blk_sm= (ss%max_i)/nproc
        if(l_blk_sm*nproc < (ss%max_i)) l_blk_sm = l_blk_sm+1
        nbegin_sm=mpime*l_blk_sm+1
        nend_sm=nbegin_sm+l_blk_sm-1
        if(nend_sm > ss%max_i) nend_sm=ss%max_i
        nsize_sm=nend_sm-nbegin_sm+1
        do jj=1,ss%max_i
           call  read_data_pw_vt_mat_lanczos(vtl_j, jj, options%prefix, .false.,is)
           call  read_data_pw_tt_mat_lanczos(ttl_j, jj, options%prefix, .false.,is)
           if(jj==1) then
              numt=ttl_j%numt
              numl=ttl_j%numl
              allocate(sl_save(numpw,numl,l_blk_sm))
              allocate(st_save(numt,numl,l_blk_sm))
           endif
           if(jj>=nbegin_sm .and. jj<=nend_sm) then
              st_save(1:numt,1:numl,jj-nbegin_sm+1)=ttl_j%tt_mat(1:numt,1:numl)
              sl_save(1:numpw,1:numl,jj-nbegin_sm+1)=vtl_j%vt_mat(1:numpw,1:numl)
           endif
           call free_memory(vtl_j)
           call free_memory(ttl_j)
        enddo
        allocate(re_e_mat_part(numt,numl,l_blk),im_e_mat_part(numt,numl,l_blk))
        allocate(pw_part_t(numpw,numl,l_blk))
     endif





 !loop on KS states
     if(options%l_big_system) then
        if(.not.options%l_list) then
           i_min_cycles=ss%i_min+i_cycles-1
           i_max_cycles=i_min_cycles
        else
           i_min_cycles=i_list(i_cycles,is)
           i_max_cycles=i_min_cycles
        endif
     else
        i_min_cycles=ss%i_min
        i_max_cycles=ss%i_max
     endif

     do ii=i_min_cycles,i_max_cycles
        write(stdout,*) 'Loop on KS:',ii, is
        FLUSH(stdout)
        


!calculates G on partial basis
!read vtl ttl
        call  read_data_pw_vt_mat_lanczos(vtl, ii, options%prefix, .false.,is)
        call  read_data_pw_tt_mat_lanczos(ttl, ii, options%prefix, .false.,is)
        if(ii==i_min_cycles.and.i_cycles==1) then
           numt=ttl%numt
           numl=ttl%numl
           allocate(re_g_mat(numl,numl,l_blk_g),im_g_mat(numl,numl,l_blk_g))
           allocate(re_g_mat_t(numl,numl,l_blk),im_g_mat_t(numl,numl,l_blk))
        endif
!multiply and  put on array
!loop on frequency
        if(.not.ss%whole_s) then
           j_min=ii
           j_max=ii
        else
           j_min=ss%i_min_whole
           j_max=ss%i_max_whole
        endif
        if( options%l_verbose) write(stdout,*) 'Doing dgemms',numl,numt,numl,numpw

        if(ss%whole_s .and. l_distribute_sm) then!if required calculates partial products
           do it=nbegin,nend
              call dgemm('N','N',numt,numl,numt,1.d0,re_e_mat_t(1,1,it-nbegin+1),numt,&
                   &ttl%tt_mat,numt,0.d0,re_e_mat_part(1,1,it-nbegin+1),numt)
              call dgemm('N','N',numt,numl,numt,1.d0,im_e_mat_t(1,1,it-nbegin+1),numt,&
                   &ttl%tt_mat,numt,0.d0,im_e_mat_part(1,1,it-nbegin+1),numt)
              call dgemm('N','N',numpw,numl,numpw,1.d0,pw_mat_t(1,1,it-nbegin+1) ,numpw, &
                   & vtl%vt_mat,numpw,0.d0,pw_part_t(1,1,it-nbegin+1),numpw)
           enddo


        endif
        FLUSH(stdout)
        do jj=j_min,j_max
           if(ss%whole_s .and. l_distribute_sm) then
              allocate(ttl_j%tt_mat(numt,numl))
              allocate(vtl_j%vt_mat(numpw,numl))
              if(jj>=nbegin_sm .and. jj<=nend_sm) then
                 vtl_j%vt_mat(1:numpw,1:numl)=sl_save(1:numpw,1:numl,jj-nbegin_sm+1)
                 ttl_j%tt_mat(1:numt,1:numl)=st_save(1:numt,1:numl,jj-nbegin_sm+1)
                 
              endif
             iproc=(jj-1)/l_blk_sm
             call mp_bcast(vtl_j%vt_mat, iproc,world_comm)
             call mp_bcast(ttl_j%tt_mat, iproc,world_comm)
           else
              call  read_data_pw_vt_mat_lanczos(vtl_j, jj, options%prefix, .false.,is)
              call  read_data_pw_tt_mat_lanczos(ttl_j, jj, options%prefix, .false.,is)
           endif
           if(.not.(ss%whole_s .and. l_distribute_sm)) then
              allocate(tmp_mat(numl,numt))
              do iw=nbegin_g,nend_g
                 if( options%l_verbose) write(stdout,*) 'Doing dgemms',numl,numt,numpw,l_blk,iw
                 FLUSH(stdout)
              
                 if(.not.l_single) then
                    call dgemm('T','N',numl,numt,numt,1.d0,ttl_j%tt_mat,numt,re_e_mat(1,1,iw-nbegin_g+1),numt,0.d0,tmp_mat,numl)
                    call dgemm('N','N',numl,numl,numt,1.d0,tmp_mat,numl,ttl%tt_mat,numt,0.d0,re_g_mat(1,1,iw-nbegin_g+1),numl)
                    call dgemm('T','N',numl,numt,numt,1.d0,ttl_j%tt_mat,numt,im_e_mat(1,1,iw-nbegin_g+1),numt,0.d0,tmp_mat,numl)
                    call dgemm('N','N',numl,numl,numt,1.d0,tmp_mat,numl,ttl%tt_mat,numt,0.d0,im_g_mat(1,1,iw-nbegin_g+1),numl)
                 else
                    do in=0,ndivt-1
                       nbegin_t=in*l_blk_t+1
                       nend_t=min(nbegin_t+l_blk_t-1,lc%numt)
                       nsize_t=nend_t-nbegin_t+1
                       if(nsize_t >= 1) then
                          e_mat_double(1:lc%numt, 1:nsize_t)=dble(re_e_mat_single(1:lc%numt, nbegin_t:nend_t, iw-nbegin_g+1))
                       
                          call dgemm('T','N',numl,nsize_t,numt,1.d0,ttl_j%tt_mat,numt,e_mat_double,numt,0.d0,&
                           &tmp_mat(1,nbegin_t),numl)
                       endif
                    enddo
                    
                    if( options%l_verbose) write(stdout,*) 'ATT1'!DEBUG
                    FLUSH(stdout)
                    call dgemm('N','N',numl,numl,numt,1.d0,tmp_mat,numl,ttl%tt_mat,numt,0.d0,&
                         &re_g_mat(1,1,iw-nbegin_g+1),numl)
                    if( options%l_verbose) write(stdout,*) 'ATT2'!DEBUG
                    FLUSH(stdout)

                    do in=0,ndivt-1
                       nbegin_t=in*l_blk_t+1
                       nend_t=min(nbegin_t+l_blk_t-1,lc%numt)
                       nsize_t=nend_t-nbegin_t+1
                       e_mat_double(1:lc%numt, 1:nsize_t)=dble(im_e_mat_single(1:lc%numt, nbegin_t:nend_t, iw-nbegin_g+1))
                       call dgemm('T','N',numl,nsize_t,numt,1.d0,ttl_j%tt_mat,numt,e_mat_double,numt,0.d0,tmp_mat(1,nbegin_t),numl)
                    enddo
                    if( options%l_verbose) write(stdout,*) 'ATT3'!DEBUG
                    FLUSH(stdout)
                    
                    call dgemm('N','N',numl,numl,numt,1.d0,tmp_mat,numl,ttl%tt_mat,numt,0.d0,im_g_mat(1,1,iw-nbegin_g+1),numl)
                 endif
     
    

              enddo
       
              deallocate(tmp_mat)
 
              write(stdout,*) 'Fourier trasform:'
              FLUSH(stdout)



!Fourier trasform
              allocate(g_tmp(numl,numl))
              allocate(g_dumm(numl,numl))
              if( options%l_verbose) write(stdout,*) 'ATT1'
              FLUSH(stdout)
           
!loop on time
              do it=0,tf%n
!each procs sums up its matrices in frequency with opportune factor
                 time=tf%times(it)
                 g_tmp(:,:)=0.d0
                 do iw=nbegin_g,nend_g
                    !factor=2.d0*dble(tf%weights_freq_g(iw)*exp((0.d0,1.d0)*tf%times(it)*tf%freqs_g_eff(iw)))/(2.d0*pi)
                    factor=2.d0*dble(tf%weights_freq_g(iw)*exp_table(it+1,iw-nbegin_g+1))/(2.d0*pi)
                    g_tmp(:,:)=g_tmp(:,:)+re_g_mat(:,:,iw-nbegin_g+1)*factor
                 enddo
#if defined(__MPI)
!the distribution of times on procs is the same of that for frequecies
                 iproc_time=it/l_blk
!all processors sums to iproc_time
                 if(iproc_time==mpime) then
                    call MPI_REDUCE(g_tmp,re_g_mat_t(1,1,it-nbegin+1),numl*numl,&
                         &MPI_DOUBLE_PRECISION,MPI_SUM,iproc_time,world_comm,ierr)
                 else
                    call MPI_REDUCE(g_tmp,g_dumm,numl*numl,MPI_DOUBLE_PRECISION,&
                         &MPI_SUM,iproc_time,world_comm,ierr)
                 endif

#else
                 re_g_mat_t(:,:,it+1)=g_tmp(:,:)

#endif
                 g_tmp(:,:)=0.d0
                 do iw=nbegin_g,nend_g
                    !factor=-2.d0*dimag(tf%weights_freq_g(iw)*exp((0.d0,1.d0)*tf%times(it)*tf%freqs_g_eff(iw)))/(2.d0*pi)
                    factor=-2.d0*dimag(tf%weights_freq_g(iw)*exp_table(it+1,iw-nbegin_g+1))/(2.d0*pi)
                    g_tmp(:,:)=g_tmp(:,:)+im_g_mat(:,:,iw-nbegin_g+1)*factor
                 enddo
#if defined(__MPI)
!the distribution of times on procs is the same of that for frequecies
                 iproc_time=it/l_blk
!all processors sums to iproc_time
                 if(iproc_time==mpime) then
                    call MPI_REDUCE(g_tmp,im_g_mat_t(1,1,it-nbegin+1),numl*numl,&
                         &MPI_DOUBLE_PRECISION,MPI_SUM,iproc_time,world_comm,ierr)
                 else
                    call MPI_REDUCE(g_tmp,g_dumm,numl*numl,MPI_DOUBLE_PRECISION,MPI_SUM,iproc_time,world_comm,ierr)
                 endif
                 
#else
                 im_g_mat_t(:,:,it+1)=g_tmp(:,:)

#endif


              enddo
              deallocate(g_tmp,g_dumm)
           else!for whole self-energy calculates directly in time domain
              do it=nbegin,nend
                 call dgemm('T','N',numl,numl,numt,1.d0,ttl_j%tt_mat,numt,re_e_mat_part(1,1,it-nbegin+1),&
                      &numt,0.d0,re_g_mat_t(1,1,it-nbegin+1),numl)
                 call dgemm('T','N',numl,numl,numt,1.d0,ttl_j%tt_mat,numt,im_e_mat_part(1,1,it-nbegin+1),&
                      &numt,0.d0,im_g_mat_t(1,1,it-nbegin+1),numl)
              enddo
           endif
           if( options%l_verbose) write(stdout,*) 'done'
           FLUSH(stdout)


   

!loop on frequency

     
           write(stdout,*) 'Products in imaginary time:'
           FLUSH(stdout)

           
           if( .not.(ss%whole_s .and.l_distribute_sm )) then
              allocate(tmp_mat(numpw,numl),re_h_mat(numpw,numpw),im_h_mat(numpw,numpw))
              do it=nbegin,nend

              
!matrix multiply
                 call dgemm('N','N',numpw,numl,numl,1.d0,vtl_j%vt_mat,numpw,re_g_mat_t(1,1,it-nbegin+1),numl,0.d0,tmp_mat,numpw)
                 call dgemm('N','T',numpw,numpw,numl,1.d0,tmp_mat,numpw, vtl%vt_mat,numpw,0.d0,re_h_mat,numpw)
           
                 call dgemm('N','N',numpw,numl,numl,1.d0,vtl_j%vt_mat,numpw,im_g_mat_t(1,1,it-nbegin+1),numl,0.d0,tmp_mat,numpw)
                 call dgemm('N','T',numpw,numpw,numl,1.d0,tmp_mat,numpw, vtl%vt_mat,numpw,0.d0,im_h_mat,numpw)

!product

                 if(ii==jj) then
                    ss%diag(ii,it+ss%n+1,is)=DDOT(numpw*numpw,re_h_mat,1,pw_mat_t(1,1,it-nbegin+1),1)+&
                         &DDOT(numpw*numpw,im_h_mat,1,pw_mat_t(1,1,it-nbegin+1),1)
                    ss%diag(ii,ss%n+1-it,is)=DDOT(numpw*numpw,re_h_mat,1,pw_mat_t(1,1,it-nbegin+1),1)-&
                         &DDOT(numpw*numpw,im_h_mat,1,pw_mat_t(1,1,it-nbegin+1),1)
        
!if required add semicore terms
                    if(options%l_semicore) then
!NOT_TO_BE_INCLUDED_START
                       allocate(tmp_vec_sc(numpw))
                       do iv_sc=1,sc%n_semicore
                          if(it==0) write(stdout,*) 'SEMICORE PRIMA',iv_sc,ii, ss%diag(ii,it+ss%n+1,is)
                          call dgemv('N',numpw,numpw,1.d0,pw_mat_t(1,1,it-nbegin+1),numpw,&
                               &sc%ppw_mat(1,iv_sc,ii),1,0.d0,tmp_vec_sc,1)
                          ss%diag(ii,it+ss%n+1,is)= ss%diag(ii,it+ss%n+1,is)- &!ATTENZIONE ERA +
                               &DDOT(numpw,tmp_vec_sc,1,sc%ppw_mat(1,iv_sc,ii),1)*exp((offset+sc%en_sc(iv_sc))*tf%times(it))
                          if(it==0) write(stdout,*) 'SEMICORE',iv_sc,ii, ss%diag(ii,it+ss%n+1,is)
                          if((offset+sc%en_sc(iv_sc))*tf%times(it)>0) write(stdout,*) 'OCIO!!!'!DEBUG ATTENZIONE
                       enddo
                       deallocate(tmp_vec_sc)
!NOT_TO_BE_INCLUDED_END
                    endif
                 endif
                 if(ss%whole_s) then
                    ss%whole(jj,ii,it+ss%n+1,is)=DDOT(numpw*numpw,re_h_mat,1,pw_mat_t(1,1,it-nbegin+1),1)+&
                         &DDOT(numpw*numpw,im_h_mat,1,pw_mat_t(1,1,it-nbegin+1),1)
                    ss%whole(jj,ii,ss%n+1-it,is)=DDOT(numpw*numpw,re_h_mat,1,pw_mat_t(1,1,it-nbegin+1),1)-&
                         &DDOT(numpw*numpw,im_h_mat,1,pw_mat_t(1,1,it-nbegin+1),1)

                 endif
              enddo
              deallocate(re_h_mat,im_h_mat)
              deallocate(tmp_mat)
           else
              allocate(tmp_mat(numl,numl))
              do it=nbegin,nend
                 call dgemm('T','N',numl,numl,numpw,1.d0,vtl_j%vt_mat,numpw,pw_part_t(1,1,it-nbegin+1),numpw,0.d0,tmp_mat,numl)
                 if(ii==jj) then
                    ss%diag(ii,it+ss%n+1,is)=DDOT(numl*numl,re_g_mat_t(1,1,it-nbegin+1),1,tmp_mat,1)+&
                         &DDOT(numl*numl,im_g_mat_t(1,1,it-nbegin+1),1,tmp_mat,1)
                    ss%diag(ii,ss%n+1-it,is)=DDOT(numl*numl,re_g_mat_t(1,1,it-nbegin+1),1,tmp_mat,1)-&
                         &DDOT(numl*numl,im_g_mat_t(1,1,it-nbegin+1),1,tmp_mat,1)
                 endif
                 
                 ss%whole(jj,ii,it+ss%n+1,is)=DDOT(numl*numl,re_g_mat_t(1,1,it-nbegin+1),1,tmp_mat,1)+&
                      &DDOT(numl*numl,im_g_mat_t(1,1,it-nbegin+1),1,tmp_mat,1)
                 ss%whole(jj,ii,ss%n+1-it,is)=DDOT(numl*numl,re_g_mat_t(1,1,it-nbegin+1),1,tmp_mat,1)-&
                      &DDOT(numl*numl,im_g_mat_t(1,1,it-nbegin+1),1,tmp_mat,1)
           
                 
              enddo
              deallocate(tmp_mat)
           endif
           if( options%l_verbose) write(stdout,*) 'done'
           FLUSH(stdout)



!mp_sum for distributing on all processors
           if(ii==jj) then
              call mp_sum(ss%diag(ii,1:2*ss%n+1,is),world_comm)
           endif
           if(ss%whole_s) then
              call mp_sum(ss%whole(jj,ii,1:2*ss%n+1,is),world_comm)
           endif

     
!FT of diagonal part of self-nergy 
    
           
      
           !deallocate(re_g_mat_t, im_g_mat_t)
       
           call free_memory(vtl_j)
           call free_memory(ttl_j)

           
        enddo !jj
    
        call free_memory(vtl)
        call free_memory(ttl)

     enddo!on i_min, i_max
 
  enddo!on i_cycles


!put factor due to FFT on imaginary axis
     ss%diag(:,:,is)=ss%diag(:,:,is)*(0.d0,1.d0)
   
     if(ss%whole_s) then
        ss%whole(:,:,:,is)=ss%whole(:,:,:,is)*(0.d0,1.d0)
     endif
     deallocate(re_g_mat,im_g_mat)
     deallocate(re_g_mat_t,im_g_mat_t)

     deallocate(pw_mat_t)
     if(.not.l_single) then
        deallocate(re_e_mat,im_e_mat)
     else
        deallocate(re_e_mat_single, im_e_mat_single, e_mat_double)
     endif
     if( ss%whole_s .and.l_distribute_sm ) then
        deallocate(st_save,sl_save)
        if(.not.l_single)  deallocate(re_e_mat_t,im_e_mat_t)
        if(.not.l_single)  deallocate(re_e_mat_part,im_e_mat_part,pw_part_t)
     endif
     deallocate(exp_table)
     call free_memory(lc)
  enddo!on spin
  call free_memory(uu)
  call free_memory(sc)
  return
end subroutine do_self_lanczos_time



subroutine solve_lanczos_single(alpha,re_e_mat,im_e_mat,lc)
!this subroutine sums to  the matrix E_{no}=<t_n|(H-alpha)^-1|t_o>

  USE kinds,            ONLY : DP
  USE basic_structures, ONLY : lanczos_chain, initialize_memory,free_memory
  USE io_global,        ONLY : stdout
  USE mp,               ONLY : mp_sum,mp_bcast
  USE mp_world,         ONLY : nproc, mpime, world_comm

  implicit none


  COMPLEX(kind=DP) :: alpha!constant for Ev+iw
  TYPE(lanczos_chain) :: lc!lanczos chain descriptor
  REAL(kind=4) :: re_e_mat(lc%numt,lc%numt) !real part of matrix to be calculated
  REAL(kind=4) :: im_e_mat(lc%numt,lc%numt) !imaginary part of matrix to be calculated
 


  INTEGER :: io,info,ii,jj
  COMPLEX(kind=DP), ALLOCATABLE :: dl(:),du(:),d(:),t(:)
  COMPLEX(kind=DP), ALLOCATABLE :: omat(:,:)
  REAL(kind=DP), ALLOCATABLE :: tmp_mat(:,:)
  INTEGER :: l_blk,nbegin,nend, iproc
  COMPLEX(kind=DP), ALLOCATABLE :: e_mat(:)



  l_blk= (lc%numt)/nproc
  if(l_blk*nproc < (lc%numt)) l_blk = l_blk+1
  nbegin=mpime*l_blk+1
  nend=nbegin+l_blk-1

  allocate(dl(lc%num_steps-1),du(lc%num_steps-1),d(lc%num_steps),t(lc%num_steps))
  re_e_mat(:,:)=0.0
  im_e_mat(:,:)=0.0
  allocate(omat(lc%numt,lc%num_steps))
  allocate(tmp_mat(lc%numt,lc%num_steps))
  allocate(e_mat(lc%numt))
!loop on o
  do io=1,lc%numt
!!set up vectors for lapack routine
     !recover matrix from processor                                                                                                                               
     tmp_mat(:,:)=0.d0
     if(io >= nbegin .and. io <= nend) then
        tmp_mat(1:lc%numt,1:lc%num_steps)=lc%o_mat(1:lc%numt,1:lc%num_steps,io-nbegin+1)
     endif
     iproc=(io-1)/l_blk
     call mp_bcast(tmp_mat(:,:), iproc, world_comm)
     omat(:,:)=dcmplx(tmp_mat(:,:),0.d0)

   

     dl(1:lc%num_steps-1)=cmplx(lc%f(1:lc%num_steps-1,io),0.d0)
     du(1:lc%num_steps-1)=cmplx(lc%f(1:lc%num_steps-1,io),0.d0)
     d(1:lc%num_steps)=cmplx(lc%d(1:lc%num_steps,io),0.d0)+alpha
     t(:)=(0.d0,0.d0)
     t(1)=(1.d0,0.d0)

!!call lapack
     call zgtsv(lc%num_steps,1,dl,d,du,t,lc%num_steps,info)
     if(info /= 0) then
        write(stdout,*) 'ZGTSV info:', info
        FLUSH(stdout)
        stop
     endif



!!calculate term
     
     call zgemm('N','N',lc%numt,1,lc%num_steps,(1.d0,0.d0),omat,lc%numt,t,lc%num_steps,(0.d0,0.d0),e_mat,lc%numt)
     re_e_mat(1:lc%numt,io)=re_e_mat(1:lc%numt,io)+real(e_mat(1:lc%numt))
     im_e_mat(1:lc%numt,io)=im_e_mat(1:lc%numt,io)+imag(e_mat(1:lc%numt))
  enddo
 
 
  do ii=1,lc%numt
     do jj=ii+1,lc%numt
        re_e_mat(ii,jj)=0.5*(re_e_mat(ii,jj)+re_e_mat(jj,ii))
        re_e_mat(jj,ii)=re_e_mat(ii,jj)
        im_e_mat(ii,jj)=0.5*(im_e_mat(ii,jj)+im_e_mat(jj,ii))
        im_e_mat(jj,ii)=im_e_mat(ii,jj)

     enddo
  enddo
 
  deallocate(dl,du,d,t)
  deallocate(omat,tmp_mat)
  deallocate(e_mat)
  return
end  subroutine solve_lanczos_single


subroutine solve_lanczos_fake_single(lc)

!this subroutine is a parallel fake routine for the solve lanczos routine

  USE kinds,            ONLY : DP
  USE basic_structures, ONLY : lanczos_chain, initialize_memory,free_memory
  USE io_global,        ONLY : stdout
  USE mp,               ONLY : mp_sum,mp_bcast
  USE mp_world,         ONLY : nproc,mpime,world_comm

  implicit none

  TYPE(lanczos_chain) :: lc!lanczos chain descriptor

  INTEGER :: l_blk,nbegin,nend, iproc
  REAL(kind=DP), ALLOCATABLE :: o_mat(:,:)
  INTEGER :: io


  allocate(o_mat(lc%numt,lc%num_steps))
  l_blk= (lc%numt)/nproc
  if(l_blk*nproc < (lc%numt)) l_blk = l_blk+1
  nbegin=mpime*l_blk+1
  nend=nbegin+l_blk-1


!loop on io
  do io=1,lc%numt
!recover matrix from processor
     o_mat(:,:)=0.d0
     if(io >= nbegin .and. io <= nend) then
        o_mat(:,:)=lc%o_mat(:,:,io-nbegin+1)
     endif
     iproc=(io-1)/l_blk
     call mp_bcast(o_mat(:,:), iproc, world_comm)


  enddo

  deallocate(o_mat)
end subroutine solve_lanczos_fake_single


