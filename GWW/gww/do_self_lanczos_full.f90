!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!




subroutine do_self_lanczos_full(ss, tf ,options,l_real_axis,energy)
!this subroutine calculates the self-energy on time using fourier transform with the lanczos scheme

  USE kinds,             ONLY : DP
  USE io_global,         ONLY : stdout, ionode, ionode_id
  USE input_gw,          ONLY : input_options
  USE basic_structures,  ONLY : v_pot,wannier_u,free_memory, initialize_memory,lanczos_chain, vt_mat_lanczos,tt_mat_lanczos,&
                                     & semicore,mat_lanczos_full,full_prods
  USE green_function,    ONLY : green, read_green, free_memory_green, initialize_green
  USE polarization,      ONLY : polaw, free_memory_polaw, read_polaw, write_polaw,invert_v_pot, initialize_polaw, &
                                  & read_polaw_global
  USE mp,                ONLY : mp_sum, mp_bcast
  USE mp_world,         ONLY : nproc,mpime,world_comm
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
  COMPLEX(kind=DP), ALLOCATABLE :: e_mat(:,:,:)
  COMPLEX(kind=DP), ALLOCATABLE :: e_mat_tmp(:,:,:)
  COMPLEX(kind=DP) :: af(1)
  REAL(kind=DP), ALLOCATABLE :: pw_mat(:,:,:),pw_mat_t(:,:,:)
  INTEGER :: numpw,numt,numl,nums
  INTEGER :: it, ii,jj
  REAL(kind=DP) :: time,factor
  COMPLEX(kind=DP) :: cfactor
  REAL(kind=DP), ALLOCATABLE :: pw_tmp(:,:),pw_dumm(:,:)
  INTEGER :: iproc_time,ierr

  TYPE(mat_lanczos_full) :: fli,flj
  TYPE(full_prods) :: fp

  COMPLEX(kind=DP), ALLOCATABLE :: g_mat(:,:,:,:), g_mat_t(:,:,:,:),tmp_mat(:,:)
  REAL(kind=DP), ALLOCATABLE :: tmp_mat1(:,:),tmp_mat2(:,:)
  REAL(kind=DP), ALLOCATABLE :: g_tmp(:,:), g_dumm(:,:), re_h_mat(:,:),im_h_mat(:,:)
  REAL(kind=DP), EXTERNAL :: DDOT
  COMPLEX(kind=DP), EXTERNAL :: zdotc
  
  LOGICAL :: l_single=.false.!if true e_mat is saved in single precision 

  REAL(kind=4), ALLOCATABLE :: re_e_mat_single(:,:,:),im_e_mat_single(:,:,:)
  REAL(kind=DP), ALLOCATABLE :: e_mat_double(:,:) 
  INTEGER :: l_blk_t, nbegin_t,nend_t, nsize_t,in
  INTEGER, PARAMETER :: ndivt=1!10
  INTEGER :: l_blk_g, nbegin_g,nend_g, nsize_g!paremter for optional dedicated frequency grid for G
  INTEGER :: j_min, j_max, is

  TYPE(semicore) :: sc
  REAL(kind=DP), ALLOCATABLE :: tmp_vec_sc(:)
  COMPLEX(kind=DP), ALLOCATABLE :: tmp_vec(:)
  INTEGER :: iv_sc

  LOGICAL, PARAMETER :: l_distribute_sm=.false.!if true the S matrices are distributed among mpi tasks instead of being read from disk
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
  INTEGER :: ipol,iv,ic

  if(options%whole_s) then
     l_single=.false.
  endif
  write(stdout,*) 'Routine do_self_lanczos_time_full'
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
     
     call initialize_memory(fp)
     if(options%n_full>0) call read_data_pw_full_prods(fp, options%prefix)

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


  
 


     allocate(e_mat(lc%numt,lc%numt,l_blk_g))

     do iw=nbegin_g,nbegin_g+l_blk_g-1
        if(iw <= tf%n_g) then
           af(1)=dcmplx(offset,-tf%freqs_g(iw))
           call solve_lanczos_complex(1,af,e_mat(1,1,iw-nbegin_g+1),lc)
        else
           call solve_lanczos_fake_complex(lc)
        endif
        
     end do

!for entire self-energy store 
!!!!!!!!!!!!!!
  

     

     call initialize_memory(fli)
     call initialize_memory(flj)


!if required read all S matrices and distribute among mpi tasks



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

        call  read_data_pw_mat_lanczos_full(fli, ii, options%prefix)
        if(ii==i_min_cycles.and.i_cycles==1) then
           nums=fli%nums
           allocate(g_mat(numpw,numpw,l_blk_g,2),g_mat_t(numpw,numpw,l_blk,2))

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
        if( options%l_verbose) write(stdout,*) 'Doing dgemms',nums,numpw


        do jj=j_min,j_max
           call  read_data_pw_mat_lanczos_full(flj, jj, options%prefix)

           
          
           allocate(tmp_mat(numpw,nums))
           do iw=nbegin_g,nend_g
              if( options%l_verbose) write(stdout,*) 'Doing dgemms',nums,numpw,l_blk,iw
              FLUSH(stdout)
              
              do ipol=1,2
                
                 call ZGEMM('N','N',numpw,nums,nums,(1.d0,0.d0),fli%f_mat(1,1,ipol),numpw,&
                      &e_mat(1,1,iw-nbegin_g+1),nums,(0.d0,0.d0),tmp_mat,numpw)
                 
                 
                 call ZGEMM('N','C',numpw,numpw,nums,(1.d0,0.d0),tmp_mat,numpw,flj%f_mat(1,1,ipol),numpw,&
                      &(0.d0,0.d0),g_mat(1,1,iw-nbegin_g+1,ipol),numpw)
                  

              enddo

           enddo
           
           deallocate(tmp_mat)
 
           write(stdout,*) 'Fourier trasform:'
           FLUSH(stdout)



!Fourier trasform
           allocate(g_tmp(numpw,numpw))
           allocate(g_dumm(numpw,numpw))
           if( options%l_verbose) write(stdout,*) 'ATT1'
           FLUSH(stdout)
           
           do ipol=1,2
!loop on time
              do it=0,tf%n
!each procs sums up its matrices in frequency with opportune factor
                 time=tf%times(it)
                 g_tmp=0.d0
                 do iw=nbegin_g,nend_g
                    factor=2.d0*dble(tf%weights_freq_g(iw)*exp_table(it+1,iw-nbegin_g+1))/(2.d0*pi)
                    g_tmp(1:numpw,1:numpw)=g_tmp(1:numpw,1:numpw)+dble(g_mat(1:numpw,1:numpw,iw-nbegin_g+1,ipol))*factor
                 enddo
#if defined(__MPI)
!the distribution of times on procs is the same of that for frequecies
                 iproc_time=it/l_blk
!all processors sums to iproc_time
                 if(iproc_time==mpime) then
                    call MPI_REDUCE(g_tmp,g_dumm,numpw*numpw,&
                         &MPI_DOUBLE_PRECISION,MPI_SUM,iproc_time,world_comm,ierr)
                    g_mat_t(1:numpw,1:numpw,it-nbegin+1,ipol)=g_dumm(1:numpw,1:numpw)
                 else
                    call MPI_REDUCE(g_tmp,g_dumm,numpw*numpw,MPI_DOUBLE_PRECISION,&
                         &MPI_SUM,iproc_time,world_comm,ierr)
                 endif

#else
                 g_mat_t(1:numpw,1:numpw,it+1,ipol)=g_tmp(1:numpw,1:numpw)

#endif
                 g_tmp=0.d0
                 do iw=nbegin_g,nend_g
                    factor=-2.d0*dimag(tf%weights_freq_g(iw)*exp_table(it+1,iw-nbegin_g+1))/(2.d0*pi)
                    g_tmp(1:numpw,1:numpw)=g_tmp(1:numpw,1:numpw)+dimag(g_mat(1:numpw,1:numpw,iw-nbegin_g+1,ipol))*factor
                 enddo
#if defined(__MPI)
!the distribution of times on procs is the same of that for frequecies
                 iproc_time=it/l_blk
!all processors sums to iproc_time
                 if(iproc_time==mpime) then
                    call MPI_REDUCE(g_tmp,g_dumm,numpw*numpw,&
                         &MPI_DOUBLE_PRECISION,MPI_SUM,iproc_time,world_comm,ierr)
                    g_mat_t(1:numpw,1:numpw,it-nbegin+1,ipol)=g_mat_t(1:numpw,1:numpw,it-nbegin+1,ipol)&
                         &+(0.d0,1.d0)*g_dumm(1:numpw,1:numpw)
                 else
                    call MPI_REDUCE(g_tmp,g_dumm,numpw*numpw,MPI_DOUBLE_PRECISION,MPI_SUM,iproc_time,world_comm,ierr)
                 endif
                 
#else
                  g_mat_t(1:numpw,1:numpw,it-nbegin+1,ipol)=g_mat_t(1:numpw,1:numpw,it-nbegin+1,ipol)+&
                       &(0.d0,1.d0)*g_tmp(1:numpw,1:numpw)

#endif


              enddo
           enddo
           deallocate(g_tmp,g_dumm)
          
           if( options%l_verbose) write(stdout,*) 'done'
           FLUSH(stdout)


   

!loop on frequency

     
           write(stdout,*) 'Products in imaginary time:'
           FLUSH(stdout)

           
          
           allocate(re_h_mat(numpw,numpw),im_h_mat(numpw,numpw))
           allocate(tmp_mat(numpw,numpw),tmp_vec(numpw))
           do it=nbegin,nend


!product

              if(ii==jj) then
                 ss%diag(ii,it+ss%n+1,is)=0.d0
                 ss%diag(ii,ss%n+1-it,is)=0.d0
                 do ipol=1,2
                    re_h_mat(1:numpw,1:numpw)=dble(g_mat_t(1:numpw,1:numpw,it-nbegin+1,ipol))
                    im_h_mat(1:numpw,1:numpw)=dimag(g_mat_t(1:numpw,1:numpw,it-nbegin+1,ipol))
                    ss%diag(ii,it+ss%n+1,is)=ss%diag(ii,it+ss%n+1,is)+DDOT(numpw*numpw,re_h_mat,1,pw_mat_t(1,1,it-nbegin+1),1)+&
                         &DDOT(numpw*numpw,im_h_mat,1,pw_mat_t(1,1,it-nbegin+1),1)
                    ss%diag(ii,ss%n+1-it,is)=ss%diag(ii,ss%n+1-it,is)+DDOT(numpw*numpw,re_h_mat,1,pw_mat_t(1,1,it-nbegin+1),1)-&
                         &DDOT(numpw*numpw,im_h_mat,1,pw_mat_t(1,1,it-nbegin+1),1)
                 enddo
                 if(options%n_full>0) then
                    tmp_mat(1:numpw,1:numpw)=pw_mat_t(1:numpw,1:numpw,it-nbegin+1)
!first valence states (t>0)
                    do iv=1,fp%numv
                       do ipol=1,2
                          if(it==0) then
                             factor=0.5d0*exp((offset+fp%ene_ks(iv))*tf%times(it))
                          else
                             factor=exp((offset+fp%ene_ks(iv))*tf%times(it))
                          endif
                       
                          call ZGEMV('N',numpw,numpw,(1.d0,0.d0),tmp_mat,numpw,&
                               &fp%gmat(1,ipol,iv,ii),1,(0.d0,0.d0),tmp_vec,1)
                    
                          ss%diag(ii,it+ss%n+1,is)= ss%diag(ii,it+ss%n+1,is)- &
                            &zdotc(numpw,fp%gmat(1,ipol,iv,ii),1,tmp_vec,1)*factor
                    
                          
                       enddo
                    enddo
!then conduction states(t<0)
                       do ic=fp%numv+1,options%n_full
                          do ipol=1,2
                             if(it==0) then
                                factor=0.5d0*exp(-(offset+fp%ene_ks(ic))*tf%times(it))
                             else
                                factor=exp(-(offset+fp%ene_ks(ic))*tf%times(it))
                             endif
                             
                             call ZGEMV('N',numpw,numpw,(1.d0,0.d0),tmp_mat,numpw,&
                                  &fp%gmat(1,ipol,ic,ii),1,(0.d0,0.d0),tmp_vec,1)
                        
                             ss%diag(ii,ss%n+1-it,is)= ss%diag(ii,ss%n+1-it,is)+ &
                                  &zdotc(numpw,tmp_vec,1,fp%gmat(1,ipol,ic,ii),1)*factor
                          enddo

                       enddo
                    endif
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
!to be done this...
                 ss%whole(jj,ii,it+ss%n+1,is)=DDOT(numpw*numpw,re_h_mat,1,pw_mat_t(1,1,it-nbegin+1),1)+&
                      &DDOT(numpw*numpw,im_h_mat,1,pw_mat_t(1,1,it-nbegin+1),1)
                 ss%whole(jj,ii,ss%n+1-it,is)=DDOT(numpw*numpw,re_h_mat,1,pw_mat_t(1,1,it-nbegin+1),1)-&
                      &DDOT(numpw*numpw,im_h_mat,1,pw_mat_t(1,1,it-nbegin+1),1)
!!!!!!!                 
              endif
           enddo
           deallocate(re_h_mat,im_h_mat)
           deallocate(tmp_mat,tmp_vec)
           
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
       
           call free_memory(flj)

           
        enddo !jj
    
        call free_memory(fli)


     enddo!on i_min, i_max
 
  enddo!on i_cycles


!put factor due to FFT on imaginary axis
     ss%diag(:,:,is)=ss%diag(:,:,is)*(0.d0,1.d0)
   
     if(ss%whole_s) then
        ss%whole(:,:,:,is)=ss%whole(:,:,:,is)*(0.d0,1.d0)
     endif
     deallocate(g_mat)
     deallocate(g_mat_t)

     deallocate(pw_mat_t)

     deallocate(e_mat)

     deallocate(exp_table)
     call free_memory(lc)
     call free_memory(fp)
  enddo!on spin
  call free_memory(uu)
  call free_memory(sc)
  return
end subroutine do_self_lanczos_full




