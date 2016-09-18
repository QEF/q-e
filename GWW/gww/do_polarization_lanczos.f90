!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

!this subroutines performs the lanczos style calculation of the polarization

subroutine calculate_compact_pola_lanczos(options,ispin)
!this subroutine calculate the compact products for the lanczos calculation
!of the polarization
   USE kinds,              ONLY : DP
   USE input_gw,           ONLY : input_options
   USE basic_structures,   ONLY : wannier_u,vt_mat_lanczos,tt_mat_lanczos,initialize_memory,free_memory
   USE lanczos,            ONLY : compact_q_lanczos,initialize_compact_q_lanczos,&
      &free_memory_compact_q_lanczos,do_compact_q_lanczos,write_compact_q_lanczos
   USE mp_world,           ONLY : nproc,mpime
   USE io_global,          ONLY : stdout 

  implicit none

  TYPE(input_options), INTENT(in) :: options
  INTEGER, INTENT(in) :: ispin!spin channel

  TYPE(wannier_u) :: uu!U matrix
  TYPE(vt_mat_lanczos) :: vtl
  TYPE(tt_mat_lanczos) :: ttl
  TYPE(compact_q_lanczos) :: cql

  INTEGER :: l_blk,nbegin,nend,iv,iu
  
  call initialize_memory(vtl)
  call initialize_memory(ttl)
  call initialize_compact_q_lanczos(cql)

  write(stdout,*) 'Routine calculate_compact_pola_lanczos'
  FLUSH(stdout)
!
!read U matrix
  call read_data_pw_u(uu,options%prefix)
!KS valence state index v distributed on processors
  l_blk= uu%nums_occ(ispin)/nproc
  if(l_blk*nproc < uu%nums_occ(ispin)) l_blk = l_blk+1
  nbegin=mpime*l_blk+1
  nend=nbegin+l_blk-1
  if(nend > uu%nums_occ(ispin)) nend = uu%nums_occ(ispin)
  do iv=nbegin,nbegin+l_blk-1
     if(iv <= uu%nums_occ(ispin))  then
!loop v' on wannier valence state
        do iu=1,uu%nums_occ(ispin)
!read in V and T
           call  read_data_pw_vt_mat_lanczos(vtl, iu, options%prefix, .true.,ispin)
           call  read_data_pw_tt_mat_lanczos(ttl, iu, options%prefix, .true.,ispin)
!calculate contribution to Q
           
           if(iu==1) then
              allocate( cql%qlm(vtl%numpw,ttl%numt))
              cql%qlm(:,:)=0.d0
           endif
!sum up with factor Uvvi
           
           call do_compact_q_lanczos(vtl,ttl,cql,dble(uu%umat(iv,iu,ispin)))
           call free_memory(vtl)
           call free_memory(ttl)
        enddo
!write Q^v on disk
        cql%ii=iv
       
        call write_compact_q_lanczos(cql)
       
       
        deallocate(cql%qlm)
     else
!put here global routines involving MPI
        do iu=1,uu%nums_occ(ispin)
             call  read_data_pw_vt_mat_lanczos(vtl, iu, options%prefix, .true.,ispin)
             call  read_data_pw_tt_mat_lanczos(ttl, iu, options%prefix, .true.,ispin)
             call free_memory(vtl)
             call free_memory(ttl)
        enddo
     endif
  enddo

  call free_memory(vtl)
  call free_memory(ttl)
  call free_memory_compact_q_lanczos(cql)
  call free_memory(uu)

  
  return
end subroutine calculate_compact_pola_lanczos

subroutine solve_lanczos(nbuf, alpha,e_mat,lc, l_verbose)
!this subroutine sums to  the matrix E_{no}=<t_n|(H-alpha)^-1|t_o>
!the matrix lc%o_mat is distributed among the processors

  USE kinds,            ONLY : DP
  USE basic_structures, ONLY : lanczos_chain, initialize_memory,free_memory 
  USE io_global,        ONLY : stdout
  USE mp,               ONLY : mp_sum, mp_bcast
  USE mp_world,         ONLY : nproc,mpime,world_comm

  implicit none

  INTEGER :: nbuf!number of matrices to treat
  COMPLEX(kind=DP) :: alpha(nbuf)!constant for Ev+iw
  TYPE(lanczos_chain) :: lc!lanczos chain descriptor
  REAL(kind=DP) :: e_mat(lc%numt,lc%numt,nbuf)!matrix to be calculated 
  LOGICAL, INTENT(in) ::  l_verbose!if true a lot of output

  INTEGER :: io,info
  COMPLEX(kind=DP), ALLOCATABLE :: dl(:),du(:),d(:),t(:)
  REAL(kind=DP), ALLOCATABLE :: tr(:,:)
  INTEGER :: l_blk,nbegin,nend
  REAL(kind=DP), ALLOCATABLE :: o_mat(:,:,:)
  INTEGER :: iv, iproc,nsize


  allocate(dl(lc%num_steps-1),du(lc%num_steps-1),d(lc%num_steps),t(lc%num_steps))
  e_mat(:,:,:)=0.d0

  l_blk= (lc%numt)/nproc
  if(l_blk*nproc < (lc%numt)) l_blk = l_blk+1
  allocate(o_mat(lc%numt,lc%num_steps,l_blk))
  allocate(tr(lc%num_steps,l_blk))
  


!loop on procs
  do iproc=0,nproc-1
     l_blk= (lc%numt)/nproc
     if(l_blk*nproc < (lc%numt)) l_blk = l_blk+1
     nbegin=iproc*l_blk+1
     nend=min(nbegin+l_blk-1,lc%numt)
     nsize=nend-nbegin+1
     if(nbegin<=lc%numt) then
        if(iproc==mpime) then
           o_mat(:,:,:)=lc%o_mat(:,:,1:nsize)
        endif
        if(l_verbose) write(stdout,*) 'proc:', iproc, lc%numt,lc%num_steps,l_blk,nsize
        if(l_verbose)FLUSH(stdout)

        call mp_bcast(o_mat, iproc, world_comm)
        if(l_verbose) write(stdout,*) 'mp_bcast done'
        if(l_verbose)FLUSH(stdout)
      
 
        do iv=1,nbuf
           do io=nbegin,nend
!!set up vectors for lapack routine
              dl(1:lc%num_steps-1)=cmplx(lc%f(1:lc%num_steps-1,io),0.d0)
              du(1:lc%num_steps-1)=cmplx(lc%f(1:lc%num_steps-1,io),0.d0)
              d(1:lc%num_steps)=cmplx(lc%d(1:lc%num_steps,io),0.d0)+alpha(iv)
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
              tr(:,io-nbegin+1)=dble(t(:))
           enddo
           if(l_verbose) write(stdout,*) 'zgtsv done'
           if(l_verbose) FLUSH(stdout)

           do io=nbegin,nend
          
              call dgemv( 'N', lc%numt,lc%num_steps,1.d0, o_mat(1,1,io-nbegin+1), lc%numt,tr(1,io-nbegin+1),1,1.d0,e_mat(1,io,iv),1)

          !   call dgemm('N','N',lc%numt,1,lc%num_steps,1.d0,o_mat(1,1,io-nbegin+1),lc%numt,tr(1,io-nbegin+1),lc%num_steps,1.d0,e_mat(1,io,iv),lc%numt)
           enddo
           if(l_verbose) write(stdout,*) 'dgemv done'
           if(l_verbose) FLUSH(stdout)

        enddo!on iv
     end if
  end do

  deallocate(tr)
  deallocate(o_mat)
  return
end  subroutine solve_lanczos

subroutine solve_lanczos_fake(lc,l_verbose)
!this subroutine is a parallel fake routine for the solve lanczos routine

  USE kinds,            ONLY : DP
  USE basic_structures, ONLY : lanczos_chain, initialize_memory,free_memory
  USE io_global,        ONLY : stdout
  USE mp,               ONLY : mp_sum,mp_bcast
  USE mp_world,         ONLY : nproc,mpime,world_comm

  implicit none

  TYPE(lanczos_chain) :: lc!lanczos chain descriptor  
  LOGICAL, INTENT(in) :: l_verbose

  INTEGER :: l_blk,nbegin,nend, iproc
  REAL(kind=DP), ALLOCATABLE :: o_mat(:,:,:)
  INTEGER :: io,nsize



  l_blk= (lc%numt)/nproc
  if(l_blk*nproc < (lc%numt)) l_blk = l_blk+1

  allocate(o_mat(lc%numt,lc%num_steps,l_blk))

  do iproc=0,nproc-1
     l_blk= (lc%numt)/nproc
     if(l_blk*nproc < (lc%numt)) l_blk = l_blk+1
     nbegin=iproc*l_blk+1
     nend=min(nbegin+l_blk-1,lc%numt)
     nsize=nend-nbegin+1
     if(nbegin<=lc%numt) then
        if(iproc==mpime) then
           o_mat(:,:,:)=lc%o_mat(:,:,1:nsize)
        endif
        if(l_verbose) write(stdout,*) 'proc:', iproc, lc%numt,lc%num_steps,l_blk,nsize
        if(l_verbose) FLUSH(stdout)
        call mp_bcast(o_mat, iproc,world_comm)
     endif
  enddo



  deallocate(o_mat)
end subroutine solve_lanczos_fake

 

subroutine do_polarization_lanczos(tf,options,ispin)
!this subroutine create polarization descriptor in imaginary frequency using the lanczos scheme
  USE kinds,              ONLY : DP
  USE input_gw,           ONLY : input_options
  USE lanczos,            ONLY : compact_q_lanczos,initialize_compact_q_lanczos,free_memory_compact_q_lanczos,&
       &read_compact_q_lanczos
  USE mp_world,           ONLY : world_comm, mpime, nproc
  USE basic_structures,   ONLY : wannier_u,lanczos_chain, initialize_memory,free_memory,&
                                 &vt_mat_lanczos,tt_mat_lanczos,partial_occ
  USE io_global,          ONLY : stdout,ionode,ionode_id
  USE times_gw,           ONLY : times_freqs
  USE polarization,       ONLY : polaw,initialize_polaw,free_memory_polaw,write_polaw,read_polaw
  USE mp,                 ONLY : mp_bcast
      USE io_files,  ONLY : prefix, tmp_dir
 

  implicit none

  INTEGER, EXTERNAL :: find_free_unit
   TYPE(times_freqs) ,  INTENT(in) :: tf!for time grid 
   TYPE(input_options), INTENT(in) :: options! for imaginary time range and number of samples   
   INTEGER,             INTENT(in) :: ispin!spin channel


   TYPE(wannier_u) :: uu
   TYPE(lanczos_chain) :: lc
   TYPE(partial_occ) :: po
   INTEGER :: l_blk,nbegin,nend,iv,iw, it,ii,jj,il
   TYPE(polaw) :: pw
  ! REAL(kind=DP), ALLOCATABLE ::  e_mat(:,:,:)
   COMPLEX(kind=DP), ALLOCATABLE :: af(:)
   TYPE(compact_q_lanczos) :: cql
   REAL(kind=DP), ALLOCATABLE :: tmp_mat(:,:)
   REAL(kind=DP), ALLOCATABLE :: tmp_mat1(:,:),tmp_mat2(:,:),tmp_mat3(:,:)
   REAL(kind=DP) :: offset

   INTEGER :: n_bufv!number of states which are bufferized
   INTEGER :: n_dim

   INTEGER :: io_tot,l_blk_io, nbegin_io,nend_io, iproc
   REAL(kind=DP), ALLOCATABLE :: cql_save(:,:,:)
   REAL(kind=DP), ALLOCATABLE :: vtl_save(:,:,:)
   REAL(kind=DP), ALLOCATABLE :: ttl_save(:,:,:)
   COMPLEX(kind=DP) :: afc
   

   TYPE(vt_mat_lanczos) :: vtl
   TYPE(tt_mat_lanczos) :: ttl
   

   INTEGER :: numpw,numl,numt!to be passed to subroutine

   LOGICAL :: l_reduce_memory=.true. !if true reduces the use of memory, only useful  for nprocs==nfreqs

   LOGICAL :: exst
   INTEGER :: iw_nbegin, iv_begin, iunrestart,off_nbegin
   LOGICAL :: l_do_restart
   REAL(kind=DP), ALLOCATABLE :: occ(:)

#if defined(__OPENMP)
   INTEGER :: omp_get_num_threads, omp_get_max_threads 
   EXTERNAL omp_set_num_threads, omp_get_num_threads, omp_get_max_threads 
#endif


   write(stdout,*) 'Routine: do_polarization_lanczos'
   FLUSH(stdout)
   call initialize_memory(vtl)
   call initialize_memory(ttl)
   call initialize_memory(po)

   n_bufv=options%n_set_pola
  ! if(l_reduce_memory) n_bufv=1  OCIO ATTENZIONE

!read uu for energies and lanczos chain
   call read_data_pw_u(uu,options%prefix)
   if(uu%nums_occ(ispin)==0) then
      call free_memory(uu)
      return
   endif

   !call free_memory(uu) !i do not need the orthonormalization matrix
!if required read partial occupancis stuff
   if(options%l_frac_occ) then
!NOT_TO_BE_INCLUDED_START
      call read_data_pw_partial_occ(po,options%prefix,ispin)
!NOT_TO_BE_INCLUDED_END
   endif
   call initialize_memory(lc)
   call read_data_pw_lanczos_chain(lc, 1, options%prefix, .true.,ispin)
   write(stdout,*) 'Lanczos dimensions', lc%numt,lc%num_steps
   do ii=1,lc%numt
!      do jj=1,lc%num_steps
!         write(stdout,*)' Lanczos elemets:', ii, jj, lc%d(jj,ii), lc%f(jj,ii)
!      enddo
      write(stdout,*)' Lanczos elemets:', ii, lc%d(1,ii), lc%f(1,ii) 
   enddo
!   allocate(e_mat(lc%numt,lc%numt,n_bufv))
   allocate(af(n_bufv),occ(n_bufv))
   occ(:)=1.d0
   call initialize_compact_q_lanczos(cql)
  
   if(uu%nums > uu%nums_occ(ispin)) then
      offset=-(uu%ene(uu%nums_occ(ispin)+1,ispin)+uu%ene(uu%nums_occ(ispin),ispin))/2.d0
   else
      offset=-uu%ene(uu%nums_occ(ispin),ispin)
   endif


!   l_blk= (2*tf%n+1)/nproc
!   if(l_blk*nproc < (2*tf%n+1)) l_blk = l_blk+1
!   nbegin=mpime*l_blk-tf%n
!   nend=nbegin+l_blk-1

   l_blk= (tf%n+1)/nproc                                                                                                   
   if(l_blk*nproc < (tf%n+1)) l_blk = l_blk+1                                                                              
   nbegin=mpime*l_blk                                                                                                   
   nend=nbegin+l_blk-1                 

!if reduce io and l_direct put compact_q_lanczos in memory

  
   if(options%l_t_wannier) then
      l_blk_io= (uu%nums_occ(ispin))/nproc
      if(l_blk_io*nproc < (uu%nums_occ(ispin))) l_blk_io = l_blk_io+1
      nbegin_io=mpime*l_blk_io+1
      nend_io=nbegin_io+l_blk_io-1
      if(nend_io > uu%nums_occ(ispin)) nend_io=uu%nums_occ(ispin)
      io_tot=nend_io-nbegin_io+1
      if(.not.l_reduce_memory) then
         do iv=1,uu%nums_occ(ispin)
            call read_compact_q_lanczos(cql, iv)
            if(iv==1)  allocate(cql_save(cql%numpw,cql%numt,l_blk_io))
            if(iv >= nbegin_io .and. iv <= nend_io) then
               cql_save(:,:,iv-nbegin_io+1)=cql%qlm(:,:)
            endif
         enddo
      else
         call read_compact_q_lanczos(cql, 1)!just for obtaining numpw,numt...poor man solution
         allocate(cql_save(1,1,1))
      endif
         
    

   else
!put all matrices vtl and ttl in memory, distributed according to valence state
      l_blk_io= (uu%nums_occ(ispin))/nproc
      if(l_blk_io*nproc < (uu%nums_occ(ispin))) l_blk_io = l_blk_io+1
      nbegin_io=mpime*l_blk_io+1
      nend_io=nbegin_io+l_blk_io-1
      if(nend_io > uu%nums_occ(ispin)) nend_io=uu%nums_occ(ispin)
      io_tot=nend_io-nbegin_io+1
      
      do iv=1,uu%nums_occ(ispin)
         call  read_data_pw_vt_mat_lanczos(vtl, iv, options%prefix, .true.,ispin)
         call  read_data_pw_tt_mat_lanczos(ttl, iv, options%prefix, .true.,ispin)
         if(iv==1)  then
            allocate(vtl_save(vtl%numpw,vtl%numl,l_blk_io))
            allocate(ttl_save(ttl%numt,ttl%numl,l_blk_io))
         endif
         if(iv >= nbegin_io .and. iv <= nend_io) then
            vtl_save(:,:,iv-nbegin_io+1)= vtl%vt_mat(:,:)
            ttl_save(:,:,iv-nbegin_io+1)= ttl%tt_mat(:,:)
         endif

         call free_memory(vtl)
         call free_memory(ttl)
      enddo

   endif

!sutup parameters for solve_lanczos
   if(options%l_t_wannier) then
      numpw=cql%numpw
      numt=cql%numt
      pw%numpw=cql%numpw
   else
      numpw=vtl%numpw
      numl=vtl%numl
      numt=ttl%numt
      pw%numpw=vtl%numpw
   endif



!loop on imaginary frequency i\omega
!check for restart
   if(ionode) then

     inquire(file=trim(tmp_dir)//trim(prefix)//'-'//'restart_polaw', exist = exst)
     if(.not. exst .or. ispin==2) then
        off_nbegin=0
        iv_begin=1
        l_do_restart=.false.
     else
        iunrestart =  find_free_unit()
        open( unit= iunrestart, file=trim(tmp_dir)//trim(prefix)//'-'//'restart_polaw', status='old')
        read(iunrestart,*) off_nbegin
        read(iunrestart,*) iv_begin
        close(iunrestart)
        if(off_nbegin<0) off_nbegin=0
        if(iv_begin<1 ) then
           iv_begin=1
           l_do_restart=.false.
        else
           iv_begin=iv_begin+n_bufv
           l_do_restart=.true.
        endif
     endif
  endif
  call mp_bcast(off_nbegin,ionode_id,world_comm)
  iw_nbegin=nbegin+off_nbegin
  call mp_bcast(iv_begin, ionode_id,world_comm)
  call mp_bcast(l_do_restart, ionode_id,world_comm)




   do iw=iw_nbegin,nbegin+l_blk-1
      if(iw <= tf%n) then

         


         write(stdout,*) 'do_polarization_lanczos1', iw, l_do_restart
         FLUSH(stdout)

!!set up polaw descriptor
!if required read polaw from disk
         if(l_do_restart.or.ispin==2) then
            call initialize_polaw(pw)
            call read_polaw(iw,pw,.false.,options%l_verbose)
            l_do_restart=.false.
         else
            call initialize_polaw(pw)
            pw%label=iw
            pw%ontime=.false.
            pw%time=tf%freqs(iw)
            pw%factor=(1.d0,0.d0)
         endif


!!loop on valence states v
         do iv=iv_begin,uu%nums_occ(ispin),n_bufv       
            write(stdout,*) 'do_polarization_lanczos iv', iv
            FLUSH(stdout)


!!!solve tridiagonal problem for -E_v+i\omega
            n_dim=0
            do il=1,n_bufv
               if((iv+il-1) <= uu%nums_occ(ispin)) then
                  af(il) = -uu%ene(iv+il-1,ispin)+dcmplx(0.d0,tf%freqs(iw))
                  if(.not.options%l_frac_occ) then
                     occ(il)=1.d0
                  else
!NOT_TO_BE_INCLUDED_START
                     occ(il)=po%f_occ(iv+il-1)
!NOT_TO_BE_INCLUDED_END
                  endif
                  n_dim=n_dim+1
               endif
            enddo
            if(options%l_verbose) write(stdout,*) 'Call solve_lanczos'
            if(options%l_verbose) FLUSH(stdout)
            if(iv==1.and.ispin/=2) then
!                if(options%l_t_wannier) then
!                   pw%numpw=cql%numpw
!                else
!                   pw%numpw=vtl%numpw
!                endif
                allocate(pw%pw(pw%numpw,pw%numpw))
                pw%pw(:,:)=0.d0
             endif
!NOT_TO_BE_INCLUDED_START
             if(iv==1 .and. options%l_frac_occ)  then
                afc=dcmplx(0.d0,tf%freqs(iw))
                call add_partial_pola(pw,po,uu,afc,options%nspin,ispin)
             endif
!NOT_TO_BE_INCLUDED_END
             call solve_lanczos_2(numpw,numt,numl,n_dim,n_bufv,af,lc, iv, uu%nums_occ(ispin),l_blk_io,&
    &cql_save,pw,n_dim,options%l_t_wannier,vtl_save,ttl_save,l_reduce_memory,uu%nspin,occ,options%l_verbose)
           
             if(options%l_verbose) write(stdout,*) 'Done'
             if(options%l_verbose) FLUSH(stdout)
         
   
#if defined(__OPENMP)
   !               write(stdout,*) 'Call dgemms, ntids=', omp_get_max_threads() 
#else
   !               write(stdout,*) 'Call dgemms'
#endif

!restart part, every 10
             if(mod(iv-iv_begin,10)==0.or.n_bufv/=1) then
                call write_polaw(pw,options%debug)
                if(ionode) then
                   iunrestart =  find_free_unit()
                   open( unit= iunrestart, file=trim(tmp_dir)//trim(prefix)//'-'//'restart_polaw', status='unknown')
                   write(iunrestart,*) iw-nbegin
                   write(iunrestart,*) iv
                   close(iunrestart)
                endif
             endif
          enddo
!!write polarization on disk
          if(ionode) then
             iunrestart =  find_free_unit()
             open( unit= iunrestart, file=trim(tmp_dir)//trim(prefix)//'-'//'restart_polaw', status='unknown')
             write(iunrestart,*) iw-nbegin
             write(iunrestart,*) iv
             close(iunrestart)
          endif

          call write_polaw(pw,options%debug)
          call free_memory_polaw(pw)
       
    else
!here just global parallel calls
    do iv=iv_begin,uu%nums_occ(ispin),n_bufv
       n_dim=0
       do il=1,n_bufv
          if((iv+il-1) <= uu%nums_occ(ispin)) then
             n_dim=n_dim+1
          endif
       enddo

       call solve_lanczos_2(numpw,numt,numl,0,n_bufv,af,lc, iv, uu%nums_occ(ispin),l_blk_io,&
      &cql_save,pw,n_dim,options%l_t_wannier,vtl_save,ttl_save,l_reduce_memory,uu%nspin,occ,.false.)
       
    enddo
    endif
    iv_begin=1         
    l_do_restart=.false.
 enddo
 if(ionode) then
    iunrestart =  find_free_unit()
    open( unit= iunrestart, file=trim(tmp_dir)//trim(prefix)//'-'//'restart_polaw', status='unknown')
    write(iunrestart,*) -1
    write(iunrestart,*) -1
    close(iunrestart)
 endif

   if(options%l_t_wannier) then
      deallocate(cql_save)
      call free_memory_compact_q_lanczos(cql)
   else
      deallocate(vtl_save,ttl_save)
   endif

  ! deallocate(e_mat)
   deallocate(af,occ)
   call free_memory(lc)
   call free_memory(uu) 
   call free_memory(vtl)
   call free_memory(ttl)
   call free_memory(po)
  return

end subroutine do_polarization_lanczos

subroutine solve_lanczos_2(numpw,numt,numl,nbuf,mbuf, alpha,lc, iv0,nbndv,&
    &lblkio,cql_save,pw,n_dim,l_t_wannier, vtl_save, ttl_save, l_reduce_memory,nspin,occ,l_verbose)
!this subroutine sums to  the matrix E_{no}=<t_n|(H-alpha)^-1|t_o>
!the matrix lc%o_mat is distributed among the processors

  USE kinds,            ONLY : DP
  USE basic_structures, ONLY : lanczos_chain, initialize_memory,free_memory 
  USE io_global,        ONLY : stdout
  USE mp,               ONLY : mp_sum, mp_bcast
  USE mp_world,         ONLY : nproc,mpime,world_comm
  USE parallel_include
  USE lanczos,          ONLY : compact_q_lanczos,initialize_compact_q_lanczos,&
                               free_memory_compact_q_lanczos, &
                               read_compact_q_lanczos

  USE polarization,     ONLY : polaw

  implicit none

  INTEGER, INTENT(in) :: numpw!dimension of polarizability basis
  INTEGER, INTENT(in) :: numt!number of t vectors
  INTEGER, INTENT(in) :: numl!number of tv vectors
  INTEGER :: nbuf!number of matrices to treat
  INTEGER :: mbuf
  COMPLEX(kind=DP) :: alpha(mbuf)!constant for Ev+iw
  TYPE(lanczos_chain) :: lc!lanczos chain descriptor
  INTEGER, INTENT(in) :: iv0 !staring valence energy
  INTEGER, INTENT(in) :: nbndv!number of valence functions
  INTEGER, INTENT(in) :: lblkio!for the dimension of the following array
  REAL(kind=DP), INTENT(in) :: cql_save(numpw,numt,lblkio)
  TYPE(polaw), INTENT(inout) :: pw
  INTEGER, INTENT(in) :: n_dim!max dimension of nbuf
  LOGICAL, INTENT(in) :: l_t_wannier !if true t vectors from products with wannier
  REAL(kind=DP) :: vtl_save(numpw,numl,lblkio)
  REAL(kind=DP) :: ttl_save(numt,numl,lblkio)
  LOGICAL, INTENT(in) :: l_reduce_memory!if true cql_save not used
  INTEGER, INTENT(in) :: nspin!spin multiplicity
  REAL(kind=DP), INTENT(in) :: occ(mbuf)!occupancies of KS valence states 
  LOGICAL, INTENT(in) :: l_verbose

  INTEGER :: io,info
  COMPLEX(kind=DP), ALLOCATABLE :: dl(:),du(:),d(:),t(:)
  REAL(kind=DP), ALLOCATABLE :: tr(:,:)
  INTEGER :: l_blk,nbegin,nend
  REAL(kind=DP), ALLOCATABLE :: e_mat_ip(:,:)
  INTEGER :: iv, iproc,nsize, iproc_cql

  INTEGER :: nbuf_ip
  COMPLEX(kind=DP), ALLOCATABLE :: alpha_ip(:)
  INTEGER, ALLOCATABLE :: displ(:), recvcounts(:)
  INTEGER :: l_blk_io, nbegin_io,nend_io
  REAL(kind=DP), ALLOCATABLE :: tmp_mat(:,:),tmp_mat1(:,:),tmp_mat2(:,:),tmp_mat3(:,:)
  REAL(kind=DP), ALLOCATABLE :: pw_ip(:,:),pw_dumm(:,:)
  REAL(kind=DP), ALLOCATABLE :: qlm_tmp(:,:,:),vtl_tmp(:,:,:),ttl_tmp(:,:,:)
  INTEGER :: ierr

  TYPE(compact_q_lanczos) :: cql

  LOGICAL :: l_qo=.true.!if true perform QO matrix multimplication first (convenient only if nprocs>nlanczos_steps)
  REAL(kind=DP), ALLOCATABLE :: qo_mat(:,:,:,:)
  LOGICAL :: l_single=.true.!if true saves qlm_tmp in single precision
  REAL(kind=4), ALLOCATABLE :: qlm_tmp_single(:,:,:)
  REAL(kind=DP) :: factor!global prefactor


  if(nspin==1) then
     factor=-4.d0
  else
     factor=-2.d0
  endif
  allocate(dl(lc%num_steps-1),du(lc%num_steps-1),d(lc%num_steps),t(lc%num_steps))
 

  ierr=0

  l_blk= (lc%numt)/nproc
  if(l_blk*nproc < (lc%numt)) l_blk = l_blk+1
  nbegin=mpime*l_blk+1
  nend=min(nbegin+l_blk-1,lc%numt)
  nsize=nend-nbegin+1


   l_blk_io= (nbndv)/nproc
   if(l_blk_io*nproc < (nbndv)) l_blk_io = l_blk_io+1
   nbegin_io=mpime*l_blk_io+1
   nend_io=nbegin_io+l_blk_io-1




  allocate(tr(lc%num_steps,l_blk))
  allocate(displ(nproc),recvcounts(nproc))
  allocate(e_mat_ip(lc%numt,l_blk))

  if(l_t_wannier) then
     if(.not.l_reduce_memory) then
        allocate(qlm_tmp(numpw,numt,n_dim))
        do iv=1,n_dim
!read compact matrix
           iproc_cql=(iv0+iv-2)/l_blk_io
           if(mpime==iproc_cql) qlm_tmp(:,:,iv) = cql_save(:,:,iv0+iv-1-nbegin_io+1)
           call mp_bcast(qlm_tmp(:,:,iv), iproc_cql,world_comm)
        enddo
     else
        call initialize_compact_q_lanczos(cql)
        if(.not.l_single) then
           allocate(qlm_tmp(numpw,numt,n_dim))
           do iv=1,n_dim
              call read_compact_q_lanczos(cql,iv0+iv-1)
              qlm_tmp(1:numpw,1:numt,iv)=cql%qlm(1:numpw,1:numt)
           enddo
        else
           allocate(qlm_tmp(numpw,numt,1))
           allocate(qlm_tmp_single(numpw,numt,n_dim))
           do iv=1,n_dim
              call read_compact_q_lanczos(cql,iv0+iv-1)
              qlm_tmp_single(1:numpw,1:numt,iv)=real(cql%qlm(1:numpw,1:numt))
           enddo
        endif
        call free_memory_compact_q_lanczos(cql)
     endif
     if(l_qo) then
         allocate(qo_mat(numpw,lc%num_steps,nsize,n_dim))
        do iv=1,n_dim
           !calculate qo_mat
           qlm_tmp(1:numpw,1:numt,1)=dble(qlm_tmp_single(1:numpw,1:numt,iv))
           do io=nbegin,nend
              if(.not.l_single) then
                 call dgemm('N','N', numpw,lc%num_steps,numt,1.d0, qlm_tmp(1,1,iv),&
   &numpw,lc%o_mat(1,1,io-nbegin+1),numt,0.d0,qo_mat(1,1,io-nbegin+1,iv),numpw)
              else
                 !qlm_tmp(1:numpw,1:numt,1)=dble(qlm_tmp_single(1:numpw,1:numt,iv))
                 call dgemm('N','N', numpw,lc%num_steps,numt,1.d0, qlm_tmp(1,1,1),numpw,&
         &lc%o_mat(1,1,io-nbegin+1),numt,0.d0,qo_mat(1,1,io-nbegin+1,iv),numpw)
              endif
           enddo
        !call dgemv( 'N', lc%numt,lc%num_steps,1.d0, lc%o_mat(1,1,io-nbegin+1), lc%numt,tr(1,io-nbegin+1),1,1.d0,e_mat_ip(1,io-nbegin+1),1)
        enddo
     endif
  else
     allocate(vtl_tmp(numpw,numl,n_dim))
     allocate(ttl_tmp(numt,numl,n_dim))
     do iv=1,n_dim
 !read compact matrix
        iproc_cql=(iv0+iv-2)/l_blk_io
        if(mpime==iproc_cql) then
           vtl_tmp(:,:,iv) = vtl_save(:,:,iv0+iv-1-nbegin_io+1)
           ttl_tmp(:,:,iv) = ttl_save(:,:,iv0+iv-1-nbegin_io+1)
        endif
        call mp_bcast(vtl_tmp(:,:,iv), iproc_cql,world_comm)
        call mp_bcast(ttl_tmp(:,:,iv), iproc_cql,world_comm)
     enddo


  endif


  do iproc=0,nproc-1
     displ(iproc+1)=iproc*l_blk
     if(displ(iproc+1)+l_blk > lc%numt) then
        recvcounts(iproc+1)=(lc%numt-displ(iproc+1))
     else
        recvcounts(iproc+1)=l_blk
     endif
  enddo
  displ(:)=displ(:)*lc%numt
  recvcounts(:)=recvcounts(:)*lc%numt


  if(l_t_wannier) then
     allocate(tmp_mat(numpw,l_blk))
  else
     allocate(tmp_mat1(numl,l_blk))
     allocate(tmp_mat2(numpw,l_blk))
     allocate(tmp_mat3(numpw,numl))
  endif
  allocate(pw_ip(numpw,numpw),pw_dumm(numpw,numpw))

!loop on procs
  do iproc=0,nproc-1

     if(iproc==mpime) nbuf_ip=nbuf
!distribute number of nbuf
!distribute af
     call mp_bcast(nbuf_ip, iproc,world_comm)
     if(nbuf_ip >0) then
        allocate(alpha_ip(nbuf_ip))
        if(iproc==mpime) alpha_ip(1:nbuf_ip)=alpha(1:nbuf_ip)
        call mp_bcast(alpha_ip,iproc,world_comm)
        !loop on freuqency
        pw_ip(:,:)=0.d0
        do iv=1,nbuf_ip
           if(l_single) then
              qlm_tmp(1:numpw,1:numt,1)=dble(qlm_tmp_single(1:numpw,1:numt,iv))
           endif

           e_mat_ip(:,:)=0.d0
!calculate part
!$OMP PARALLEL SHARED(tr,lc,alpha_ip,iv) PRIVATE(dl,du,d,t,info,io)
!          allocate(dl(lc%num_steps-1),du(lc%num_steps-1),d(lc%num_steps),t(lc%num_steps))
!$OMP DO
           do io=nbegin,nend
!!set up vectors for lapack routine
              dl(1:lc%num_steps-1)=cmplx(lc%f(1:lc%num_steps-1,io),0.d0)
              du(1:lc%num_steps-1)=cmplx(lc%f(1:lc%num_steps-1,io),0.d0)
              d(1:lc%num_steps)=cmplx(lc%d(1:lc%num_steps,io),0.d0)+alpha_ip(iv)
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
              tr(1:lc%num_steps,io-nbegin+1)=dble(t(1:lc%num_steps))
           enddo
!$OMP END DO 
!           deallocate(dl,du,d,t)
!$OMP END PARALLEL
           if(l_verbose) write(stdout,*) 'zgtsv done'
           if(l_verbose) FLUSH(stdout)

           if(.not.l_qo) then
              do io=nbegin,nend
                 
                 call dgemv( 'N', lc%numt,lc%num_steps,1.d0, lc%o_mat(1,1,io-nbegin+1), &
   &lc%numt,tr(1,io-nbegin+1),1,1.d0,e_mat_ip(1,io-nbegin+1),1)
              enddo
           endif
           if(l_verbose) write(stdout,*) 'dgemv done',nsize
           if(l_verbose) FLUSH(stdout)

           if(nbegin <= lc%numt) then
              if(l_t_wannier) then
                 if(.not.l_qo) then
                    if(.not.l_single) then
                       call dgemm('N','N',numpw,nsize,numt,1.d0,qlm_tmp(1,1,iv),numpw,e_mat_ip,lc%numt,0.d0,tmp_mat,numpw)
                    else
                       call dgemm('N','N',numpw,nsize,numt,1.d0,qlm_tmp(1,1,1),numpw,e_mat_ip,lc%numt,0.d0,tmp_mat,numpw)
                    endif
                    if(l_verbose) write(stdout,*) 'dgemm1 done'
                    if(l_verbose) FLUSH(stdout)
                 else
                    do io=nbegin,nend
                       call dgemv( 'N', numpw,lc%num_steps,1.d0, qo_mat(1,1,io-nbegin+1,iv), &
   &numpw,tr(1,io-nbegin+1),1,0.d0,tmp_mat(1,io-nbegin+1),1)
                    enddo

                 endif
                 
                 if(.not.l_single) then
                    call dgemm('N','T',numpw,numpw,nsize,factor*occ(iv),tmp_mat,numpw,qlm_tmp(1,nbegin,iv),numpw,1.d0,pw_ip,numpw)
                 else
                    call dgemm('N','T',numpw,numpw,nsize,factor*occ(iv),tmp_mat,numpw,qlm_tmp(1,nbegin,1),numpw,1.d0,pw_ip,numpw)
                 endif
                 if(l_verbose) write(stdout,*) 'dgemm2 done', pw%numpw, numpw,ierr
                 if(l_verbose) FLUSH(stdout)
              else
                 call dgemm('T','N',numl,nsize,numt,1.d0,ttl_tmp(1,1,iv),numt,e_mat_ip,numt,0.d0,tmp_mat1,numl)
                 if(l_verbose) write(stdout,*) 'dgemm1 done'
                 if(l_verbose) FLUSH(stdout)
                 call dgemm('N','N',numpw,nsize,numl,1.d0,vtl_tmp(1,1,iv),numpw,tmp_mat1,numl,0.d0,tmp_mat2,numpw)
                 if(l_verbose) write(stdout,*) 'dgemm2 done'
                 if(l_verbose) FLUSH(stdout)
                 call dgemm('N','N', numpw,numl,nsize,1.d0,tmp_mat2,numpw,ttl_tmp(nbegin,1,iv),numt,0.d0,tmp_mat3,numpw)
                 if(l_verbose)write(stdout,*) 'dgemm3 done'
                 if(l_verbose) FLUSH(stdout)
                 call dgemm('N','T', numpw,numpw,numl,factor*occ(iv),tmp_mat3,numpw,vtl_tmp(1,1,iv),numpw,1.d0,pw_ip,numpw)

                 

    !              call dgemm('T','N', ttl%numl,ttl%numt,ttl%numt,1.d0,ttl%tt_mat,ttl%numt,e_mat(1,1,il),lc%numt,0.d0,tmp_mat1,ttl%numl)
    !              call dgemm('N','N', pw%numpw,ttl%numt,ttl%numl,1.d0,vtl%vt_mat,vtl%numpw,tmp_mat1,ttl%numl,0.d0,tmp_mat2,pw%numpw)
    !              call dgemm('N','N', pw%numpw,ttl%numl,ttl%numt,1.d0,tmp_mat2,pw%numpw,ttl%tt_mat,ttl%numt,0.d0,tmp_mat3,pw%numpw)
    !              call dgemm('N','T', pw%numpw,pw%numpw,ttl%numl,-4.d0,tmp_mat3,pw%numpw,vtl%vt_mat,pw%numpw,1.d0,pw%pw,pw%numpw)

              endif
           endif



!gather_collect result
!           call mpi_gatherv(e_mat_ip, lc%numt*nsize,MPI_DOUBLE_PRECISION,e_mat(1,1,iv),recvcounts,displ, MPI_DOUBLE_PRECISION,iproc,world_comm)
        enddo

#if defined(__MPI)
        if(iproc==mpime) then
           pw_ip(:,:)=pw_ip(:,:)+pw%pw(:,:)
           CALL MPI_REDUCE(pw_ip, pw%pw, pw%numpw*pw%numpw,MPI_DOUBLE_PRECISION,MPI_SUM, iproc,world_comm,ierr)
        else
           !ATTENTION pw%numpw could not be initialized in THIS CASE
           CALL MPI_REDUCE(pw_ip, pw_dumm, numpw*numpw,MPI_DOUBLE_PRECISION,MPI_SUM, iproc,world_comm,ierr)
        endif

        if(l_verbose) write(stdout,*) 'mpi_reduce done'
        if(l_verbose) FLUSH(stdout)
#else
        pw%pw(:,:)=pw%pw(:,:)+pw_ip(:,:)
#endif


        deallocate(alpha_ip)
     endif
  enddo

  if(l_qo) deallocate(qo_mat)
  deallocate(displ,recvcounts)
  deallocate(e_mat_ip)
  deallocate(dl,du,d,t)
  deallocate(tr)
  if(l_t_wannier) then
     deallocate(tmp_mat)
!     if(l_reduce_memory) call free_memory_compact_q_lanczos(cql)
  else
     deallocate(tmp_mat1,tmp_mat2,tmp_mat3)
  endif
  deallocate(pw_ip,pw_dumm)
  if(l_t_wannier) then
     deallocate(qlm_tmp)
     if(l_single) deallocate(qlm_tmp_single)
  else
     deallocate(vtl_tmp, ttl_tmp)
  endif
  return
end  subroutine solve_lanczos_2



subroutine add_partial_pola(pw,po,uu,af,nspin,ispin)
!NOT_TO_BE_INCLUDED_START
!this subroutine adds the contribution to the irreducible polarizability due to 
!fractionally occupied states

  USE kinds,              ONLY : DP
  USE basic_structures,   ONLY : partial_occ,wannier_u
  USE io_global,          ONLY : stdout,ionode,ionode_id
  USE polarization,       ONLY : polaw

  implicit none

  TYPE(polaw), INTENT(inout) :: pw!polarizability matrix to be modified
  TYPE(partial_occ), INTENT(in) ::  po!fractionary occupaction terms
  TYPE(wannier_u), INTENT(in) :: uu!for KS energies
  COMPLEX(kind=DP),INTENT(in) :: af!!complex frequency 
  INTEGER, INTENT(in) :: nspin!spin multiplicity
  INTEGER, INTENT(in) :: ispin!spin channel


  REAL(kind=DP) :: fact
  INTEGER :: iv,jv

  do iv=po%nums_occ_min+1,po%nums_occ
     do jv=1,iv-1
        fact=-(2.d0/dble(nspin))*(po%f_occ(jv)-po%f_occ(iv))/dble(uu%ene(iv,ispin)-uu%ene(jv,ispin)+af)
!invertiti per energie
        call dgemm('N','T',pw%numpw,pw%numpw,1,fact,po%ppp_mat(1,jv,iv),po%numpw,&
             &po%ppp_mat(1,jv,iv),po%numpw,1.d0,pw%pw,pw%numpw)
     enddo
  enddo

  return
!NOT_TO_BE_INCLUDED_END
end subroutine add_partial_pola
