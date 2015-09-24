!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

!this subroutines performs the lanczos style calculation of the self-energy

subroutine calculate_compact_self_lanczos(options)
!this subroutine calculate the compact products for the lanczos calculation
!of the self-energy
   USE kinds,              ONLY : DP
   USE input_gw,           ONLY : input_options
   USE basic_structures,   ONLY : wannier_u,vt_mat_lanczos,tt_mat_lanczos,initialize_memory,free_memory
   USE lanczos,            ONLY : compact_q_lanczos,initialize_compact_q_lanczos,&
      &free_memory_compact_q_lanczos,do_compact_q_lanczos,write_compact_q_lanczos
   USE mp_world,          ONLY : nproc,mpime
   USE io_global,          ONLY : ionode


  implicit none

  TYPE(input_options), INTENT(in) :: options

  TYPE(wannier_u) :: uu!U matrix
  TYPE(vt_mat_lanczos) :: vtl
  TYPE(tt_mat_lanczos) :: ttl
  TYPE(compact_q_lanczos) :: cql

  INTEGER :: l_blk,nbegin,nend,iv,iu
  
  call initialize_memory(vtl)
  call initialize_memory(ttl)
  call initialize_compact_q_lanczos(cql)

!
!read U matrix
  call read_data_pw_u(uu,options%prefix)


  do iu=1,options%max_i
!read in V and T
     call  read_data_pw_vt_mat_lanczos(vtl, iu, options%prefix, .false.,1)
     call  read_data_pw_tt_mat_lanczos(ttl, iu, options%prefix, .false.,1)
!calculate contribution to Q
     allocate( cql%qlm(vtl%numpw,ttl%numt))
     cql%qlm(:,:)=0.d0
!sum up with factor Uvvi
     call do_compact_q_lanczos(vtl,ttl,cql,1.d0)
!write Q^v on disk   
     cql%ii=iu
     if(ionode) call write_compact_q_lanczos(cql)
     deallocate(cql%qlm)
  enddo


  call free_memory(vtl)
  call free_memory(ttl)
  call free_memory_compact_q_lanczos(cql)
  call free_memory(uu)

  return
end subroutine calculate_compact_self_lanczos


subroutine solve_lanczos_complex(nbuf,alpha,e_mat,lc)
!this subroutine sums to  the matrix E_{no}=<t_n|(H-alpha)^-1|t_o> 

  USE kinds,            ONLY : DP
  USE basic_structures, ONLY : lanczos_chain, initialize_memory,free_memory
  USE io_global,        ONLY : stdout
  USE mp,               ONLY : mp_sum,mp_bcast
  USE mp_world,         ONLY : nproc, mpime, world_comm

  implicit none

  INTEGER :: nbuf!number of matrices to treat
  COMPLEX(kind=DP) :: alpha(nbuf)!constant for Ev+iw 
  TYPE(lanczos_chain) :: lc!lanczos chain descriptor   
  COMPLEX(kind=DP) :: e_mat(lc%numt,lc%numt,nbuf)!matrix to be calculated 


  INTEGER :: io,info,ii,jj,il
  COMPLEX(kind=DP), ALLOCATABLE :: dl(:),du(:),d(:),t(:)
  COMPLEX(kind=DP), ALLOCATABLE :: omat(:,:)
  REAL(kind=DP), ALLOCATABLE :: tmp_mat(:,:)
  INTEGER :: l_blk,nbegin,nend, iproc

  l_blk= (lc%numt)/nproc
  if(l_blk*nproc < (lc%numt)) l_blk = l_blk+1
  nbegin=mpime*l_blk+1
  nend=nbegin+l_blk-1

  allocate(dl(lc%num_steps-1),du(lc%num_steps-1),d(lc%num_steps),t(lc%num_steps))
   e_mat(:,:,:)=(0.d0,0.d0)
  allocate(omat(lc%numt,lc%num_steps))
  allocate(tmp_mat(lc%numt,lc%num_steps))
!loop on o
  do io=1,lc%numt
!!set up vectors for lapack routine 
     !recover matrix from processor                                                                                                                                        
     tmp_mat(:,:)=0.d0
     if(io >= nbegin .and. io <= nend) then
        tmp_mat(:,:)=lc%o_mat(:,:,io-nbegin+1)
     endif
     !call mp_sum(tmp_mat(:,:),world_comm)
     iproc=(io-1)/l_blk
     call mp_bcast(tmp_mat(:,:), iproc, world_comm)
     omat(:,:)=dcmplx(tmp_mat(:,:),0.d0)

     do il=1,nbuf

        dl(1:lc%num_steps-1)=cmplx(lc%f(1:lc%num_steps-1,io),0.d0)
        du(1:lc%num_steps-1)=cmplx(lc%f(1:lc%num_steps-1,io),0.d0)
        d(1:lc%num_steps)=cmplx(lc%d(1:lc%num_steps,io),0.d0)+alpha(il)
        t(:)=(0.d0,0.d0)
        t(1)=(1.d0,0.d0)

!!call lapack
        call zgtsv(lc%num_steps,1,dl,d,du,t,lc%num_steps,info)
        if(info /= 0) then
           write(stdout,*) 'ZGTSV info:', info
           FLUSH(stdout)
           stop
        endif
     !ATTENZIONE
   !  if(io==1) then
   !     write(stdout,*) 'Solve lanczos first',t(1),t(lc%num_steps)
   !     write(stdout,*) 'omat', lc%o_mat(1,1,io), lc%o_mat(1,lc%num_steps,io),&
   !          & lc%o_mat(lc%num_steps,1,io), lc%o_mat(lc%num_steps,lc%num_steps,io)
   !  endif
   !  if(io==lc%numt) then
   !     write(stdout,*) 'Solve lanczos last',t(1),t(lc%num_steps)
   !      write(stdout,*) 'omat', lc%o_mat(1,1,io), lc%o_mat(1,lc%num_steps,io),&
   !          & lc%o_mat(lc%num_steps,1,io), lc%o_mat(lc%num_steps,lc%num_steps,io)
   !  endif

!!calculate term


     !omat(:,:)=dcmplx(lc%o_mat(:,:,io),0.d0)
        call zgemm('N','N',lc%numt,1,lc%num_steps,(1.d0,0.d0),omat,lc%numt,t,lc%num_steps,(1.d0,0.d0),e_mat(:,io,il),lc%numt)
     enddo
  enddo
  do il=1,nbuf
     do ii=1,lc%numt
        do jj=ii+1,lc%numt
           e_mat(ii,jj,il)=0.5d0*(e_mat(ii,jj,il)+e_mat(jj,ii,il))
           e_mat(jj,ii,il)=e_mat(ii,jj,il)
        enddo
     enddo
  enddo
  deallocate(dl,du,d,t)
  deallocate(omat,tmp_mat)
  return
end  subroutine solve_lanczos_complex




!routine do_self_lanczos
subroutine do_self_lanczos(ss, tf ,options)
!this subroutine calculte the self-energy on the frequency fit grid using the lanczos scheme

  USE kinds,             ONLY : DP
  USE io_global,         ONLY : stdout, ionode, ionode_id
  USE input_gw,          ONLY : input_options
  USE basic_structures,  ONLY : v_pot,wannier_u,free_memory, initialize_memory,lanczos_chain, vt_mat_lanczos,tt_mat_lanczos
  USE green_function,    ONLY : green, read_green, free_memory_green, initialize_green
  USE polarization,      ONLY : polaw, free_memory_polaw, read_polaw, write_polaw,invert_v_pot, initialize_polaw, &
                                  & read_polaw_global
  USE mp,                ONLY : mp_sum, mp_bcast
  USE mp_world,          ONLY : nproc,mpime,world_comm
  USE times_gw,          ONLY : times_freqs
  USE self_energy_storage, ONLY : self_storage,write_self_storage_ondisk,free_memory_self_storage
  USE lanczos
  USE constants,          ONLY : tpi
  USE start_end ! debug

  implicit none
  
  TYPE(times_freqs), INTENT(in) :: tf!for time frequency grids 
  TYPE(input_options) :: options
  TYPE(self_storage)  :: ss

  TYPE(compact_q_lanczos), ALLOCATABLE :: cql(:)
  TYPE(vt_mat_lanczos), ALLOCATABLE :: vtl(:)
  TYPE(tt_mat_lanczos), ALLOCATABLE :: ttl(:)

  TYPE(wannier_u) :: uu
  TYPE(lanczos_matrix) :: lm

!1)
!loop on iw (parallelized)
!calculate screened polarization
!write back on polaw


  TYPE(v_pot) :: vp,vpi
  TYPE(polaw) :: ww!dressed polarization 
  TYPE(lanczos_chain) :: lc
  INTEGER :: l_blk, nbegin,nend,  l_blk_freq, nbegin_freq,nend_freq
  INTEGER :: iw_tot,l_blk_iw, nbegin_iw,nend_iw
  REAL(kind=DP), ALLOCATABLE:: wtemp(:,:), vtemp(:)
  INTEGER :: iw,ii, jw, il, kw
  COMPLEX(kind=DP), ALLOCATABLE :: af(:)
  REAL(kind=DP), ALLOCATABLE :: re_e_mat(:,:),im_e_mat(:,:) ,tmp_mat(:,:), re_c_mat(:,:,:,:), im_c_mat(:,:,:,:)
  REAL(kind=DP), ALLOCATABLE ::  tmp_mat1(:,:),tmp_mat2(:,:)
  INTEGER :: n_buf,iw_min,iw_max,iiw,ix,iy,it,ip,ixx
  REAL(kind=DP), EXTERNAL :: ddot
  REAL(kind=DP) :: offset
  COMPLEX(kind=DP), ALLOCATABLE :: gw_tab(:,:), exp_tab(:,:,:), exp_iw(:,:),exp_ix(:,:)
  LOGICAL :: l_conv!if true do just the convolution
  LOGICAL :: l_half=.true.
  REAL(kind=DP), ALLOCATABLE :: re_c_tmp(:,:), im_c_tmp(:,:)
  INTEGER, ALLOCATABLE :: iw_proc(:) 
  INTEGER :: n_buf_ii, ii_buf, ii_begin,ii_end,ipos

  LOGICAL :: l_direct=.true.!if true direct way of performing products for each state
  INTEGER :: numpw

  COMPLEX(kind=DP), ALLOCATABLE :: e_mat_tmp(:,:,:)
  INTEGER :: n_dim
  INTEGER :: iproc
  COMPLEX(kind=DP) :: csca
  

  REAL(kind=DP), ALLOCATABLE :: polaw_store(:,:,:) !to reduce io
  INTEGER :: io_tot,l_blk_io, nbegin_io,nend_io

  COMPLEX(kind=DP), ALLOCATABLE :: cterm(:)
  REAL(kind=DP) :: inv_epsi,v_head
  INTEGER :: iw1,jw1
  INTEGER :: n_max
  REAL(kind=DP), ALLOCATABLE :: diag_tmp(:)
  LOGICAL :: l_iw_min
  REAL(kind=DP), ALLOCATABLE :: mat_tmp(:,:)

  INTEGER :: i, j

!calculate offset
!read in DFT energies
                                                             
  call read_data_pw_u(uu,options%prefix)


  
  if(uu%nums > uu%nums_occ(1)) then
     offset=-(uu%ene(uu%nums_occ(1)+1,1)+uu%ene(uu%nums_occ(1),1))/2.d0
  else
     offset=-uu%ene(uu%nums_occ(1),1)
  endif
  
 

  call free_memory(uu)


   
!initialize self_storage ss

  if(options%l_lanczos_conv) then
     l_conv=.true.
  else
     l_conv=.false.
  endif

  ss%ontime=.false.
  ss%max_i=options%max_i
  ss%i_min=options%i_min
  ss%i_max=options%i_max
  ss%n=tf%n_grid_fit
  ss%tau=options%tau
  ss%whole_s=options%whole_s
  ss%n_grid_fit=tf%n_grid_fit
  

  if(ss%whole_s) then
   allocate(ss%whole(ss%max_i,ss%max_i,2*ss%n+1,1))
   ss%whole(:,:,:,:)=(0.d0,0.d0)
   allocate(ss%whole_freq_fit(ss%max_i,ss%max_i,2*ss%n+1,1))
   ss%whole_freq_fit(:,:,:,:)=(0.d0,0.d0)
   nullify(ss%diag)
   nullify(ss%diag_freq_fit)
  else
   allocate(ss%diag(ss%max_i,2*ss%n+1,1))
   ss%diag(:,:,:)=(0.d0,0.d0)
   nullify(ss%whole)
   allocate(ss%diag_freq_fit(ss%max_i,2*ss%n+1,1))
   ss%diag_freq_fit(:,:,:)=(0.d0,0.d0)
   nullify(ss%whole_freq_fit)
  endif
!for compatibility
  allocate(ss%ene_remainder(ss%max_i,1))
  ss%ene_remainder(:,:)=0.d0
  



  if(ss%whole_s) then
     write(stdout,*) 'Whole s not implemented YET'
     stop
  endif




  
   nullify(vp%vmat)
   nullify(vpi%vmat)
   call initialize_polaw(ww)
   
    allocate(diag_tmp(options%n_set_ii))

!2) create compact_self_lanczos
   if(.not.l_direct) call calculate_compact_self_lanczos(options)

!3)
!read in lanczos chain for self-energy
   call initialize_memory(lc)
   call read_data_pw_lanczos_chain(lc, 1, options%prefix, .false.,1)
   write(stdout,*) 'Lanczos dimensions', lc%numt,lc%num_steps
   write(stdout,*) 'Lanczos first diagonal', lc%d(1,1),lc%d(lc%num_steps,1)
   write(stdout,*) 'Lanczos last diagonal', lc%d(1,lc%numt),lc%d(lc%num_steps,lc%numt)
   write(stdout,*) 'Lanczos first tri', lc%f(1,1),lc%f(lc%num_steps-1,1)
   write(stdout,*) 'Lanczos last tri', lc%f(1,lc%numt),lc%f(lc%num_steps-1,lc%numt)

   FLUSH(stdout)
   call initialize_lanczos_matrix(lm)
   lm%numt=lc%numt
   allocate(lm%e_mat(lm%numt,lm%numt))
   allocate(re_e_mat(lc%numt,lc%numt))
   allocate(im_e_mat(lc%numt,lc%numt))

!if required allocates table for products and exponents
   if(.not.l_conv) then
      allocate(gw_tab(-tf%n:tf%n,-tf%n:tf%n))
      l_blk_freq= (2*ss%n+1)/nproc
      if(l_blk_freq*nproc < (2*ss%n+1)) l_blk_freq = l_blk_freq+1
      nbegin_freq=mpime*l_blk_freq-ss%n
      nend_freq=nbegin_freq+l_blk_freq-1
      allocate(exp_tab(-tf%n:tf%n,-tf%n:tf%n,l_blk_freq))
      exp_tab(:,:,:)=(0.d0,0.d0)
     ! allocate(exp_iw(-tf%n:tf%n,-tf%n:tf%n))
     ! do iw=-tf%n,tf%n
     !    do it=-tf%n,tf%n
     !       exp_iw(it,iw)=exp((0.d0,1.d0)*tf%freqs(iw)*tf%times(it))
     !    enddo
     ! enddo
     ! allocate(exp_ix(-tf%n:tf%n,l_blk_freq))
     ! do ix=nbegin_freq,nbegin_freq+l_blk_freq-1
     !    if(ix<=ss%n) then
     !       do it=-tf%n,tf%n
     !          exp_ix(it,ix-nbegin_freq+1)=exp((0.d0,-1.d0)*tf%freqs_fit(ix)*tf%times(it))
     !       enddo
     !    endif
     ! enddo
      do ix=nbegin_freq,nbegin_freq+l_blk_freq-1
         if(ix<=ss%n) then
            write(stdout,*) 'Exp table:', ix
            FLUSH(stdout)
            do jw=-tf%n,tf%n
               do iw=-tf%n,tf%n
                  do it=-tf%n,tf%n
                     exp_tab(iw,jw,ix-nbegin_freq+1)=exp_tab(iw,jw,ix-nbegin_freq+1)+tf%weights_freq(iw)*&
                          &tf%weights_freq(jw)*tf%weights_time(it)*&
                          &exp((0.d0,1.d0)*(tf%freqs(iw)+tf%freqs(jw)-tf%freqs_fit(ix))*tf%times(it))
                        !  *exp_iw(it,iw)*exp_iw(it,jw)*exp_ix(it,ix-nbegin_freq+1)
                  enddo
               enddo
            enddo
         endif
      enddo
      !deallocate(exp_iw,exp_ix)
   endif


   if(.not.l_direct) then
      allocate(cql(options%n_set_ii))
      do ii=1,options%n_set_ii
         call initialize_compact_q_lanczos(cql(ii))
      enddo
   else
      allocate(vtl(options%n_set_ii))
      allocate(ttl(options%n_set_ii))
      do ii=1,options%n_set_ii
         call initialize_memory(vtl(ii))
         call initialize_memory(ttl(ii))
      enddo
   endif
!!calculate (H-iw)^-1 and write on disk (parallel)
   l_blk= (2*tf%n+1)/nproc
   if(l_blk*nproc < (2*tf%n+1)) l_blk = l_blk+1
   nbegin=mpime*l_blk-tf%n
   nend=nbegin+l_blk-1
   if(nend > tf%n) nend=tf%n
   allocate(e_mat_tmp(lm%numt,lm%numt,options%n_set_self))
   allocate(af(options%n_set_self))
   do iw=nbegin,nbegin+l_blk-1,options%n_set_self
      n_dim=0
      do il=1,options%n_set_self
         if((iw+il-1) <= tf%n) then
            n_dim=n_dim+1
            af(il) = dcmplx(offset,-tf%freqs(iw+il-1))
         endif
      enddo
 
     if(n_dim > 0) then
         write(stdout,*) 'Doing solve_lanczos_complex',iw
         FLUSH(stdout)
         call solve_lanczos_complex(n_dim,af,e_mat_tmp,lc)
         write(stdout,*) 'Done'
         FLUSH(stdout)
         do il=1,n_dim
            lm%iw=iw+il-1
            lm%e_mat(:,:)=e_mat_tmp(:,:,il)
            call  write_lanczos_matrix(lm)
         enddo
      else
         call solve_lanczos_fake_complex(lc)
      endif



   enddo
   deallocate(e_mat_tmp,af)


!if reduced I/O read all polaw in memory
   if(options%l_reduce_io) then
      if(.not. tf%grid_freq == 5) then
         n_max=tf%n/2
      else
         n_max=(tf%n-(1+2*tf%second_grid_i)*tf%second_grid_n)/2+&
                 &(1+2*tf%second_grid_i)*tf%second_grid_n
      endif

      l_blk_io= (n_max+1)/nproc
      if(l_blk_io*nproc < (n_max+1)) l_blk_io = l_blk_io+1
      nbegin_io=mpime*l_blk_io
      nend_io=nbegin_io+l_blk_io-1
      if(nend_io > n_max) nend_io=n_max
      io_tot=nend_io-nbegin_io+1
      do jw=nbegin_io,nend_io
         call read_polaw(jw, ww,options%debug,options%l_verbose)
         if(jw==nbegin_io) then
            allocate(polaw_store(ww%numpw,ww%numpw,io_tot))
         endif
         polaw_store(:,:,jw-nbegin_io+1)=ww%pw(:,:)
      enddo
!distribute ww%numpw
      if(.not.ionode) ww%numpw=0
      !call mp_sum(ww%numpw,world_comm)
      call mp_bcast(ww%numpw, ionode_id,world_comm)
      if(.not.ionode) ww%factor=(0.d0,0.d0)
      !call mp_sum(ww%factor,world_comm)
      call mp_bcast(ww%factor, ionode_id,world_comm)
         
   endif
!!loop on i KS states
!!loop/module n_set  on iw'                         

   n_buf_ii=(options%i_max-options%i_min+1)/options%n_set_ii
   if(n_buf_ii*options%n_set_ii < (options%i_max-options%i_min+1)) n_buf_ii=n_buf_ii+1
   allocate(cterm(options%n_set_ii))
   
   do ii_buf=1,n_buf_ii

      ii_begin=options%i_min+(ii_buf-1)*options%n_set_ii
      ii_end=min(ii_begin+options%n_set_ii-1, options%i_max)
      

!      do ii=options%i_min,options%i_max
      do ii=ii_begin,ii_end
         write(stdout,*) 'KS STATE', ii
         FLUSH(stdout)

!!read in Q^i
         write(stdout,*) 'Call read_compact_q_lanczos'
         FLUSH(stdout)
         if(.not.l_direct) then
            call read_compact_q_lanczos(cql(ii-ii_begin+1), ii)
         else
            call  read_data_pw_vt_mat_lanczos(vtl(ii-ii_begin+1), ii, options%prefix, .false.,1)
            call  read_data_pw_tt_mat_lanczos(ttl(ii-ii_begin+1), ii, options%prefix, .false.,1)
         endif
         write(stdout,*) 'Done'
         FLUSH(stdout)
      enddo

     
!!loop/module n_set  on iw' 
      n_buf=(2*tf%n+1)/options%n_set
      if(n_buf*options%n_set < 2*tf%n+1) n_buf=n_buf+1
    
      do iiw=1,n_buf
!!build C matrix
         iw_min=-tf%n+options%n_set*(iiw-1)
         iw_max=min(iw_min+options%n_set-1,tf%n)
         write(stdout,*) 'Solve'
         FLUSH(stdout)
         iw_tot=iw_max-iw_min+1
         l_blk_iw= (iw_tot)/nproc
         if(l_blk_iw*nproc < (iw_tot)) l_blk_iw = l_blk_iw+1
         nbegin_iw=mpime*l_blk_iw+iw_min
         nend_iw=nbegin_iw+l_blk_iw-1
         if(nend_iw > iw_max) nend_iw = iw_max
         
              
         if(.not.l_direct) then
            numpw=cql(1)%numpw
         else
            numpw=vtl(1)%numpw
         endif
         allocate(re_c_mat(numpw,numpw,options%n_set_ii,l_blk_iw))
         allocate(im_c_mat(numpw,numpw,options%n_set_ii,l_blk_iw))
         re_c_mat(:,:,:,:)=0.d0
         im_c_mat(:,:,:,:)=0.d0
         if(.not.l_conv) then
            allocate(re_c_tmp(numpw,numpw))
            allocate(im_c_tmp(numpw,numpw))
         endif
         !allocate(iw_proc(iw_tot))
         !iw_proc(:)=-1
         !do ip=0,nproc-1
         !   iw_proc(ip*l_blk_iw+1:min(ip*l_blk_iw+l_blk_iw,iw_tot))=ip
         !enddo

         do iw=nbegin_iw,nend_iw
             
!!build C matrices
!!!calculate G matrix on {t} basis  
!!!solve tridiagonal problem for +i\omega  
            !af = dcmplx(offset,-tf%freqs(iw))
            !call solve_lanczos_complex(af,lm%e_mat,lc)
            write(stdout,*) 'Reading lanczos matrix',iw
            FLUSH(stdout)
            call read_lanczos_matrix(lm,iw)
            write(stdout,*) 'Done'
            FLUSH(stdout)
            
!!calculate C^i_lm matrix    
               re_e_mat(:,:)=dble(lm%e_mat(:,:))
               im_e_mat(:,:)=dimag(lm%e_mat(:,:))
               if(.not.l_direct) then
                  allocate(tmp_mat(cql(1)%numpw,cql(1)%numt))
                  write(stdout,*) 'Doing dgemms'
                  write(stdout,*) 're_e_mat', re_e_mat(1,1),re_e_mat(1,cql(1)%numt),&
   &re_e_mat(cql(1)%numt,1),re_e_mat(cql(1)%numt,cql(1)%numt)
                 
                 
                  
                  do ii=ii_begin,ii_end
                     ipos=ii-ii_begin+1
                     FLUSH(stdout)
                     call dgemm('N','N',cql(ipos)%numpw,cql(ipos)%numt,cql(ipos)%numt,&
   &1.d0,cql(ipos)%qlm,cql(ipos)%numpw,re_e_mat,lc%numt,0.d0,tmp_mat,cql(ipos)%numpw)
                     call dgemm('N','T',cql(ipos)%numpw,cql(ipos)%numpw,cql(ipos)%numt,&
   &1.d0,tmp_mat,cql(ipos)%numpw,cql(ipos)%qlm,cql(ipos)%numpw,0.d0,re_c_mat(1,1,ipos,iw-nbegin_iw+1),cql(ipos)%numpw)
                     call dgemm('N','N',cql(ipos)%numpw,cql(ipos)%numt,cql(ipos)%numt,&
   &1.d0,cql(ipos)%qlm,cql(ipos)%numpw,im_e_mat,lc%numt,0.d0,tmp_mat,cql(ipos)%numpw)
                     call dgemm('N','T',cql(ipos)%numpw,cql(ipos)%numpw,cql(ipos)%numt,&
   &1.d0,tmp_mat,cql(ipos)%numpw,cql(ipos)%qlm,cql(ipos)%numpw,0.d0,im_c_mat(1,1,ipos,iw-nbegin_iw+1),cql(ipos)%numpw)
                     

                  enddo
                  write(stdout,*) 'Done'
                  FLUSH(stdout)
                  deallocate(tmp_mat)

               else
                  write(stdout,*) 'Doing dgemms',ttl(1)%numl,ttl(1)%numt,vtl(1)%numl,vtl(1)%numpw
                  FLUSH(stdout)
                  allocate(tmp_mat(ttl(1)%numl,ttl(1)%numt))
                  allocate(tmp_mat1(ttl(1)%numl,ttl(1)%numl))
                  allocate(tmp_mat2(numpw,ttl(1)%numl))
                  do ii=ii_begin,ii_end
                     ipos=ii-ii_begin+1
                     write(stdout,*) 'ATTENZIONE1', ttl(ipos)%numl,ttl(ipos)%numt,vtl(ipos)%numl,vtl(ipos)%numpw
                     FLUSH(stdout)

                     call dgemm('T','N',ttl(ipos)%numl,ttl(ipos)%numt,ttl(ipos)%numt,1.d0,ttl(ipos)%tt_mat,ttl(ipos)%numt,&
                          &re_e_mat,lc%numt,0.d0,tmp_mat,ttl(ipos)%numl)

                     write(stdout,*) 'ATTENZIONE2', ii
                     FLUSH(stdout)

                     call dgemm('N','N',ttl(ipos)%numl,ttl(ipos)%numl,ttl(ipos)%numt,1.d0,tmp_mat,ttl(ipos)%numl,&
                          ttl(ipos)%tt_mat,ttl(ipos)%numt,0.d0,tmp_mat1,ttl(ipos)%numl)

                     write(stdout,*) 'ATTENZIONE3', ii, vtl(ipos)%numpw,vtl(ipos)%numl,vtl(ipos)%numl
                     write(stdout,*) 'ATTENZIONE3', size(vtl(ipos)%vt_mat(:,1)), size(vtl(ipos)%vt_mat(1,:)), ttl(1)%numl, numpw
                     FLUSH(stdout)
                     
                     call dgemm('N','N',vtl(ipos)%numpw,vtl(ipos)%numl,vtl(ipos)%numl,1.d0,vtl(ipos)%vt_mat,vtl(ipos)%numpw,&
                          &tmp_mat1,vtl(ipos)%numl,0.d0,tmp_mat2,numpw)

                     write(stdout,*) 'ATTENZIONE4', ii
                     FLUSH(stdout)

                   
                     call dgemm('N','T',vtl(ipos)%numpw,vtl(ipos)%numpw,vtl(ipos)%numl,1.d0,tmp_mat2,numpw,vtl(ipos)%vt_mat,&
                          vtl(ipos)%numpw,0.d0,re_c_mat(:,:,ipos,iw-nbegin_iw+1),numpw)

                     write(stdout,*) 'ATTENZIONE5', ii
                     FLUSH(stdout)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

                     call dgemm('T','N',ttl(ipos)%numl,ttl(ipos)%numt,ttl(ipos)%numt,1.d0,ttl(ipos)%tt_mat,ttl(ipos)%numt,&
                          &im_e_mat,lc%numt,0.d0,tmp_mat,ttl(ipos)%numl)
                     
                     write(stdout,*) 'ATTENZIONE6', ii
                     FLUSH(stdout)


                     call dgemm('N','N',ttl(ipos)%numl,ttl(ipos)%numl,ttl(ipos)%numt,1.d0,tmp_mat,ttl(ipos)%numl,&
                          ttl(ipos)%tt_mat,ttl(ipos)%numt,0.d0,tmp_mat1,ttl(ipos)%numl)

                     write(stdout,*) 'ATTENZIONE7', ii
                     FLUSH(stdout)


                     call dgemm('N','N',vtl(ipos)%numpw,vtl(ipos)%numl,vtl(ipos)%numl,1.d0,vtl(ipos)%vt_mat,vtl(ipos)%numpw,&
                          &tmp_mat1,vtl(ipos)%numl,0.d0,tmp_mat2,numpw)

                     write(stdout,*) 'ATTENZIONE8', ii
                     FLUSH(stdout)


                  
                     call dgemm('N','T',vtl(ipos)%numpw,vtl(ipos)%numpw,vtl(ipos)%numl,1.d0,tmp_mat2,numpw,vtl(ipos)%vt_mat,&
                          vtl(ipos)%numpw,0.d0,im_c_mat(:,:,ipos,iw-nbegin_iw+1),numpw)
                     

                     write(stdout,*) 'ATTENZIONE9', ii
                     FLUSH(stdout)

                   
                 enddo
                   write(stdout,*) 'Done'
                  FLUSH(stdout)
                  deallocate(tmp_mat,tmp_mat1,tmp_mat2)
               endif
            enddo
            write(stdout,*) 'Done'
            FLUSH(stdout)

!!!loop on w' parallelized
            if(.not.l_half) then
               l_blk= (tf%n+1)/nproc
               if(l_blk*nproc < (tf%n+1)) l_blk = l_blk+1
               nbegin=mpime*l_blk
               nend=nbegin+l_blk-1
               if(nend > tf%n) nend=tf%n
            else
               if(.not. tf%grid_freq == 5) then
                  l_blk= (tf%n/2+1)/nproc
                  if(l_blk*nproc < (tf%n/2+1)) l_blk = l_blk+1
                  nbegin=mpime*l_blk
                  nend=nbegin+l_blk-1
                  if(nend > tf%n/2) nend=tf%n/2
               else
                  n_max=tf%n/2+1+(1+2*tf%second_grid_i)*tf%second_grid_n
                  l_blk= (n_max)/nproc
                  if(l_blk*nproc < (n_max)) l_blk = l_blk+1
                  nbegin=mpime*l_blk
                  nend=nbegin+l_blk-1
                  if(nend > n_max) nend=n_max
               endif
            endif
            if(.not.l_conv) gw_tab(:,:)=(0.d0,0.d0)
            if(.not. tf%grid_freq == 5) then
               n_max=tf%n/2
            else
                 n_max=(tf%n-(1+2*tf%second_grid_i)*tf%second_grid_n)/2+&
                 &(1+2*tf%second_grid_i)*tf%second_grid_n
            endif
            do jw=0,n_max
!         do jw=nbegin,nbegin+l_blk-1
               if(jw <= tf%n) then
                  write(stdout,*) 'FREQ', jw
                  FLUSH(stdout)
                  if(.not.options%l_reduce_io) then
                     call read_polaw(jw, ww,options%debug,options%l_verbose)
                  else
                     call free_memory_polaw(ww)
                     allocate(ww%pw(ww%numpw,ww%numpw))
                     if(jw>=nbegin_io .and. jw <= nend_io) then
                        !ww%pw(:,:)=polaw_store(:,:,jw-nbegin_io+1)
                        call dcopy( SIZE(ww%pw), polaw_store(1,1,jw-nbegin_io+1), 1, ww%pw(1,1), 1 )
                     else
                        !ww%pw(:,:)=0.d0
                        do j = 1, SIZE( ww%pw, 2 )
                          do i = 1, SIZE( ww%pw, 1 )
                            ww%pw(i,j)=0.d0
                          end do
                        end do
                     endif
                     iproc=jw/l_blk_io
                     do kw=1,ww%numpw
                        !call mp_sum(ww%pw(:,kw),world_comm)
                         call mp_bcast(ww%pw(:,kw),iproc,world_comm)
                     enddo
                     !call mp_bcast(ww%pw(:,:),iproc,world_comm)
                  endif

                  if(l_conv) then
!loop on iw
                     if(tf%grid_freq == 3) then
                        l_iw_min=.false.
                        do iw=nbegin_iw,nend_iw
                           ix=iw+jw
                           ixx=iw-jw
                           if((ix >= -tf%n .and. ix <= tf%n).or.(ixx >= -tf%n .and. ixx <= tf%n)) then
                              if(.not.l_iw_min) then
                                 l_iw_min=.true.
                                 iw_min=iw
                              endif
                              iw_max=iw
                           endif
                        enddo
                        if(l_iw_min) then
                           allocate(mat_tmp(options%n_set_ii,l_blk_iw))
!the weights are all the same for this kind of grid
                           call dgemv('T',ww%numpw*ww%numpw,options%n_set_ii*(iw_max-iw_min+1),&
        &tf%weights_freq(1),re_c_mat(1,1,1,iw_min-nbegin_iw+1),&
                                   &ww%numpw*ww%numpw,ww%pw,1,0.d0,mat_tmp,1)
                           do iw=iw_min,iw_max
                              ix=iw+jw
                              ixx=iw-jw
                              if(ix >= -tf%n .and. ix <= tf%n) then
                                 ss%diag(ii_begin:ii_end,ix+tf%n+1,1)=ss%diag(ii_begin:ii_end,ix+tf%n+1,1)+&
        &dcmplx(mat_tmp(:,iw-iw_min+1),0.d0)*ww%factor
                              endif
                              if(jw/=0) then
                                 if(ixx >= -tf%n .and. ixx <= tf%n) then
                                    ss%diag(ii_begin:ii_end,ixx+tf%n+1,1)=ss%diag(ii_begin:ii_end,ixx+tf%n+1,1)+&
        &dcmplx(mat_tmp(:,iw-iw_min+1),0.d0)*ww%factor
                                 endif
                              endif
                           enddo
                           call dgemv('T',ww%numpw*ww%numpw,options%n_set_ii*(iw_max-iw_min+1),tf%weights_freq(1),&
       &im_c_mat(1,1,1,iw_min-nbegin_iw+1),ww%numpw*ww%numpw,ww%pw,1,0.d0,mat_tmp,1)

                           do iw=iw_min,iw_max
                              ix=iw+jw
                              ixx=iw-jw
                              if(ix >= -tf%n .and. ix <= tf%n) then
                                 ss%diag(ii_begin:ii_end,ix+tf%n+1,1)=ss%diag(ii_begin:ii_end,ix+tf%n+1,1)+&
       &dcmplx(0.d0,mat_tmp(:,iw-iw_min+1))*ww%factor
                              endif
                              if(jw/=0) then
                                 if(ixx >= -tf%n .and. ixx <= tf%n) then
                                    ss%diag(ii_begin:ii_end,ixx+tf%n+1,1)=ss%diag(ii_begin:ii_end,ixx+tf%n+1,1)+&
       &dcmplx(0.d0,mat_tmp(:,iw-iw_min+1))*ww%factor
                                 endif
                              endif
                           enddo

                           deallocate(mat_tmp)
                        endif
                     else if (tf%grid_freq == 5) then
                        l_iw_min=.false.
                        do iw=nbegin_iw,nend_iw
                           iw1=tf%whois_freq(iw)
                           jw1=tf%whois_freq(jw)

                           ix=iw1+jw1
                           jw1=tf%whois_freq(-jw)
                           ixx=iw1+jw1

                           if((ix >= -ss%n .and. ix <= ss%n).or.(ixx >= -ss%n .and. ixx <= ss%n)) then
                              if(.not.l_iw_min) then
                                 l_iw_min=.true.
                                 iw_min=iw
                              endif
                              iw_max=iw
                           endif
                        enddo
                        if(l_iw_min) then
                           allocate(mat_tmp(options%n_set_ii,l_blk_iw))
!the weights are all the same for this kind of grid     
                              call dgemv('T',ww%numpw*ww%numpw,options%n_set_ii*(iw_max-iw_min+1),&
        &1.d0,re_c_mat(1,1,1,iw_min-nbegin_iw+1),&
                                   &ww%numpw*ww%numpw,ww%pw,1,0.d0,mat_tmp,1)
                              
                              do iw=nbegin_iw,nend_iw
                                 iw1=tf%whois_freq(iw)
                                 jw1=tf%whois_freq(jw)
                                 ix=iw1+jw1
                                 if(ix >= -ss%n .and. ix <= ss%n) then
                                    ss%diag(ii_begin:ii_end,ix+ss%n+1,1)=ss%diag(ii_begin:ii_end,ix+ss%n+1,1)+&
    &dcmplx(mat_tmp(:,iw-iw_min+1),0.d0)*tf%weights_freq(jw)*ww%factor*tf%relative_weight(iw)*tf%relative_weight(jw)
                                 endif
                                 jw1=tf%whois_freq(-jw)
                                 ix=iw1+jw1

                                 if(ix >= -ss%n .and. ix <= ss%n) then
                                    ss%diag(ii_begin:ii_end,ix+ss%n+1,1)=ss%diag(ii_begin:ii_end,ix+ss%n+1,1)+&
   &dcmplx(mat_tmp(:,iw-iw_min+1),0.d0)*tf%weights_freq(-jw)*ww%factor*tf%relative_weight(iw)*tf%relative_weight(-jw)
               
                                 endif

                              enddo
                              call dgemv('T',ww%numpw*ww%numpw,options%n_set_ii*(iw_max-iw_min+1),&
     &1.d0,im_c_mat(1,1,1,iw_min-nbegin_iw+1),ww%numpw*ww%numpw,ww%pw,1,0.d0,mat_tmp,1)

                              do iw=nbegin_iw,nend_iw
                                 iw1=tf%whois_freq(iw)
                                 jw1=tf%whois_freq(jw)
                                 ix=iw1+jw1
                                 if(ix >= -ss%n .and. ix <= ss%n) then
                                    ss%diag(ii_begin:ii_end,ix+ss%n+1,1)=ss%diag(ii_begin:ii_end,ix+ss%n+1,1)+&
    &dcmplx(0.d0,mat_tmp(:,iw-iw_min+1))*tf%weights_freq(jw)*ww%factor*tf%relative_weight(iw)*tf%relative_weight(jw)
                                 endif
                                 jw1=tf%whois_freq(-jw)
                                 ix=iw1+jw1

                                 if(ix >= -ss%n .and. ix <= ss%n) then
                                    ss%diag(ii_begin:ii_end,ix+ss%n+1,1)=ss%diag(ii_begin:ii_end,ix+ss%n+1,1)+&
   &dcmplx(0.d0,mat_tmp(:,iw-iw_min+1))*tf%weights_freq(-jw)*ww%factor*tf%relative_weight(iw)*tf%relative_weight(-jw)

                                 endif

                              enddo


                              deallocate(mat_tmp)
                           endif
                     else
                  

                     endif
                  else
!here does a explicit convolution
                     do iw=nbegin_iw,nend_iw
!                  do iw=iw_min,iw_max
!distribute re_c_mat  and im_c_mat
                     !if(iw_proc(iw-iw_min+1)==mpime) then
                     !   re_c_tmp(:,:)=re_c_mat(:,:,iw-nbegin_iw+1)
                     !   im_c_tmp(:,:)=im_c_mat(:,:,iw-nbegin_iw+1)
                     !endif
                     !call
                     !mp_bcast(re_c_tmp(:,:),iw_proc(iw-iw_min+1),world_comm)
                     !call
                     !mp_bcast(im_c_tmp(:,:),iw_proc(iw-iw_min+1),world_comm)

                        re_c_tmp(:,:)=re_c_mat(:,:,ipos,iw-nbegin_iw+1)
                        im_c_tmp(:,:)=im_c_mat(:,:,ipos,iw-nbegin_iw+1)!ATTENZIONE TUTTO DA FARE



                        gw_tab(iw,jw)=gw_tab(iw,jw)+&
                             &ddot(ww%numpw*ww%numpw,re_c_tmp(:,1),1,ww%pw(:,1),1)*ww%factor
                        gw_tab(iw,jw)=gw_tab(iw,jw)+(0.d0,1.d0)*&
                             &ddot(ww%numpw*ww%numpw,im_c_tmp(:,1),1,ww%pw,1)*ww%factor
                        if(jw/=0)then
                           gw_tab(iw,-jw)=gw_tab(iw,-jw)+&
                                &ddot(ww%numpw*ww%numpw,re_c_tmp(:,1),1,ww%pw(:,1),1)*ww%factor
                           gw_tab(iw,-jw)=gw_tab(iw,-jw)+(0.d0,1.d0)*&
                                &ddot(ww%numpw*ww%numpw,im_c_tmp(:,1),1,ww%pw,1)*ww%factor
                        endif
                     enddo
                  endif
               else
!here global io and parallel routines
     
               endif
            enddo
            if(.not.l_conv) then
               call mp_sum(gw_tab(:,:),world_comm)
               do ix=nbegin_freq,nbegin_freq+l_blk_freq-1
                  if(ix<=ss%n) then
                     do jw=-tf%n,tf%n
                        do iw=-tf%n,tf%n
                           ss%diag(ii,ix+ss%n+1,1)=ss%diag(ii,ix+ss%n+1,1)+gw_tab(iw,jw)*exp_tab(iw,jw,ix-nbegin_freq+1)
                        enddo
                     enddo
                  endif
               enddo
            endif
            deallocate(re_c_mat,im_c_mat)
            if(.not.l_conv) deallocate(re_c_tmp,im_c_tmp)
            !deallocate(iw_proc)
         enddo
         
         if(.not.l_direct) then
            do ii=1,options%n_set_ii
               call free_memory_compact_q_lanczos(cql(ii))
            enddo
         else
            do ii=1,options%n_set_ii
               call free_memory(vtl(ii))
               call free_memory(ttl(ii))
            enddo
         endif
         call free_memory_polaw(ww)
      enddo
      deallocate(cterm)
      
      if(options%l_reduce_io) then
         if(nend_io >= nbegin_io) deallocate(polaw_store)
         


      endif
         

      if(.not.l_direct) then
          deallocate(cql)
       else
          deallocate(vtl,ttl)
       endif

!set global factor
   call mp_sum(ss%diag(:,:,:),world_comm)
   if(l_conv) then
      ss%diag(:,:,:)=ss%diag(:,:,:)*(1.d0/tpi)
   else
      ss%diag(:,:,:)=ss%diag(:,:,:)*(1.d0/tpi**2.d0)
   endif
   ss%diag_freq_fit(:,:,:)=ss%diag(:,:,:)

   call free_memory(lc)
   deallocate(re_e_mat, im_e_mat)
   call free_memory_lanczos_matrix(lm)
   if(.not.l_conv) then
      deallocate(gw_tab,exp_tab)
   endif

   deallocate(diag_tmp)
  return
end  subroutine do_self_lanczos

subroutine solve_lanczos_fake_complex(lc)

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
   !  call mp_sum(o_mat(:,:),world_comm)!this should be much faster than mp_bcast
     iproc=(io-1)/l_blk
     call mp_bcast(o_mat(:,:), iproc,world_comm)


  enddo

  deallocate(o_mat)
end subroutine solve_lanczos_fake_complex

