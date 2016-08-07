!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

!routines for the calculation of the self-energy
!lanczos-style
!ONLY FOR NORMCONSERVING PSEUDOS !!!!!

subroutine self_basis_lanczos(n_set,nstates,numpw, nsteps,ispin,lfull,nfull)
!this subroutine calculates the basis for every v 
!the minimal orthonormal basis for the Psi_i(r)*(v(r,r')*w^P'_i(r')) products

   USE io_global,            ONLY : stdout, ionode, ionode_id
   USE io_files,             ONLY : prefix, tmp_dir, diropn
   USE kinds,    ONLY : DP
   USE wannier_gw
   USE gvect
   USE constants, ONLY : e2, pi, tpi, fpi
   USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2
   USE wvfct,     ONLY : npwx, npw, nbnd
   USE gvecw,     ONLY : ecutwfc
   USE wavefunctions_module, ONLY : evc, psic
   USE mp, ONLY : mp_sum, mp_barrier, mp_bcast
   USE mp_world, ONLY : world_comm, mpime, nproc
   USE mp_pools, ONLY : intra_pool_comm
   USE gvecs,              ONLY : nls, nlsm, doublegrid
   USE fft_custom_gwl
   USE mp_wave, ONLY : mergewf,splitwf
   USE fft_base,             ONLY : dfftp, dffts
   USE fft_interfaces,       ONLY : fwfft, invfft
  
   implicit none

   INTEGER, EXTERNAL :: find_free_unit

   INTEGER, INTENT(in)  :: n_set  !defines the number of states to be read from disk at the same tim\e
   INTEGER, INTENT(in) :: nstates!number of orthonormal states to retain
   INTEGER, INTENT(in) :: numpw!dimension of polarization basis
   INTEGER, INTENT(in) :: nsteps!number of lanczos steps
   INTEGER, INTENT(in) :: ispin!spin channel 1,2
   LOGICAL, INTENT(in) :: lfull!if true prepare terms for full-relativistic calculation
   INTEGER, INTENT(in) :: nfull!in case of full-relativistic calculation the number of KS states (*1/2) to be treated in G
                               !as proper two dimensional spinors


   INTEGER :: iv,iw,ig,ii,jj
   REAL(kind=DP), ALLOCATABLE :: wv_real(:),tmp_r(:),tmp_r2(:)
   COMPLEX(kind=DP), ALLOCATABLE :: tmp_g(:), wp_prod(:,:,:)
   INTEGER :: iungprod,iunrprod, iungresult,iuntmat
   LOGICAL :: exst
   REAL(kind=DP), ALLOCATABLE :: omat(:,:),omat_hold(:,:)

   REAL(kind=DP), ALLOCATABLE :: eigen(:),work(:)
   INTEGER :: lwork,info,liwork
   COMPLEX(kind=DP), ALLOCATABLE :: wp_g(:,:)!product terms in g wfc grid
   REAL(kind=DP), ALLOCATABLE :: t_mat(:,:),t_mat_hold(:,:),t_mat_hold2(:,:)
   CHARACTER(4) :: nfile
   REAL(kind=DP), ALLOCATABLE :: fac(:)
   REAL(kind=DP) :: qq,exxdiv
   LOGICAL :: l_reduce_io=.true.!MUST BE THE SAME AS IN POLA_LANCZOS
   COMPLEX(kind=DP), ALLOCATABLE :: p_basis(:,:)!polarizability basis
   LOGICAL :: l_dsyevr=.true.!if true uses dsyevr
   REAL(kind=DP), ALLOCATABLE :: vectors(:,:)!for dsyevr
   INTEGER, ALLOCATABLE :: iwork(:), ifail(:)
   INTEGER, ALLOCATABLE :: isuppz(:)
   INTEGER :: n_found
   INTEGER, PARAMETER :: n_int=20
   REAL(kind=DP) :: qx,qy,qz
   INTEGER :: ix,iy,iz,n_int_loc
   REAL(kind=DP) :: qq_fact

   COMPLEX(kind=DP), ALLOCATABLE :: evc_t(:,:),p_basis_t(:,:)
   COMPLEX(kind=DP), ALLOCATABLE :: evc_g(:)
   COMPLEX(kind=DP), ALLOCATABLE :: wp_g_t(:,:)

   REAL(kind=DP), ALLOCATABLE :: p_basis_r(:,:)

   INTEGER :: ivv,nbuf
   INTEGER, PARAMETER :: offset=0!ATTENZIONE RESTART it is to used a reduced polarizability basis it should be 0

   INTEGER :: first_state!attenzione should be 1  RESTART
   INTEGER :: last_state!ATTENZIONE RESTART should be num_nbnds
   INTEGER :: il,iu
   REAL(kind=DP) :: vl,vu
   REAL(kind=DP), ALLOCATABLE :: t_eigen_hold(:)
   REAL(kind=DP), ALLOCATABLE :: fumat(:,:)


   TYPE(fft_cus) :: fc

   write(stdout,*) 'Routine self_basis_lanczos'
   FLUSH(stdout)


   if(s_first_state==0) then
      first_state=1
   else
      first_state=s_first_state
   endif

   if(s_last_state==0) then
      last_state=num_nbnds
   else
      last_state=s_last_state
   endif


!   first_state=1!1
!   last_state=num_nbnds

 !calculate V(G)                                                                                                             !
   
  if(.not.l_truncated_coulomb) then
      exxdiv=0.d0!exx_divergence_new()
   else
      exxdiv = 0.d0
   endif

   allocate(fac(max_ngm))
 

  if(l_truncated_coulomb) then
     do ig=1,max_ngm
    
         qq = g(1,ig)**2.d0 + g(2,ig)**2.d0 + g(3,ig)**2.d0
         
      
         if (qq > 1.d-8) then
            fac(ig)=(e2*fpi/(tpiba2*qq))*(1.d0-dcos(dsqrt(qq)*truncation_radius*tpiba))
         else
            fac(ig)=e2*fpi*(truncation_radius**2.d0/2.d0)
         endif
      enddo
     
  
      fac(:)=fac(:)/omega
   else
      fac(:)=0.d0
      fac(1:npw)=vg_q(1:npw)
   endif

   if(l_verbose) write(stdout,*) 'Call initialize_fft_custom'
   fc%ecutt=ecutwfc
   fc%dual_t=dual_vs
   FLUSH(stdout)
   call initialize_fft_custom(fc)
   !allocate(evc_g(fc%ngmt_g))

   

   allocate(wv_real(fc%nrxxt))

!read w^P'_i on file on real space
!open product of wanniers filed


   iungprod = find_free_unit()
   CALL diropn( iungprod, 'wiwjwfc_red', max_ngm*2, exst )

   if(.not.l_reduce_io) then
      iunrprod = find_free_unit()
      CALL diropn( iunrprod, 'wiwjwfc_red_r', dfftp%nnr, exst )
   else
      allocate(p_basis(max_ngm,numpw))
       do iw=1,numpw
          call davcio(p_basis(:,iw),max_ngm*2,iungprod,iw,-1)
      enddo
   endif

   iungresult = find_free_unit()
   if(ispin==1) then
      CALL diropn( iungresult, 'vw_lanczos_ss',npw*2, exst)
   else
      CALL diropn( iungresult, 'vw_lanczos_2ss',npw*2, exst)
   endif
   if(.not.l_reduce_io) then
      allocate(tmp_g(max_ngm),tmp_r(dfftp%nnr))
      do iw=1,numpw
         call davcio(tmp_g,max_ngm*2,iungprod,iw,-1)
         !trasform to r-space
         psic(:)=(0.d0,0.d0)
         do ig=1,max_ngm
            psic(nl(ig))=tmp_g(ig)*fac(ig)
            psic(nlm(ig))=CONJG(tmp_g(ig))*fac(ig)
         enddo
         CALL invfft ('Dense', psic, dfftp)
         tmp_r(:)=dble(psic(:))
         call davcio(tmp_r,dfftp%nnr,iunrprod,iw,1)
      enddo
      deallocate(tmp_g,tmp_r)
   endif
   close(iungprod)
   if(.not.l_reduce_io) close(iunrprod)



!now polarizability basis  times v are put on the ordering of the redueced grid, if required
  allocate(p_basis_t(fc%npwt,numpw))
  
  do ii=1,numpw
     p_basis(1:npw,ii)=fac(1:npw)*p_basis(1:npw,ii)
  enddo
  if(fc%dual_t==4.d0) then
     p_basis_t(:,:)=p_basis(:,:)
  else
     call reorderwfp_col(numpw,npw,fc%npwt,p_basis(1,1),p_basis_t(1,1), npw,fc%npwt, &
      & ig_l2g,fc%ig_l2gt,fc%ngmt_g,mpime, nproc,intra_pool_comm )
  endif
    
     !trasform to real space
  allocate(p_basis_r(fc%nrxxt,numpw))
  do ii=1,numpw,2
     psic(:)=(0.d0,0.d0)
     if(ii==numpw) then
        psic(fc%nlt(1:fc%npwt))  = p_basis_t(1:fc%npwt,ii)
        psic(fc%nltm(1:fc%npwt)) = CONJG( p_basis_t(1:fc%npwt,ii) )
     else
        psic(fc%nlt(1:fc%npwt))=p_basis_t(1:fc%npwt,ii)+(0.d0,1.d0)*p_basis_t(1:fc%npwt,ii+1)
        psic(fc%nltm(1:fc%npwt)) = CONJG( p_basis_t(1:fc%npwt,ii) )+(0.d0,1.d0)*CONJG( p_basis_t(1:fc%npwt,ii+1) )
     endif
     CALL cft3t( fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, 2 )
     p_basis_r(1:fc%nrxxt,ii)= DBLE(psic(1:fc%nrxxt))
     if(ii/=numpw) p_basis_r(1:fc%nrxxt,ii+1)= DIMAG(psic(1:fc%nrxxt))
     
  enddo

  

!now valence wavefunctions are put on the ordering of the reduced grid
  allocate(evc_t(fc%npwt,num_nbnds))
  if(fc%dual_t==4.d0) then
     evc_t(:,1:num_nbnds)=evc(:,1:num_nbnds)
  else
     call reorderwfp_col(num_nbnds,npw,fc%npwt,evc(1,1),evc_t(1,1), npw,fc%npwt, &
           & ig_l2g,fc%ig_l2gt,fc%ngmt_g,mpime, nproc,intra_pool_comm )
     !do iv=1,num_nbnds
     !   call mergewf(evc(:,iv),evc_g,npw,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
     !   call splitwf(evc_t(:,iv),evc_g,fc%npwt,fc%ig_l2gt,mpime,nproc,ionode_id,intra_pool_comm)
     !enddo
  endif






   if(l_verbose) write(stdout,*) 'self_basis_lanczos 1'
   FLUSH(stdout)
   nbuf=min(5,nproc) 
   allocate(wp_prod(fc%npwt,numpw,nbuf))

   allocate(tmp_r(fc%nrxxt),tmp_r2(fc%nrxxt))
   allocate(omat(numpw,numpw),omat_hold(numpw,numpw))
   allocate(t_mat(numpw,nstates), t_mat_hold(numpw,nstates), t_mat_hold2(numpw,nstates))
   allocate(wp_g(npw,nstates))
   allocate(wp_g_t(fc%npwt,nstates))
   allocate(t_eigen_hold(nstates))
   
!loop on kohn-sham states
   do ivv=first_state,last_state,nbuf
      call start_clock('sl_loop')
      do iv=ivv,min(ivv+nbuf-1,last_state)
!put iv on real space
      
         psic(:)=(0.d0,0.d0)
         psic(fc%nlt(1:fc%npwt))  = evc_t(1:fc%npwt,iv)
         psic(fc%nltm(1:fc%npwt)) = CONJG( evc_t(1:fc%npwt,iv) )
         CALL cft3t( fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, 2 )
         wv_real(1:fc%nrxxt)= DBLE(psic(1:fc%nrxxt))
         
!!loop on products of wanniers

         if(.not.l_reduce_io) then
            iunrprod = find_free_unit()
            CALL diropn( iunrprod, 'wiwjwfc_red_r', dfftp%nnr, exst )
         endif
         if(l_verbose) write(stdout,*) 'do fft'
         FLUSH(stdout)

         do ii=offset+1,numpw,2
!!read n_set w^P'_i from disk
            if(.not.l_reduce_io) then
               call davcio(tmp_r,dfftp%nnr,iunrprod,ii,-1)
            endif
            
            tmp_r(1:fc%nrxxt)=p_basis_r(1:fc%nrxxt,ii)*wv_real(1:fc%nrxxt)
            if(ii/=numpw) then
               tmp_r2(1:fc%nrxxt)=p_basis_r(1:fc%nrxxt,ii+1)*wv_real(1:fc%nrxxt)
            else
               tmp_r2(1:fc%nrxxt)=0.d0
            endif
!!form products with w_v and trasfrom in G space
            psic(1:fc%nrxxt)=dcmplx(tmp_r(1:fc%nrxxt),tmp_r2(1:fc%nrxxt))
            CALL cft3t( fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, -2 )
            if(ii==numpw) then
               wp_prod(1:fc%npwt, ii,iv-ivv+1) = psic(fc%nlt(1:fc%npwt))
            else
               wp_prod(1:fc%npwt, ii,iv-ivv+1)= 0.5d0*(psic(fc%nlt(1:fc%npwt))+conjg( psic(fc%nltm(1:fc%npwt))))
               wp_prod(1:fc%npwt, ii+1,iv-ivv+1)= (0.d0,-0.5d0)*(psic(fc%nlt(1:fc%npwt)) - conjg(psic(fc%nltm(1:fc%npwt))))
            endif
        

      
!!form products with w_v and trasfrom in G space
      
        

         enddo

!if required project off first nfull states
         if(lfull .and. nfull > 1 ) then
            allocate(fumat(nfull,numpw))
            call DGEMM('T','N',nfull,numpw,2*fc%npwt,2.d0,evc_t,2*fc%npwt,wp_prod(1,1,iv-ivv+1),2*fc%npwt,0.d0,fumat,nfull)
            if(fc%gstart_t==2) then
               do ii=1,nfull
                  do jj=1,numpw
                     fumat(ii,jj)=fumat(ii,jj)-dble(conjg(evc_t(1,ii))*wp_prod(1,jj,iv-ivv+1))
                  enddo
               enddo
            endif
            call mp_sum(fumat,world_comm)
            call DGEMM('N','N',2*fc%npwt,numpw,nfull,-1.d0,evc_t,2*fc%npwt,fumat,nfull,1.d0,wp_prod(1,1,iv-ivv+1),2*fc%npwt)

            deallocate(fumat)
         endif


         if(l_verbose) write(stdout,*) 'calculate omat'
         FLUSH(stdout)

         if(.not.l_reduce_io) close(iunrprod)
!!calculate overlap matrix
         call dgemm('T','N',numpw-offset,numpw-offset,2*fc%npwt,2.d0,wp_prod(1,offset+1,iv-ivv+1),&
     &2*fc%npwt,wp_prod(1,offset+1,iv-ivv+1),2*fc%npwt,0.d0,omat,numpw)
         if(fc%gstart_t==2) then
            do ii=1,numpw-offset
               do jj=1,numpw-offset
                  omat(jj,ii)=omat(jj,ii)-dble(conjg(wp_prod(1,offset+jj,iv-ivv+1))*wp_prod(1,offset+ii,iv-ivv+1))
               enddo
            enddo
         endif
         do ii=1,numpw-offset
            call mp_sum(omat(1:numpw-offset,ii),world_comm)
         enddo

         if(iv-ivv==mpime) then
            omat_hold(:,:)=omat(:,:)
         endif
      enddo

!!
!!solve eigen/values vector problem 
!!
      

     if(l_verbose) write(stdout,*) 'solve eig'
     FLUSH(stdout)
     do iv=ivv,min(ivv+nbuf-1,last_state)
        if(iv-ivv==mpime) then
           if(.not.l_dsyevr) then
              allocate(eigen(numpw-offset))
              allocate(work(1))
              call  DSYEV( 'V', 'U', numpw-offset, omat_hold, numpw, eigen, work, -1, info )
              lwork=work(1)
              deallocate(work)
              allocate(work(lwork))
              call  DSYEV( 'V', 'U', numpw-offset, omat_hold, numpw, eigen, work, lwork, info )
              deallocate(work)
              if(info/=0) then
                 write(stdout,*) 'ROUTINE self_basis_lanczos, INFO:', info
                 stop
              endif
         
         !do iw=1,numpw
         !   write(stdout,*) 'EIGEN:',iv,iw, eigen(iw)
         !enddo
         !FLUSH(stdout)
           else
              allocate(eigen(numpw-offset))
              allocate(vectors(numpw-offset,nstates))
              allocate(isuppz(2*nstates))
              allocate(work(1),iwork(1))
              omat_hold(1:numpw,1:numpw)=omat_hold(1:numpw,1:numpw)*1.d10!to cope with DSYEVR instabilities
              call DSYEVR('V','I','U',numpw-offset,omat_hold,numpw,0.d0,0.d0,&
         &numpw-offset-nstates+1,numpw-offset,0.d0,n_found,eigen,&
                   & vectors,numpw-offset,isuppz,work, -1,iwork,-1, info)
              lwork=work(1)
              liwork=iwork(1)
              deallocate(work,iwork)
              allocate(work(lwork))
              allocate(iwork(liwork))
              vl=0.d0
              vu=0.d0
              il=numpw-offset-nstates+1
              iu=numpw-offset
              call DSYEVR('V','I','U',numpw-offset,omat_hold,numpw,vl,vu,il,iu,0.d0,n_found,eigen,&
                   & vectors,numpw-offset,isuppz,work,lwork,iwork,liwork, info)
              eigen(1:numpw)=eigen(1:numpw)*1.d-10!to cope with DSYEVR instabilities
              if(info/=0) then
                 write(stdout,*) 'ROUTINE pola_lanczos DSYEVR, INFO:', info
                 stop
              endif
              deallocate(isuppz)
              deallocate(work,iwork)
              do iw=1,nstates, nstates-1
                 write(stdout,*) 'EIGEN S LOCAL:',iv,iw, eigen(iw)
              enddo
              FLUSH(stdout)

           endif
           t_mat_hold(:,:)=0.d0
           t_mat_hold2(:,:)=0.d0
           if(.not.l_dsyevr) then
              do ii=1,nstates
                 do jj=1,numpw-offset
                    t_mat_hold(jj+offset,ii)=omat_hold(jj,numpw-offset-ii+1)*(dsqrt(eigen(numpw-offset-ii+1)))
                 enddo
                 t_eigen_hold(ii)=eigen(numpw-ii+1)
              enddo
           else
              do ii=1,nstates
                 do jj=1,numpw-offset
                    t_mat_hold(jj+offset,ii)=vectors(jj,ii)*(dsqrt(eigen(ii)))
                 enddo
                 t_eigen_hold(ii)=eigen(ii)
              enddo
           endif
           if(.not.l_dsyevr) then
              do ii=1,nstates
                 t_mat_hold2(offset+1:numpw,ii)=omat_hold(1:numpw-offset,numpw-offset-ii+1)*(1.d0/dsqrt(eigen(numpw-offset-ii+1)))
              enddo
           else
              do ii=1,nstates
                 t_mat_hold2(offset+1:numpw,ii)=vectors(1:numpw-offset,ii)*(1.d0/dsqrt(eigen(ii)))
              enddo
           endif
           deallocate(eigen)
           if(l_dsyevr) deallocate(vectors)

        endif
     enddo
!!find liner dependent products

!!find transformation matrix and write on disk
!

      if(l_verbose) write(stdout,*) 'write on file'
      FLUSH(stdout)
      allocate(eigen(nstates))
      do iv=ivv,min(ivv+nbuf-1,last_state)
        if(iv-ivv == mpime) then
           t_mat(:,:)=t_mat_hold(:,:)
           eigen(1:nstates)=t_eigen_hold(1:nstates)
        endif
        call mp_bcast(t_mat,iv-ivv,world_comm)
        call mp_bcast(eigen(1:nstates),iv-ivv,world_comm)
        if(ionode) then
           iuntmat = find_free_unit()
           write(nfile,'(4i1)') iv/1000,mod(iv,1000)/100,mod(iv,100)/10,mod(iv,10)
           if(ispin==1) then
              open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.s_mat_lanczos'//nfile, status='unknown',form='unformatted')
           else
              open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.s_mat_lanczos2'//nfile, status='unknown',form='unformatted')
           endif
           write(iuntmat) iv
           write(iuntmat) num_nbndv(1)!for compatibility with polarization file
           write(iuntmat) numpw
           write(iuntmat) nstates
           do ii=1,nstates
              write(iuntmat) t_mat(1:numpw,ii)
           enddo
           close(iuntmat)
        endif

!write on disk file with eigen values                                                  
        if(ionode) then
           iuntmat = find_free_unit()
           write(nfile,'(4i1)') iv/1000,mod(iv,1000)/100,mod(iv,100)/10,mod(iv,10)
           if(ispin==1) then
              open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.s_eig_lanczos'//nfile, status='unknown',form='unformatted')
           else
              open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.s_eig_lanczos2'//nfile, status='unknown',form='unformatted')
           endif
           write(iuntmat) nstates
           write(iuntmat) eigen(1:nstates)
           close(iuntmat)
        endif



!!find liner dependent products
     
        if(iv-ivv == mpime) then
           t_mat(:,:)=t_mat_hold2(:,:)
        endif
        call mp_bcast(t_mat,iv-ivv,world_comm)

        call dgemm('N','N',2*fc%npwt,nstates,numpw,1.d0,wp_prod(1,1,iv-ivv+1),2*fc%npwt,t_mat,numpw,0.d0,wp_g_t,2*fc%npwt)


!put the correct order

        if(l_verbose) write(stdout,*) 'do merge split',iv,ivv
        FLUSH(stdout)

        if(fc%dual_t==4.d0) then
           wp_g(:,:)=wp_g_t(:,:)
        else
           call reorderwfp_col(nstates,fc%npwt,npw,wp_g_t(1,1),wp_g(1,1),fc%npwt,npw, &
                & fc%ig_l2gt,ig_l2g,fc%ngmt_g,mpime, nproc,intra_pool_comm )
           !do ii=1,nstates
           !   call mergewf(wp_g_t(:,ii),evc_g,fc%npwt,fc%ig_l2gt,mpime,nproc,ionode_id,intra_pool_comm)
           !   call splitwf(wp_g(:,ii),evc_g,npw,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
           !enddo
        endif

        if(l_verbose) write(stdout,*) 'do davcio'
        FLUSH(stdout)


  !!write on disk
        do ii=1,nstates
           call davcio(wp_g(:,ii),npw*2,iungresult,ii+(iv-first_state)*nstates,1)
        enddo

     

     enddo
     deallocate(eigen)
     call stop_clock('sl_loop')
  enddo
  deallocate(tmp_r,tmp_r2,omat,omat_hold,p_basis_r)
  deallocate(t_mat,t_mat_hold,t_mat_hold2)
  deallocate(wp_g,wp_g_t)
  close(iungresult)
  deallocate(wv_real,wp_prod)
  if(l_verbose) write(stdout,*) 'Exiting self_basis_lanczos'
  FLUSH(stdout)
  deallocate(p_basis_t)
  if(l_verbose) write(stdout,*) 'Call deallocate_fft_custom'
  FLUSH(stdout)
  call deallocate_fft_custom(fc)
  deallocate(t_eigen_hold)
  if(l_reduce_io) deallocate(p_basis)
  return
 end subroutine self_basis_lanczos


subroutine global_self_lanczos(nstates,nstates_eff,threshold,nglobal,nsteps,numpw,g_threshold,ispin,l_eigen,istate,lfull)
!this subroutine from the orthonormal basis at each ks-state i 
!construct a global basis for the lanczos calculation of the
!self-energy
  USE io_global,            ONLY : stdout, ionode, ionode_id
  USE io_files,             ONLY : prefix, tmp_dir, diropn
  USE kinds,    ONLY : DP
  USE wannier_gw, ONLY : num_nbnds,max_ngm,l_truncated_coulomb,truncation_radius,&
                               &num_nbndv,l_pmatrix,vg_q,s_first_state,s_last_state, l_verbose,&
                               &l_contour,l_big_system,l_list,n_list,i_list,optimal_options
  USE gvect
  USE wvfct,    ONLY : npwx, npw, nbnd
  USE mp, ONLY : mp_sum, mp_barrier, mp_bcast
  USE mp_world, ONLY : world_comm, mpime,nproc
  USE wavefunctions_module, ONLY : evc, psic
  USE gvect
  USE gvecs,              ONLY : nls, nlsm, doublegrid
  USE klist,    ONLY : igk_k
  USE constants, ONLY : e2, pi, tpi, fpi
  USE cell_base, ONLY : at, alat, tpiba, omega, tpiba2
  USE wvfct,     ONLY : et 
  USE fft_base,             ONLY : dfftp, dffts
  USE fft_interfaces,       ONLY : fwfft, invfft
  
  implicit none
  
  INTEGER, EXTERNAL :: find_free_unit

  INTEGER, INTENT(in)      :: nstates!number of orthonormal states for each v
  INTEGER, INTENT(in)      :: nstates_eff!number of orthonormal states for each v
  REAL(kind=DP),INTENT(in) :: threshold!threshold for orthonormalization algorithm
  INTEGER, INTENT(out)     :: nglobal!total number of final orthonormal states
  INTEGER, INTENT(in)      :: nsteps!number of lanczos steps                 
  INTEGER, INTENT(in)      :: numpw!number of wannier products for testing
  REAL(kind=DP), INTENT(in)      :: g_threshold!threshold for eigen values of trial green function
  INTEGER, INTENT(in)      :: ispin!spin channel 1,2
  LOGICAL, INTENT(in)      :: l_eigen!if true partial t states are scaled with the corresponding eigenvalue
  INTEGER, INTENT(in)      :: istate!considered state for l_big_systems case
  LOGICAL, INTENT(in)      :: lfull!if true writes on disk global s vectors on charge R grid for further processing


  INTEGER :: iunv,iuntmat
  LOGICAL :: exst
  INTEGER :: ii,jj,iv,ic,is
  COMPLEX(kind=DP), ALLOCATABLE :: old_basis(:,:), new_basis(:,:),v_basis(:,:)
  INTEGER :: nglobal_old
  REAL(kind=DP), ALLOCATABLE :: t_mat(:,:)
  CHARACTER(4) :: nfile
!for test:
  REAL(kind=DP) :: sca,sca1
  INTEGER :: iungprod,ig,iw
  REAL(kind=DP), ALLOCATABLE :: wv_real(:),tmp_r(:)
  COMPLEX(kind=DP), ALLOCATABLE :: tmp_g(:),wp_prod(:)
  REAL(kind=DP), ALLOCATABLE :: fac(:)
  REAL(kind=DP) :: qq,exxdiv

!For trial Green's function   
  LOGICAL :: l_greent=.false.!if true calculate the trial green function
  REAL(kind=DP), ALLOCATABLE :: o_t_psi(:,:)
  REAL(kind=DP), ALLOCATABLE :: gtrail(:,:)
  REAL(kind=DP) :: offset
  REAL(kind=DP), ALLOCATABLE :: eigen(:),work(:)
  INTEGER :: lwork,info
  LOGICAL :: l_test=.false.!ATTENZIONE
  REAL(kind=DP) :: proj_tot

  INTEGER :: nbuffer,ndelta!for avoiding nested allocation/deallocation cycles
  LOGICAL :: l_update_memory
  INTEGER, PARAMETER :: offset_s=0!RESTART ATTENZIONE THEN PUT 0!!!!!!!
  INTEGER :: iuns
  INTEGER :: ip
  LOGICAL, PARAMETER :: l_restart = .false. !if true do a restart RESTART ATTENZIONE
  INTEGER :: first_state, last_state! ATTENZIONE should be 1 and num_nbnd

  REAL(kind=DP), ALLOCATABLE :: s_eigen(:)
  INTEGER :: idumm

  INTEGER :: iunc
  TYPE(optimal_options) :: options
  LOGICAL, PARAMETER :: l_reortho=.true.
 
  nbuffer=6*numpw
  ndelta=numpw




  if(.not.l_big_system) then
     if(s_first_state==0) then
        first_state=1
     else
        first_state=s_first_state
     endif
     
     if(s_last_state==0) then
        last_state=num_nbnds
     else
        last_state=s_last_state
     endif
  else
     first_state=istate
     last_state=istate
  endif



 ! first_state=1
 ! last_state=num_nbnds


 !calculate V(G)             
  
  allocate(s_eigen(nstates))
  if(l_eigen) then
     if(ionode) then
        iv=first_state
        iuntmat = find_free_unit()
        write(nfile,'(4i1)') iv/1000,mod(iv,1000)/100,mod(iv,100)/10,mod(iv,10)
        if(ispin==1) then
           open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.s_eig_lanczos'//nfile, status='old',form='unformatted')
        else
           open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.s_eig_lanczos2'//nfile, status='old',form='unformatted')
        endif
        read(iuntmat) idumm
        read(iuntmat) s_eigen(1:nstates)
        close(iuntmat)
     endif
     call mp_bcast(s_eigen, ionode_id,world_comm)
  endif
  

  allocate(fac(max_ngm))

   if(l_truncated_coulomb) then
      do ig=1,max_ngm
         qq = g(1,ig)**2.d0 + g(2,ig)**2.d0 + g(3,ig)**2.d0
         if (qq > 1.d-8) then
            fac(ig)=(e2*fpi/(tpiba2*qq))*(1.d0-dcos(dsqrt(qq)*truncation_radius*tpiba))
         else
            fac(ig)=e2*fpi*(truncation_radius**2.d0/2.d0)
         endif
      end do
      fac(:)=fac(:)/omega
   else
      fac(:)=0.d0
      fac(1:npw)=vg_q(1:npw)
   endif

   if(.not.l_restart) then

!set first basis from first valence state
  allocate(old_basis(npw,nbuffer))
  iunv = find_free_unit()
  if(ispin==1) then
     CALL diropn( iunv, 'vw_lanczos_ss',npw*2, exst)
  else
     CALL diropn( iunv, 'vw_lanczos_2ss',npw*2, exst)
  endif

  if(.not.l_eigen) then
     do ii=1,nstates_eff
        if(.not.l_big_system) then
           !call davcio(old_basis(:,ii),npw*2,iunv,ii+offset_s+(first_state-1)*(nstates+offset_s),-1)
           call davcio(old_basis(:,ii),npw*2,iunv,ii+offset_s,-1)
        else
           call davcio(old_basis(:,ii),npw*2,iunv,ii+offset_s+(istate-s_first_state)*(nstates+offset_s),-1)
        endif
     enddo
     nglobal=nstates_eff
  else
     nglobal=1
     if(.not.l_big_system) then
        !call davcio(old_basis(:,nglobal),npw*2,iunv,1+offset_s+(first_state-1)*(nstates+offset_s),-1)
        call davcio(old_basis(:,nglobal),npw*2,iunv,1+offset_s,-1)
     else
        call davcio(old_basis(:,nglobal),npw*2,iunv,1+offset_s+(istate-s_first_state)*(nstates+offset_s),-1)
     endif
     do ii=2,nstates_eff
        if(s_eigen(ii) > threshold) then
           nglobal=nglobal+1
           if(.not.l_big_system) then
              !call davcio(old_basis(:,nglobal),npw*2,iunv,ii+offset_s+(first_state-1)*(nstates+offset_s),-1)
              call davcio(old_basis(:,nglobal),npw*2,iunv,ii+offset_s,-1)
           else
              call davcio(old_basis(:,nglobal),npw*2,iunv,ii+offset_s+(istate-s_first_state)*(nstates+offset_s),-1)
           endif
        endif
     enddo
  endif

!loop on valence states (from 2nd)
  allocate(v_basis(npw,nstates_eff))
  allocate(new_basis(npw,nbuffer))
  do iv=first_state+1,last_state

!!read from disk
     do ii=1,nstates_eff
        if(.not.l_big_system) then
           !call davcio(v_basis(:,ii),npw*2,iunv,ii+offset_s+(iv-1)*(nstates+offset_s),-1)
           call davcio(v_basis(:,ii),npw*2,iunv,ii+offset_s+(iv-first_state)*(nstates+offset_s),-1)
        else
           call davcio(v_basis(:,ii),npw*2,iunv,ii+offset_s+(iv-s_first_state)*(nstates+offset_s),-1)
        endif
     enddo

     if(l_eigen) then
        if(ionode) then
           iuntmat = find_free_unit()
           write(nfile,'(4i1)') iv/1000,mod(iv,1000)/100,mod(iv,100)/10,mod(iv,10)
           if(ispin==1) then
              open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.s_eig_lanczos'//nfile, status='old',form='unformatted')
           else
              open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.s_eig_lanczos2'//nfile, status='old',form='unformatted')
           endif
           read(iuntmat) idumm
           read(iuntmat) s_eigen(1:nstates)
           close(iuntmat)
        endif
        call mp_bcast(s_eigen, ionode_id,world_comm)
     endif



     if(nglobal+nstates_eff >nbuffer) then
        deallocate(new_basis)
        allocate(new_basis(npw,nbuffer+ndelta))
        l_update_memory=.true.
     else
        l_update_memory=.false.
     endif

!     allocate(new_basis(npw,nglobal+nstates_eff))

!!calculate basis
     nglobal_old=nglobal
     if(.not.l_pmatrix) then
        !call orthonormalize_two_manifolds( old_basis, nglobal_old,v_basis, nstates, threshold, new_basis, nglobal)
        call orthonormalize_two_manifolds_prj( old_basis, nglobal_old,v_basis, nstates_eff, threshold, new_basis, nglobal,&
             &l_eigen,s_eigen)
     else
        call orthonormalize_two_manifolds_scalapack(old_basis, nglobal_old,v_basis, nstates_eff, threshold, new_basis, nglobal)
     endif

!!set arrays for next iteration
     if(l_update_memory) then
        deallocate(old_basis)
        allocate(old_basis(npw,nbuffer+ndelta))
        nbuffer=nbuffer+ndelta
     endif
     old_basis(:,1:nglobal)=new_basis(:,1:nglobal)

     !deallocate(old_basis)
     !allocate(old_basis(npw,nglobal))
     !old_basis(:,1:nglobal)=new_basis(:,1:nglobal)
     !deallocate(new_basis)

  enddo
  deallocate(new_basis)
  if(l_greent) then
!NOT_TO_BE_INCLUDED_START
!if required turn the basis to a more diagonal form
!a)calculate overlaps <t_i|\Psi_'>
     if(l_verbose) write(stdout,*) 'GREENT 1'
     FLUSH(stdout)
     allocate(o_t_psi(nglobal,nbnd))
     call dgemm('T','N',nglobal,nbnd,2*npw,2.d0,old_basis,2*npw,evc,2*npwx,0.d0,o_t_psi,nglobal)
     if(gstart==2) then
         do ii=1,nbnd
            do jj=1,nglobal
               o_t_psi(jj,ii)=o_t_psi(jj,ii)-dble(conjg(old_basis(1,jj))*evc(1,ii))
            enddo
         enddo
      endif
      call mp_sum(o_t_psi(:,:),world_comm)
      offset=(et(num_nbndv(1)+1,1)-et(num_nbndv(1),1))/2.d0
      if(l_verbose) write(stdout,*) 'GREENT 2'
     FLUSH(stdout)

!b)calculate matrix G^T_ij=<t_i|Psi_i'><Psi_i'|t_j>/(E_i'+offset)
      allocate(gtrail(nglobal,nglobal))
!mettere DGEMM!!
      gtrail(:,:)=0.d0
      do jj=1,nglobal
         do ii=1,nglobal
            do is=1,nbnd
               gtrail(ii,jj)=gtrail(ii,jj)+o_t_psi(ii,is)*o_t_psi(jj,is)/(et(is,1)-offset)
            enddo
         enddo
      enddo
      if(l_verbose) write(stdout,*) 'GREENT 3'
     FLUSH(stdout)

     deallocate(o_t_psi)
!c)diagonalize
     allocate(eigen(nglobal))
     if(ionode) then
        allocate(work(1))
        call  DSYEV( 'V', 'U', nglobal, gtrail, nglobal, eigen, work, -1, info )
        lwork=work(1)
        deallocate(work)
        allocate(work(lwork))
        call  DSYEV( 'V', 'U', nglobal, gtrail, nglobal, eigen, work, lwork, info )
        deallocate(work)
        if(info/=0) then
           write(stdout,*) 'ROUTINE self_lanczos, INFO:', info
           stop
        endif
     else
        gtrail(:,:)=0.d0
        eigen(:)=0.d0
     endif
     if(l_verbose) write(stdout,*) 'GREENT 4'
     FLUSH(stdout)
     
     do ii=1,nglobal
        !call mp_bcast(gtrail(:,ii), ionode_id,world_comm)
        call mp_sum(gtrail(:,ii),world_comm)
     enddo
     !call mp_bcast(eigen(:), ionode_id,world_comm)
     call mp_sum(eigen(:),world_comm)
     
     do ii=1,nglobal
        if(l_verbose) write(stdout,*) 'EIGEN GTRAIL:',ii, eigen(ii)
     enddo
     FLUSH(stdout)

!d)calculate t_i' eigen states of G^T
   allocate(new_basis(npw,nglobal))
   call dgemm('N','N',npw*2,nglobal,nglobal,1.d0,old_basis,2*npw,gtrail,nglobal,0.d0,new_basis,2*npw)

!e)take eigen states corresponding to eigevalues large than threshold
   nglobal_old=nglobal
   nglobal=0
   do ii=1,nglobal_old
      if(abs(eigen(ii)) >= g_threshold) then
         nglobal=nglobal+1
         old_basis(:,nglobal)=new_basis(:,ii)
      endif
   enddo
 !  old_basis(:,:)=new_basis(:,:)
   write(stdout,*) 'NUMBER T STATES:',nglobal
   FLUSH(stdout)
     
   
   deallocate(eigen)
   deallocate(gtrail,new_basis)
!NOT_TO_BE_INCLUDED_END
endif

  write(stdout,*) 'TOTAL NUMBER OF GLOBAL S VECTORS: ', nglobal

!re-orthonormalized them it should be important in large systems
  if(l_reortho.and.  .not.l_big_system) then
     if(l_verbose) write(stdout,*) 'Call optimal driver'
     FLUSH(stdout)
     options%l_complete=.true.
     options%idiago=0
     call optimal_driver(nglobal,old_basis,npw,options,idumm, info)
  endif



!call lanczos chain routine
  if(.not.l_big_system) then
    call lanczos_state(old_basis, nglobal, 1, nsteps,1,ispin)
 else
    call lanczos_state(old_basis, nglobal, 1, nsteps,istate-s_first_state+1,ispin)
 endif


!if required calculate overlaps with KS states
!NOT_TO_BE_INCLUDED_START
 if(l_contour) then
    if(.not.l_big_system) then
       call  contour_terms(nglobal,old_basis,ispin,1)
    else
       call  contour_terms(nglobal,old_basis,ispin,istate)
    endif
 endif
!NOT_TO_BE_INCLUDED_END


  !do ip=0,nproc-1
  !   if(mpime==ip) then
        iuns = find_free_unit()
        CALL diropn( iuns, 's_vectors',npw*2, exst)
        do ii=1,nglobal
           call davcio(old_basis(:,ii),npw*2,iuns,ii,1)
        enddo
        close(iuns)
  !   endif
  !   call mp_barrier
 ! enddo
     
  else
     if(l_verbose) write(stdout,*) 'Doing restart iuns'
     FLUSH(stdout)

     nglobal=721!RESTART here put te value from the ouputfile

     write(stdout,*) 'Total number of s vectors:', nglobal
     FLUSH(stdout)

     allocate(old_basis(npw,nglobal))

     iunv = find_free_unit()
     CALL diropn( iunv, 'vw_lanczos_ss',npw*2, exst)
     allocate(v_basis(npw,nstates_eff))

     do ip=0,nproc-1
        if(mpime==ip) then
           iuns = find_free_unit()
           CALL diropn( iuns, 's_vectors',npw*2, exst)
           do ii=1,nglobal
              call davcio(old_basis(:,ii),npw*2,iuns,ii,-1)
           enddo
           close(iuns)
        endif
        call mp_barrier( world_comm )
     enddo



  endif
  
  !if required writes global s vectors on disk
  if(lfull) then
       iunc = find_free_unit()
       if(ispin==1) then
          CALL diropn( iunc, 'sreal2full',dffts%nnr, exst)
       else
          CALL diropn( iunc, 'sreal2full2_',dffts%nnr, exst)
       endif
       do ii=1,nglobal,2
          psic(1:dffts%nnr)=(0.d0,0.d0)
          do ig=1,npw
             if(ii<nglobal) then
                psic(nls(ig))=old_basis(ig,ii)+(0.d0,1.d0)*old_basis(ig,ii+1)
                psic(nlsm(ig))=conjg(old_basis(ig,ii))+(0.d0,1.d0)*conjg(old_basis(ig,ii+1))
             else
                psic(nls(ig))=old_basis(ig,ii)
                psic(nlsm(ig))=conjg(old_basis(ig,ii))
             endif
          enddo
          CALL invfft ('Wave', psic, dffts)
          call davcio(dble(psic(1:dffts%nnr)),dffts%nnr,iunc,ii,1)
          if(ii<nglobal) call davcio(dimag(psic(1:dffts%nnr)),dffts%nnr,iunc,ii+1,1)
       enddo
       close(iunc)

  endif




!calculate matrix element and write on disk  
  allocate(t_mat(nglobal,nstates_eff))
  do iv=first_state,last_state!RESTART ATTENZIONE it was 1!!!!!
     write(nfile,'(4i1)') iv/1000,mod(iv,1000)/100,mod(iv,100)/10,mod(iv,10)
     if(ionode) then
        iuntmat=find_free_unit()
        if(ispin==1) then
           open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.st_mat_lanczos'//nfile, status='unknown',form='unformatted')
        else
           open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.st_mat_lanczos2'//nfile, status='unknown',form='unformatted')
        endif
     endif
!     do ip=0,nproc-1,3

           do ii=1,nstates_eff
              if(.not.l_big_system) then
                 !call davcio(v_basis(:,ii),npw*2,iunv,ii+offset_s+(iv-1)*(nstates+offset_s),-1)
                 call davcio(v_basis(:,ii),npw*2,iunv,ii+offset_s+(iv-first_state)*(nstates+offset_s),-1)
              else
                 call davcio(v_basis(:,ii),npw*2,iunv,ii+offset_s+(iv-s_first_state)*(nstates+offset_s),-1)
              endif
!              call davcio(v_basis(:,ii),npw*2,iunv,ii+offset_s+(iv-first_state)*(nstates+offset_s),-1)
           enddo

 !    enddo
     call dgemm('T','N',nglobal,nstates_eff,2*npw,2.d0,old_basis,2*npw,v_basis,2*npw,0.d0,t_mat,nglobal)
     if(gstart==2) then
        do ii=1,nstates_eff
           do jj=1,nglobal
              t_mat(jj,ii)=t_mat(jj,ii)-dble(conjg(old_basis(1,jj))*v_basis(1,ii))
           enddo
        enddo
     endif
     call mp_sum(t_mat(:,:),world_comm)
     if(ionode) then
        write(iuntmat) nglobal
        write(iuntmat) nstates_eff
        write(iuntmat) iv
        do ii=1,nstates_eff
           write(iuntmat) t_mat(1:nglobal,ii)
        enddo
        close(iuntmat)
     endif
  enddo

  close(iunv)
!THE FOLLOWING PART IS FOR TESTING POURPOSES

!test that the basis {t_i} is orthonormal
  deallocate(t_mat)

  if(l_test) then
     write(stdout,*) 'TEST1'
     FLUSH(stdout)
     allocate(t_mat(nglobal,nglobal))
     write(stdout,*) 'TEST2'
     FLUSH(stdout)

     call dgemm('T','N',nglobal,nglobal,2*npw,2.d0,old_basis,2*npw,old_basis,2*npw,0.d0,t_mat,nglobal)
     if(gstart==2) then
        do ii=1,nglobal
           do jj=1,nglobal
              t_mat(jj,ii)=t_mat(jj,ii)-dble(conjg(old_basis(1,jj))*old_basis(1,ii))
           enddo
        enddo
     endif
     write(stdout,*) 'TEST3'
     FLUSH(stdout)

     call mp_sum(t_mat(:,:),world_comm)
!!write diagonal terms
     do ii=1,nglobal
        sca=0.d0
        do jj=1,nglobal
           if(ii/=jj) sca=sca+abs(t_mat(ii,jj))
        enddo
        write(stdout,*) 'Diagonal',ii,t_mat(ii,ii),sca
        FLUSH(stdout)
     enddo
     deallocate(t_mat)
     allocate(t_mat(nglobal,1))
!test for representability
     iungprod = find_free_unit()
     CALL diropn( iungprod, 'wiwjwfc_red', max_ngm*2, exst )
     allocate(wv_real(dfftp%nnr),tmp_r(dfftp%nnr),tmp_g(ngm),wp_prod(npw))
     proj_tot=0.d0
     do iv=1,num_nbnds-1
!put iv on real space
        psic(:)=(0.d0,0.d0)
        psic(nls (igk_k(1:npw,1)))  = evc(1:npw,iv)
        psic(nlsm(igk_k(1:npw,1))) = CONJG( evc(1:npw,iv) )
        CALL invfft ('Wave', psic, dffts)
        wv_real(:)= DBLE(psic(:))
        
      
!loop on wannier_products
        do iw=1,numpw
           call davcio(tmp_g,max_ngm*2,iungprod,iw,-1)
           !trasform to r-space
           psic(:)=(0.d0,0.d0)
           do ig=1,max_ngm
              psic(nl(ig))=tmp_g(ig)*fac(ig)
              psic(nlm(ig))=CONJG(tmp_g(ig))*fac(ig)
           enddo
           CALL invfft ('Dense', psic, dfftp)
           tmp_r(:)=dble(psic(:))
!!form products with w_v and trasfrom in G space
           psic(:)=cmplx(tmp_r(:)*wv_real(:),0.d0)
           CALL fwfft ('Wave', psic, dffts)
           wp_prod(1:npw) = psic(nls(igk_k(1:npw,1)))

!!do scalar product
           call dgemm('T','N',nglobal,1,2*npw,2.d0,old_basis,2*npw,wp_prod,2*npw,0.d0,t_mat,nglobal)
           if(gstart==2) then
              do jj=1,nglobal
                 t_mat(jj,1)=t_mat(jj,1)-dble(conjg(old_basis(1,jj))*wp_prod(1))
              enddo
           endif
           call mp_sum(t_mat(:,1),world_comm)
           sca=0.d0
           do ii=1,nglobal
              sca=sca+t_mat(ii,1)**2.d0
           end do
!calculate norm
           sca1=0.d0
           do ig=1,npw
              sca1=sca1+2.d0*dble(conjg(wp_prod(ig))*wp_prod(ig))
           enddo
           if(gstart==2) sca1=sca1-dble(conjg(wp_prod(1))*wp_prod(1))
           call mp_sum(sca1,world_comm)
         
           write(stdout,*) 'Projection',iv,iw,sca/sca1,sca1
           proj_tot=proj_tot+sca/sca1
!         do ii=1,nglobal,50
!            write(stdout,*) 'Q terms',iv,iw,ii, t_mat(ii,1)
!         enddo
           FLUSH(stdout)
        enddo
     enddo
     write(stdout,*) 'Average projection', proj_tot/dble(numpw*2)
     FLUSH(stdout)
     deallocate(t_mat)
     deallocate(wv_real,tmp_g,tmp_r,wp_prod)
     close(iungprod)
     deallocate(fac)
!END OF TESTING PART
  endif

  deallocate(old_basis)
  deallocate(v_basis)
  deallocate(s_eigen)

  return
end subroutine global_self_lanczos


subroutine self_basis_lanczos_real(n_set,nstates,numpw, nsteps,ispin)
!NOT_TO_BE_INCLUDED_START
!this subroutine calculates the basis for every v 
!the minimal orthonormal basis for the Psi_i(r)*(v(r,r')*w^P'_i(r')) products

   USE io_global,            ONLY : stdout, ionode, ionode_id
   USE io_files,             ONLY : prefix, tmp_dir, diropn
   USE kinds,    ONLY : DP
   USE wannier_gw
   USE gvect
   USE constants, ONLY : e2, pi, tpi, fpi
   USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2
   USE wvfct,     ONLY : npwx, npw, nbnd
   USE gvecw,     ONLY : ecutwfc
   USE wavefunctions_module, ONLY : evc, psic
   USE mp, ONLY : mp_sum, mp_barrier, mp_bcast
   USE mp_world, ONLY : world_comm, mpime, nproc
   USE mp_pools, ONLY : intra_pool_comm
   USE gvecs,              ONLY : nls, nlsm, doublegrid
   !USE exx, ONLY : exx_divergence_new, yukawa
   USE fft_custom_gwl
   USE mp_wave, ONLY : mergewf,splitwf
   USE fft_base,             ONLY : dfftp, dffts
   USE fft_interfaces,       ONLY : fwfft, invfft
  
   implicit none

   INTEGER, EXTERNAL :: find_free_unit

   INTEGER, INTENT(in)  :: n_set  !defines the number of states to be read from disk at the same tim\e
   INTEGER, INTENT(in) :: nstates!number of orthonormal states to retain
   INTEGER, INTENT(in) :: numpw!dimension of polarization basis
   INTEGER, INTENT(in) :: nsteps!number of lanczos steps
   INTEGER, INTENT(in) :: ispin!spin channel 1,2


   INTEGER :: iv,iw,ig,ii,jj
   REAL(kind=DP), ALLOCATABLE :: wv_real(:),tmp_r(:),tmp_r2(:)
   COMPLEX(kind=DP), ALLOCATABLE :: tmp_g(:)
   REAL(kind=DP), ALLOCATABLE :: wp_prod(:,:,:)
   INTEGER :: iungprod,iunrprod, iungresult,iuntmat
   LOGICAL :: exst
   REAL(kind=DP), ALLOCATABLE :: omat(:,:),omat_hold(:,:)

   REAL(kind=DP), ALLOCATABLE :: eigen(:),work(:)
   INTEGER :: lwork,info,liwork
   COMPLEX(kind=DP), ALLOCATABLE :: wp_g(:,:)!product terms in g wfc grid
   REAL(kind=DP), ALLOCATABLE :: t_mat(:,:),t_mat_hold(:,:),t_mat_hold2(:,:)
   CHARACTER(4) :: nfile
   REAL(kind=DP), ALLOCATABLE :: fac(:)
   REAL(kind=DP) :: qq,exxdiv
    COMPLEX(kind=DP), ALLOCATABLE :: p_basis(:,:)!polarizability basis
   LOGICAL :: l_dsyevr=.true.!if true uses dsyevr
   REAL(kind=DP), ALLOCATABLE :: vectors(:,:)!for dsyevr
   INTEGER, ALLOCATABLE :: iwork(:), ifail(:)
   INTEGER, ALLOCATABLE :: isuppz(:)
   INTEGER :: n_found
   INTEGER, PARAMETER :: n_int=20
   REAL(kind=DP) :: qx,qy,qz
   INTEGER :: ix,iy,iz,n_int_loc
   REAL(kind=DP) :: qq_fact

   COMPLEX(kind=DP), ALLOCATABLE :: evc_t(:,:),p_basis_t(:,:)
   COMPLEX(kind=DP), ALLOCATABLE :: evc_g(:), wp_g_t2(:,:),wp_g_t3(:)
   REAL(kind=DP), ALLOCATABLE :: wp_g_t(:,:)
   
   REAL(kind=DP), ALLOCATABLE :: p_basis_r(:,:)

   INTEGER :: ivv,nbuf
   INTEGER, PARAMETER :: offset=0!ATTENZIONE RESTART it is to used a reduced polarizability basis it should be 0

   INTEGER :: first_state!attenzione should be 1  RESTART
   INTEGER :: last_state!ATTENZIONE RESTART should be num_nbnds
   INTEGER :: il,iu
   REAL(kind=DP) :: vl,vu
   REAL(kind=DP), ALLOCATABLE :: t_eigen_hold(:)

   TYPE(fft_cus) :: fc

   write(stdout,*) 'Routine self_basis_lanczos'
   FLUSH(stdout)


   if(s_first_state==0) then
      first_state=1
   else
      first_state=s_first_state
   endif

   if(s_last_state==0) then
      last_state=num_nbnds
   else
      last_state=s_last_state
   endif


!   first_state=1!1
!   last_state=num_nbnds

 !calculate V(G)                                                                                                             !
   
  if(.not.l_truncated_coulomb) then
      exxdiv=0.d0!exx_divergence_new()
   else
      exxdiv = 0.d0
   endif

   allocate(fac(max_ngm))
 

  if(l_truncated_coulomb) then
     do ig=1,max_ngm
    
         qq = g(1,ig)**2.d0 + g(2,ig)**2.d0 + g(3,ig)**2.d0
         
      
         if (qq > 1.d-8) then
            fac(ig)=(e2*fpi/(tpiba2*qq))*(1.d0-dcos(dsqrt(qq)*truncation_radius*tpiba))
         else
            fac(ig)=e2*fpi*(truncation_radius**2.d0/2.d0)
         endif
      enddo
     
  
      fac(:)=fac(:)/omega
   else
      fac(:)=0.d0
      fac(1:npw)=vg_q(1:npw)
   endif

   if(l_verbose) write(stdout,*) 'Call initialize_fft_custom'
   fc%ecutt=ecutwfc
   fc%dual_t=dual_vs
   FLUSH(stdout)
   call initialize_fft_custom(fc)
   allocate(evc_g(fc%ngmt_g*nproc))


   allocate(wv_real(fc%nrxxt))

!read w^P'_i on file on real space
!open product of wanniers filed


   iungprod = find_free_unit()
   CALL diropn( iungprod, 'wiwjwfc_red', max_ngm*2, exst )


   allocate(p_basis(max_ngm,numpw))
   do iw=1,numpw
      call davcio(p_basis(:,iw),max_ngm*2,iungprod,iw,-1)
   enddo
   

   iungresult = find_free_unit()
   if(ispin==1) then
      CALL diropn( iungresult, 'vw_lanczos_ss',npw*2, exst)
   else
      CALL diropn( iungresult, 'vw_lanczos_2ss',npw*2, exst)
   endif

   close(iungprod)
 



!now polarizability basis  times v are put on the ordering of the redueced grid, if required
  allocate(p_basis_t(fc%npwt,numpw))
  
  do ii=1,numpw
     p_basis(1:npw,ii)=fac(1:npw)*p_basis(1:npw,ii)
  enddo
  if(fc%dual_t==4.d0) then
     p_basis_t(:,:)=p_basis(:,:)
  else
     call reorderwfp (numpw,npw, fc%npwt,p_basis(:,:),p_basis_t(:,:), &
          &npw,fc%npwt, ig_l2g,fc%ig_l2gt, fc%ngmt_g , mpime, nproc,ionode_id, intra_pool_comm )

     !do ii=1,numpw
     !   call mergewf(p_basis(:,ii),evc_g,npw,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
     !   call splitwf(p_basis_t(:,ii),evc_g,fc%npwt,fc%ig_l2gt,mpime,nproc,ionode_id,intra_pool_comm)
     !enddo
     !trasform to real space
     allocate(p_basis_r(fc%nrxxt,numpw))
     do ii=1,numpw,2
        psic(:)=(0.d0,0.d0)
        if(ii==numpw) then
           psic(fc%nlt(1:fc%npwt))  = p_basis_t(1:fc%npwt,ii)
           psic(fc%nltm(1:fc%npwt)) = CONJG( p_basis_t(1:fc%npwt,ii) )
        else
           psic(fc%nlt(1:fc%npwt))=p_basis_t(1:fc%npwt,ii)+(0.d0,1.d0)*p_basis_t(1:fc%npwt,ii+1)
           psic(fc%nltm(1:fc%npwt)) = CONJG( p_basis_t(1:fc%npwt,ii) )+(0.d0,1.d0)*CONJG( p_basis_t(1:fc%npwt,ii+1) )
        endif
        CALL cft3t( fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, 2 )
        p_basis_r(1:fc%nrxxt,ii)= DBLE(psic(1:fc%nrxxt))
        if(ii/=numpw) p_basis_r(1:fc%nrxxt,ii+1)= DIMAG(psic(1:fc%nrxxt))

     enddo

  endif

!now valence wavefunctions are put on the ordering of the reduced grid
  allocate(evc_t(fc%npwt,num_nbnds))
  if(fc%dual_t==4.d0) then
     evc_t(:,1:num_nbnds)=evc(:,1:num_nbnds)
  else
     call reorderwfp (num_nbnds,npw, fc%npwt,evc(:,:),evc_t(:,:), &
          &npw,fc%npwt, ig_l2g,fc%ig_l2gt, fc%ngmt_g , mpime, nproc,ionode_id, intra_pool_comm )

     !do iv=1,num_nbnds
     !   call mergewf(evc(:,iv),evc_g,npw,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
     !   call splitwf(evc_t(:,iv),evc_g,fc%npwt,fc%ig_l2gt,mpime,nproc,ionode_id,intra_pool_comm)
     !enddo
  endif






   if(l_verbose) write(stdout,*) 'self_basis_lanczos 1'
   FLUSH(stdout)
   !nbuf=min(5,nproc) ATTENZIONE iera cussi
   nbuf=nproc
   allocate(wp_prod(fc%nrxxt,numpw,nbuf))

   allocate(tmp_r(fc%nrxxt),tmp_r2(fc%nrxxt))
   allocate(omat(numpw,numpw),omat_hold(numpw,numpw))
   allocate(t_mat(numpw,nstates), t_mat_hold(numpw,nstates), t_mat_hold2(numpw,nstates))
   allocate(wp_g(npw,nstates))
   allocate(wp_g_t(fc%nrxxt,nstates))
   allocate(t_eigen_hold(nstates))
   allocate(wp_g_t3(npw*nproc))

!loop on kohn-sham states
   do ivv=first_state,last_state,nbuf
      call start_clock('sl_loop')
      do iv=ivv,min(ivv+nbuf-1,last_state)
!put iv on real space
      
         psic(:)=(0.d0,0.d0)
         psic(fc%nlt(1:fc%npwt))  = evc_t(1:fc%npwt,iv)
         psic(fc%nltm(1:fc%npwt)) = CONJG( evc_t(1:fc%npwt,iv) )
         CALL cft3t( fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, 2 )
         wv_real(1:fc%nrxxt)= DBLE(psic(1:fc%nrxxt))
         
!!loop on products of wanniers

        
         if(l_verbose) write(stdout,*) 'do fft'
         FLUSH(stdout)

         do ii=offset+1,numpw
            wp_prod(1:fc%nrxxt, ii,iv-ivv+1)=p_basis_r(1:fc%nrxxt,ii)*wv_real(1:fc%nrxxt)
         enddo


         if(l_verbose) write(stdout,*) 'calculate omat'
         FLUSH(stdout)

         
!!calculate overlap matrix
         call start_clock('sl_dgemm')
         call dgemm('T','N',numpw-offset,numpw-offset,fc%nrxxt,1.d0,wp_prod(1,offset+1,iv-ivv+1),&
     &fc%nrxxt,wp_prod(1,offset+1,iv-ivv+1),fc%nrxxt,0.d0,omat,numpw)
         call stop_clock('sl_dgemm')
         do ii=1,numpw-offset
            call mp_sum(omat(1:numpw-offset,ii),world_comm)
            omat(1:numpw-offset,ii)=omat(1:numpw-offset,ii)/dble(fc%nr1t*fc%nr2t*fc%nr3t)
         enddo

         if(iv-ivv==mpime) then
            omat_hold(:,:)=omat(:,:)
         endif
      enddo

!!
!!solve eigen/values vector problem 
!!
      

     if(l_verbose) write(stdout,*) 'solve eig'
     FLUSH(stdout)
     do iv=ivv,min(ivv+nbuf-1,last_state)
        if(iv-ivv==mpime) then
           call start_clock('sl_dsyevX')
           if(.not.l_dsyevr) then
              allocate(eigen(numpw-offset))
              allocate(work(1))
              call  DSYEV( 'V', 'U', numpw-offset, omat_hold, numpw, eigen, work, -1, info )
              lwork=work(1)
              deallocate(work)
              allocate(work(lwork))
              call  DSYEV( 'V', 'U', numpw-offset, omat_hold, numpw, eigen, work, lwork, info )
              deallocate(work)
              if(info/=0) then
                 write(stdout,*) 'ROUTINE self_basis_lanczos, INFO:', info
                 stop
              endif
         
         !do iw=1,numpw
         !   write(stdout,*) 'EIGEN:',iv,iw, eigen(iw)
         !enddo
         !FLUSH(stdout)
           else
              allocate(eigen(numpw-offset))
              allocate(vectors(numpw-offset,nstates))
              allocate(isuppz(2*nstates))
              allocate(work(1),iwork(1))
              call DSYEVR('V','I','U',numpw-offset,omat_hold,numpw,0.d0,0.d0,&
         &numpw-offset-nstates+1,numpw-offset,0.d0,n_found,eigen,&
                   & vectors,numpw-offset,isuppz,work, -1,iwork,-1, info)
              lwork=work(1)
              liwork=iwork(1)
              deallocate(work,iwork)
              allocate(work(lwork))
              allocate(iwork(liwork))
              vl=0.d0
              vu=0.d0
              il=numpw-offset-nstates+1
              iu=numpw-offset
              call DSYEVR('V','I','U',numpw-offset,omat_hold,numpw,vl,vu,il,iu,0.d0,n_found,eigen,&
                   & vectors,numpw-offset,isuppz,work,lwork,iwork,liwork, info)
              if(info/=0) then
                 write(stdout,*) 'ROUTINE pola_lanczos DSYEVR, INFO:', info
                 stop
              endif
              deallocate(isuppz)
              deallocate(work,iwork)
              !do iw=1,nstates
              !write(stdout,*) 'EIGEN:',iv,iw, eigen(iw)
              !enddo
              !FLUSH(stdout)

           endif
           call stop_clock('sl_dsyevX')
           t_mat_hold(:,:)=0.d0
           t_mat_hold2(:,:)=0.d0
           if(.not.l_dsyevr) then
              do ii=1,nstates
                 do jj=1,numpw-offset
                    t_mat_hold(jj+offset,ii)=omat_hold(jj,numpw-offset-ii+1)*(dsqrt(eigen(numpw-offset-ii+1)))
                 enddo
                 t_eigen_hold(ii)=eigen(numpw-ii+1)
              enddo
           else
              do ii=1,nstates
                 do jj=1,numpw-offset
                    t_mat_hold(jj+offset,ii)=vectors(jj,ii)*(dsqrt(eigen(ii)))
                 enddo
                 t_eigen_hold(ii)=eigen(ii)
              enddo
           endif
           if(.not.l_dsyevr) then
              do ii=1,nstates
                 t_mat_hold2(offset+1:numpw,ii)=omat_hold(1:numpw-offset,numpw-offset-ii+1)*(1.d0/dsqrt(eigen(numpw-offset-ii+1)))
              enddo
           else
              do ii=1,nstates
                 t_mat_hold2(offset+1:numpw,ii)=vectors(1:numpw-offset,ii)*(1.d0/dsqrt(eigen(ii)))
              enddo
           endif
           deallocate(eigen)
           if(l_dsyevr) deallocate(vectors)

        endif
     enddo
!!find liner dependent products

!!find transformation matrix and write on disk
!

      if(l_verbose) write(stdout,*) 'write on file'
      FLUSH(stdout)
      allocate(eigen(nstates))
      do iv=ivv,min(ivv+nbuf-1,last_state)
        if(iv-ivv == mpime) then
           t_mat(:,:)=t_mat_hold(:,:)
           eigen(1:nstates)=t_eigen_hold(1:nstates)
        endif
        call start_clock('sl_mpbcast')
        call mp_bcast(t_mat,iv-ivv,world_comm)
        call mp_bcast(eigen(1:nstates),iv-ivv,world_comm)
        call stop_clock('sl_mpbcast')
        if(ionode) then
           iuntmat = find_free_unit()
           write(nfile,'(4i1)') iv/1000,mod(iv,1000)/100,mod(iv,100)/10,mod(iv,10)
           if(ispin==1) then
              open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.s_mat_lanczos'//nfile, status='unknown',form='unformatted')
           else
              open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.s_mat_lanczos2'//nfile, status='unknown',form='unformatted')
           endif
           write(iuntmat) iv
           write(iuntmat) num_nbndv(1)!for compatibility with polarization file
           write(iuntmat) numpw
           write(iuntmat) nstates
           do ii=1,nstates
              write(iuntmat) t_mat(1:numpw,ii)
           enddo
           close(iuntmat)
        endif

!write on disk file with eigen values                                                  
        if(ionode) then
           iuntmat = find_free_unit()
           write(nfile,'(4i1)') iv/1000,mod(iv,1000)/100,mod(iv,100)/10,mod(iv,10)
           if(ispin==1) then
              open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.s_eig_lanczos'//nfile, status='unknown',form='unformatted')
           else
              open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.s_eig_lanczos2'//nfile, status='unknown',form='unformatted')
           endif
           write(iuntmat) nstates
           write(iuntmat) eigen(1:nstates)
           close(iuntmat)
        endif



!!find liner dependent products
     
        if(iv-ivv == mpime) then
           t_mat(:,:)=t_mat_hold2(:,:)
        endif
        call start_clock('sl_mpbcast')
        call mp_bcast(t_mat,iv-ivv,world_comm)
        call stop_clock('sl_mpbcast')

        call start_clock('sl_dgemm')
        call dgemm('N','N',fc%nrxxt,nstates,numpw,1.d0,wp_prod(1,1,iv-ivv+1),fc%nrxxt,t_mat,numpw,0.d0,wp_g_t,fc%nrxxt)
        call stop_clock('sl_dgemm')

!put the correct order
        psic=(0.d0,0.d0)
        write(stdout,*) 'do merge split',iv,ivv
        FLUSH(stdout)

        call start_clock('sl_merge')
    
        allocate(wp_g_t2(fc%npwt,nstates))
        do ii=1,nstates,2
           if(ii==nstates) then
              psic(1:fc%nrxxt)=cmplx(wp_g_t(1:fc%nrxxt,ii),0.d0)
           else
              psic(1:fc%nrxxt)=cmplx(wp_g_t(1:fc%nrxxt,ii),wp_g_t(1:fc%nrxxt,ii+1))
           endif
           CALL cft3t( fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, -2 )
           if(ii==nstates) then
              wp_g_t2(1:fc%npwt,ii) = psic(fc%nlt(1:fc%npwt))
           else
              wp_g_t2(1:fc%npwt, ii)= 0.5d0*(psic(fc%nlt(1:fc%npwt))+conjg( psic(fc%nltm(1:fc%npwt))))
              wp_g_t2(1:fc%npwt, ii+1)= (0.d0,-0.5d0)*(psic(fc%nlt(1:fc%npwt)) - conjg(psic(fc%nltm(1:fc%npwt))))

           endif
        enddo

        if(fc%dual_t==4.d0) then
           wp_g(1:npw,1:nstates)=wp_g_t2(1:fc%npwt,1:nstates)
        else
           call reorderwfp (nstates,fc%npwt, npw,wp_g_t2,wp_g, &
                &fc%npwt,npw, fc%ig_l2gt,ig_l2g, fc%ngmt_g , mpime, nproc,ionode_id, intra_pool_comm )

          ! call mergewfv(min(ii+nproc-1,nstates)-ii+1,fc%ngmt_g,wp_g_t2,&
          !      &evc_g,fc%npwt,fc%ig_l2gt,mpime,nproc,ionode_id,intra_pool_comm)
          ! call splitwfv(min(ii+nproc-1,nstates)-ii+1,fc%ngmt_g,wp_g_t3,&
          !      &evc_g,npw,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
                        

        endif
          
          
          call stop_clock('sl_merge')
          deallocate(wp_g_t2)
        if(l_verbose) write(stdout,*) 'do davcio'
        FLUSH(stdout)


  !!write on disk
        do ii=1,nstates
           call davcio(wp_g(:,ii),npw*2,iungresult,ii+(iv-first_state)*nstates,1)
        enddo

     

     enddo
     deallocate(eigen)
     call stop_clock('sl_loop')
  enddo
  deallocate(tmp_r,tmp_r2,omat,omat_hold,p_basis_r)
  deallocate(t_mat,t_mat_hold,t_mat_hold2)
  deallocate(wp_g,wp_g_t)
  close(iungresult)
  deallocate(wv_real,wp_prod)
  if(l_verbose) write(stdout,*) 'Exiting self_basis_lanczos'
  FLUSH(stdout)
  deallocate(p_basis_t,evc_g)
  if(l_verbose) write(stdout,*) 'Call deallocate_fft_custom'
  FLUSH(stdout)
  call deallocate_fft_custom(fc)
  deallocate(t_eigen_hold)
  deallocate(p_basis)
  return
!NOT_TO_BE_INCLUDED_END
end subroutine self_basis_lanczos_real




      SUBROUTINE mergewfv ( nn,lda,pw, pwt, ngwl, ig_l2g, mpime, nproc, root, comm )

! ... This subroutine merges the pieces of a wave functions (pw) splitted across 
! ... processors into a total wave function (pwt) containing al the components
! ... in a pre-defined order (the same as if only one processor is used)

      USE kinds
      USE parallel_include

      IMPLICIT NONE

      COMPLEX(DP), intent(in) :: PW(ngwl*nproc)
      COMPLEX(DP), intent(out) :: PWT(lda*nproc)
      INTEGER, intent(in) :: nn !number of wfcs
      INTEGER, intent(in) :: lda!leading dimension of pwt
      INTEGER, INTENT(IN) :: mpime     ! index of the calling processor ( starting from 0 )
      INTEGER, INTENT(IN) :: nproc     ! number of processors
      INTEGER, INTENT(IN) :: root      ! root processor ( the one that should receive the data )
      INTEGER, INTENT(IN) :: comm    ! communicator
      INTEGER, INTENT(IN) :: ig_l2g(ngwl)
      INTEGER, INTENT(IN) :: ngwl

      INTEGER, ALLOCATABLE :: ig_ip(:)
      COMPLEX(DP), ALLOCATABLE :: pw_ip(:)

      INTEGER :: ierr, i, ip, ngw_ip, ngw_lmax, itmp, igwx, gid, J

#if defined __MPI
      INTEGER :: istatus(MPI_STATUS_SIZE)
#endif

!
! ... Subroutine Body
!

      igwx = MAXVAL( ig_l2g(1:ngwl) )

#if defined __MPI

      gid = comm

! ... Get local and global wavefunction dimensions
      CALL MPI_ALLREDUCE( ngwl, ngw_lmax, 1, MPI_INTEGER, MPI_MAX, gid, IERR )
      CALL MPI_ALLREDUCE( igwx, itmp, 1, MPI_INTEGER, MPI_MAX, gid, IERR )
      igwx = itmp

#endif

      IF( igwx > SIZE( pwt ) ) &
        CALL errore(' mergewf ',' wrong size for pwt ',SIZE(pwt) )

#if defined __MPI

      DO ip = 1, nproc

        IF( (ip-1) /= root ) THEN

! ...     In turn each processors send to root the wave components and their indexes in the 
! ...     global array
          IF ( mpime == (ip-1) ) THEN
            CALL MPI_SEND( ig_l2g, ngwl, MPI_INTEGER, ROOT, IP, gid, IERR )
            CALL MPI_SEND( pw(1), ngwl*nn, MPI_DOUBLE_COMPLEX, ROOT, IP+NPROC, gid, IERR )
          END IF
          IF ( mpime == root) THEN
            ALLOCATE(ig_ip(ngw_lmax))
            ALLOCATE(pw_ip(ngw_lmax*nn))
            CALL MPI_RECV( ig_ip, ngw_lmax, MPI_INTEGER, (ip-1), IP, gid, istatus, IERR )
            CALL MPI_GET_COUNT( istatus, MPI_INTEGER, ngw_ip, ierr )
            CALL MPI_RECV( pw_ip, ngw_lmax*nn, MPI_DOUBLE_COMPLEX, (ip-1), IP+NPROC, gid, istatus, IERR )
           ! CALL MPI_GET_COUNT( istatus, MPI_DOUBLE_COMPLEX, ngw_ip, ierr )
           ! ngw_ip=ngw_ip/nn
            DO J=1,nn
               DO I = 1, ngw_ip
                  PWT(ig_ip(i)+lda*(J-1)) = pw_ip(i+(ngw_ip)*(J-1))
               END DO
            END DO
            DEALLOCATE(ig_ip)
            DEALLOCATE(pw_ip)
          END IF

        ELSE

          IF(mpime == root) THEN
             DO J=1,nn
                DO I = 1, ngwl
                   PWT(ig_l2g(i)+lda*(J-1)) = pw(i+ngwl*(J-1))
                END DO
            END DO
          END IF

        END IF

        CALL MPI_BARRIER( gid, IERR )

      END DO

#else

      do J=1,nn
         DO I = 1, ngwl
            ! WRITE( stdout,*) 'MW ', ig_l2g(i), i
            PWT( ig_l2g(i)+lda*(J-1) ) = pw(i+ngwl*(J-1))
         END DO
      END do

#endif

      RETURN
    END SUBROUTINE mergewfv

!=----------------------------------------------------------------------------=!


      SUBROUTINE splitwfv ( nn, lda, pw, pwt, ngwl, ig_l2g, mpime, nproc, root, comm )

! ... This subroutine splits a total wave function (pwt) containing al the components
! ... in a pre-defined order (the same as if only one processor is used), across 
! ... processors (pw).

      USE kinds
      USE parallel_include
      IMPLICIT NONE

      COMPLEX(DP), INTENT(OUT) :: PW(ngwl*nproc)
      COMPLEX(DP), INTENT(IN) :: PWT(lda*nproc)
      INTEGER, intent(in) :: nn !number of wfcs                       
      INTEGER, intent(in) :: lda!leading dimension of pwt                
      INTEGER, INTENT(IN) :: mpime, nproc, root
      INTEGER, INTENT(IN) :: comm    ! communicator
      INTEGER, INTENT(IN) :: ig_l2g(ngwl)
      INTEGER, INTENT(IN) :: ngwl

      INTEGER, ALLOCATABLE :: ig_ip(:)
      COMPLEX(DP), ALLOCATABLE :: pw_ip(:)

      INTEGER ierr, i, ngw_ip, ip, ngw_lmax, gid, igwx, itmp, J

#if defined __MPI
      integer istatus(MPI_STATUS_SIZE)
#endif

!
! ... Subroutine Body
!

      igwx = MAXVAL( ig_l2g(1:ngwl) )

#if defined __MPI

      gid = comm

! ... Get local and global wavefunction dimensions
      CALL MPI_ALLREDUCE(ngwl, ngw_lmax, 1, MPI_INTEGER, MPI_MAX, gid, IERR )
      CALL MPI_ALLREDUCE(igwx, itmp    , 1, MPI_INTEGER, MPI_MAX, gid, IERR )
      igwx = itmp

#endif

      IF( igwx > SIZE( pwt ) ) &
        CALL errore(' splitwf ',' wrong size for pwt ',SIZE(pwt) )

#if defined __MPI

      DO ip = 1, nproc
! ...   In turn each processor send to root the the indexes of its wavefunction conponents
! ...   Root receive the indexes and send the componens of the wavefunction read from the disk (pwt)
        IF ( (ip-1) /= root ) THEN
          IF ( mpime == (ip-1) ) THEN
            CALL MPI_SEND( ig_l2g, ngwl, MPI_INTEGER, ROOT, IP, gid,IERR)
            CALL MPI_RECV( pw(1), ngwl*nn, MPI_DOUBLE_COMPLEX, ROOT, IP+NPROC, gid, istatus, IERR )
          END IF
          IF ( mpime == root ) THEN
            ALLOCATE(ig_ip(ngw_lmax))
            ALLOCATE(pw_ip(ngw_lmax*nn))
            CALL MPI_RECV( ig_ip, ngw_lmax, MPI_INTEGER, (ip-1), IP, gid, istatus, IERR )
            CALL MPI_GET_COUNT(istatus, MPI_INTEGER, ngw_ip, ierr)
            DO J=1,nn
               DO i = 1, ngw_ip
                  pw_ip(i+ngw_ip*(J-1)) = PWT(ig_ip(i)+lda*(J-1))
               END DO
            END DO
            CALL MPI_SEND( pw_ip, ngw_ip*nn, MPI_DOUBLE_COMPLEX, (ip-1), IP+NPROC, gid, IERR )
            DEALLOCATE(ig_ip)
            DEALLOCATE(pw_ip)
          END IF
        ELSE
          IF ( mpime == root ) THEN
             DO J=1,nn
                DO i = 1, ngwl
                   pw(i+ngwl*(J-1)) = PWT(ig_l2g(i)+lda*(J-1)) 
                END DO
             END DO
          END IF
        END IF
        CALL MPI_BARRIER(gid, IERR)
      END DO

#elif ! defined __MPI

      DO J=1,nn
         DO I = 1, ngwl
            pw(i+ngwl*(J-1)) = pwt( ig_l2g(i)+lda*(J-1) )
         END DO
      END DO
#else

      CALL errore(' SPLITWF ',' no communication protocol ',0)

#endif

      RETURN
      END SUBROUTINE splitwfv

!=----------------------------------------------------------------------------=!

