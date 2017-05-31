!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

!routines for the calculation of the polarization
!lanczos-style
!ONLY FOR NORMCONSERVING PSEUDOS !!!!!

subroutine pola_basis_lanczos(n_set,nstates,numpw, nsteps,ispin)
!this subroutine calculates the basis for every v 
!the minimal orthonormal basis for the w_v(r)*w^P'_i(r) products

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
   INTEGER, INTENT(in) :: ispin! spin channel

   INTEGER :: iv,iw,ig,ii,jj
   REAL(kind=DP), ALLOCATABLE :: wv_real(:),tmp_r(:),tmp_r2(:)
   COMPLEX(kind=DP), ALLOCATABLE :: tmp_g(:), wp_prod(:,:,:)
   INTEGER :: iungprod,iunrprod, iungresult,iuntmat
   LOGICAL :: exst
   REAL(kind=DP), ALLOCATABLE :: omat(:,:),omat_hold(:,:)

   REAL(kind=DP), ALLOCATABLE :: eigen(:),work(:)
   INTEGER :: lwork,info,liwork
   COMPLEX(kind=DP), ALLOCATABLE :: wp_g(:,:)!product terms in g wfc grid
   COMPLEX(kind=DP), ALLOCATABLE :: wp_g_t(:,:)!
   REAL(kind=DP), ALLOCATABLE :: t_mat(:,:),t_mat_hold(:,:), t_mat_hold2(:,:)
   CHARACTER(4) :: nfile
   LOGICAL :: l_reduce_io=.true.!if true reduces io
   COMPLEX(kind=DP), ALLOCATABLE :: p_basis(:,:)!polarizability basis
   LOGICAL :: l_dsyevr=.true.!if true uses dsyevr 
   REAL(kind=DP), ALLOCATABLE :: vectors(:,:)!for dsyevr
   INTEGER, ALLOCATABLE :: iwork(:), ifail(:)
   INTEGER, ALLOCATABLE :: isuppz(:)
   INTEGER :: n_found

   LOGICAL :: l_fft_custom=.false.!if true uses custom fft grid
   COMPLEX(kind=DP), ALLOCATABLE :: evc_t(:,:),p_basis_t(:,:)
   COMPLEX(kind=DP), ALLOCATABLE :: evc_g(:)
   LOGICAL :: l_sumrule=.false.!if true imposes the sum rule over the norm of Pc|\Phi_\mu\Psi_v> for each of them
   REAL(kind=DP), ALLOCATABLE :: norms(:)
   REAL(kind=DP) :: norm_t, c_norm,norm
   REAL(kind=DP), ALLOCATABLE :: p_basis_r(:,:) !polarizabilty basis in real custom space

   INTEGER :: ivv,nbuf
   REAL(kind=DP) :: vl,vu
   INTEGER :: il,iu
   REAL(kind=DP), ALLOCATABLE :: t_eigen_hold(:)

   TYPE(fft_cus) :: fc

   write(stdout,*) 'Routine pola_basis_lanczos'
   FLUSH(stdout)

   fc%ecutt=ecutwfc
   fc%dual_t=dual_vt

   if(l_verbose) write(stdout,*) 'Call initialize_fft_custom'
   FLUSH(stdout)
   call initialize_fft_custom(fc)
   allocate(evc_g(fc%ngmt_g))


   allocate(wv_real(fc%nrxxt))


   allocate(norms(numpw))

!read w^P'_i on file on real space
!open product of wanniers filed


   iungprod = find_free_unit()
   CALL diropn( iungprod, 'wiwjwfc_red', max_ngm*2, exst )

   if(.not.l_reduce_io) then
      iunrprod = find_free_unit()
      CALL diropn( iunrprod, 'wiwjwfc_red_r', dfftp%nnr, exst )
   endif

   iungresult = find_free_unit()
   CALL diropn( iungresult, 'vw_lanczos',npw*2, exst)

   if(.not.l_reduce_io) then
      allocate(tmp_g(max_ngm),tmp_r(dfftp%nnr))
      do iw=1,numpw
         call davcio(tmp_g,max_ngm*2,iungprod,iw,-1)
         !trasform to r-space
         psic(:)=(0.d0,0.d0)
         do ig=1,max_ngm
            psic(nl(ig))=tmp_g(ig)
            psic(nlm(ig))=CONJG(tmp_g(ig))
         enddo
         CALL invfft ('Dense', psic, dfftp)
         tmp_r(:)=dble(psic(:))
         call davcio(tmp_r,dfftp%nnr,iunrprod,iw,1)
      enddo
      deallocate(tmp_g,tmp_r)
      close(iunrprod)
   else
!read polarizability basis functions
      allocate(p_basis(max_ngm,numpw))
       do iw=1,numpw
          call davcio(p_basis(:,iw),max_ngm*2,iungprod,iw,-1)
      enddo
   endif
   close(iungprod)
  

   if(l_verbose) write(stdout,*) 'pola_basis_lanczos 1'
   FLUSH(stdout)
   
!now polarizability basis  are put on the ordering of the redueced grid, if required
  allocate(p_basis_t(fc%npwt,numpw))
  if(fc%dual_t==4.d0) then
     p_basis_t(:,:)=p_basis(:,:)
  else
     call reorderwfp_col(numpw,npw,fc%npwt,p_basis(1,1),p_basis_t(1,1), npw,fc%npwt, &
             & ig_l2g,fc%ig_l2gt,fc%ngmt_g,mpime, nproc,intra_pool_comm )
     !do ii=1,numpw
     !   call mergewf(p_basis(:,ii),evc_g,npw,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
     !   call splitwf(p_basis_t(:,ii),evc_g,fc%npwt,fc%ig_l2gt,mpime,nproc,ionode_id,intra_pool_comm)
     !enddo
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
     CALL cft3t(fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, 2 )
     p_basis_r(1:fc%nrxxt,ii)= DBLE(psic(1:fc%nrxxt))
     if(ii/=numpw) p_basis_r(1:fc%nrxxt,ii+1)= DIMAG(psic(1:fc%nrxxt))
     
  enddo
  


!now valence wavefunctions are put on the ordering of the reduced grid
  allocate(evc_t(fc%npwt,num_nbndv(ispin)))
  if(fc%dual_t==4.d0) then
     evc_t(:,1:num_nbndv(ispin))=evc(:,1:num_nbndv(ispin))
  else
     call reorderwfp_col(num_nbndv(ispin),npw,fc%npwt,evc(1,1),evc_t(1,1), npw,fc%npwt, &
           & ig_l2g,fc%ig_l2gt,fc%ngmt_g,mpime, nproc,intra_pool_comm )
  endif


!loop on v

  allocate(tmp_r(fc%nrxxt),tmp_r2(fc%nrxxt))
  allocate(omat(numpw,numpw),omat_hold(numpw,numpw))
  allocate(t_mat(numpw,nstates), t_mat_hold(numpw,nstates), t_mat_hold2(numpw,nstates))
  allocate(wp_g(npw,nstates))
  allocate(wp_g_t(fc%npwt,nstates))
  allocate(t_eigen_hold(nstates))

  nbuf=min(5,nproc)
  allocate(wp_prod(fc%npwt,numpw,nbuf))

   do ivv=1,num_nbndv(ispin),nbuf
!put iv on real space
      
      do iv=ivv,min(ivv+nbuf-1,num_nbndv(ispin))
         psic(:)=(0.d0,0.d0)
         psic(fc%nlt(1:fc%npwt))  = evc_t(1:fc%npwt,iv)
         psic(fc%nltm(1:fc%npwt)) = CONJG( evc_t(1:fc%npwt,iv) )
         CALL cft3t(fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, 2 )
         wv_real(1:fc%nrxxt)= DBLE(psic(1:fc%nrxxt))
   

!!loop on products of wanniers
      
         if(.not.l_reduce_io) then
            iunrprod = find_free_unit()
            CALL diropn( iunrprod, 'wiwjwfc_red_r', dfftp%nnr, exst )
         endif
    !  allocate(tmp_r(fc%nrxxt))

         if(l_verbose) write(stdout,*) 'do fft'
         FLUSH(stdout)

         do ii=1,numpw,2
!!read n_set w^P'_i from disk
            if(.not.l_reduce_io) then
               call davcio(tmp_r,dfftp%nnr,iunrprod,ii,-1)
               write(stdout,*) 'ERROR l_reduce_io must be true'
               FLUSH(stdout)
               stop
            endif

           tmp_r(1:fc%nrxxt)=p_basis_r(1:fc%nrxxt,ii)*wv_real(1:fc%nrxxt)
          
            if(ii/=numpw) then
               tmp_r2(1:fc%nrxxt)=p_basis_r(1:fc%nrxxt,ii+1)*wv_real(1:fc%nrxxt) 
            else
               tmp_r2(1:fc%nrxxt)=0.d0
            endif

           
!!form products with w_v and trasfrom in G space
            psic(1:fc%nrxxt)=dcmplx(tmp_r(1:fc%nrxxt),tmp_r2(1:fc%nrxxt))
            CALL cft3t(fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, -2 )
          
            if(ii==numpw) then
               wp_prod(1:fc%npwt, ii,iv-ivv+1) = psic(fc%nlt(1:fc%npwt))
            else
               wp_prod(1:fc%npwt, ii,iv-ivv+1)= 0.5d0*(psic(fc%nlt(1:fc%npwt))+conjg( psic(fc%nltm(1:fc%npwt))))
               wp_prod(1:fc%npwt, ii+1,iv-ivv+1)= (0.d0,-0.5d0)*(psic(fc%nlt(1:fc%npwt)) - conjg(psic(fc%nltm(1:fc%npwt))))
            endif
            
         enddo
         if(l_verbose) write(stdout,*) 'do pc_operator'
         FLUSH(stdout)
        call pc_operator_t_m(numpw,wp_prod(1,1,iv-ivv+1),evc_t,ispin, fc) 

         if(l_verbose) write(stdout,*) 'calculate omat'
         FLUSH(stdout)


    
         if(.not.l_reduce_io) close(iunrprod)
!!calculate overlap matrix
    
         call dgemm('T','N',numpw,numpw,2*fc%npwt,2.d0,wp_prod(1,1,iv-ivv+1),2*fc%npwt,&
              &wp_prod(1,1,iv-ivv+1),2*fc%npwt,0.d0,omat,numpw)
         if(fc%gstart_t==2) then
            do ii=1,numpw
               do jj=1,numpw
                  omat(jj,ii)=omat(jj,ii)-dble(conjg(wp_prod(1,jj,iv-ivv+1))*wp_prod(1,ii,iv-ivv+1))
               enddo
            enddo
         endif
         do ii=1,numpw
            call mp_sum(omat(1:numpw,ii), world_comm)
         enddo
!set up norms
!         do ii=1,numpw
!            norms(ii)=omat(ii,ii)
!         enddo


         if(iv-ivv==mpime) then
            omat_hold(:,:)=omat(:,:)
         endif
      enddo
  
!!
!!solve eigen/values vector problem 
!!
      
     if(l_verbose)  write(stdout,*) 'solve eig'
     FLUSH(stdout)

    
      FLUSH(stdout)
      do iv=ivv,min(ivv+nbuf-1,num_nbndv(ispin))
         if(l_verbose) write(stdout,*) 'solve eig', iv
         FLUSH(stdout)
         if(iv-ivv==mpime) then
            
            if(.not.l_dsyevr) then
               allocate(eigen(numpw))
               allocate(work(1))
               call  DSYEV( 'V', 'U', numpw, omat_hold, numpw, eigen, work, -1, info )
               lwork=work(1)
               deallocate(work)
               allocate(work(lwork))
               call  DSYEV( 'V', 'U', numpw, omat_hold, numpw, eigen, work, lwork, info )
               deallocate(work)
               if(info/=0) then
                  write(stdout,*) 'ROUTINE pola_basis_lanczos, INFO:', info
                  stop
               endif
            

             !  do iw=1,numpw
             !     write(stdout,*) 'EIGEN:',iv,iw, eigen(iw)
             !  enddo
             !  FLUSH(stdout)
            else
               if(l_verbose) write(stdout,*) 'ATT1'
               FLUSH(stdout)
               allocate(eigen(numpw))
               allocate(vectors(numpw,nstates))
               allocate(isuppz(2*nstates))
               allocate(work(1),iwork(1))
               if(l_verbose)  write(stdout,*) 'ATT2'
                FLUSH(stdout)
               call DSYEVR('V','I','U',numpw,omat_hold,numpw,0.d0,0.d0,numpw-nstates+1,numpw,0.d0,n_found,eigen,&
                    & vectors,numpw,isuppz,work, -1,iwork,-1, info)
               lwork=work(1)
               liwork=iwork(1)
               deallocate(work,iwork)
               allocate(work(lwork))
               allocate(iwork(liwork))
               if(l_verbose)  write(stdout,*) 'ATT3',numpw,nstates,size(omat_hold(:,1)),size(omat_hold(1,:)),lwork,liwork
               FLUSH(stdout)
               vl=0.d0
               vu=0.d0
               il=numpw-nstates+1
               iu=numpw
               call DSYEVR('V','I','U',numpw,omat_hold,numpw,vl,vu,il,iu,0.d0,n_found,eigen,&
                    & vectors,numpw,isuppz,work,lwork,iwork,liwork, info)
               if(info/=0) then
                  write(stdout,*) 'ROUTINE pola_lanczos DSYEVR, INFO:', info
                  stop
               endif
               if(l_verbose)  write(stdout,*) 'ATT4'
               FLUSH(stdout)
               deallocate(isuppz)
               deallocate(work,iwork)
        
               do iw=1,nstates,nstates-1
                  write(stdout,*) 'EIGEN T LOCAL:',iv,iw, eigen(iw)
               enddo
               FLUSH(stdout)

            endif
      
      
!!find transformation matrix and write on disk
!
        if(l_verbose) write(stdout,*) 'pola_basis_lanczos t_mat'
        FLUSH(stdout)


      

        if(.not.l_dsyevr) then
           do ii=1,nstates
              do jj=1,numpw
                 t_mat_hold(jj,ii)=omat_hold(jj,numpw-ii+1)*(dsqrt(eigen(numpw-ii+1)))
              enddo
              t_eigen_hold(ii)=eigen(numpw-ii+1)
           enddo
        else
           do ii=1,nstates
              do jj=1,numpw
                 t_mat_hold(jj,ii)=vectors(jj,ii)*(dsqrt(eigen(ii)))
              enddo
              t_eigen_hold(ii)=eigen(ii)
           enddo
        endif

        !!find liner dependent products

        if(.not.l_dsyevr) then
           do ii=1,nstates
              t_mat_hold2(:,ii)=omat_hold(:,numpw-ii+1)*(1.d0/dsqrt(eigen(numpw-ii+1)))
           enddo
        else
           do ii=1,nstates
              t_mat_hold2(:,ii)=vectors(:,ii)*(1.d0/dsqrt(eigen(ii)))
           enddo
        endif
     
        deallocate(eigen)
        if(l_dsyevr) deallocate(vectors)

     endif
     

    enddo
    allocate(eigen(nstates))
     do iv=ivv,min(ivv+nbuf-1,num_nbndv(ispin))
        if(iv-ivv == mpime) then
           t_mat(:,:)=t_mat_hold(:,:)
           eigen(1:nstates)=t_eigen_hold(1:nstates)
        endif
        call mp_bcast(t_mat,iv-ivv,world_comm)
        call mp_bcast(eigen(1:nstates),iv-ivv,world_comm)
!if required imposes sum rule

!        if(l_sumrule) then
!           norm_t=0.d0
!           do jj=1,numpw
!              do ii=1,nstates
!                 norm_t=norm_t+t_mat(jj,ii)**2.d0
!              enddo
!           enddo
!           norm=0.d0
!           do jj=1,numpw
!              norm=norm+norms(jj)
!           enddo
!           c_norm=dsqrt(norm/norm_t)
!           write(stdout,*) 'Sum rule:',c_norm
!           t_mat(:,:)=t_mat(:,:)*c_norm
!        endif

        if(ionode) then
           iuntmat = find_free_unit()
           write(nfile,'(4i1)') iv/1000,mod(iv,1000)/100,mod(iv,100)/10,mod(iv,10)
           if(ispin==1) then
              open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.p_mat_lanczos'//nfile, status='unknown',form='unformatted')
           else
              open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.p_mat_lanczos2'//nfile, status='unknown',form='unformatted')
           endif
           write(iuntmat) iv
           write(iuntmat) num_nbndv(ispin)
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
              open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.p_eig_lanczos'//nfile, status='unknown',form='unformatted')
           else
              open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.p_eig_lanczos2'//nfile, status='unknown',form='unformatted')
           endif
           write(iuntmat) nstates
           write(iuntmat) eigen(1:nstates)
           close(iuntmat)
        endif


        if(l_verbose) write(stdout,*) 'pola_basis update wp_g'
        FLUSH(stdout)


!!find liner dependent products

        if(iv-ivv == mpime) then
           t_mat(:,:)=t_mat_hold2(:,:)
        endif
        call mp_bcast(t_mat,iv-ivv,world_comm)



        if(l_verbose) write(stdout,*) 'pola_basis update wp_g dgemm'
        FLUSH(stdout)

        call dgemm('N','N',2*fc%npwt,nstates,numpw,1.d0,wp_prod(1,1,iv-ivv+1),2*fc%npwt,t_mat,numpw,0.d0,wp_g_t,2*fc%npwt)

        write(stdout,*) 'pola_basis update merge-split',iv,ivv
        FLUSH(stdout)

!put the correct order

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
        if(l_verbose) write(stdout,*) 'pola_basis update davcio',iv
        FLUSH(stdout)


  !!write on disk
        do ii=1,nstates
           call davcio(wp_g(:,ii),npw*2,iungresult,ii+(iv-1)*nstates,1)
        enddo

        if(l_verbose) write(stdout,*) 'pola_basis update done'
        FLUSH(stdout)


       
     enddo
     deallocate(eigen)
  enddo
  deallocate(t_mat,t_mat_hold,t_mat_hold2)
  deallocate(omat,omat_hold,tmp_r,tmp_r2,p_basis_r)
  deallocate(wp_g,wp_g_t)
  close(iungresult)
  deallocate(norms)
  deallocate(wv_real,wp_prod)
  deallocate(t_eigen_hold)
  
  if(l_verbose) write(stdout,*) 'Exiting pola_basis_lanczos'
  FLUSH(stdout)
   
  if(l_reduce_io) deallocate(p_basis)
  deallocate(p_basis_t,evc_t)
  
  if(l_verbose) write(stdout,*) 'Call deallocate_fft_custom'
  FLUSH(stdout)
  deallocate(evc_g)
  call deallocate_fft_custom(fc)
   

   return
 end subroutine pola_basis_lanczos


subroutine pc_operator(state,ispin,l_cond)
!this operator project the wavefunction state on the conduction 
!subspace, the valence wavefunction are in evc
!ONLY FOR GAMMA POINT NOW!!!!
   USE io_global,            ONLY : stdout
   USE kinds,    ONLY : DP
   USE gvect
   USE wvfct,    ONLY : npwx, npw, nbnd
   USE wavefunctions_module, ONLY : evc, psic
   USE mp, ONLY : mp_sum, mp_barrier, mp_bcast
   USE mp_world, ONLY : world_comm
   USE wannier_gw, ONLY : num_nbndv,num_nbnds
   
  implicit none

  COMPLEX(kind=DP), INTENT(inout) :: state(npw)!state to be projected
  INTEGER, INTENT(in) :: ispin!spin channel
  LOGICAL :: l_cond!if true project out alson conduction states till num_nbnds

  INTEGER :: iv,ig
  REAL(kind=DP), ALLOCATABLE :: prod(:)

  if(.not.l_cond) then
     if(num_nbndv(ispin)==0) return
     allocate(prod(num_nbndv(ispin)))
     call dgemm('T','N', num_nbndv(ispin),1,2*npw,2.d0,evc,2*npwx,state,2*npw,&
          & 0.d0,prod,num_nbndv(ispin))
     do iv=1,num_nbndv(ispin)
        if(gstart==2) prod(iv)=prod(iv)-dble(conjg(evc(1,iv))*state(1))
     enddo
     call mp_sum(prod(:), world_comm)
     call dgemm('N','N',2*npw,1,num_nbndv(ispin),-1.d0,evc,2*npwx,prod,&
          &num_nbndv(ispin),1.d0,state,2*npw)
  else
     allocate(prod(num_nbnds))
     call dgemm('T','N', num_nbnds,1,2*npw,2.d0,evc,2*npwx,state,2*npw,&
          & 0.d0,prod,num_nbnds)
     do iv=1,num_nbnds
        if(gstart==2) prod(iv)=prod(iv)-dble(conjg(evc(1,iv))*state(1))
     enddo
     call mp_sum(prod(:), world_comm)
     call dgemm('N','N',2*npw,1,num_nbnds,-1.d0,evc,2*npwx,prod,&
          &num_nbnds,1.d0,state,2*npw)

  endif
 
  deallocate(prod)
  return
end subroutine pc_operator

subroutine pc_operator_t(state,evc_t,ispin, fc)
!this operator project the wavefunction state on the conduction 
!subspace, the valence wavefunction are in evc
!ONLY FOR GAMMA POINT NOW!!!!
   USE io_global,            ONLY : stdout
   USE kinds,    ONLY : DP
   USE gvect
   USE wvfct,    ONLY : npwx, npw, nbnd
   USE wavefunctions_module, ONLY : evc, psic
   USE mp, ONLY : mp_sum, mp_barrier, mp_bcast
   USE mp_world, ONLY : world_comm
   USE wannier_gw, ONLY : num_nbndv
   USE fft_custom_gwl
   USE fft_base,             ONLY : dfftp, dffts
   USE fft_interfaces,       ONLY : fwfft, invfft


  implicit none

  TYPE(fft_cus), INTENT(in) :: fc
  COMPLEX(kind=DP), INTENT(inout) :: state(fc%npwt)!state to be projected
  ! COMPLEX(kind=DP), INTENT(inout) :: evc_t(fc%npwt,num_nbndv(ispin))
  ! above syntax not accepted by all compilers
  COMPLEX(kind=DP), INTENT(inout) :: evc_t(fc%npwt,*)!valence states
  INTEGER, INTENT(in) :: ispin!spin channel

  INTEGER :: iv,ig
  REAL(kind=DP), ALLOCATABLE :: prod(:)

  allocate(prod(num_nbndv(ispin)))
   call dgemm('T','N', num_nbndv(ispin),1,2*fc%npwt,2.d0,evc_t,2*fc%npwt,state,2*fc%npwt,&
        & 0.d0,prod,num_nbndv(ispin))
  do iv=1,num_nbndv(ispin)
     if(fc%gstart_t==2) prod(iv)=prod(iv)-dble(conjg(evc_t(1,iv))*state(1))
  enddo
  call mp_sum(prod(:), world_comm)
  call dgemm('N','N',2*fc%npwt,1,num_nbndv(ispin),-1.d0,evc_t,2*fc%npwt,prod,&
       &num_nbndv(ispin),1.d0,state,2*fc%npwt)

 
  deallocate(prod)
  return
end subroutine pc_operator_t 


subroutine lanczos_state(zstates, nstates, itype, nsteps,istate,ispin)
!this subroutine perform nsteps collective lanczos iterations
!on orthonormal zstates state
!GAMMA POINT ONLY!!!
 
  USE io_global,            ONLY : stdout, ionode, ionode_id
  USE io_files,             ONLY : prefix, tmp_dir
  USE kinds,    ONLY : DP
  USE wannier_gw
  USE gvect
  USE constants, ONLY : e2, pi, tpi, fpi
  USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2
  USE wvfct,    ONLY : g2kin, npwx, npw, nbnd
  USE wavefunctions_module, ONLY : evc, psic
  USE mp, ONLY : mp_sum, mp_barrier, mp_bcast
  USE mp_world, ONLY : mpime, nproc, world_comm
  USE gvecs,              ONLY : nls, nlsm, doublegrid
  USE g_psi_mod,            ONLY : h_diag, s_diag
  USE becmod,           ONLY : becp,allocate_bec_type,deallocate_bec_type
  USE uspp,                 ONLY : vkb, nkb, okvan
  USE klist,                ONLY : xk,igk_k
  USE noncollin_module,     ONLY : noncolin, npol
  USE fft_base,             ONLY : dfftp, dffts
  USE fft_interfaces,       ONLY : fwfft, invfft
   


  implicit none

  INTEGER, EXTERNAL :: find_free_unit

  COMPLEX(kind=DP), INTENT(in) :: zstates(npw,nstates)!states for starting lanczos chains
  INTEGER, INTENT(in) :: nstates!number of states
  INTEGER, INTENT(in) :: itype!matrices to be saved: 0 for polarization; 1 for self-energy; 2 for other uses
  INTEGER, INTENT(in) :: nsteps!number of Lanczos iteration to be performed
  INTEGER, INTENT(in) :: istate!corresponding KS state(for labelling output files)
  INTEGER, INTENT(in) :: ispin!spin channel 1,2
  


  COMPLEX(kind=DP), ALLOCATABLE :: psi_1(:,:),psi_2(:,:),psi_3(:,:)
  COMPLEX(kind=DP), ALLOCATABLE :: u_0(:,:),u_1(:,:)
  REAL(kind=DP), ALLOCATABLE :: alpha(:),beta(:),gamma(:), n_1(:)
  REAL(kind=DP), ALLOCATABLE :: d(:,:),f(:,:),c(:)
  REAL(kind=DP), ALLOCATABLE :: omat(:,:)!overlap with intermediate lanczos states to be saved on disk for each iteration!
  CHARACTER(4) :: nfile
  INTEGER :: is,ig,ii,jj,it
  INTEGER :: iunlan
  LOGICAL :: omat_div = .false.
  INTEGER :: l_blk,nbegin,nend,nsize
  INTEGER :: ip,nbegin_ip, nend_ip, nsize_ip
  REAL(kind=DP), ALLOCATABLE :: omat_tot(:,:,:)
  INTEGER :: iv

  
 
!  allocate(alpha(nstates),beta(nstates),gamma(nstates),n_1(nstates))
!  allocate(d(nsteps,nstates),f(nsteps,nstates),c(nstates))

!if omat_div is true distribute memory on MPI tasks, but it is slower
  if(nstates < 1600) then
     omat_div=.false.
  else
     omat_div=.true.
  endif

  


  call start_clock('lanczos_state')

!to decide whether to calculate the overlap matrix column by column or distributed on processors
 
  if(omat_div) then
     l_blk= (nstates)/nproc
     if(l_blk*nproc < (nstates)) l_blk = l_blk+1
     nbegin=mpime*l_blk+1
     nend=nbegin+l_blk-1
     if(nend > nstates) nend=nstates
     nsize=nend-nbegin+1
  
     allocate(omat(nstates,l_blk))
     allocate(omat_tot(nstates,l_blk,nsteps))
  else
     l_blk=nstates
     allocate(omat(nstates,nstates))
  endif
  
  allocate(psi_1(npw,l_blk),psi_2(npw,l_blk),psi_3(npw,l_blk))
  allocate(u_0(npw,l_blk),u_1(npw,l_blk))
  allocate(alpha(l_blk),beta(l_blk),gamma(l_blk),n_1(l_blk))
  allocate(d(nsteps,nstates),f(nsteps,nstates),c(l_blk))
  d(:,:)=0.d0
  f(:,:)=0.d0

!loop on l_blk
  ALLOCATE( h_diag( npwx,npol ) )
  ALLOCATE( s_diag( npwx,npol ) )

  call allocate_bec_type ( nkb, l_blk, becp)

  IF ( nkb > 0 )  CALL init_us_2( npw, igk_k(1,1), xk(1,1), vkb )
  g2kin(1:npw) = ( (g(1,igk_k(1:npw,1)) )**2 + &
       ( g(2,igk_k(1:npw,1)) )**2 + &
       ( g(3,igk_k(1:npw,1)) )**2 ) * tpiba2

  ip=0
  do iv=1,nstates,l_blk
     nbegin_ip=iv
     nend_ip=min(nbegin_ip+l_blk-1,nstates)
     nsize_ip=nend_ip-nbegin_ip+1



!first step
     psi_1(1:npw,1:nsize_ip)=zstates(1:npw,nbegin_ip:nend_ip)

!for h_psi allocations are required

     
     write(stdout,*) 'lanczos_state:', istate,ispin
     FLUSH(stdout)

     
!calculate H|\phi_i>
     call h_psi( npw, npw, nsize_ip,psi_1(1,1), u_0 )
     if(l_selfconsistent) call h_psi_self( npw, npw, nsize_ip,psi_1(1,1), u_0 )
     if(l_scissor) call h_psi_scissor(ispin, npw, npw, nsize_ip,psi_1(1,1), u_0 )
!calculate n_1
     n_1(1:nsize_ip)=0.d0
     do is=1,nsize_ip
        do ig=1,npw
           n_1(is)=n_1(is)+2.d0*dble(conjg(u_0(ig,is))*u_0(ig,is))
        enddo
        if(gstart==2) n_1(is)=n_1(is)-dble(conjg(u_0(1,is))*u_0(1,is))
     enddo
     call mp_sum(n_1(1:nsize_ip), world_comm)
     n_1(1:nsize_ip)=dsqrt(n_1(1:nsize_ip))

!calculate alpha
     alpha(1:nsize_ip)=0.d0
     do is=1,nsize_ip
        do ig=1,npw
           alpha(is)=alpha(is)+2.d0*dble(conjg(psi_1(ig,is))*u_0(ig,is))
        enddo
        if(gstart==2) alpha(is)=alpha(is)-dble(conjg(psi_1(1,is))*u_0(1,is))
     enddo
     call mp_sum(alpha(1:nsize_ip), world_comm)
     alpha(1:nsize_ip)=alpha(1:nsize_ip)/n_1(1:nsize_ip)

!calculate psi_2 and beta
     do is=1,nsize_ip
        psi_2(:,is)=u_0(:,is)/n_1(is)-alpha(is)*psi_1(:,is)
     enddo
     beta(1:nsize_ip)=0.d0
     do is=1,nsize_ip
        do ig=1,npw
           beta(is)=beta(is)+2.d0*dble(conjg(psi_2(ig,is))*psi_2(ig,is))
        enddo
        if(gstart==2) beta(is)=beta(is)-dble(conjg(psi_2(1,is))*psi_2(1,is))
     enddo
     call mp_sum(beta(1:nsize_ip), world_comm)
     beta(1:nsize_ip)=dsqrt(beta(1:nsize_ip))
     do is=1,nsize_ip
        psi_2(:,is)=psi_2(:,is)/beta(is)
     enddo

!calculate d

     do is=1,nsize_ip
        do ig=1,npw
           d(1,is+nbegin_ip-1)=d(1,is+nbegin_ip-1)+2.d0*dble(conjg(psi_1(ig,is))*u_0(ig,is))
        enddo
        if(gstart==2) d(1,is+nbegin_ip-1)=d(1,is+nbegin_ip-1)-dble(conjg(psi_1(1,is))*u_0(1,is))
     enddo
     call mp_sum(d(1,nbegin_ip:nend_ip), world_comm)

     if(l_verbose) write(stdout,*) 'Lanczos Diagonal 1', d(1,nbegin_ip:nend_ip)
     FLUSH(stdout)
  
!calculate f


     do is=1,nsize_ip
        do ig=1,npw
           f(1,is+nbegin_ip-1)=f(1,is+nbegin_ip-1)+2.d0*dble(conjg(psi_2(ig,is))*u_0(ig,is))
        enddo
        if(gstart==2) f(1,is+nbegin_ip-1)=f(1,is+nbegin_ip-1)-dble(conjg(psi_2(1,is))*u_0(1,is))
     enddo
     call mp_sum(f(1,nbegin_ip:nend_ip), world_comm)

     if(l_verbose) write(stdout,*) 'ATTENZIONE1'
     FLUSH(stdout)

!calculate overlaps and write on output file
     if(omat_div) then
        call dgemm('T','N',nstates,nsize_ip,2*npw,2.d0,zstates,2*npw,psi_1(1,1),2*npw,0.d0,omat,nstates)
        if(gstart==2) then
           do ii=1,nsize_ip
              do jj=1,nstates
                 omat(jj,ii)=omat(jj,ii)-dble(conjg(zstates(1,jj))*psi_1(1,ii))
              enddo
           enddo
        endif
        do ii=1,nsize_ip
           call mp_sum(omat(:,ii), world_comm)
        enddo
        
        if(ip==mpime) omat_tot(1:nstates,1:nsize_ip,1)=omat(1:nstates,1:nsize_ip)
    
     else
        omat(:,:)=0.d0
        call dgemm('T','N',nstates,nstates,2*npw,2.d0,zstates,2*npw,psi_1,2*npw,0.d0,omat,nstates)
        if(gstart==2) then
           do ii=1,nstates
              do jj=1,nstates
                 omat(jj,ii)=omat(jj,ii)-dble(conjg(zstates(1,jj))*psi_1(1,ii))
              enddo
           enddo
        endif
        if(l_verbose) write(stdout,*) 'ATTENZIONE2'
        FLUSH(stdout)

        do ii=1,nstates
           call mp_sum(omat(:,ii), world_comm)
        enddo

        if(ionode) then
           iunlan=find_free_unit()
           write(nfile,'(4i1)') istate/1000,mod(istate,1000)/100,mod(istate,100)/10,mod(istate,10)
           if(ispin==1) then
              if(itype==0) then
                 open( unit= iunlan, file=trim(tmp_dir)//trim(prefix)//'.p_iter_lanczos',&
                      &status='unknown',form='unformatted')
              else if(itype==1) then
                 open( unit= iunlan, file=trim(tmp_dir)//trim(prefix)//'.s_iter_lanczos'//'_'//nfile,&
                      &status='unknown',form='unformatted')
              else
                 open( unit= iunlan, file=trim(tmp_dir)//trim(prefix)//'.o_iter_lanczos',&
                      &status='unknown',form='unformatted')
              endif
           else 
              if(itype==0) then
                 open( unit= iunlan, file=trim(tmp_dir)//trim(prefix)//'.p_iter_lanczos2', &
                      &status='unknown',form='unformatted')
              else if(itype==1) then
                 open( unit= iunlan, file=trim(tmp_dir)//trim(prefix)//'.s_iter_lanczos2'//'_'//nfile, &
                      &status='unknown',form='unformatted')
              else
                 open( unit= iunlan, file=trim(tmp_dir)//trim(prefix)//'.o_iter_lanczos2', &
                      &status='unknown',form='unformatted')
              endif
           endif
           write(iunlan) nstates
           write(iunlan) istate
           write(iunlan) nsteps
           do is=1,nstates
              write(iunlan) omat(1:nstates,is)
           enddo
        endif
     endif
  
     if(l_verbose) write(stdout,*) 'ATTENZIONE3'
     FLUSH(stdout)


!calculate second overlap
     if(omat_div) then
    
        
        call dgemm('T','N',nstates,nsize_ip,2*npw,2.d0,zstates,2*npw,psi_2(1,1),2*npw,0.d0,omat,nstates)
        if(gstart==2) then
           do ii=1,nsize_ip
              do jj=1,nstates
                 omat(jj,ii)=omat(jj,ii)-dble(conjg(zstates(1,jj))*psi_2(1,ii))
              enddo
           enddo
        endif
        do ii=1,nsize_ip
           call mp_sum(omat(:,ii), world_comm)
        enddo
           
        if(ip==mpime) omat_tot(1:nstates,1:nsize_ip,2)=omat(1:nstates,1:nsize_ip)
     else
        omat(:,:)=0.d0
        call dgemm('T','N',nstates,nstates,2*npw,2.d0,zstates,2*npw,psi_2,2*npw,0.d0,omat,nstates)
        if(gstart==2) then
           do ii=1,nstates
              do jj=1,nstates
                 omat(jj,ii)=omat(jj,ii)-dble(conjg(zstates(1,jj))*psi_2(1,ii))
              enddo
           enddo
        endif
        if(l_verbose) write(stdout,*) 'ATTENZIONE4'
        FLUSH(stdout)

        do ii=1,nstates
           call mp_sum(omat(:,ii), world_comm)
        enddo
        if(ionode) then
           do is=1,nstates
              write(iunlan) omat(1:nstates,is)
           enddo
        endif
     endif
     
!do iterate
     do it=2,nsteps
        if(l_verbose) write(stdout,*) 'lanczos h_psi'
        FLUSH(stdout)

!calculate H|\phi_i+1>
        call h_psi( npw, npw, nsize_ip,psi_2(1,1), u_1 )
        if(l_selfconsistent) call h_psi_self( npw, npw, nsize_ip,psi_2(1,1), u_1 )
        if(l_scissor) call h_psi_scissor( ispin,npw, npw, nsize_ip,psi_2(1,1), u_1 )

        if(l_verbose) write(stdout,*) 'lanczos alfa beta gamma'
        FLUSH(stdout)




!calculate n_1
        n_1(1:nsize_ip)=0.d0
        do is=1,nsize_ip
           do ig=1,npw
              n_1(is)=n_1(is)+2.d0*dble(conjg(u_1(ig,is))*u_1(ig,is))
           enddo
           if(gstart==2) n_1(is)=n_1(is)-dble(conjg(u_1(1,is))*u_1(1,is))
        enddo
        call mp_sum(n_1(1:nsize_ip), world_comm)
        n_1(1:nsize_ip)=dsqrt(n_1(1:nsize_ip))

!calculate alpha
        alpha(1:nsize_ip)=0.d0
        do is=1,nsize_ip
           do ig=1,npw
              alpha(is)=alpha(is)+2.d0*dble(conjg(psi_1(ig,is))*u_1(ig,is))
           enddo
           if(gstart==2) alpha(is)=alpha(is)-dble(conjg(psi_1(1,is))*u_1(1,is))
        enddo
        call mp_sum(alpha(1:nsize_ip), world_comm)
        alpha(1:nsize_ip)=alpha(1:nsize_ip)/n_1(1:nsize_ip)

!calculate beta
        beta(1:nsize_ip)=0.d0
        do is=1,nsize_ip
           do ig=1,npw
              beta(is)=beta(is)+2.d0*dble(conjg(psi_2(ig,is))*u_1(ig,is))
           enddo
           if(gstart==2) beta(is)=beta(is)-dble(conjg(psi_2(1,is))*u_1(1,is))
        enddo
        call mp_sum(beta(1:nsize_ip), world_comm)
        beta(1:nsize_ip)=beta(1:nsize_ip)/n_1(1:nsize_ip)

!calculate psi_3 and gamma
        do is=1,nsize_ip
           psi_3(:,is)=u_1(:,is)/n_1(is)-alpha(is)*psi_1(:,is)-beta(is)*psi_2(:,is)
        enddo
        gamma(1:nsize_ip)=0.d0
        do is=1,nsize_ip
           do ig=1,npw
              gamma(is)=gamma(is)+2.d0*dble(conjg(psi_3(ig,is))*psi_3(ig,is))
           enddo
           if(gstart==2) gamma(is)=gamma(is)-dble(conjg(psi_3(1,is))*psi_3(1,is))
        enddo
        call mp_sum(gamma(1:nsize_ip), world_comm)
        gamma(1:nsize_ip)=dsqrt(gamma(1:nsize_ip))
        do is=1,nsize_ip
           psi_3(:,is)=psi_3(:,is)/gamma(is)
        enddo

        if(l_verbose) write(stdout,*) 'lanczos d f omat'
        FLUSH(stdout)


!calculate d
        do is=1,nsize_ip
           do ig=1,npw
              d(it,is+nbegin_ip-1)=d(it,is+nbegin_ip-1)+2.d0*dble(conjg(psi_2(ig,is))*u_1(ig,is))
           enddo
           if(gstart==2) d(it,is+nbegin_ip-1)=d(it,is+nbegin_ip-1)-dble(conjg(psi_2(1,is))*u_1(1,is))
        enddo
        call mp_sum(d(it,nbegin_ip:nend_ip), world_comm)

!calculate f
        do is=1,nsize_ip
           do ig=1,npw
              f(it,is+nbegin_ip-1)=f(it,is+nbegin_ip-1)+2.d0*dble(conjg(psi_3(ig,is))*u_1(ig,is))
           enddo
           if(gstart==2) f(it,is+nbegin_ip-1)=f(it,is+nbegin_ip-1)-dble(conjg(psi_3(1,is))*u_1(1,is))
        enddo
        call mp_sum(f(it,nbegin_ip:nend_ip), world_comm)

!calculate overlap
        if(it /=nsteps) then
           if(omat_div) then
              call dgemm('T','N',nstates,nsize_ip,2*npw,2.d0,zstates,2*npw,psi_3(1,1),2*npw,0.d0,omat,nstates)
              if(gstart==2) then
                 do ii=1,nsize_ip
                    do jj=1,nstates
                       omat(jj,ii)=omat(jj,ii)-dble(conjg(zstates(1,jj))*psi_3(1,ii))
                    enddo
                 enddo
              endif
              do ii=1,nsize_ip
                 call mp_sum(omat(:,ii), world_comm)
              enddo
              if(ip==mpime) omat_tot(1:nstates,1:nsize_ip,it+1)=omat(1:nstates,1:nsize_ip)

           else
              omat(:,:)=0.d0
              call dgemm('T','N',nstates,nstates,2*npw,2.d0,zstates,2*npw,psi_3,2*npw,0.d0,omat,nstates)
              if(gstart==2) then
                 do ii=1,nstates
                    do jj=1,nstates
                       omat(jj,ii)=omat(jj,ii)-dble(conjg(zstates(1,jj))*psi_3(1,ii))
                    enddo
                 enddo
              endif
              do ii=1,nstates
                 call mp_sum(omat(:,ii), world_comm)
              enddo
              if(ionode) then
                 do is=1,nstates
                    write(iunlan) omat(1:nstates,is)
                 enddo
              endif
           endif
        endif
!update arrays
        psi_1(1:npw,1:nsize_ip)=psi_2(1:npw,1:nsize_ip)
        psi_2(1:npw,1:nsize_ip)=psi_3(1:npw,1:nsize_ip)
        u_0(1:npw,1:nsize_ip)=u_1(1:npw,1:nsize_ip)
        
     enddo
     ip=ip+1
  enddo

!if omat is distribute here writes onfile

  if(omat_div) then
      if(ionode) then
         iunlan=find_free_unit()
         write(nfile,'(4i1)') istate/1000,mod(istate,1000)/100,mod(istate,100)/10,mod(istate,10)
         if(ispin==1) then
            if(itype==0) then
               open( unit= iunlan, file=trim(tmp_dir)//trim(prefix)//'.p_iter_lanczos', &
                    &status='unknown',form='unformatted')
            else
               open( unit= iunlan, file=trim(tmp_dir)//trim(prefix)//'.s_iter_lanczos'//'_'//nfile, &
                    &status='unknown',form='unformatted')
            endif
         else
            if(itype==0) then
               open( unit= iunlan, file=trim(tmp_dir)//trim(prefix)//'.p_iter_lanczos2', &
                    &status='unknown',form='unformatted')
            else
               open( unit= iunlan, file=trim(tmp_dir)//trim(prefix)//'.s_iter_lanczos2'//'_'//nfile, &
                    &status='unknown',form='unformatted')
            endif
         endif
         write(iunlan) nstates
         write(iunlan) istate
         write(iunlan) nsteps
      endif
      do it=1,nsteps
         do ip=0,nproc-1
            nbegin_ip=ip*l_blk+1
            nend_ip=min(nbegin_ip+l_blk-1,nstates)
            nsize_ip=nend_ip-nbegin_ip+1
            if(nsize_ip >=1) then
               if(mpime==ip) omat(1:nstates,1:nsize_ip)=omat_tot(1:nstates,1:nsize_ip,it)
               do is=1,nsize_ip
                  call mp_bcast(omat(1:nstates,is),ip,world_comm)
                  if(ionode) write(iunlan) omat(1:nstates,is)
               enddo
            endif
         enddo
      enddo

   endif
!write tridiagonal matrix on disk
  if(ionode) then
     do is=1,nstates
        write(iunlan) d(1:nsteps,is)
     enddo
     do is=1,nstates
        write(iunlan) f(1:nsteps,is)
     enddo
  endif

  if(ionode) close(iunlan)

  deallocate(psi_1,psi_2,psi_3)
  deallocate(u_0,u_1)
  deallocate(alpha,beta,gamma,n_1)
  deallocate(f,d,omat,c)

  deallocate(h_diag,s_diag)
  call deallocate_bec_type(becp)

  call stop_clock('lanczos_state')
  return
end subroutine lanczos_state


subroutine orthonormalize_two_manifolds( state1, n1,state2, n2, threshold, state_out, n_out)
!this subroutine form am orthormal basis set from 2 manifold (with orthonormal basis sets)
!ONLY FOR NORM_CONSERVING CASE

   USE io_global,            ONLY : stdout, ionode, ionode_id
   USE kinds,    ONLY : DP
   USE gvect
   USE wvfct,    ONLY : npwx, npw, nbnd
   USE mp, ONLY : mp_sum, mp_barrier, mp_bcast
   USE mp_world, ONLY : world_comm
   USE fft_base,             ONLY : dfftp, dffts
   USE fft_interfaces,       ONLY : fwfft, invfft
   USE wannier_gw, ONLY : l_verbose
   
  implicit none

  COMPLEX(kind=DP), INTENT(in)  :: state1(npw,n1)!1st orthonormal basis 
  INTEGER, INTENT(in)           :: n1!number of 1st basis elements
  COMPLEX(kind=DP), INTENT(in)  :: state2(npw,n2)!2nd orthonormal basis
  INTEGER, INTENT(in)           :: n2!number of 2nd basis elements
  REAL(kind=DP), INTENT(in)     :: threshold!threshold for orthonormality
  COMPLEX(kind=DP), INTENT(out) :: state_out(npw,n1+n2)!output basis set
  INTEGER, INTENT(out)          :: n_out!number of output states

  INTEGER :: ii,jj
  REAL(kind=DP), ALLOCATABLE :: omat(:,:),tmp_mat(:,:) 
  REAL(kind=DP), ALLOCATABLE :: eigen(:),work(:)
  INTEGER :: lwork,info,liwork
  REAL(kind=DP), ALLOCATABLE :: omat1(:,:),omat2(:,:)
 
!buid overlap matrix
  if(l_verbose) write(stdout,*) 'orthonormalize dgemm'
  FLUSH(stdout)

  allocate(omat(n1+n2,n1+n2))
  omat(:,:)=0.d0
  do ii=1,n1+n2
     omat(ii,ii)=1.d0
  enddo
  allocate(tmp_mat(n1,n2))
  call dgemm('T','N',n1,n2,2*npw,2.d0,state1,2*npw,state2,2*npw,0.d0,tmp_mat,n1)
  if(gstart==2) then
     do ii=1,n2
        do jj=1,n1
           tmp_mat(jj,ii)=tmp_mat(jj,ii)-dble(conjg(state1(1,jj))*state2(1,ii))
        enddo
     enddo
  endif
  if(l_verbose) write(stdout,*) 'orthonormalize mp_sum'
  FLUSH(stdout)

  do ii=1,n2
     call mp_sum(tmp_mat(:,ii), world_comm)
  enddo
  if(l_verbose) write(stdout,*) 'orthonormalize copy array'
  FLUSH(stdout)

  omat(1:n1,n1+1:n1+n2)=tmp_mat(1:n1,1:n2)
  deallocate(tmp_mat)

!diagonalize

   allocate(eigen(n1+n2))
   if(l_verbose) write(stdout,*) 'orthonormalize dsyev'
   FLUSH(stdout)
  
   if(ionode) then
      allocate(work(1))
      call  DSYEV( 'V', 'U', n1+n2, omat, n1+n2, eigen, work, -1, info )
      lwork=work(1)
      deallocate(work)
      allocate(work(lwork))
      call  DSYEV( 'V', 'U', n1+n2, omat, n1+n2, eigen, work, lwork, info )
      deallocate(work)
      if(info/=0) then
         write(stdout,*) 'ROUTINE orthonormalize_two_manifolds, INFO:', info
         stop
      endif
   else
      eigen(:)=0.d0
      omat(:,:)=0.d0
   endif
   if(l_verbose) write(stdout,*) 'orthonormalize mp_bcast now mp_sum'
   FLUSH(stdout)

   do ii=1,n1+n2
      !call mp_bcast(omat(:,ii), ionode_id,world_comm)
      call mp_sum(omat(:,ii), world_comm)
   enddo
   !call mp_bcast(eigen(:), ionode_id,world_comm)
   call mp_sum(eigen(:), world_comm)

   do ii=1,n1+n2
      if(l_verbose) write(stdout,*) 'EIGEN:',ii, eigen(ii)
   enddo
   FLUSH(stdout)
   if(l_verbose) write(stdout,*) 'orthonormalize copy'
  FLUSH(stdout)

!construct orthonormal basis set 
 !  state_out(:,:)=(0.d0,0.d0)

   n_out=0
   do ii=1,n1+n2
      if(eigen(ii) >= threshold) then
         n_out=n_out+1
      endif
   enddo
   allocate(omat1(n1,n_out),omat2(n2, n_out))

   do ii=1,n_out
      omat1(1:n1,ii)=omat(1:n1,n1+n2-n_out+ii)/dsqrt(eigen(n1+n2-n_out+ii))
   enddo
   do ii=1,n_out
      omat2(1:n2,ii)=omat(n1+1:n1+n2,n1+n2-n_out+ii)/dsqrt(eigen(n1+n2-n_out+ii))
   enddo

   call dgemm('N','N',2*npw,n_out,n1,1.d0,state1,2*npw,omat1,n1,0.d0,state_out,2*npw)
   call dgemm('N','N',2*npw,n_out,n2,1.d0,state2,2*npw,omat2,n2,1.d0,state_out,2*npw)

   deallocate(omat1,omat2)

!  n_out=0  
!   do ii=1,n1+n2
!      if(eigen(ii) >= threshold) then
!         n_out=n_out+1
!         do jj=1,n1
!            state_out(:,n_out)=state_out(:,n_out)+omat(jj,ii)*state1(:,jj)/dsqrt(eigen(ii))
!         enddo
!         do jj=1,n2
!            state_out(:,n_out)=state_out(:,n_out)+omat(jj+n1,ii)*state2(:,jj)/dsqrt(eigen(ii))
!         enddo
!      endif
!   enddo

   write(stdout,*) 'orthonormalize_two_manifolds: basis dimension:', n_out
   FLUSH(stdout)

  deallocate (omat)
  return
end subroutine orthonormalize_two_manifolds


subroutine global_pola_lanczos(nstates,nstates_eff,threshold,nglobal,nsteps,numpw,ispin,l_eigen)
!this subroutine from the orthonormal basis at each v 
!construct a global basis for the lanczos calculation of the
!polarization
  USE io_global,            ONLY : stdout, ionode, ionode_id
  USE io_files,             ONLY : prefix, tmp_dir, diropn
  USE kinds,    ONLY : DP
  USE wannier_gw, ONLY : num_nbndv,max_ngm,l_pmatrix
  USE gvect
  USE wvfct,    ONLY : npwx, npw, nbnd
  USE mp, ONLY : mp_sum, mp_barrier, mp_bcast
  USE mp_world, ONLY : world_comm
  USE wavefunctions_module, ONLY : evc, psic
  USE gvect
  USE gvecs,              ONLY : nls, nlsm, doublegrid
  USE fft_base,             ONLY : dfftp, dffts
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE wannier_gw, ONLY : l_verbose
  USE klist, ONLY : igk_k

  implicit none
  
  INTEGER, EXTERNAL :: find_free_unit
  
  INTEGER, INTENT(in)      :: nstates!number of orthonormal states for each v
  INTEGER, INTENT(in)      :: nstates_eff!effective number of orthonormal states for each v
  REAL(kind=DP),INTENT(in) :: threshold!threshold for orthonormalization algorithm
  INTEGER, INTENT(out)     :: nglobal!total number of final orthonormal states
  INTEGER, INTENT(in)      :: nsteps!number of lanczos steps                 
  INTEGER, INTENT(in)      :: numpw!number of wannier products for testing
  INTEGER, INTENT(in)      :: ispin!spin channel 1,2
  LOGICAL, INTENT(in)      :: l_eigen!if true partial t states are scaled with the corresponding eigenvalue

  INTEGER :: iunv,iuntmat
  LOGICAL :: exst
  INTEGER :: ii,jj,iv,ic
  COMPLEX(kind=DP), ALLOCATABLE :: old_basis(:,:), new_basis(:,:),v_basis(:,:)
  INTEGER :: nglobal_old
  REAL(kind=DP), ALLOCATABLE :: t_mat(:,:)
  CHARACTER(4) :: nfile
!for test:
  REAL(kind=DP) :: sca,sca1
  INTEGER :: iungprod,ig,iw
  REAL(kind=DP), ALLOCATABLE :: wv_real(:),tmp_r(:)
  COMPLEX(kind=DP), ALLOCATABLE :: tmp_g(:),wp_prod(:)

  LOGICAL :: l_test=.false.
  REAL(kind=DP)::proj_tot

  INTEGER :: nbuffer,ndelta!for avoiding nested allocation/deallocation cycles
  LOGICAL :: l_update_memory

  INTEGER, PARAMETER :: offset=0!ATTENZIONE THEN PUT 0!!!!!!!
  REAL(kind=DP), ALLOCATABLE :: eigen(:)
  INTEGER :: idumm

  if(num_nbndv(ispin) == 0) return

  nbuffer=6*numpw
  ndelta=numpw

!set first basis from first valence state
!if required read eigenvectors too
  allocate(eigen(nstates))
  if(l_eigen) then
     if(ionode) then
        iv=1
        iuntmat = find_free_unit()
        write(nfile,'(4i1)') iv/1000,mod(iv,1000)/100,mod(iv,100)/10,mod(iv,10)
        if(ispin==1) then
           open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.p_eig_lanczos'//nfile, &
                &status='old',form='unformatted')
        else
           open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.p_eig_lanczos2'//nfile, &
                &status='old',form='unformatted')
        endif
        read(iuntmat) idumm
        read(iuntmat) eigen(1:nstates)
        close(iuntmat)
     endif
     call mp_bcast(eigen, ionode_id,world_comm)
  endif

  allocate(old_basis(npw,nbuffer))
  iunv = find_free_unit()
  CALL diropn( iunv, 'vw_lanczos',npw*2, exst)
  if(.not.l_eigen) then
     do ii=1,nstates_eff
        call davcio(old_basis(:,ii),npw*2,iunv,ii+offset,-1)
     enddo
     nglobal=nstates_eff
  else
     nglobal=1
     call davcio(old_basis(:,nglobal),npw*2,iunv,1+offset,-1)
     do ii=2,nstates_eff
        if(eigen(ii) > threshold) then
           nglobal=nglobal+1
           call davcio(old_basis(:,nglobal),npw*2,iunv,ii+offset,-1)
        endif
     enddo
  endif
!loop on valence states (from 2nd)
  allocate(v_basis(npw,nstates_eff))
  allocate(new_basis(npw,nbuffer))
 

  do iv=2,num_nbndv(ispin)

!!read from disk
     do ii=1,nstates_eff
        call davcio(v_basis(:,ii),npw*2,iunv,ii+offset+(iv-1)*(nstates+offset),-1)
     enddo


     if(l_eigen) then
        if(ionode) then
           iuntmat = find_free_unit()
           write(nfile,'(4i1)') iv/1000,mod(iv,1000)/100,mod(iv,100)/10,mod(iv,10)
           if(ispin==1) then
              open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.p_eig_lanczos'//nfile, &
                   &status='old',form='unformatted')
           else
              open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.p_eig_lanczos2'//nfile, &
                   &status='old',form='unformatted')
           endif
           read(iuntmat) idumm
           read(iuntmat) eigen(1:nstates)
           close(iuntmat)
        endif
        call mp_bcast(eigen, ionode_id,world_comm)

     endif



     if(nglobal+nstates_eff >nbuffer) then
        deallocate(new_basis)
        allocate(new_basis(npw,nbuffer+ndelta))
        l_update_memory=.true.
     else
        l_update_memory=.false.
     endif

!!calculate basis
     nglobal_old=nglobal
     if(l_verbose) write(stdout,*) 'Call orthonormalize_two_manifolds'
     FLUSH(stdout)
     if(.not.l_pmatrix) then
        !call orthonormalize_two_manifolds( old_basis, nglobal_old,v_basis, nstates, threshold, new_basis, nglobal)
        call orthonormalize_two_manifolds_prj( old_basis, nglobal_old,v_basis, nstates_eff, threshold, new_basis, nglobal,&
                        l_eigen,eigen)
     else
        call orthonormalize_two_manifolds_scalapack(old_basis, nglobal_old,v_basis, nstates_eff, threshold, new_basis, nglobal) 
     endif
     if(l_verbose) write(stdout,*) 'Done orthonormalize_two_manifolds',ispin
     FLUSH(stdout)


!!set arrays for next iteration
     if(l_update_memory) then
        deallocate(old_basis)
        allocate(old_basis(npw,nbuffer+ndelta))
        nbuffer=nbuffer+ndelta
     endif
     old_basis(:,1:nglobal)=new_basis(:,1:nglobal)
     !deallocate(new_basis)

  enddo
  deallocate(new_basis)
  write(stdout,*) 'TOTAL NUMBER OF GLOBAL T VECTORS: ', nglobal
!call lanczos chain routine
  call lanczos_state(old_basis, nglobal, 0, nsteps,1,ispin)

!calculate matrix element and write on disk  
  allocate(t_mat(nglobal,nstates_eff))
  do iv=1,num_nbndv(ispin)
     write(nfile,'(4i1)') iv/1000,mod(iv,1000)/100,mod(iv,100)/10,mod(iv,10)
     if(ionode) then
        iuntmat=find_free_unit()
        if(ispin==1) then
           open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.pt_mat_lanczos'//nfile, &
                &status='unknown',form='unformatted')
        else
           open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.pt_mat_lanczos2'//nfile, &
                &status='unknown',form='unformatted')
        endif
     endif
     do ii=1,nstates_eff
        call davcio(v_basis(:,ii),npw*2,iunv,ii+offset+(iv-1)*(nstates+offset),-1)
     enddo
     call dgemm('T','N',nglobal,nstates_eff,2*npw,2.d0,old_basis,2*npw,v_basis,2*npw,0.d0,t_mat,nglobal)
     if(gstart==2) then
        do ii=1,nstates_eff
           do jj=1,nglobal
              t_mat(jj,ii)=t_mat(jj,ii)-dble(conjg(old_basis(1,jj))*v_basis(1,ii))
           enddo
        enddo
     endif
     call mp_sum(t_mat(:,:), world_comm)
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

     call mp_sum(t_mat(:,:), world_comm)
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
     do iv=1,num_nbndv(ispin),num_nbndv(ispin)-1!ATTENZIONE
!put iv on real space
        psic(:)=(0.d0,0.d0)
        psic(nls(igk_k(1:npw,1)))  = evc(1:npw,iv)
        psic(nlsm(igk_k(1:npw,1))) = CONJG( evc(1:npw,iv) )
        CALL invfft ('Wave', psic, dffts)
        wv_real(:)= DBLE(psic(:))
!loop on wannier_products
        do iw=1,numpw
           call davcio(tmp_g,max_ngm*2,iungprod,iw,-1)
     !trasform to r-space
           psic(:)=(0.d0,0.d0)
           do ig=1,max_ngm
              psic(nl(ig))=tmp_g(ig)
              psic(nlm(ig))=CONJG(tmp_g(ig))
           enddo
           CALL invfft ('Dense', psic, dfftp)
           tmp_r(:)=dble(psic(:))
!!form products with w_v and trasfrom in G space
           psic(:)=cmplx(tmp_r(:)*wv_real(:),0.d0)
           CALL fwfft ('Wave', psic, dffts)
           wp_prod(1:npw) = psic(nls(igk_k(1:npw,1)))
!!project on conduction subspace
           call pc_operator(wp_prod(:),ispin, .false.)
!!do scalar product
           call dgemm('T','N',nglobal,1,2*npw,2.d0,old_basis,2*npw,wp_prod,2*npw,0.d0,t_mat,nglobal)
           if(gstart==2) then
              do jj=1,nglobal
                 t_mat(jj,1)=t_mat(jj,1)-dble(conjg(old_basis(1,jj))*wp_prod(1))
              enddo
           endif
           call mp_sum(t_mat(:,1), world_comm)
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
           call mp_sum(sca1, world_comm)
         
           write(stdout,*) 'Projection',iv,iw,sca/sca1
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
!END OF TESTING PART
  endif
  deallocate(old_basis)
  deallocate(v_basis)
  deallocate(eigen)

  return
end subroutine global_pola_lanczos

subroutine orthonormalize_two_manifolds_scalapack( state1, n1,state2, n2, threshold, state_out, n_out)
!this subroutine form am orthormal basis set from 2 manifold (with orthonormal basis sets)
!ONLY FOR NORM_CONSERVING CASE

   USE io_global,            ONLY : stdout, ionode, ionode_id
   USE kinds,    ONLY : DP
   USE gvect
   USE wvfct,    ONLY : npwx, npw, nbnd
   USE mp, ONLY : mp_sum, mp_barrier, mp_bcast
   USE mp_world, ONLY : world_comm
   USE wannier_gw, ONLY : p_mpime,p_nproc, npcol, nprow,icontxt,myrow,mycol
   
  implicit none

  COMPLEX(kind=DP), INTENT(in)  :: state1(npw,n1)!1st orthonormal basis 
  INTEGER, INTENT(in)           :: n1!number of 1st basis elements
  COMPLEX(kind=DP), INTENT(in)  :: state2(npw,n2)!2nd orthonormal basis
  INTEGER, INTENT(in)           :: n2!number of 2nd basis elements
  REAL(kind=DP), INTENT(in)     :: threshold!threshold for orthonormality
  COMPLEX(kind=DP), INTENT(out) :: state_out(npw,n1+n2)!output basis set
  INTEGER, INTENT(out)          :: n_out!number of output states

  INTEGER :: ii,jj
  REAL(kind=DP), ALLOCATABLE :: omat(:,:),tmp_mat(:,:) 
  REAL(kind=DP), ALLOCATABLE :: eigen(:),work(:)
  INTEGER :: lwork,info,liwork
  INTEGER, ALLOCATABLE :: iwork(:)
  INTEGER :: n, n_r,n_c,n_dimr,n_dimc,n1_r,n1_dimr,n2_c,n2_dimc
  INTEGER :: icrow, iccol, iproc,ilrow,ilcol
  INTEGER, EXTERNAL :: indxg2l,indxg2p
  REAL(kind=DP), EXTERNAL :: ddot
  REAL(kind=DP) :: sca
  INTEGER :: desc_a(9),desc_b(9)
 
#if defined(__SCALAPACK)
!buid overlap matrix
  n=n1+n2
  n_r=ceiling(real(n)/real(max(nprow,npcol)))
  n_c=ceiling(real(n)/real(max(nprow,npcol)))
  n_dimr=ceiling (real(n)/real(n_r*nprow))*n_r
  n_dimc=ceiling (real(n)/real(n_c*npcol))*n_c

  n1_r=ceiling(real(n1)/real(max(nprow,npcol)))
  n1_dimr=ceiling (real(n1)/real(n1_r*nprow))*n1_r

  n2_c=ceiling(real(n2)/real(max(nprow,npcol)))
  n2_dimc=ceiling (real(n2)/real(n2_c*npcol))*n2_c

  allocate(omat(n_dimr,n_dimc))
  omat(:,:)=0.d0

  do ii=1,n
     icrow = indxg2p(ii,n_r,0,0,nprow)
     iccol = indxg2p(ii,n_c,0,0,npcol)
     iproc=icrow*npcol+iccol  
     if(myrow==icrow .and. mycol==iccol) then
        ilrow=indxg2l(ii,n_r,0,0,nprow)
        ilcol=indxg2l(ii,n_c,0,0,npcol)
        omat(ilrow,ilcol)=1.d0
     endif
  enddo

  
  write(stdout,*) 'orthonormalize para1'
  FLUSH(stdout)

  do ii=1,n1
     do jj=1,n2
        sca=2.d0*ddot(2*npw,state1(:,ii),1,state2(:,jj),1) 
        if(gstart==2) sca=sca-dble(conjg(state1(1,ii))*state2(1,jj))
        call mp_sum(sca, world_comm)
        icrow = indxg2p(ii,n_r,0,0,nprow)
        iccol = indxg2p(jj+n1,n_c,0,0,npcol)
        iproc=icrow*npcol+iccol
        if(myrow==icrow .and. mycol==iccol) then
           ilrow=indxg2l(ii,n_r,0,0,nprow)
           ilcol=indxg2l(jj+n1,n_c,0,0,npcol)
           omat(ilrow,ilcol)=sca
        endif
     enddo
  enddo

  allocate(tmp_mat(n_dimr,n_dimc))

! A = omat                                                                                                                    
  desc_a(1)=1
  desc_a(2)=icontxt
  desc_a(3)=n
  desc_a(4)=n
  desc_a(5)=n_r
  desc_a(6)=n_c
  desc_a(7)=0
  desc_a(8)=0
  desc_a(9)=n_dimr

!B = tmp_mat

  desc_b(1)=1
  desc_b(2)=icontxt
  desc_b(3)=n
  desc_b(4)=n
  desc_b(5)=n_r
  desc_b(6)=n_c
  desc_b(7)=0
  desc_b(8)=0
  desc_b(9)=n_dimr


!diagonalize 

  allocate(work(1))
  allocate(eigen(n))

  write(stdout,*) 'orthonormalize para2'
  FLUSH(stdout)

  call pdsyev('V','U',n,omat,1,1,desc_a,eigen,tmp_mat,1,1,desc_b,work,-1,info)
  lwork=work(1)
  deallocate(work)
  allocate(work(lwork))
  call pdsyev('V','U',n,omat,1,1,desc_a,eigen,tmp_mat,1,1,desc_b,work,lwork,info)
  deallocate(work)
  if(info/=0) then
     write(stdout,*) 'ROUTINE orthonormalize_two_manifolds_scalapack, INFO:', info
     stop
  endif
   
  write(stdout,*) 'orthonormalize para3'
  FLUSH(stdout)

   do ii=1,n,n
      write(stdout,*) 'EIGEN:',ii, eigen(ii)
   enddo
   FLUSH(stdout)


   state_out(:,:)=(0.d0,0.d0)
   n_out=0
   do ii=1,n
      if(eigen(ii) >= threshold) then
         n_out=n_out+1
         do jj=1,n1
            icrow = indxg2p(jj,n_r,0,0,nprow)
            iccol = indxg2p(ii,n_c,0,0,npcol)
            iproc=icrow*npcol+iccol
            if(myrow==icrow .and. mycol==iccol) then
               ilrow=indxg2l(jj,n_r,0,0,nprow)
               ilcol=indxg2l(ii,n_c,0,0,npcol)
               sca=tmp_mat(ilrow,ilcol)
            endif
            call mp_bcast(sca, iproc,world_comm)
            state_out(:,n_out)=state_out(:,n_out)+sca*state1(:,jj)/dsqrt(eigen(ii))

         enddo
         do jj=1,n2
            icrow = indxg2p(jj+n1,n_r,0,0,nprow)
            iccol = indxg2p(ii,n_c,0,0,npcol)
            iproc=icrow*npcol+iccol
            if(myrow==icrow .and. mycol==iccol) then
               ilrow=indxg2l(jj+n1,n_r,0,0,nprow)
               ilcol=indxg2l(ii,n_c,0,0,npcol)
               sca=tmp_mat(ilrow,ilcol)
            endif
            call mp_bcast(sca, iproc,world_comm)
            state_out(:,n_out)=state_out(:,n_out)+sca*state2(:,jj)/dsqrt(eigen(ii))
         enddo
      endif
   enddo
   write(stdout,*) 'orthonormalize para4'
   FLUSH(stdout)

   write(stdout,*) 'orthonormalize_two_manifolds: basis dimension:', n_out
   FLUSH(stdout)

  deallocate (omat,tmp_mat,eigen)
#endif
  return
end subroutine orthonormalize_two_manifolds_scalapack

subroutine orthonormalize_two_manifolds_prj( state1, n1,state2, n2, threshold, state_out, n_out,l_w,weight)
!this subroutine form am orthormal basis set from 2 manifold (with orthonormal basis sets)
!ONLY FOR NORM_CONSERVING CASE
!first projects out of the second manifold the first one
!the orthonormalizes the second manifold

   USE io_global,            ONLY : stdout, ionode, ionode_id
   USE kinds,    ONLY : DP
   USE gvect
   USE wvfct,    ONLY : npwx, npw, nbnd
   USE mp, ONLY : mp_sum, mp_barrier, mp_bcast
   USE mp_world, ONLY : world_comm
   USE wannier_gw, ONLY : l_verbose

  implicit none

  COMPLEX(kind=DP), INTENT(in)  :: state1(npw,n1)!1st orthonormal basis
  INTEGER, INTENT(in)           :: n1!number of 1st basis elements
  COMPLEX(kind=DP), INTENT(inout)  :: state2(npw,n2)!2nd orthonormal basis
  INTEGER, INTENT(in)           :: n2!number of 2nd basis elements
  REAL(kind=DP), INTENT(in)     :: threshold!threshold for orthonormality
  COMPLEX(kind=DP), INTENT(out) :: state_out(npw,n1+n2)!output basis set
  INTEGER, INTENT(out)          :: n_out!number of output states
  LOGICAL, INTENT(in)           :: l_w!if true considere the weigth in eig
  REAL(kind=DP), INTENT(in)     :: weight(n2)!weigths


  INTEGER :: ii,jj
  REAL(kind=DP), ALLOCATABLE :: omat(:,:),tmp_mat(:,:)
  REAL(kind=DP), ALLOCATABLE :: eigen(:),work(:)
  INTEGER :: lwork,info,liwork
  REAL(kind=DP), ALLOCATABLE :: omat1(:,:),omat2(:,:)



!buid overlap matrix
  if(l_verbose) write(stdout,*) 'orthonormalize dgemm'
  FLUSH(stdout)

  allocate(tmp_mat(n1,n2))
  call dgemm('T','N',n1,n2,2*npw,2.d0,state1,2*npw,state2,2*npw,0.d0,tmp_mat,n1)
  if(gstart==2) then
     do ii=1,n2
        do jj=1,n1
           tmp_mat(jj,ii)=tmp_mat(jj,ii)-dble(conjg(state1(1,jj))*state2(1,ii))
        enddo
     enddo
  endif

  if(l_verbose) write(stdout,*) 'orthonormalize mp_sum'
  FLUSH(stdout)

  do ii=1,n2
     call mp_sum(tmp_mat(:,ii), world_comm)
  enddo

  call dgemm('N','N',2*npw, n2,n1,-1.d0,state1,2*npw,tmp_mat,n1,1.d0,state2,2*npw)

  deallocate(tmp_mat)

  if(l_w) then
     do ii=1,n2
        state2(1:npw,ii)=state2(1:npw,ii)*weight(ii)
     enddo
  endif


  allocate(omat(n2,n2))
 
  if(gstart==2) state2(1,1:n2)=dcmplx(dble(state2(1,1:n2)),0.d0)
  call dgemm('T','N',n2,n2,2*npw,2.d0,state2,2*npw,state2,2*npw,0.d0,omat,n2)
  if(gstart==2) then
     do ii=1,n2
        do jj=1,n2
           omat(jj,ii)=omat(jj,ii)-dble(conjg(state2(1,jj))*state2(1,ii))
        enddo
     enddo
  endif

  if(l_verbose) write(stdout,*) 'orthonormalize mp_sum'
  FLUSH(stdout)

  do ii=1,n2
     call mp_sum(omat(:,ii), world_comm)
  enddo

 
!diagonalize

   allocate(eigen(n2))
   if(l_verbose) write(stdout,*) 'orthonormalize dsyev'
   FLUSH(stdout)

   if(ionode) then
      allocate(work(1))
      if(l_w) omat(1:n2,1:n2)=omat(1:n2,1:n2)/weight(1)!to avoid numerical instabilities in DSYEV
      call  DSYEV( 'V', 'U', n2, omat,n2, eigen, work, -1, info )
      lwork=work(1)
      deallocate(work)
      allocate(work(lwork))
      call  DSYEV( 'V', 'U', n2, omat, n2, eigen, work, lwork, info )
      deallocate(work)
      if(info/=0) then
         write(stdout,*) 'ROUTINE orthonormalize_two_manifolds, INFO:', info
         stop
      endif
      if(l_w) eigen(1:n2)=eigen(1:n2)*weight(1)!to avoid numerical instabilities in DSYEV
   else
      eigen(:)=0.d0
      omat(:,:)=0.d0
   endif
   if(l_verbose) write(stdout,*) 'orthonormalize mp_bcast now mp_sum'
   FLUSH(stdout)

   do ii=1,n2
      !call mp_bcast(omat(:,ii), ionode_id,world_comm)
      call mp_sum(omat(:,ii), world_comm)
   enddo
   !call mp_bcast(eigen(:), ionode_id,world_comm)
   call mp_sum(eigen(:), world_comm)

!   do ii=1,n2
   do ii=1,n2,n2-1
      write(stdout,*) 'EIGEN GLOBAL:',ii, eigen(ii)
   enddo
   FLUSH(stdout)
   if(l_verbose) write(stdout,*) 'orthonormalize copy'
  FLUSH(stdout)

!construct orthonormal basis set
 !  state_out(:,:)=(0.d0,0.d0)

   n_out=0
   do ii=1,n2
      if(eigen(ii) >= threshold) then
         n_out=n_out+1
      endif
   enddo
   

  
   do ii=n2-n_out+1,n2
      omat(1:n2,ii)=omat(1:n2,ii)/dsqrt(eigen(ii))
   enddo

   
   call dgemm('N','N',2*npw,n_out,n2,1.d0,state2,2*npw,omat(:,n2-n_out+1:n2),n2,0.d0,state_out,2*npw)

   state_out(:,n_out+1:n_out+n1)=state1(:,1:n1)
   n_out=n_out+n1

   write(stdout,*) 'orthonormalize_two_manifolds: basis dimension:', n_out
   FLUSH(stdout)

  deallocate (omat,eigen)
  return
end subroutine orthonormalize_two_manifolds_prj

subroutine pc_operator_test(state)
!this operator project the wavefunction state on the conduction
!subspace, the valence wavefunction are in evc
!ONLY FOR GAMMA POINT NOW!!!!
   USE io_global,            ONLY : stdout
   USE kinds,    ONLY : DP
   USE gvect
   USE wvfct,    ONLY : npwx, npw, nbnd
   USE wavefunctions_module, ONLY : evc, psic
   USE mp, ONLY : mp_sum, mp_barrier, mp_bcast
   USE mp_world, ONLY : world_comm
   USE wannier_gw, ONLY : num_nbndv

  implicit none

  COMPLEX(kind=DP), INTENT(inout) :: state(npw)!state to be projected
  INTEGER :: iv,ig
  REAL(kind=DP), ALLOCATABLE :: prod(:)

  allocate(prod(nbnd-num_nbndv(1)))
  prod(:)=0.d0
  call dgemm('T','N', nbnd-num_nbndv(1),1,2*npw,2.d0,evc(:,num_nbndv(1)+1:nbnd),2*npwx,state,2*npw,0.d0,prod,nbnd-num_nbndv(1))
  do iv=num_nbndv(1)+1,nbnd
     if(gstart==2) prod(iv-num_nbndv(1))=prod(iv-num_nbndv(1))-dble(conjg(evc(1,iv))*state(1))
  enddo
  call mp_sum(prod(:), world_comm)
  call dgemm('N','N',2*npw,1,nbnd-num_nbndv(1),1.d0,evc(:,num_nbndv(1)+1:nbnd),2*npwx,prod,nbnd-num_nbndv(1),0.d0,state,2*npw)


  deallocate(prod)
  return
end subroutine pc_operator_test

subroutine pc_operator_t_m(numpw,state,evc_t,ispin,fc)
!this operator project the wavefunction state on the conduction
!subspace, the valence wavefunction are in evc
!it works on an arry of states
!ONLY FOR GAMMA POINT NOW!!!!
   USE io_global,            ONLY : stdout
   USE kinds,    ONLY : DP
   USE gvect
   USE wvfct,    ONLY : npwx, npw, nbnd
   USE wavefunctions_module, ONLY : evc, psic
   USE mp, ONLY : mp_sum, mp_barrier, mp_bcast
   USE mp_world, ONLY : world_comm
   USE wannier_gw, ONLY : num_nbndv
   USE fft_custom_gwl

  implicit none

  TYPE(fft_cus), INTENT(in) :: fc
  INTEGER, INTENT(in) :: numpw!number of vectors
  INTEGER, INTENT(in) :: ispin!spin channel
  COMPLEX(kind=DP), INTENT(inout) :: state(fc%npwt,numpw)!state to be projected
  ! COMPLEX(kind=DP), INTENT(inout) :: evc_t(fc%npwt,num_nbndv(ispin))
  ! above syntax not accepted by all compilers
  COMPLEX(kind=DP), INTENT(inout) :: evc_t(fc%npwt,*)!valence states

  INTEGER :: ii,iv,ig
  REAL(kind=DP), ALLOCATABLE :: prod(:,:)

  allocate(prod(num_nbndv(ispin),numpw))
  call dgemm('T','N', num_nbndv(ispin),numpw,2*fc%npwt,2.d0,evc_t,2*fc%npwt,state,2*fc%npwt,&
       & 0.d0,prod,num_nbndv(ispin))
  

 if(fc%gstart_t==2) then
    do ii=1,numpw
       do iv=1,num_nbndv(ispin)
          prod(iv,ii)=prod(iv,ii)-dble(conjg(evc_t(1,iv))*state(1,ii))
       enddo
    enddo
 endif
  do ii=1,numpw
     call mp_sum(prod(:,ii), world_comm)
  enddo
  call dgemm('N','N',2*fc%npwt,numpw,num_nbndv(ispin),-1.d0,evc_t,2*fc%npwt,prod,&
       &num_nbndv(ispin),1.d0,state,2*fc%npwt)
  

  deallocate(prod)
  return
end subroutine pc_operator_t_m

subroutine pc_operator_t_r(numpw,state,evc_r,ispin,fc)
!NOT_TO_BE_INCLUDED_START
!this operator project the wavefunction state on the conduction                                                                      
!subspace, the valence wavefunction are in evc

   USE io_global,            ONLY : stdout
   USE kinds,    ONLY : DP
   USE gvect
   USE wvfct,    ONLY : npwx, npw, nbnd
   USE wavefunctions_module, ONLY : evc, psic
   USE mp, ONLY : mp_sum, mp_barrier, mp_bcast
   USE mp_world, ONLY : world_comm
   USE wannier_gw, ONLY : num_nbndv
   USE fft_custom_gwl

  implicit none

  TYPE(fft_cus), INTENT(in) :: fc
  INTEGER, INTENT(in) :: numpw!number of vectors
  INTEGER, INTENT(in) :: ispin!spin channel
  REAL(kind=DP), INTENT(inout) :: state(fc%nrxxt,numpw)!state to be projected
  ! REAL(kind=DP), INTENT(inout) :: evc_r(fc%nrxxt,num_nbndv(ispin))
  ! above syntax not accepted by all compilers
  REAL(kind=DP), INTENT(inout) :: evc_r(fc%nrxxt,*)!valence states
  INTEGER :: ii,iv,ig
  REAL(kind=DP), ALLOCATABLE :: prod(:,:)

  allocate(prod(num_nbndv(ispin),numpw))
  call dgemm('T','N', num_nbndv(ispin),numpw,fc%nrxxt,1.d0,evc_r,fc%nrxxt,state,fc%nrxxt,&
       & 0.d0,prod,num_nbndv(ispin))

  do ii=1,numpw
     call mp_sum(prod(:,ii), world_comm)
     prod(:,ii)=prod(:,ii)/dble(fc%nr1t*fc%nr2t*fc%nr3t)
  enddo
  call dgemm('N','N',fc%nrxxt,numpw,num_nbndv(ispin),-1.d0,evc_r,fc%nrxxt,prod,&
       &num_nbndv(ispin),1.d0,state,fc%nrxxt)


  deallocate(prod)
  return
!NOT_TO_BE_INCLUDED_END
end subroutine pc_operator_t_r





subroutine  h_psi_self( lda, n, m, psi, hpsi )
!NOT_TO_BE_INCLUDED_START
!add to hpsi part dur to self-consistent GW calculation

  ! ... input:
  ! ...    lda   leading dimension of arrays psi, spsi, hpsi
  ! ...    n     true dimension of psi, spsi, hpsi
  ! ...    m     number of states psi
  ! ...    psi
  !
  ! ... output:
  ! ...    hpsi  H*psi
  !
  USE kinds,    ONLY : DP
  USE gvect,    ONLY : gstart
  USE wvfct,    ONLY : npwx, npw, nbnd,et
  USE wavefunctions_module, ONLY : evc
  USE wannier_gw, ONLY : n_gw_states, ene_gw, delta_self
  USE mp, ONLY : mp_sum
  USE mp_world, ONLY : world_comm
  !
  implicit none

  INTEGER, INTENT(IN)     :: lda, n, m
  COMPLEX(kind=DP), INTENT(IN)  :: psi(lda,m) 
  COMPLEX(kind=DP), INTENT(OUT) :: hpsi(lda,m)   

  INTEGER :: ii,jj
  REAL(kind=DP), ALLOCATABLE :: prod(:,:)

!apply \Delta1
  hpsi(1:n,1:m)=hpsi(1:n,1:m)+delta_self*psi(1:n,1:m)


!apply \Sum_i (e^GW_i-e^DFT_i-Delta)|\psi_i><\psi_i|

  allocate(prod(n_gw_states,m))
  prod(:,:)=0.d0
  call dgemm('T','N', n_gw_states,m,2*npw,2.d0,evc,2*npwx,psi,2*lda,0.d0,prod,n_gw_states)
  do ii=1,n_gw_states
     do jj=1,m
        if(gstart==2) prod(ii,jj)=prod(ii,jj)-dble(conjg(evc(1,ii))*psi(1,jj))
     enddo
  enddo
  call mp_sum(prod,world_comm)

  do jj=1,m
     do ii=1,n_gw_states
        prod(ii,jj)=prod(ii,jj)*(ene_gw(ii,1)-et(ii,1)-delta_self)
     enddo
  enddo
  call dgemm('N','N',2*npw,m,n_gw_states,1.d0,evc,2*npwx,prod,n_gw_states,1.d0,hpsi,2*lda)


  deallocate(prod)
  return
!NOT_TO_BE_INCLUDED_END
end subroutine h_psi_self


subroutine  h_psi_scissor( ispin,lda, n, m, psi, hpsi )
!NOT_TO_BE_INCLUDED_START                                    
!add to hpsi part dur to self-consistent GW calculation
  ! ... input:          
  ! ...    lda   leading dimension of arrays psi, spsi, hpsi
  ! ...    n     true dimension of psi, spsi, hpsi
  ! ...    m     number of states psi 
  ! ...    psi
  ! ... output:
  ! ...    hpsi  H*psi

  USE kinds,    ONLY : DP
  USE gvect,    ONLY : gstart
  USE wvfct,    ONLY : npwx, npw, nbnd,et
  USE wavefunctions_module, ONLY : evc
  USE wannier_gw, ONLY : num_nbndv,scissor
  USE mp, ONLY : mp_sum
  USE mp_world, ONLY : world_comm
  USE constants, ONLY : rytoev

  implicit none

  INTEGER, INTENT(in)  :: ispin!spin channel
  INTEGER, INTENT(IN)     :: lda, n, m
  COMPLEX(kind=DP), INTENT(IN)  :: psi(lda,m) 
  COMPLEX(kind=DP), INTENT(OUT) :: hpsi(lda,m)   

  INTEGER :: ii,jj
  REAL(kind=DP), ALLOCATABLE :: prod(:,:)



  allocate(prod(num_nbndv(ispin),m))
  prod=0.d0
  call dgemm('T','N', num_nbndv(ispin),m,2*npw,2.d0,evc,2*npwx,psi,2*lda,0.d0,prod,num_nbndv(ispin))
  do ii=1,num_nbndv(ispin)
     do jj=1,m
        if(gstart==2) prod(ii,jj)=prod(ii,jj)-dble(conjg(evc(1,ii))*psi(1,jj))
     enddo
  enddo
  call mp_sum(prod,world_comm)

  do jj=1,m
     do ii=1,num_nbndv(ispin)
        prod(ii,jj)=prod(ii,jj)*(scissor(1)-scissor(2))/rytoev
     enddo
  enddo
  call dgemm('N','N',2*npw,m,num_nbndv(ispin),1.d0,evc,2*npwx,prod,num_nbndv(ispin),&
       &1.d0+scissor(2)/rytoev,hpsi,2*lda)


  deallocate(prod)
  return
!NOT_TO_BE_INCLUDED_END
end subroutine h_psi_scissor


subroutine pola_basis_lanczos_real(n_set,nstates,numpw, nsteps,ispin)
!NOT_TO_BE_INCLUDED_START
!this subroutine calculates the basis for every v 
!the minimal orthonormal basis for the w_v(r)*w^P'_i(r) products

   USE io_global,            ONLY : stdout, ionode, ionode_id
   USE io_files,             ONLY : prefix, tmp_dir, diropn
   USE kinds,    ONLY : DP
   USE wannier_gw
   USE gvect
   USE constants, ONLY : e2, pi, tpi, fpi
   USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2
   USE wvfct,    ONLY : npwx, npw, nbnd
   USE gvecw,     ONLY : ecutwfc
   USE wavefunctions_module, ONLY : evc, psic
   USE mp, ONLY : mp_sum, mp_barrier, mp_bcast
   USE mp_pools, ONLY : intra_pool_comm
   USE mp_world, ONLY : world_comm, mpime, nproc
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
   INTEGER, INTENT(in) :: ispin! spin channel

   INTEGER :: iv,iw,ig,ii,jj,ir
   REAL(kind=DP), ALLOCATABLE :: wv_real(:),tmp_r(:),tmp_r2(:)
   COMPLEX(kind=DP), ALLOCATABLE :: tmp_g(:)
   REAL(kind=DP), ALLOCATABLE :: wp_prod(:,:,:)
   INTEGER :: iungprod,iunrprod, iungresult,iuntmat
   LOGICAL :: exst
   REAL(kind=DP), ALLOCATABLE :: omat(:,:),omat_hold(:,:)

   REAL(kind=DP), ALLOCATABLE :: eigen(:),work(:)
   INTEGER :: lwork,info,liwork
   COMPLEX(kind=DP), ALLOCATABLE :: wp_g(:,:),wp_g_t2(:,:)!product terms in g wfc grid
   REAL(kind=DP), ALLOCATABLE :: wp_g_t(:,:)!
   REAL(kind=DP), ALLOCATABLE :: t_mat(:,:),t_mat_hold(:,:), t_mat_hold2(:,:)
   CHARACTER(4) :: nfile
    COMPLEX(kind=DP), ALLOCATABLE :: p_basis(:,:)!polarizability basis
   LOGICAL :: l_dsyevr=.true.!if true uses dsyevr 
   REAL(kind=DP), ALLOCATABLE :: vectors(:,:)!for dsyevr
   INTEGER, ALLOCATABLE :: iwork(:), ifail(:)
   INTEGER, ALLOCATABLE :: isuppz(:)
   INTEGER :: n_found

   LOGICAL :: l_fft_custom=.false.!if true uses custom fft grid
   COMPLEX(kind=DP), ALLOCATABLE :: evc_t(:,:),p_basis_t(:,:)
   REAL(kind=DP), ALLOCATABLE :: evc_r(:,:)
   COMPLEX(kind=DP), ALLOCATABLE :: evc_g(:)
   LOGICAL :: l_sumrule=.false.!if true imposes the sum rule over the norm of Pc|\Phi_\mu\Psi_v> for each of them
   REAL(kind=DP), ALLOCATABLE :: norms(:)
   REAL(kind=DP) :: norm_t, c_norm,norm
   REAL(kind=DP), ALLOCATABLE :: p_basis_r(:,:) !polarizabilty basis in real custom space

   INTEGER :: ivv,nbuf
   REAL(kind=DP) :: vl,vu
   INTEGER :: il,iu
   REAL(kind=DP), ALLOCATABLE :: t_eigen_hold(:)
   REAL(kind=DP) :: sca

   TYPE(fft_cus) :: fc

   write(stdout,*) 'Routine pola_basis_lanczos_real'
   FLUSH(stdout)

   fc%ecutt=ecutwfc
   fc%dual_t=dual_vt

   if(l_verbose) write(stdout,*) 'Call initialize_fft_custom'
   FLUSH(stdout)
   call initialize_fft_custom(fc)
!   allocate(evc_g(fc%ngmt_g))
!   allocate(wv_real(dfftp%nnr))
!   allocate(wp_prod(npw,numpw))

   allocate(wv_real(fc%nrxxt))
!   allocate(wp_g_t2(fc%npwt))

   allocate(norms(numpw))

!read w^P'_i on file on real space
!open product of wanniers filed


   iungprod = find_free_unit()
   CALL diropn( iungprod, 'wiwjwfc_red', max_ngm*2, exst )



   iungresult = find_free_unit()
   CALL diropn( iungresult, 'vw_lanczos',npw*2, exst)

  
!read polarizability basis functions
   allocate(p_basis(max_ngm,numpw))
   do iw=1,numpw
      call davcio(p_basis(:,iw),max_ngm*2,iungprod,iw,-1)
   enddo
   
   close(iungprod)
  

   if(l_verbose) write(stdout,*) 'pola_basis_lanczos 1'
   FLUSH(stdout)
   
!now polarizability basis  are put on the ordering of the redueced grid, if required
  allocate(p_basis_t(fc%npwt,numpw))
!  if(fc%dual_t==4.d0) then
!     p_basis_t(:,:)=p_basis(:,:)
  !  else
  call reorderwfp (numpw,npw, fc%npwt,p_basis(:,:),p_basis_t(:,:), &
              &npw,fc%npwt, ig_l2g,fc%ig_l2gt, fc%ngmt_g , mpime, nproc,ionode_id, intra_pool_comm )

    ! do ii=1,numpw
    !    call mergewf(p_basis(:,ii),evc_g,npw,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
    !    call splitwf(p_basis_t(:,ii),evc_g,fc%npwt,fc%ig_l2gt,mpime,nproc,ionode_id,intra_pool_comm)
    ! enddo
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
        CALL cft3t(fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, 2 )
        p_basis_r(1:fc%nrxxt,ii)= DBLE(psic(1:fc%nrxxt))
        if(ii/=numpw) p_basis_r(1:fc%nrxxt,ii+1)= DIMAG(psic(1:fc%nrxxt))

     enddo
!  endif


!now valence wavefunctions are put on the ordering of the reduced grid
  allocate(evc_t(fc%npwt,num_nbndv(ispin)))
  allocate(evc_r(fc%nrxxt,num_nbndv(ispin)))
  if(fc%dual_t==4.d0) then
     evc_t(:,1:num_nbndv(ispin))=evc(:,1:num_nbndv(ispin))
  else
     call reorderwfp (num_nbndv(ispin),npw, fc%npwt,evc(:,:),evc_t(:,:), &
              &npw,fc%npwt, ig_l2g,fc%ig_l2gt, fc%ngmt_g , mpime, nproc,ionode_id, intra_pool_comm )

   !   do iv=1,num_nbndv(ispin)
   !      call mergewf(evc(:,iv),evc_g,npw,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
   !      call splitwf(evc_t(:,iv),evc_g,fc%npwt,fc%ig_l2gt,mpime,nproc,ionode_id,intra_pool_comm)
   !   enddo
  endif

  evc_r=0.d0
   do iv=1,num_nbndv(ispin)
      psic(:)=(0.d0,0.d0)
      psic(fc%nlt(1:fc%npwt))  = evc_t(1:fc%npwt,iv)
      psic(fc%nltm(1:fc%npwt)) = CONJG( evc_t(1:fc%npwt,iv) )
      CALL cft3t(fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, 2 )
      evc_r(1:fc%nrxxt,iv)= DBLE(psic(1:fc%nrxxt))
   enddo


!loop on v

  allocate(tmp_r(fc%nrxxt),tmp_r2(fc%nrxxt))
  allocate(omat(numpw,numpw),omat_hold(numpw,numpw))
  allocate(t_mat(numpw,nstates), t_mat_hold(numpw,nstates), t_mat_hold2(numpw,nstates))
  allocate(wp_g(npw,nstates))
  allocate(wp_g_t(fc%nrxxt,nstates))
  allocate(t_eigen_hold(nstates))

!set the number of products to be distributed
  nbuf=min(12,nproc)
  allocate(wp_prod(fc%nrxxt,numpw,nbuf))
  wp_prod=0.d0
  wp_g_t=0.d0
  
   do ivv=1,num_nbndv(ispin),nbuf
!put iv on real space
      
      do iv=ivv,min(ivv+nbuf-1,num_nbndv(ispin))
   
         wv_real(1:fc%nrxxt)= evc_r(1:fc%nrxxt,iv)
   

!!loop on products of wanniers
      
         
    !  allocate(tmp_r(fc%nrxxt))

         if(l_verbose) write(stdout,*) 'do fft',fc%nrxxt, numpw,iv-ivv+1
         FLUSH(stdout)
         do ii=1,numpw
            wp_prod(1:fc%nrxxt, ii,iv-ivv+1)=p_basis_r(1:fc%nrxxt,ii)*wv_real(1:fc%nrxxt) 
         enddo

       
         call pc_operator_t_r(numpw,wp_prod(1,1,iv-ivv+1),evc_r,ispin, fc)

         if(l_verbose) write(stdout,*) 'calculate omat'
         FLUSH(stdout)


    
       !!calculate overlap matrix
    
         call dgemm('T','N',numpw,numpw,fc%nrxxt,1.d0,wp_prod(1,1,iv-ivv+1),fc%nrxxt,&
              &wp_prod(1,1,iv-ivv+1),fc%nrxxt,0.d0,omat,numpw)
         do ii=1,numpw
            call mp_sum(omat(1:numpw,ii),world_comm)
            omat(1:numpw,ii)=omat(1:numpw,ii)/dble(fc%nr1t*fc%nr2t*fc%nr3t)
         enddo
         
!set up norms
!         do ii=1,numpw
!            norms(ii)=omat(ii,ii)
!         enddo


         if(iv-ivv==mpime) then
            omat_hold(:,:)=omat(:,:)
         endif
      enddo
  
!!
!!solve eigen/values vector problem 
!!
      
     if(l_verbose)  write(stdout,*) 'solve eig'
     FLUSH(stdout)

    
      FLUSH(stdout)
      do iv=ivv,min(ivv+nbuf-1,num_nbndv(ispin))
         if(l_verbose) write(stdout,*) 'solve eig', iv
         FLUSH(stdout)
         if(iv-ivv==mpime) then
            
            if(.not.l_dsyevr) then
               allocate(eigen(numpw))
               allocate(work(1))
               call  DSYEV( 'V', 'U', numpw, omat_hold, numpw, eigen, work, -1, info )
               lwork=work(1)
               deallocate(work)
               allocate(work(lwork))
               call  DSYEV( 'V', 'U', numpw, omat_hold, numpw, eigen, work, lwork, info )
               deallocate(work)
               if(info/=0) then
                  write(stdout,*) 'ROUTINE pola_basis_lanczos, INFO:', info
                  stop
               endif
            

             !  do iw=1,numpw
             !     write(stdout,*) 'EIGEN:',iv,iw, eigen(iw)
             !  enddo
             !  FLUSH(stdout)
            else
               if(l_verbose) write(stdout,*) 'ATT1'
               FLUSH(stdout)
               allocate(eigen(numpw))
               allocate(vectors(numpw,nstates))
               allocate(isuppz(2*nstates))
               allocate(work(1),iwork(1))
               if(l_verbose)  write(stdout,*) 'ATT2'
                FLUSH(stdout)
               call DSYEVR('V','I','U',numpw,omat_hold,numpw,0.d0,0.d0,numpw-nstates+1,numpw,0.d0,n_found,eigen,&
                    & vectors,numpw,isuppz,work, -1,iwork,-1, info)
               lwork=work(1)
               liwork=iwork(1)
               deallocate(work,iwork)
               allocate(work(lwork))
               allocate(iwork(liwork))
               if(l_verbose)  write(stdout,*) 'ATT3',numpw,nstates,size(omat_hold(:,1)),size(omat_hold(1,:)),lwork,liwork
               FLUSH(stdout)
               vl=0.d0
               vu=0.d0
               il=numpw-nstates+1
               iu=numpw
               call DSYEVR('V','I','U',numpw,omat_hold,numpw,vl,vu,il,iu,0.d0,n_found,eigen,&
                    & vectors,numpw,isuppz,work,lwork,iwork,liwork, info)
               if(info/=0) then
                  write(stdout,*) 'ROUTINE pola_lanczos DSYEVR, INFO:', info
                  stop
               endif
               if(l_verbose)  write(stdout,*) 'ATT4'
               FLUSH(stdout)
               deallocate(isuppz)
               deallocate(work,iwork)
        
               do iw=1,nstates
                  write(stdout,*) 'EIGEN:',iv,iw, eigen(iw)
               enddo
               FLUSH(stdout)

            endif
      
      
!!find transformation matrix and write on disk
!
        if(l_verbose) write(stdout,*) 'pola_basis_lanczos t_mat'
        FLUSH(stdout)


      

        if(.not.l_dsyevr) then
           do ii=1,nstates
              do jj=1,numpw
                 t_mat_hold(jj,ii)=omat_hold(jj,numpw-ii+1)*(dsqrt(eigen(numpw-ii+1)))
              enddo
              t_eigen_hold(ii)=eigen(numpw-ii+1)
           enddo
        else
           do ii=1,nstates
              do jj=1,numpw
                 t_mat_hold(jj,ii)=vectors(jj,ii)*(dsqrt(eigen(ii)))
              enddo
              t_eigen_hold(ii)=eigen(ii)
           enddo
        endif

        !!find liner dependent products

        if(.not.l_dsyevr) then
           do ii=1,nstates
              t_mat_hold2(:,ii)=omat_hold(:,numpw-ii+1)*(1.d0/dsqrt(eigen(numpw-ii+1)))
           enddo
        else
           do ii=1,nstates
              t_mat_hold2(:,ii)=vectors(:,ii)*(1.d0/dsqrt(eigen(ii)))
           enddo
        endif
     
        deallocate(eigen)
        if(l_dsyevr) deallocate(vectors)

     endif
     

    enddo
    allocate(eigen(nstates))
     do iv=ivv,min(ivv+nbuf-1,num_nbndv(ispin))
        if(iv-ivv == mpime) then
           t_mat(:,:)=t_mat_hold(:,:)
           eigen(1:nstates)=t_eigen_hold(1:nstates)
        endif
        call mp_bcast(t_mat,iv-ivv,world_comm)
        call mp_bcast(eigen(1:nstates),iv-ivv,world_comm)
!if required imposes sum rule

!        if(l_sumrule) then
!           norm_t=0.d0
!           do jj=1,numpw
!              do ii=1,nstates
!                 norm_t=norm_t+t_mat(jj,ii)**2.d0
!              enddo
!           enddo
!           norm=0.d0
!           do jj=1,numpw
!              norm=norm+norms(jj)
!           enddo
!           c_norm=dsqrt(norm/norm_t)
!           write(stdout,*) 'Sum rule:',c_norm
!           t_mat(:,:)=t_mat(:,:)*c_norm
!        endif

        if(ionode) then
           iuntmat = find_free_unit()
           write(nfile,'(4i1)') iv/1000,mod(iv,1000)/100,mod(iv,100)/10,mod(iv,10)
           if(ispin==1) then
              open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.p_mat_lanczos'//nfile, &
                   &status='unknown',form='unformatted')
           else
              open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.p_mat_lanczos2'//nfile, &
                   &status='unknown',form='unformatted')
           endif
           write(iuntmat) iv
           write(iuntmat) num_nbndv(ispin)
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
              open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.p_eig_lanczos'//nfile, &
                   &status='unknown',form='unformatted')
           else
              open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.p_eig_lanczos2'//nfile, &
                   &status='unknown',form='unformatted')
           endif
           write(iuntmat) nstates
           write(iuntmat) eigen(1:nstates)
           close(iuntmat)
        endif


        if(l_verbose) write(stdout,*) 'pola_basis update wp_g'
        FLUSH(stdout)


!!find liner dependent products

        if(iv-ivv == mpime) then
           t_mat(:,:)=t_mat_hold2(:,:)
        endif
        call mp_bcast(t_mat,iv-ivv,world_comm)



        if(l_verbose) write(stdout,*) 'pola_basis update wp_g dgemm'
        FLUSH(stdout)

        call dgemm('N','N',fc%nrxxt,nstates,numpw,1.d0,wp_prod(1,1,iv-ivv+1),fc%nrxxt,t_mat,numpw,0.d0,wp_g_t,fc%nrxxt)

        write(stdout,*) 'pola_basis update merge-split',iv,ivv
        FLUSH(stdout)

!put the correct order
       
        psic=0.d0
        allocate(wp_g_t2(fc%npwt,nstates)) 
        do ii=1,nstates,2
           if(ii==nstates) then
              psic(1:fc%nrxxt)=wp_g_t(1:fc%nrxxt,ii)
           else
              psic(1:fc%nrxxt)=cmplx(wp_g_t(1:fc%nrxxt,ii),wp_g_t(1:fc%nrxxt,ii+1))
           endif
           CALL cft3t( fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, -2 )
           if(ii==nstates) then
              wp_g_t2(1:fc%npwt,ii) = psic(fc%nlt(1:fc%npwt))              
!project on conduction manifold
              call  pc_operator_t(wp_g_t2(:,ii),evc_t,ispin,fc)
           else
              wp_g_t2(1:fc%npwt, ii)= 0.5d0*(psic(fc%nlt(1:fc%npwt))+conjg( psic(fc%nltm(1:fc%npwt))))
              wp_g_t2(1:fc%npwt, ii+1)= (0.d0,-0.5d0)*(psic(fc%nlt(1:fc%npwt)) - conjg(psic(fc%nltm(1:fc%npwt))))
              call  pc_operator_t(wp_g_t2(:,ii),evc_t,ispin,fc)
              call  pc_operator_t(wp_g_t2(:,ii+1),evc_t,ispin,fc)
           endif

          
        enddo


        if(fc%dual_t==4.d0) then
           wp_g(1:npw,1:nstates)=wp_g_t2(1:fc%npwt,1:nstates)
        else
           call reorderwfp (nstates,fc%npwt, npw,wp_g_t2,wp_g, &
                &fc%npwt,npw, fc%ig_l2gt,ig_l2g, fc%ngmt_g , mpime, nproc,ionode_id, intra_pool_comm )
           ! call mergewf(wp_g_t2(:),evc_g,fc%npwt,fc%ig_l2gt,mpime,nproc,ionode_id,intra_pool_comm)
           ! call splitwf(wp_g(:,ii),evc_g,npw,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
    
        endif
             


        deallocate(wp_g_t2)
        if(l_verbose) write(stdout,*) 'pola_basis update davcio',iv
        FLUSH(stdout)


  !!write on disk
        do ii=1,nstates
           call davcio(wp_g(:,ii),npw*2,iungresult,ii+(iv-1)*nstates,1)
        enddo

        if(l_verbose) write(stdout,*) 'pola_basis update done'
        FLUSH(stdout)


       
     enddo
     deallocate(eigen)
  enddo
  deallocate(evc_r)
  deallocate(t_mat,t_mat_hold,t_mat_hold2)
  deallocate(omat,omat_hold,tmp_r,tmp_r2,p_basis_r)
  deallocate(wp_g,wp_g_t)
  close(iungresult)
  deallocate(norms)
  deallocate(wv_real,wp_prod)
  deallocate(t_eigen_hold)
  
  if(l_verbose) write(stdout,*) 'Exiting pola_basis_lanczos'
  FLUSH(stdout)
   
  deallocate(p_basis)
  deallocate(p_basis_t,evc_t)
  
  if(l_verbose) write(stdout,*) 'Call deallocate_fft_custom'
  FLUSH(stdout)
  !deallocate(evc_g)
  call deallocate_fft_custom(fc)
   

   return
!NOT_TO_BE_INCLUDED_END
 end subroutine pola_basis_lanczos_real

