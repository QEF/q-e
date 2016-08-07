!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

!routines for the set up of terms related with the calòculation
!of the polarizability for partially occupied states


subroutine pola_partial(numpw,ispin)
!NOT_TO_BE_INCLUDED_START
!this subroutine calculates the basis for every v                                                                
!the minimal orthonormal basis for the w_v(r)*w^P'_i(r) products                                                                                
!the occupancies are from  variable wg and ste to range 0,1            

   USE io_global,            ONLY : stdout, ionode, ionode_id
   USE io_files,             ONLY : prefix, tmp_dir, diropn
   USE kinds,    ONLY : DP
   USE wannier_gw
   USE gvect
   USE constants, ONLY : e2, pi, tpi, fpi
   USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2
   USE wvfct,     ONLY : npwx, npw, nbnd, wg
   USE gvecw,     ONLY : ecutwfc
   USE wavefunctions_module, ONLY : evc, psic
   USE mp, ONLY : mp_sum, mp_barrier, mp_bcast
   USE mp_world, ONLY : mpime,nproc, world_comm
   USE mp_pools, ONLY : intra_pool_comm
   USE gvecs,              ONLY : nls, nlsm, doublegrid
   USE fft_custom_gwl
   USE mp_wave, ONLY : mergewf,splitwf
   USE fft_base,             ONLY : dfftp, dffts
   USE fft_interfaces,       ONLY : fwfft, invfft
   USE lsda_mod,             ONLY : lsda, nspin
     


   implicit none

   INTEGER, EXTERNAL :: find_free_unit
   INTEGER, INTENT(in) :: numpw!dimension of polarization basis                                                                               
   INTEGER, INTENT(in) :: ispin! spin channel 

   COMPLEX(kind=DP), ALLOCATABLE :: evc_g(:),p_basis(:,:),p_basis_t(:,:),evc_t(:,:)
   INTEGER :: iw,ii,iv,jv
   INTEGER :: iungprod,iun
   LOGICAL :: exst
   REAL(kind=DP), ALLOCATABLE :: evc_r(:,:), prod_r(:), psi_psi_phi(:)
   COMPLEX(kind=DP), ALLOCATABLE :: prod_g(:)

   TYPE(fft_cus) :: fc

   fc%ecutt=ecutwfc
   fc%dual_t=dual_vt
   
   write(stdout,*) 'Call initialize_fft_custom'
   FLUSH(stdout)
   call initialize_fft_custom(fc)
   allocate(evc_g(fc%ngmt_g))

   iungprod = find_free_unit()
   CALL diropn( iungprod, 'wiwjwfc_red', max_ngm*2, exst )
!read polarizability basis functions                                                                                    
   allocate(p_basis(max_ngm,numpw))
   do iw=1,numpw
      call davcio(p_basis(:,iw),max_ngm*2,iungprod,iw,-1)
   enddo

   close(iungprod)

!now polarizability basis  are put on the ordering of the redueced grid, if required                                                                        
   allocate(p_basis_t(fc%npwt,numpw))
   if(fc%dual_t==4.d0) then
      p_basis_t(1:npw,1:numpw)=p_basis(1:npw,1:numpw)
   else
      do ii=1,numpw
         call mergewf(p_basis(:,ii),evc_g,npw,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
         call splitwf(p_basis_t(:,ii),evc_g,fc%npwt,fc%ig_l2gt,mpime,nproc,ionode_id,intra_pool_comm)
      enddo
   endif


   allocate(evc_t(fc%npwt,num_nbndv(ispin)))
   if(fc%dual_t==4.d0) then
      evc_t(1:npw,1:num_nbndv(ispin))=evc(1:npw,1:num_nbndv(ispin))
   else
      do ii=1,num_nbndv(ispin)
         call mergewf(evc(:,ii),evc_g,npw,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
         call splitwf(evc_t(:,ii),evc_g,fc%npwt,fc%ig_l2gt,mpime,nproc,ionode_id,intra_pool_comm)
      enddo
   endif
!put KS wfcs in real space
   allocate(evc_r(fc%nrxxt,num_nbndv(ispin)))
   do iv=1,num_nbndv(ispin),2
      psic(:)=(0.d0,0.d0)
      if(iv==num_nbndv(ispin)) then
         psic(fc%nlt(1:fc%npwt))  = evc_t(1:fc%npwt,iv)
         psic(fc%nltm(1:fc%npwt)) = CONJG( evc_t(1:fc%npwt,iv) )
      else
         psic(fc%nlt(1:fc%npwt))=evc_t(1:fc%npwt,iv)+(0.d0,1.d0)*evc_t(1:fc%npwt,iv+1)
         psic(fc%nltm(1:fc%npwt)) = CONJG( evc_t(1:fc%npwt,iv) )+(0.d0,1.d0)*CONJG( evc_t(1:fc%npwt,iv+1) )
        endif
        CALL cft3t( fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, 2 )
        evc_r(1:fc%nrxxt,iv)= DBLE(psic(1:fc%nrxxt))
        if(iv/=num_nbndv(ispin)) evc_r(1:fc%nrxxt,iv+1)= DIMAG(psic(1:fc%nrxxt))

     enddo
!open file and write occupancies
     if(ionode) then
        iun=find_free_unit()
        if(ispin==1) then
           open( unit= iun, file=trim(tmp_dir)//trim(prefix)//'.occ_mat', status='unknown',form='unformatted')
        else
           open( unit= iun, file=trim(tmp_dir)//trim(prefix)//'.occ_mat2', status='unknown',form='unformatted')
        endif
        write(iun) num_nbndv_min(ispin)
        write(iun) num_nbndv(ispin)
        write(iun) numpw
        if(nspin==1) then
           write(iun) wg(1:num_nbndv(ispin),ispin)/2.d0
        else
           write(iun) wg(1:num_nbndv(ispin),ispin)
        endif
     endif




!loop on partially occupied states
     allocate(prod_r(fc%nrxxt))
     allocate(prod_g(fc%npwt),psi_psi_phi(numpw))
     do iv=num_nbndv_min(ispin)+1,num_nbndv(ispin)
!loop on all states
        do jv=1,num_nbndv(ispin)
           prod_r(1:fc%nrxxt)=evc_r(1:fc%nrxxt,iv)*evc_r(1:fc%nrxxt,jv)
           psic(1:fc%nrxxt)=dcmplx(prod_r(1:fc%nrxxt),0.d0)
           CALL cft3t( fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, -2 )
           prod_g(1:fc%npwt) = psic(fc%nlt(1:fc%npwt))
           

!do product with polarizability basis
           call dgemm('T','N',numpw,1,2*fc%npwt,2.0,p_basis_t,2*fc%npwt,prod_g,2*fc%npwt,0.d0,psi_psi_phi,numpw)
           if(fc%gstart_t==2) then
              do ii=1,numpw
                 psi_psi_phi(ii)=psi_psi_phi(ii)-dble(conjg(p_basis_t(1,ii))*prod_g(1))
              enddo
           endif
           call mp_sum(psi_psi_phi,world_comm)
           if(ionode) write(iun) psi_psi_phi(1:numpw)
        enddo
     enddo
!write
 
     


     if(ionode) close(iun)
     write(stdout,*) 'Call deallocate_fft_custom'
     FLUSH(stdout)
     deallocate(evc_g,p_basis,p_basis_t)
     deallocate(evc_t,evc_r)
     deallocate(prod_r,prod_g)
     deallocate(psi_psi_phi)
     call deallocate_fft_custom(fc)

     return
!NOT_TO_BE_INCLUDED_END
end  subroutine pola_partial
