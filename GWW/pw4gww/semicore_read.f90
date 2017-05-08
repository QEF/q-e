!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

!this routine calculate the terms  \int dr \psi_i(r)\_psi_v(sc)(r) (v\Phi_\mu)(r)
!where \psi_i(r)\_psi_v(sc)(r)   are read from the disk                                                          

  subroutine semicore_read(num_nbnds,numpw, ispin)
!NOT_TO_BE_INCLUDED_START
     USE io_global,             ONLY : stdout, ionode,ionode_id
     USE io_files,              ONLY : diropn,prefix, tmp_dir
     use pwcom
     USE wavefunctions_module,  ONLY : evc
     USE kinds,                 ONLY : DP
     USE gvect,                 ONLY : g, ig_l2g, nl, nlm, gstart
     USE gvecs,                 ONLY : nls, nlsm
     USE mp,           ONLY : mp_sum, mp_barrier, mp_bcast
     USE mp_global,             ONLY : inter_pool_comm, intra_pool_comm
     USE mp_wave, ONLY : mergewf,splitwf
     USE mp_world, ONLY : mpime, nproc, world_comm
     USE fft_base,             ONLY : dfftp, dffts
     USE fft_interfaces,       ONLY : fwfft, invfft
     USE wavefunctions_module, ONLY : psic
     USE wvfct,    ONLY : et
     USE lsda_mod,             ONLY : nspin
     USE wannier_gw,           ONLY : max_ngm, l_truncated_coulomb,vg_q,truncation_radius
     USE constants, ONLY : e2, pi, tpi, fpi
     USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2


    implicit none

    INTEGER, EXTERNAL :: find_free_unit
    INTEGER, INTENT(in) :: num_nbnds!total KS states considered
    INTEGER, INTENT(in) :: numpw!dimension of polarizability basis
    INTEGER, INTENT(in) :: ispin!spin channel

    INTEGER :: num_nbnds_sc,n_semicore,npwx_g_sc
    INTEGER :: iv, iun, ii,jj,iunw,ii_max  
    INTEGER :: ig
    REAL(kind=DP), ALLOCATABLE :: et_sc(:)
    COMPLEX(kind=DP), ALLOCATABLE :: tmp_g(:),tmp_wfc(:)
    REAL(kind=DP), ALLOCATABLE :: pp_sc(:,:,:),tmp_r(:),prods(:)
    INTEGER :: iungprod,iw
    LOGICAL :: exst
    REAL(kind=DP), ALLOCATABLE :: fac(:)
    REAL(kind=DP) :: qq
    COMPLEX(kind=DP), ALLOCATABLE :: psi_all(:,:)!for all-electron (with semicore) wfcs 
    REAL(kind=DP), ALLOCATABLE :: o_mat(:,:)
    REAL(kind=DP) :: sca
    INTEGER, ALLOCATABLE :: order(:)

    allocate(tmp_wfc(npwx))
    allocate(order(num_nbnds))


    allocate(fac(max_ngm))
    allocate(psi_all(npwx,num_nbnds))
    allocate(o_mat(num_nbnds,num_nbnds))


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

!open files and allocate
       
    if(ionode) then
       iun =  find_free_unit()
       if(ispin==1) then
          open( unit= iun, file=trim(tmp_dir)//trim(prefix)//'.sc_states', status='old',form='unformatted')
       else
          open( unit= iun, file=trim(tmp_dir)//trim(prefix)//'.sc_states2', status='old',form='unformatted')
       endif
       read(iun) num_nbnds_sc
       read(iun) n_semicore
       read(iun) npwx_g_sc
    endif
    call mp_bcast( num_nbnds_sc, ionode_id,world_comm)
    call mp_bcast(n_semicore, ionode_id,world_comm)
    call mp_bcast(npwx_g_sc,ionode_id,world_comm)
    allocate(pp_sc(dfftp%nnr,n_semicore,num_nbnds))

    allocate(tmp_g(npwx_g_sc))
    allocate(et_sc(num_nbnds_sc))
    if(ionode)  read(iun) et_sc(1:num_nbnds_sc)
    call mp_bcast(et_sc, ionode_id,world_comm)
    write(stdout,*) 'NUM. SEMICORE:', n_semicore
    write(stdout,*) 'NUM. BANDS SC:',num_nbnds_sc
    write(stdout,*) 'NUM. BANDS:',num_nbnds
    write(stdout,*) 'NPWX_G_SC:',npwx_g_sc
    write(stdout,*) 'ET_SC:',et_sc(1:num_nbnds_sc)
!write header on file for each spin channel
    iunw =  find_free_unit()
    if(ispin==1) then
       open( unit= iunw, file=trim(tmp_dir)//trim(prefix)//'.sc_gvphi', status='unknown',form='unformatted')
    else
       open( unit= iunw, file=trim(tmp_dir)//trim(prefix)//'.sc_gvphi2', status='unknown',form='unformatted')
    endif
    write(iunw) n_semicore
    write(iunw) et_sc(1:n_semicore)
    write(iunw) num_nbnds
    write(iunw) numpw

!read in all semicore products and split to charge grid
    do ii=1,num_nbnds_sc-n_semicore!as they are written
       do iv=1,n_semicore
          write(stdout,*) 'Reading state:', ii,iv
          FLUSH(stdout)
          if(ionode) read(iun) tmp_g(1:npwx_g_sc)
          call splitwf(tmp_wfc,tmp_g,npw,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
!check for consistency
          if(gstart==2) write(stdout,*) 'It should be zero:', tmp_wfc(1)

             
!trasform them to R grid

          psic(:)=(0.d0,0.d0)
          psic(nls(igk_k(1:npw,1)))  = tmp_wfc(1:npw)
          psic(nlsm(igk_k(1:npw,1))) = CONJG( tmp_wfc(1:npw) )
          CALL invfft ('Wave', psic, dffts)
          
          if(ii<num_nbnds)  pp_sc(1:dfftp%nnr,iv,ii)=dble(psic(1:dfftp%nnr))


       enddo
    enddo

!now read KS states
    do ii=1,num_nbnds_sc-n_semicore
       write(stdout,*) 'Reading KS all-electron state:', ii
       if(ionode) read(iun) tmp_g(1:npwx_g_sc)
       call splitwf(psi_all(:,ii),tmp_g,npw,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
    enddo


!in case put to 0 remaining states
    if(num_nbnds > (num_nbnds_sc-n_semicore)) then
       pp_sc(1:dfftp%nnr,1:n_semicore,num_nbnds_sc-n_semicore+1:num_nbnds)=0.d0
    endif

    if(num_nbnds > (num_nbnds_sc-n_semicore)) then
       psi_all(1:npwx,num_nbnds_sc-n_semicore+1:num_nbnds)=(0.d0,0.d0)
    endif


    deallocate(tmp_g)
    close(iun)
!check all-electrons states

    call dgemm('T','N',num_nbnds,num_nbnds,2*npw,2.d0,psi_all,2*npwx,psi_all,2*npwx,0.d0,o_mat,num_nbnds)
    if(gstart==2) then
       do ii=1,num_nbnds
          do jj=1,num_nbnds
             o_mat(ii,jj)=o_mat(ii,jj)-psi_all(1,ii)*psi_all(1,jj)
          enddo
       enddo
    endif
    call mp_sum(o_mat,world_comm)
    do ii=1,num_nbnds
       write(stdout,*) 'Orthonormality:',ii,o_mat(1:num_nbnds,ii)
    enddo




!do products of all electron with pseudo states and determines the correspondance

    call dgemm('T','N',num_nbnds,num_nbnds,2*npw,2.d0,psi_all,2*npwx,evc,2*npwx,0.d0,o_mat,num_nbnds)
    if(gstart==2) then
       do ii=1,num_nbnds
          do jj=1,num_nbnds
             o_mat(ii,jj)=o_mat(ii,jj)-psi_all(1,ii)*evc(1,jj)
          enddo
       enddo
    endif
    call mp_sum(o_mat,world_comm)

    do ii=1,num_nbnds
       sca=0.d0
       ii_max=0
       do  jj=1,num_nbnds
          if(abs(o_mat(jj,ii))>sca) then
             ii_max=jj
             sca=abs(o_mat(jj,ii))
          endif
       enddo
       write(stdout,*) 'KS state:',ii,'corresponds to AE state:',ii_max,sca
       order(ii)=ii_max
    enddo

       
    do ii=1,num_nbnds
       sca=0.d0
       ii_max=0
       do  jj=1,num_nbnds
          if(abs(o_mat(ii,jj))>sca) then
             ii_max=jj
             sca=abs(o_mat(ii,jj))
          endif
       enddo
       write(stdout,*) 'AE state',ii,'corresponds to pseudo  state:',ii_max,sca
    enddo


!open polarizability basis

    allocate(tmp_g(max_ngm),tmp_r(dfftp%nnr),prods(n_semicore))
    iungprod = find_free_unit()
    CALL diropn( iungprod, 'wiwjwfc_red', max_ngm*2, exst )
!loop on pol basis vectors
    do iw=1,numpw
       call davcio(tmp_g,max_ngm*2,iungprod,iw,-1)
!trasform to r-space 
       psic(:)=(0.d0,0.d0)
       do ig=1,max_ngm
          psic(nl(ig))=tmp_g(ig)*fac(ig)
          psic(nlm(ig))=CONJG(tmp_g(ig))*fac(ig)
       enddo
       CALL invfft ('Dense', psic, dfftp)
       tmp_r(1:dfftp%nnr)=dble(psic(1:dfftp%nnr))

       do ii=1,num_nbnds
!!do products
          call dgemv ('T',dfftp%nnr,n_semicore,1.d0,pp_sc(1,1,order(ii)),dfftp%nnr,tmp_r,1,0.d0,prods, 1)

          call mp_sum(prods(1:n_semicore),world_comm)
          prods(1:n_semicore)=prods(1:n_semicore)/dble(dfftp%nr1*dfftp%nr2*dfftp%nr3)
!! write on disk
          write(iunw) prods(1:n_semicore)
       enddo
    enddo
    close(iungprod)
    deallocate(tmp_g,tmp_r,prods)
    close(iunw)
    deallocate(pp_sc)


    deallocate(tmp_wfc)
    deallocate(psi_all)
    deallocate(o_mat)
    deallocate(order)
    return
!NOT_TO_BE_INCLUDED_END
  end subroutine semicore_read
