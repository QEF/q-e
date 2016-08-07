!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

!this routine calculate the terms  \psi_i(r)\_psi_v(sc)(r)
!and write them on disk on global G grid

  subroutine semicore(n_semicore, num_nbnds,ispin)
!NOT_TO_BE_INCLUDED_START
     USE io_global,             ONLY : stdout, ionode,ionode_id
     USE io_files,              ONLY : diropn,prefix,tmp_dir
     use pwcom
     USE wavefunctions_module,  ONLY : evc
     USE kinds,                 ONLY : DP
     USE gvect,                 ONLY : ig_l2g
     USE mp,           ONLY : mp_sum, mp_barrier, mp_bcast
     USE mp_wave,      ONLY : mergewf,splitwf
     USE mp_pools,     ONLY : intra_pool_comm, inter_pool_comm, intra_pool_comm
     USE mp_world,     ONLY : world_comm, mpime, nproc
     USE fft_base,             ONLY : dfftp, dffts
     USE fft_interfaces,       ONLY : fwfft, invfft
     USE wavefunctions_module, ONLY : psic
     USE wvfct,    ONLY : et

    implicit none

    INTEGER, EXTERNAL :: find_free_unit
    INTEGER, INTENT(in) :: n_semicore!number of semicore states
    INTEGER, INTENT(in) :: num_nbnds!total KS states considered
    INTEGER, INTENT(in) :: ispin!spin channel

    REAL(kind=DP), ALLOCATABLE :: psi_sc(:,:)
    COMPLEX(kind=DP), ALLOCATABLE :: prod(:), prod_g(:,:), prod_g_tot(:)
    INTEGER :: iv, iun, ii
    INTEGER :: npwx_g


!fft trasform semicore states to R space
    allocate(psi_sc(dfftp%nnr,n_semicore))
    allocate(prod(dfftp%nnr), prod_g(npw,2),prod_g_tot(ngm_g))
    
    do iv=1,n_semicore,2
       psic(:)=(0.d0,0.d0)
       if(iv<n_semicore) then
          psic(nls(igk_k(1:npw,1)))  = evc(1:npw,iv) + &
               ( 0.D0, 1.D0 ) * evc(1:npw,iv+1)
          psic(nlsm(igk_k(1:npw,1))) = CONJG( evc(1:npw,iv) - &
                    ( 0.D0, 1.D0 ) * evc(1:npw,iv+1) )

       else
          psic(nls(igk_k(1:npw,1)))  = evc(1:npw,iv)
          psic(nlsm(igk_k(1:npw,1))) = CONJG( evc(1:npw,iv) )
       endif
       CALL invfft ('Wave', psic, dffts)
       psi_sc(1:dfftp%nnr,iv)=dble(psic(1:dfftp%nnr))
       if(iv< n_semicore) psi_sc(1:dfftp%nnr,iv+1)=dimag(psic(1:dfftp%nnr))

    enddo


!write header of file with KS energies of all states (just in case..)
    
    npwx_g=npwx
    call mp_sum(npwx_g, world_comm)

    if(ionode) then
       iun =  find_free_unit()
       if(ispin==1) then
          open( unit= iun, file=trim(tmp_dir)//trim(prefix)//'.sc_states', status='unknown',form='unformatted')
       else
          open( unit= iun, file=trim(tmp_dir)//trim(prefix)//'.sc_states2', status='unknown',form='unformatted')
       endif
       write(iun) num_nbnds
       write(iun) n_semicore
       write(iun) npwx_g
       write(iun) et(1:num_nbnds,ispin)
    endif
    

!loop (double) on KS states
    do ii=n_semicore+1,num_nbnds

!fft

       psic(:)=(0.d0,0.d0)
       psic(nls(igk_k(1:npw,1)))  = evc(1:npw,ii)
       psic(nlsm(igk_k(1:npw,1))) = CONJG( evc(1:npw,ii) )

       CALL invfft ('Wave', psic, dffts)

!calculate product
       do iv=1,n_semicore,2
          if(iv<n_semicore) then
             prod(1:dfftp%nnr)=dcmplx(dble(psic(1:dfftp%nnr))*psi_sc(1:dfftp%nnr,iv),&
                  &dble(psic(1:dfftp%nnr))*psi_sc(1:dfftp%nnr,iv+1))
          else
             prod(1:dfftp%nnr)=dble(psic(1:dfftp%nnr))
             prod(1:dfftp%nnr)=prod(1:dfftp%nnr)*psi_sc(1:dfftp%nnr,iv)
          endif
!fft back

          CALL fwfft ('Dense', prod, dfftp)
          if(iv==n_semicore) then
             prod_g(1:npw,1)=prod(nl(1:npw))
             if(gstart==2) then
                write(stdout,*) 'Putting to zero:', iv,ii, prod_g(1,1)
                prod_g(1,1)=(0.d0,0.d0)
             endif
          else
             prod_g(1:npw, 1)= 0.5d0*(prod(nl(1:npw))+conjg( prod(nlm(1:npw))))
             prod_g(1:npw, 2)= (0.d0,-0.5d0)*(prod(nl(1:npw)) - conjg(prod(nlm(1:npw))))
             if(gstart==2) then
                write(stdout,*)'Putting to zero:', iv,ii, prod_g(1,1)
                write(stdout,*)'Putting to zero:', iv+1,ii, prod_g(1,2)
                prod_g(1,1:2)=(0.d0,0.d0)
             endif
          endif
                    
!merge

          prod_g_tot(:)=(0.d0,0.d0)
          call mergewf(prod_g(:,1),prod_g_tot,npw,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
!put it on disk 
          if(ionode) then
             write(iun) prod_g_tot(1:npwx_g)
          endif
          if(iv<n_semicore) then

!merge
             prod_g_tot(:)=(0.d0,0.d0)
             call mergewf(prod_g(:,2),prod_g_tot,npw,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
!put it on disk
             if(ionode) then
             write(iun) prod_g_tot(1:npwx_g)
          endif


          endif
       enddo
 



    enddo


!now put on disk KS states on global G 
    do ii=n_semicore+1,num_nbnds
       prod_g_tot(:)=(0.d0,0.d0)
       call mergewf(evc(:,ii),prod_g_tot,npw,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
!put it on disk
       if(ionode) then
          write(iun) prod_g_tot(1:npwx_g)
       endif
    enddo
    

    deallocate(psi_sc,prod_g_tot)
    close(iun)
    return
!NOT_TO_BE_INCLUDED_END
  end subroutine semicore
