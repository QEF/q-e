!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

!this subroutine writes polarizability basis in real space on charge grid on disk
  subroutine write_pola_basis(numpw, itask)

    USE io_global,            ONLY : stdout, ionode, ionode_id
   USE io_files,             ONLY : prefix, tmp_dir, diropn
   USE kinds,    ONLY : DP
   USE wannier_gw
   USE gvect
   USE constants, ONLY : e2, pi, tpi, fpi
   USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2
   USE wvfct,    ONLY : npwx, npw, nbnd
   USE wavefunctions_module, ONLY : evc, psic
   USE mp, ONLY : mp_sum, mp_barrier, mp_bcast
   USE mp_pools, ONLY : intra_pool_comm
   USE mp_world, ONLY : mpime,nproc
   USE gvecs,              ONLY : nls, nlsm, doublegrid
   USE mp_wave, ONLY : mergewf,splitwf
   USE fft_base,             ONLY : dfftp, dffts
   USE fft_interfaces,       ONLY : fwfft, invfft
  
    
    implicit none

    INTEGER, INTENT(in) :: numpw!dimesion of polarizability basis
    INTEGER, INTENT(in) :: itask!0 writes (v\Phi)  1 writes (\Phi)

    COMPLEX(kind=DP), allocatable :: p_basis(:,:)
    INTEGER, external :: find_free_unit
    INTEGER :: iungprod,iw,ii,iun,ig
    LOGICAL :: exst
    REAL(kind=DP), allocatable :: p_basis_r(:,:)
    REAL(kind=DP), ALLOCATABLE :: fac(:)
    REAL(kind=DP) :: qq

    iungprod = find_free_unit()
    CALL diropn( iungprod, 'wiwjwfc_red', max_ngm*2, exst )
    
    
    allocate(p_basis(max_ngm,numpw))
    allocate(p_basis_r(dfftp%nnr,2))
    allocate(fac(npw))


    if(itask==0) then
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
    else
       fac=1.d0
    endif


    do iw=1,numpw
       call davcio(p_basis(:,iw),max_ngm*2,iungprod,iw,-1)
       do ig=1,npw
          p_basis(ig,iw)=p_basis(ig,iw)*fac(ig)
       enddo
    enddo
    close(iungprod)
    
    iun=find_free_unit()
    if(itask==0) then
       CALL diropn( iun, 'basis2full',dfftp%nnr, exst )
    else
       CALL diropn( iun, 'basis2simple',dfftp%nnr, exst )
    endif

    ii=0
    do iw=1,numpw,2
       psic=0.d0!(1:dfftp%nnr)=(0.d0,0.d0)
       if(iw<numpw) then
          psic(nls(1:npw))  = p_basis(1:npw,iw) + ( 0.D0, 1.D0 ) * p_basis(1:npw,iw+1)
          psic(nlsm(1:npw))  = conjg(p_basis(1:npw,iw) - ( 0.D0, 1.D0 ) * p_basis(1:npw,iw+1))
       else
          psic(nls(1:npw))  = p_basis(1:npw,iw) 
          psic(nlsm(1:npw))  = conjg(p_basis(1:npw,iw))
       endif
       CALL invfft ('Wave', psic, dffts)
       p_basis_r(1:dfftp%nnr,1)=dble(psic(1:dfftp%nnr))
       if(iw<numpw)  p_basis_r(1:dfftp%nnr,2)=dimag(psic(1:dfftp%nnr))
       ii=ii+1
       call davcio(p_basis_r(1,1),dfftp%nnr,iun,ii,1)
       if(iw<numpw) then
          ii=ii+1
          call davcio(p_basis_r(1,2),dfftp%nnr,iun,ii,1)
       endif
    enddo



    close(iun)
    deallocate(p_basis)
    deallocate(p_basis_r)
    deallocate(fac)
    return
  end subroutine write_pola_basis
