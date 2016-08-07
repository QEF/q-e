!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

!this routines calculates the terms <\psi_i|s_\alpha> where s_alpha is the global
!s_basis and write them on disk
!the KS states are taken on evc

  subroutine contour_terms(n_s,s_basis,ispin,istate)
!NOT_TO_BE_INCLUDED_START
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE io_files,             ONLY : prefix, tmp_dir, nwordwfc,iunwfc
    USE kinds,    ONLY : DP
    USE wannier_gw, ONLY : num_nbnds,num_nbndv,s_first_state,s_last_state, l_verbose
    USE wvfct,    ONLY : npwx, npw, nbnd
    USE mp, ONLY : mp_sum, mp_barrier, mp_bcast
    USE mp_world, ONLY : mpime,nproc,world_comm
    USE wavefunctions_module, ONLY : evc
    USE gvect,    ONLY : gstart
     
    implicit none

    INTEGER, EXTERNAL :: find_free_unit
    INTEGER, INTENT(in) :: n_s!number of global s vectors
    COMPLEX(kind=DP), INTENT(in) :: s_basis( npw,n_s)!s vectors
    INTEGER, INTENT(in) :: ispin!spin channel
    INTEGER, INTENT(in) :: istate!KS states relative to global s vectors for big_system option

    REAL(kind=DP), ALLOCATABLE :: c_mat(:,:)
    INTEGER :: ii,jj,iun,ig
    CHARACTER(4) :: nfile


    allocate(c_mat(n_s,num_nbnds))

    call dgemm('T','N',n_s,num_nbnds,2*npw,2.d0,s_basis,2*npw,evc,2*npwx,0.d0,c_mat,n_s)
       
    if(gstart==2) then
       do ii=1,n_s
          do jj=1,num_nbnds
             c_mat(ii,jj)= c_mat(ii,jj)-dble(s_basis(1,ii)*conjg(evc(1,jj)))
          enddo
       enddo
    endif

!DEBUG
       

    call mp_sum(c_mat,world_comm)

    if(ionode) then
       iun= find_free_unit()
        write(nfile,'(4i1)') istate/1000,mod(istate,1000)/100,mod(istate,100)/10,mod(istate,10)
       if(ispin==1) then
          open( unit= iun, file=trim(tmp_dir)//trim(prefix)//'.s_contour'//nfile , status='unknown',form='unformatted')
       else
          open( unit= iun, file=trim(tmp_dir)//trim(prefix)//'.s_contour2'//nfile , status='unknown',form='unformatted')
       endif
       write(iun) num_nbnds
       write(iun) n_s
       do jj=1,num_nbnds
          write(iun) c_mat(1:n_s,jj)!GIUSTO CUSSI
          !write(iun) c_mat(1:n_s,4)!ATTENZIONE DEBUG
       enddo
       close(iun)
    endif


    deallocate(c_mat)
    return

!NOT_TO_BE_INCLUDED_END
  end subroutine contour_terms
