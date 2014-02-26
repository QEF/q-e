!                                                                                                                                     
! Copyright (C) 2001-2013 Quantum ESPRESSO group                                                                                      
! This file is distributed under the terms of the                                                                                     
! GNU General Public License. See the file `License'                                                                                  
! in the root directory of the present distribution,                                                                                  
! or http://www.gnu.org/copyleft/gpl.txt .                                                                                            
!                                                                                                                                     
!                             
  subroutine write_wannier_matrix(e_xc,e_h, ispin)
!this subroutine writes the inverse transfromation matrix from KS eigenstates
!to ML wanniers on file, to be read by GWW code
!the INVERSE matrix is calculated here



  

  USE kinds, ONLY : DP
  USE wannier_gw, ONLY : u_trans, num_nbndv, l_selfconsistent,ene_gw,delta_self,n_gw_states
  USE wvfct,    ONLY : et,nbnd
  USE io_global, ONLY : stdout
  USE io_files, ONLY : prefix, tmp_dir
  USE lsda_mod,    ONLY : nspin

  implicit none

  INTEGER, EXTERNAL :: find_free_unit

  REAL(kind=DP) :: e_xc(nbnd,nspin)!exchange and correlation energies
  REAL(kind=DP) :: e_h(nbnd,nspin)!hartree energies
  INTEGER, INTENT(in) :: ispin!spin channel

  COMPLEX(kind=DP) :: sca

  INTEGER :: iunu, iw,jw
  INTEGER :: ivpt(nbnd), info
  COMPLEX(kind=DP) :: cdet(2),det
  COMPLEX(kind=DP), ALLOCATABLE :: cdwork(:)

  REAL(kind=DP), ALLOCATABLE :: et_new(:)
  INTEGER :: is
  
  do iw=1,nbnd
     do jw=iw,nbnd
        sca=u_trans(iw,jw,ispin)
        u_trans(iw,jw,ispin)=conjg(u_trans(jw,iw,ispin))
        u_trans(jw,iw,ispin)=conjg(sca)
     enddo
  enddo
  

  iunu = find_free_unit()
  
  open(unit=iunu,file=trim(tmp_dir)//trim(prefix)//'.wannier',status='unknown',form='unformatted')

  write(iunu) nspin
  write(iunu) nbnd

  do is=1,nspin
     write(iunu) num_nbndv(is)

     if(.not.l_selfconsistent) then
        write(iunu) et(1:nbnd,is)
     else
        allocate(et_new(nbnd))
        et_new(1:n_gw_states)=ene_gw(1:n_gw_states,is)
        if(nbnd>n_gw_states) et_new(n_gw_states+1:nbnd)=et(n_gw_states+1:nbnd,is)+delta_self
        write(iunu) et_new(1:nbnd)
        deallocate(et_new)
     endif
     if(l_selfconsistent) e_xc(:,is)=0.d0
     write(iunu) e_xc(1:nbnd,is)
     write(iunu) e_h(1:nbnd,is)


     do iw=1,nbnd
        write(iunu) u_trans(1:nbnd,iw,is)
     enddo
  enddo
  close(iunu)


  return
  end subroutine
  

 subroutine read_wannier_matrix
!this read the inverse transfromation matrix from KS eigenstates
!to ML wanniers on file, to be read by GWW code
!the INVERSE matrix is calculated here





  USE kinds, ONLY : DP
  USE wannier_gw, ONLY : u_trans, num_nbndv
  USE wvfct,    ONLY : et,nbnd
  USE io_global, ONLY : stdout,ionode,ionode_id
  USE io_files, ONLY : prefix, tmp_dir
  USE mp, ONLY : mp_bcast
  USE mp_world, ONLY : world_comm
  USE lsda_mod, ONLY :nspin

  implicit none

  INTEGER, EXTERNAL :: find_free_unit

  INTEGER :: iunu, iw, is
  INTEGER :: idumm
  REAL(kind=DP), ALLOCATABLE :: rdummv(:)
  
  allocate(rdummv(nbnd))


  if(ionode) then
     iunu = find_free_unit()
     open(unit=iunu,file=trim(tmp_dir)//trim(prefix)//'.wannier',status='old',form='unformatted')
  
     read(iunu) idumm
     read(iunu) idumm
  endif
  do is=1,nspin
     if(ionode) then
        read(iunu) idumm
        read(iunu) rdummv(1:nbnd)
        read(iunu) rdummv(1:nbnd)
        read(iunu) rdummv(1:nbnd)
     endif
 

     do iw=1,nbnd
        if(ionode) read(iunu) u_trans(1:nbnd,iw,is)
        call mp_bcast(u_trans(1:nbnd,iw,is),ionode_id,world_comm)       
     enddo
  enddo
  if(ionode) close(iunu)
  
  deallocate(rdummv)


  return
end subroutine read_wannier_matrix
