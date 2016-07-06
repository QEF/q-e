 subroutine read_wannier_matrix
!this read the inverse transfromation matrix from KS eigenstates
!to ML wanniers on file, to be read by GWW code
!the INVERSE matrix is calculated here

  USE kinds, ONLY : DP
!  USE wannier_gw, ONLY : u_trans, num_nbndv
  USE wvfct,    ONLY : et,nbnd
  USE io_global, ONLY : stdout,ionode,ionode_id
!  USE io_files, ONLY : find_free_unit, prefix
  USE io_files, ONLY : prefix, tmp_dir
  USE mp, ONLY : mp_bcast
  USE mp_world,             ONLY : world_comm
  USE lsda_mod, ONLY :nspin
  USE bse_basic_structures, ONLY : u_trans
  use bse_wannier, ONLY:num_nbndv
 
  implicit none
  INTEGER, EXTERNAL :: find_free_unit



  INTEGER :: iunu, iw, is
  INTEGER :: idumm
  REAL(kind=DP), ALLOCATABLE :: rdummv(:)
  
  call start_clock('read_wannier_matrix')
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
        call mp_bcast(u_trans(1:nbnd,iw,is),ionode_id, world_comm)       
     enddo
!DEBUG
!     u_trans=0.d0
!     do iw=1,nbnd
!        u_trans(iw,iw,is)=1.d0
!     enddo

  enddo
  if(ionode) close(iunu)
  
  deallocate(rdummv)


  call stop_clock('read_wannier_matrix')
  return
end subroutine read_wannier_matrix
