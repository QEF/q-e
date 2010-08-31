! FOR GWW
!
! Author: P. Umari
!
  subroutine write_wannier_matrix_c
!this subroutine writes the inverse transfromation matrix from KS eigenstates
!to ML wanniers on file, to be read by GWW code
!the INVERSE matrix is calculated here

! #ifdef __GWW



  USE kinds,      ONLY : DP
  USE wannier_gw, ONLY : u_trans, num_nbndv, lnonorthogonal, num_nbndc_set, nbnd_normal
  USE wvfct,      ONLY : nbnd,et
  USE io_global,  ONLY : stdout
  USE io_files,   ONLY : find_free_unit, prefix

  implicit none

  COMPLEX(kind=DP) :: sca
  INTEGER :: iunu, iw,jw
  INTEGER :: ivpt(num_nbndc_set), info
  COMPLEX(kind=DP) :: cdet(2),det
  COMPLEX(kind=DP), ALLOCATABLE :: cdwork(:)
  COMPLEX(kind=DP), ALLOCATABLE :: u_trans_c(:,:)

  allocate(u_trans_c(num_nbndc_set,num_nbndc_set))

  do iw=1,num_nbndc_set
     do jw=1,num_nbndc_set
        u_trans_c(iw,jw)=u_trans(iw+num_nbndv,jw+num_nbndv)
     enddo
  enddo

  if(.not.lnonorthogonal) then
     do iw=1,num_nbndc_set
        do jw=iw,num_nbndc_set
           sca=u_trans_c(iw,jw)
           u_trans_c(iw,jw)=conjg(u_trans_c(jw,iw))
           u_trans_c(jw,iw)=conjg(sca)
        enddo
     enddo
  else
     allocate(cdwork(num_nbndc_set))
     CALL zgefa(u_trans_c,num_nbndc_set,num_nbndc_set,ivpt,info)
     CALL errore('write_wannier_matrix','error in zgefa',abs(info))
     CALL zgedi(u_trans_c,num_nbndc_set,num_nbndc_set,ivpt,cdet,cdwork,11)
     det=cdet(1)*10.d0**cdet(2)
     write(stdout,*) 'DETERMINANT OF A MATRIX:', det
     deallocate(cdwork)
  endif

  iunu = find_free_unit()

  open(unit=iunu,file=trim(prefix)//'.wannier_prim',status='unknown',form='unformatted')


  write(iunu) num_nbndc_set
  write(iunu) num_nbndv
  write(iunu) nbnd_normal


  do iw=1,num_nbndc_set
     write(iunu) u_trans_c(1:num_nbndc_set,iw)
  enddo

  close(iunu)

  deallocate(u_trans_c)

! #endif
  return
end subroutine write_wannier_matrix_c


