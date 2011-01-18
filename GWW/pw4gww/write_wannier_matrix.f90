! FOR GWW
!
! Author: P. Umari
! Modified by G. Stenuit and L. Martin-Samos
!
  subroutine write_wannier_matrix(e_xc,e_h)
!this subroutine writes the inverse transfromation matrix from KS eigenstates
!to ML wanniers on file, to be read by GWW code
!the INVERSE matrix is calculated here

! #ifdef __GWW



  USE kinds, ONLY : DP
  USE wannier_gw, ONLY : u_trans, num_nbndv, lnonorthogonal,nbnd_normal
  USE wvfct,    ONLY : et
  USE io_global, ONLY : stdout
  USE io_files, ONLY : find_free_unit, prefix

  implicit none

  REAL(kind=DP) :: e_xc(nbnd_normal)!exchange and correlation energies
  REAL(kind=DP) :: e_h(nbnd_normal)!hartree energies
  COMPLEX(kind=DP) :: sca

  INTEGER :: iunu, iw,jw
  INTEGER :: ivpt(nbnd_normal), info
  COMPLEX(kind=DP) :: cdet(2),det
  COMPLEX(kind=DP), ALLOCATABLE :: cdwork(:)


  write(stdout,*) "nbnb_normal:",nbnd_normal
  write(stdout,*) "ubound e_h e_xc", ubound(e_h), ubound(e_xc)
  write(stdout,*) "ubound u_trans", ubound(u_trans,1), ubound(u_trans,2)
  call flush_unit(stdout)

  if(.not.lnonorthogonal) then
    do iw=1,nbnd_normal
      do jw=iw,nbnd_normal
        sca=u_trans(iw,jw)
        u_trans(iw,jw)=conjg(u_trans(jw,iw))
        u_trans(jw,iw)=conjg(sca)
      enddo
    enddo
  else
    allocate(cdwork(nbnd_normal))
    write(stdout,*) "before zgefa"
    call flush_unit(stdout)
    CALL zgefa(u_trans,nbnd_normal,nbnd_normal,ivpt,info)
    CALL errore('write_wannier_matrix','error in zgefa',abs(info))
    write(stdout,*) "before zgedi"
    call flush_unit(stdout)
    CALL zgedi(u_trans,nbnd_normal,nbnd_normal,ivpt,cdet,cdwork,11)
    det=cdet(1)*10.d0**cdet(2)
    write(stdout,*) 'DETERMINANT OF A MATRIX:', det
    deallocate(cdwork)
  endif

  iunu = find_free_unit()

  open(unit=iunu,file=trim(prefix)//'.wannier',status='unknown',form='unformatted')

  write(iunu) nbnd_normal
  write(iunu) num_nbndv

  write(iunu) et(1:nbnd_normal,1)
  write(iunu) e_xc(1:nbnd_normal)
  write(iunu) e_h(1:nbnd_normal)


  do iw=1,nbnd_normal
     write(iunu) u_trans(1:nbnd_normal,iw)

  enddo

  close(iunu)

! #endif
  return
  end subroutine


