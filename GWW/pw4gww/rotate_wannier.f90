! FOR GWW
!
! Author: P. Umari
! Modified by G. Stenuit and L. Martin-Samos
!
!-----------------------------------------
subroutine rotate_wannier( rot_u,ispin)
!----------------------------------------
!
! this routine reads the wavefunctions from iun_wannier
! (GAMMA-ONLY CALCULATIONS) and rotate the wavefunctions
! according to rot_u
! only ispin states used (not implemented ye
! ONLY -NORMCONSERVING

!#ifdef __GWW

  USE kinds,    ONLY : DP
  USE io_global, ONLY : stdout
  USE us
  USE wvfct,    ONLY : igk, g2kin, npwx, npw, nbnd, nbndx, ecutwfc
  USE gvect
  USE basis
  USE klist
  USE constants, ONLY : e2, pi, tpi, fpi
  USE io_files, ONLY: nwordwfc, iunwfc,prefix, find_free_unit, diropn
  USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2
! same stuff investigation
  USE wavefunctions_module, ONLY: evc

  implicit none

  INTEGER, INTENT(in) :: ispin!+1 or -1
!  INTEGER, INTENT(in) :: iun_wannier !units for reading wfc
  REAL(kind=DP), INTENT(in) :: rot_u(nbnd,nbnd)

! modified sym as rotate_wannier_gamma
  integer :: iun_wannier

  COMPLEX(kind=DP), ALLOCATABLE  :: evc0(:,:)!reads wavefunctions here
  COMPLEX(kind=DP), ALLOCATABLE  :: evc1(:,:)!reads wavefunctions here
  integer i,j,k,ig
  INTEGER :: igk0(npwx)
  INTEGER :: npw0
  REAL(kind=dp) :: g2kin_bp(npwx)
  REAL(kind=dp) :: add
  COMPLEX(kind=DP) :: sca

  INTEGER :: ik

  LOGICAL :: exst


  allocate( evc0(npwx,nbnd))
  allocate( evc1(npwx,nbnd))

 !reads wfcs from iun_wannier

  iun_wannier = find_free_unit()
  call diropn(iun_wannier,"wfc_w",2*nwordwfc,exst)

  do ik=1,nks

    CALL gk_sort(xk(1,ik),ngm,g,ecutwfc/tpiba2, &
                &    npw0,igk0,g2kin_bp)
    CALL davcio(evc,nwordwfc,iunwfc,ik,-1)

! maybe not necessary as evc is still in mem
!  CALL davcio(evc0,nwordwfc,iun_wannier,1,-1)
! added just association evc0 with evc
    evc0(:,:) = evc(:,:)

    evc1=(0.d0,0.d0)


!rotate
    do i=1,nbnd
      do j=1,nbnd
         do ig=1,npw0
           evc1(ig,i)=evc1(ig,i)+rot_u(j,i)*evc0(ig,j)
         enddo
      enddo
    enddo


!check for debug


! do i=1,nbnd
!      do j=1,nbnd
!         sca=(0.d0,0.d0)
!         do ig=1,npw0
!           sca=sca+conjg(evc1(ig,i))*evc1(ig,j)
!         enddo
!         write(*,*) 'rotata_wannier_check :', i,j, sca
!      enddo
! enddo

!write back on file


! this is strange to be investigated
!  evc(1:npw0,1:nbnd)=evc1(1:npw0,1:nbnd)

!-------------------------------------------------------
! writing wannierized evc --> evc1 in file
  write(stdout,*) 'writing wannier wfcs on file'!ATTENZIONE

      CALL davcio(evc1,nwordwfc,iun_wannier,ik,1)
  enddo
  close(iun_wannier)
! -----------------------------------------------------


  DEALLOCATE(evc0)
  DEALLOCATE(evc1)

!#endif

  return

end subroutine rotate_wannier


!-----------------------------------------
subroutine rotate_wannier_gamma( rot_u,ispin, itrasp )
!----------------------------------------
!
! (GAMMA-ONLY CALCULATIONS) and rotate the wavefunctions
! according to rot_u
! only ispin states used (not implemented ye
! ONLY -NORMCONSERVING

!#ifdef  __GWW

  USE kinds,    ONLY : DP
  USE io_global, ONLY : stdout
  USE us
  USE wvfct,    ONLY : igk, g2kin, npwx, npw, nbndx
  USE gvect
  USE basis
  USE klist
  USE constants, ONLY : e2, pi, tpi, fpi
  USE io_files, ONLY: nwordwfc, prefix, find_free_unit, diropn
  USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2
  USE wavefunctions_module, ONLY: evc
  USE wannier_gw,  ONLY : nbnd_normal

  implicit none

  INTEGER, INTENT(in) :: ispin!+1 or -1
  REAL(kind=DP), INTENT(in) :: rot_u(nbnd_normal,nbnd_normal)
  INTEGER, INTENT(in) :: itrasp!if 1 takes U^T

! add lecture unit of wannierized evc --> evc1
  integer :: iun_wannier
  logical :: exst
!----------------------------------------------
! add evc1
  COMPLEX(kind=DP), ALLOCATABLE  :: evc1(:,:)!reads wavefunctions here
!--------------------------------------
  REAL(kind=DP), ALLOCATABLE  :: evc0(:,:),evc_re(:,:),evc_im(:,:)!reads wavefunctions here
  integer i,j,k,ig


  allocate( evc0(npw,nbnd_normal))
  allocate( evc_re(npw,nbnd_normal))
  allocate( evc_im(npw,nbnd_normal))
!-----------------------------------------
! add allocation
  allocate( evc1(npw,nbnd_normal))
!----------------------------------------






!now real part
  if(itrasp/=1) then
     evc0(:,:)=dble(evc(:,:))
     call dgemm('N','N',npw,nbnd_normal,nbnd_normal,1.d0,evc0,npw,rot_u,nbnd_normal,0.d0,evc_re,npw)
!now imaginary part
     evc0(:,:)=dimag(evc(:,:))
     call dgemm('N','N',npw,nbnd_normal,nbnd_normal,1.d0,evc0,npw,rot_u,nbnd_normal,0.d0,evc_im,npw)
  else
     evc0(:,:)=dble(evc(:,:))
     call dgemm('N','T',npw,nbnd_normal,nbnd_normal,1.d0,evc0,npw,rot_u,nbnd_normal,0.d0,evc_re,npw)
     !now imaginary part
     evc0(:,:)=dimag(evc(:,:))
     call dgemm('N','T',npw,nbnd_normal,nbnd_normal,1.d0,evc0,npw,rot_u,nbnd_normal,0.d0,evc_im,npw)
  endif


!  do i=1,nbnd
!     do ig=1,npw
!        evc(ig,i)=dcmplx(evc_re(ig,i),evc_im(ig,i))
!     enddo
!  enddo

! ------------------------------------------------
! commenting this evc should not be modified
!  evc(:,1:nbnd_normal)=dcmplx(evc_re(:,1:nbnd_normal),evc_im(:,1:nbnd_normal))

  evc1(:,1:nbnd_normal)=dcmplx(evc_re(:,1:nbnd_normal),evc_im(:,1:nbnd_normal))

!-------------------------------------------------------
! added davcio because it wasn't present in the original rotate_wannier_gamma
! maybe needed
  write(stdout,*) 'writing wannier wfcs on file'!ATTENZIONE

  iun_wannier = find_free_unit()

  call diropn(iun_wannier,"wfc_w",2*nwordwfc, exst)
  CALL davcio(evc1,nwordwfc,iun_wannier,1,1)
  close(iun_wannier)
! -----------------------------------------------------

!rotate
!  do i=1,nbnd
!     do j=1,nbnd
!        do ig=1,npw
!           evc(ig,i)=evc(ig,i)+rot_u(j,i)*evc0(ig,j)
!        enddo
!     enddo
!  enddo




  DEALLOCATE(evc0)
! add deallocation of evc1
  DEALLOCATE(evc1)
  deallocate(evc_re,evc_im)


!#endif

  return

end subroutine rotate_wannier_gamma

