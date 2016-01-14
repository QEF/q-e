!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

!-----------------------------------------
subroutine rotate_wannier( rot_u,ispin, iun_wannier)
!----------------------------------------
!
! this routine reads the wavefunctions from iun_wannier
! (GAMMA-ONLY CALCULATIONS) and rotate the wavefunctions
! according to rot_u
! only ispin states used (not implemented ye
! ONLY -NORMCONSERVING


  USE kinds,    ONLY : DP
  USE us
  USE wvfct,    ONLY : npwx, npw, nbnd
  USE gvecw,    ONLY : gcutw
  USE gvect
  USE basis
  USE klist
  USE constants, ONLY : e2, pi, tpi, fpi
  USE io_files, ONLY: nwordwfc
  USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2
  USE wavefunctions_module, ONLY: evc 
 
  implicit none

  INTEGER, INTENT(in) :: ispin!+1 or -1
  INTEGER, INTENT(in) :: iun_wannier !units for reading wfc
  REAL(kind=DP), INTENT(in) :: rot_u(nbnd,nbnd)



  COMPLEX(kind=DP), ALLOCATABLE  :: evc0(:,:)!reads wavefunctions here
  COMPLEX(kind=DP), ALLOCATABLE  :: evc1(:,:)!reads wavefunctions here
  integer i,j,k,ig
  INTEGER :: igk0(npwx)
  INTEGER :: npw0
  REAL(kind=dp) :: g2kin_bp(npwx)
  REAL(kind=dp) :: add
  COMPLEX(kind=DP) :: sca
  

  allocate( evc0(npwx,nbnd))
  allocate( evc1(npwx,nbnd))
  
 !reads wfcs from iun_wannier

  CALL gk_sort(xk(1,1),ngm,g,gcutw,npw0,igk0,g2kin_bp) 
  CALL davcio(evc0,2*nwordwfc,iun_wannier,1,-1)

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

  evc(1:npw0,1:nbnd)=evc1(1:npw0,1:nbnd)

  write(*,*) 'writing wannier wfcs on file'!ATTENZIONE
  CALL davcio(evc1,2*nwordwfc,iun_wannier,1,1)
  
  DEALLOCATE(evc0)
  DEALLOCATE(evc1)


  return

end subroutine rotate_wannier


!-----------------------------------------
subroutine rotate_wannier_gamma( rot_u,ispin, itrasp)
!----------------------------------------
!
! (GAMMA-ONLY CALCULATIONS) and rotate the wavefunctions
! according to rot_u
! only ispin states used (not implemented ye
! ONLY -NORMCONSERVING


  USE kinds,    ONLY : DP
  USE us
  USE wvfct,    ONLY : npwx, npw,nbnd
  USE gvect
  USE basis
  USE klist
  USE constants, ONLY : e2, pi, tpi, fpi
  USE io_files, ONLY: nwordwfc
  USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2
  USE wavefunctions_module, ONLY: evc
 

  implicit none

  INTEGER, INTENT(in) :: ispin!+1 or -1
  REAL(kind=DP), INTENT(in) :: rot_u(nbnd,nbnd)
  INTEGER, INTENT(in) :: itrasp!if 1 takes U^T



  REAL(kind=DP), ALLOCATABLE  :: evc0(:,:),evc_re(:,:),evc_im(:,:)!reads wavefunctions here
  integer i,j,k,ig


  allocate( evc0(npw,nbnd))
  allocate( evc_re(npw,nbnd))
  allocate( evc_im(npw,nbnd))


  
  


!now real part
  if(itrasp/=1) then
     evc0(:,:)=dble(evc(:,:))
     call dgemm('N','N',npw,nbnd,nbnd,1.d0,evc0,npw,rot_u,nbnd,0.d0,evc_re,npw)
!now imaginary part
     evc0(:,:)=dimag(evc(:,:))
     call dgemm('N','N',npw,nbnd,nbnd,1.d0,evc0,npw,rot_u,nbnd,0.d0,evc_im,npw)
  else
     evc0(:,:)=dble(evc(:,:))
     call dgemm('N','T',npw,nbnd,nbnd,1.d0,evc0,npw,rot_u,nbnd,0.d0,evc_re,npw)
     !now imaginary part                                                                                                
     evc0(:,:)=dimag(evc(:,:))
     call dgemm('N','T',npw,nbnd,nbnd,1.d0,evc0,npw,rot_u,nbnd,0.d0,evc_im,npw)
  endif

  
!  do i=1,nbnd
!     do ig=1,npw
!        evc(ig,i)=dcmplx(evc_re(ig,i),evc_im(ig,i))
!     enddo
!  enddo


  evc(:,1:nbnd)=dcmplx(evc_re(:,1:nbnd),evc_im(:,1:nbnd))


!rotate
!  do i=1,nbnd
!     do j=1,nbnd
!        do ig=1,npw
!           evc(ig,i)=evc(ig,i)+rot_u(j,i)*evc0(ig,j)
!        enddo
!     enddo
!  enddo




  DEALLOCATE(evc0)
  deallocate(evc_re,evc_im)


  return

end subroutine rotate_wannier_gamma

