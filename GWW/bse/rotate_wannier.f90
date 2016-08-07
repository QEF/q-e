!-----------------------------------------
subroutine rotate_wannier_gamma_bse( rot_u,a_in,a_out,ispin, itrasp)
!----------------------------------------
!
! (GAMMA-ONLY CALCULATIONS) and rotate the wavefunctions
! according to rot_u
! only ispin states used (not implemented ye
! ONLY -NORMCONSERVING


  USE kinds,    ONLY : DP
  USE us
  USE wvfct,    ONLY : g2kin, npwx, npw, nbndx,nbnd
  USE gvect
  USE basis
  USE klist
  USE constants, ONLY : e2, pi, tpi, fpi
  USE io_files, ONLY: nwordwfc
  USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2
!  USE wavefunctions_module, ONLY: evc
  use exciton
  USE io_global, ONLY : stdout 


  implicit none

  INTEGER, INTENT(in) :: ispin!+1 or -1
   type(exc):: a_in, a_out
  REAL(kind=DP) :: rot_u(a_in%numb_v,a_in%numb_v)
  INTEGER, INTENT(in) :: itrasp!if 1 takes U^T
 



  REAL(kind=DP), ALLOCATABLE  :: evc0(:,:),evc_re(:,:),evc_im(:,:)!reads wavefunctions here
  integer i,j,k,ig
  logical debug

  call start_clock('rotate_wannier_gamma_bse')
  debug=.false.

  allocate( evc0(npw,a_in%numb_v))
  allocate( evc_re(npw,a_in%numb_v))
  allocate( evc_im(npw,a_in%numb_v))

  if(debug) then 
     write(stdout,*) 'rotate wannier #1'
  endif

!now real part
  if(itrasp/=1) then
     evc0(1:a_in%npw,1:a_in%numb_v)=dble(a_in%a(1:a_in%npw,1:a_in%numb_v))
     call dgemm('N','N',npw,a_in%numb_v,a_in%numb_v,1.d0,evc0,npw,rot_u,a_in%numb_v,0.d0,evc_re,npw)
!now imaginary part
     evc0(1:a_in%npw,1:a_in%numb_v)=dimag(a_in%a(1:a_in%npw,1:a_in%numb_v))
     call dgemm('N','N',npw,a_in%numb_v,a_in%numb_v,1.d0,evc0,npw,rot_u,a_in%numb_v,0.d0,evc_im,npw)
  else
     evc0(1:a_in%npw,1:a_in%numb_v)=dble(a_in%a(1:a_in%npw,1:a_in%numb_v))
     call dgemm('N','T',npw,a_in%numb_v,a_in%numb_v,1.d0,evc0,npw,rot_u,a_in%numb_v,0.d0,evc_re,npw)
     !now imaginary part                                                                                                
     evc0(1:a_in%npw,1:a_in%numb_v)=dimag(a_in%a(1:a_in%npw,1:a_in%numb_v))
     call dgemm('N','T',npw,a_in%numb_v,a_in%numb_v,1.d0,evc0,npw,rot_u,a_in%numb_v,0.d0,evc_im,npw)
  endif

  
!  do i=1,nbnd
!     do ig=1,npw
!        evc(ig,i)=dcmplx(evc_re(ig,i),evc_im(ig,i))
!     enddo
!  enddo


 a_out%a(1:a_in%npw,1:a_in%numb_v)=dcmplx(evc_re(1:a_in%npw,1:a_in%numb_v),evc_im(1:a_in%npw,1:a_in%numb_v))
  


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


  call stop_clock('rotate_wannier_gamma_bse')
  return

end subroutine rotate_wannier_gamma_bse

