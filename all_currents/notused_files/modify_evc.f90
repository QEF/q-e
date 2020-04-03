subroutine modify_evc(evc_mod,axes)
!
USE wvfct,                ONLY : npw, npwx,nbnd
USE uspp,                 ONLY : nkb
USE ions_base,            ONLY : nat
USE wavefunctions_module, ONLY :psic,evc
USE fft_interfaces,       ONLY : fwfft, invfft
USE fft_base, ONLY:dffts
USE mp_world, ONLY:mpime
USE gvecs,              ONLY : nls, nlsm
USE cell_base,  ONLY:at,alat
USE zero_mod
use ions_base, only :tau,ityp
USE uspp_param,       ONLY : nh
USE atom, ONlY :rgrid
!
implicit none
!
integer, intent(in):: axes
complex(DP), intent(out):: evc_mod(npwx,nbnd,nat)
!
!list contiene per ogni proiettore nkb l'atomo attorno al quale è centrata.
real(DP) :: evc_r(1:dffts%nnr)
real(DP) ::u(3),u_x(3),u_y(3),u_z(3),u_mid(3),modul
integer ::iqq,n_x,n_y,n_z,ix,iy,iz,na,jkb,ijkb,ikb,iat,ibnd
real :: cutoff
!
!
do ibnd=1,nbnd
   psic=0.d0
   psic(nls(1:npw))=evc_uno(1:npw,ibnd)
   psic(nlsm(1:npw))=CONJG(evc(1:npw,ibnd))
   call invfft ('Wave', psic, dffts)  
   evc_r(1:dffts%nnr)=psic(1:dffts%nnr)
   do iat=1,nat
      do iz=1,dffts%nr3p(mpime+1)     
         do iy=1,dffts%nr2           
            do ix=1,dffts%nr1                 
               iqq=(iz-1)*(dffts%nr1x*dffts%nr2x)+(iy-1)*dffts%nr1+ix
               u_x(1:3)=real(ix-1)/real(dffts%nr1)*at(1:3,1)*alat
               u_y(1:3)=real(iy-1)/real(dffts%nr2)*at(1:3,2)*alat
               u_z(1:3)=real(iz+nr3s_start-1-1)/real(dffts%nr3)*at(1:3,3)*alat
               u(1:3) = u_x(1:3)+u_y(1:3)+u_z(1:3)                   
               call mid(u(1:3)-alat*tau(1:3,iat),u_mid(1:3))
               evc_r(iqq)=evc_r(iqq)*u_mid(axes)
            end do
         end do
      end do
      psic=0.d0
      psic(1:dffts%nnr)=evc_r(1:dffts%nnr)
      call fwfft ('Wave', psic, dffts)
      evc_mod(1:npw,ibnd,iat)=psic(nls(1:npw))
    end do
end do
!
end subroutine modify_evc
 
