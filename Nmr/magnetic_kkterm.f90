!*******************************************************************

subroutine magnetic_kkterm(ik,outvv, outhh)

!*******************************************************************

  use kinds,                only: dp
  use wvfct,               only: npw, npwx, nbnd
  use wavefunctions_module, only: evc
  use klist,                only: wk
  use control_ph,           only: alpha_pv

  implicit none

!  input
  integer, intent(in)      :: ik         ! kpoint

! output

  complex(kind=dp), intent(out) :: outvv(3,3), outhh(3,3)

! local

  integer :: p0,p1,             &  !direction
       ibnd,i
  
  complex(kind=dp), allocatable :: dev(:,:), dpsi(:,:) 
  complex(kind=dp) :: ZDOTC
  

  allocate (dev(npwx,nbnd))
  allocate (dpsi(npwx,nbnd))

! no shift here
  alpha_pv = 0.d0


  do p0=1,3
    
     dev=(0.d0,0.d0)
     call take_nloc_k_kq(ik, ik, evc, p0, dev )
     
     call solve_cg(dev,ik, dpsi)

     do p1=1,3
        dev=(0.d0,0.d0)
        call take_nloc_k_kq(ik, ik, dpsi, p1, dev)
        do ibnd=1,nbnd
           outvv(p0,p1) = outvv(p0,p1) + &
                ZDOTC(npw,evc(1,ibnd),1,dev(1,ibnd),1) &
                * wk(ik) / real(nbnd, dp)
        enddo

        dev=(0.d0,0.d0)
        call grad(ik, dpsi, p1, dev)
        do ibnd=1,nbnd
           outhh(p0,p1) = outhh(p0,p1) + &
                ZDOTC(npw,evc(1,ibnd),1,dev(1,ibnd),1) &
                * wk(ik) / real(nbnd, dp)
        enddo
        
     enddo
  enddo
  
  deallocate(dpsi)
  deallocate(dev)
  

end subroutine magnetic_kkterm
