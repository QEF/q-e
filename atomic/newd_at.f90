!
!-------------------------------------------------------------------------
       subroutine newd_at
!-------------------------------------------------------------------------
!
!     this ruotine computes the new D coeeficients
!
!
use ld1inc

       integer :: &
             ib,jb,n,is,nst

       real(kind=dp) :: &
            int_0_inf_dr, &   ! the integral function
            gi(ndm)          ! the gi function

!
!    screening the D coefficients
!
       if (pseudotype.eq.3) then
          do ib=1,nbeta
             do jb=1,ib
                if (lls(ib).eq.lls(jb).and.abs(jjs(ib)-jjs(jb)).lt.1d-7) then
                   nst=(lls(ib)+1)*2
                   do is=1,nspin
                      do n=1,ikk(ib)
                         gi(n)=qvan(n,ib,jb)*vpstot(n,is)
                      enddo
                      ddd(ib,jb,is)= bmat(ib,jb) &
                              + int_0_inf_dr(gi,r,r2,dx,ikk(ib),nst)
                      ddd(jb,ib,is)=ddd(ib,jb,is)
                   enddo
                endif
             enddo
          enddo
       endif

       return
       end
