!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

subroutine qqberry2( gqq,gqqm, ipol)

!  this subroutine computes the array gqq and gqqm
!  gqq=int_dr qq(r)exp(iGr)=<Beta_r|exp(iGr)|Beta_r'>
!  gqqm=int_dr qq(r)exp(-iGr)=<Beta_r|exp(-iGr)|Beta_r'>

!   gqq output: as defined above

  use smallbox_grid_dimensions, only: nr1b, nr2b, nr3b, &
            nr1bx, nr2bx, nr3bx, nnrb => nnrbx
  use uspp_param, only: lqmax, nqlc, kkbeta, nbeta, nh, nhm
  use uspp, only: indv, lpx, lpl, ap,nhtolm
  use atom, only: r, rab
  use core
  use gvecw, only: ngw
  use reciprocal_vectors, only: mill_l
  use parameters
  use  constants
  use cvan, only: oldvan, nvb
  use  ions_base
  use ions_base, only : nas => nax
  use cell_base, only: a1, a2, a3
  use reciprocal_vectors, only: ng0 => gstart, gx, g
  use mp, only: mp_sum
  
  implicit none

  complex(8) gqq(nhm,nhm,nas,nsp)
  complex(8) gqqm(nhm,nhm,nas,nsp)
  real(8) gmes
  integer ipol, lx

! local variables

  integer ig, is, iv, jv, i, istart, il,l,ir, igi,ia
  real(8), allocatable:: fint(:),jl(:)
  real(8), allocatable:: qrl(:,:,:), qradb2(:,:,:,:) 
  real(8) c, xg
  complex(8) qgbs,sig
  integer ivs, jvs, ivl, jvl, lp, ijv
  real(8), allocatable:: ylm(:,:)

  lx = lqmax

  allocate( fint( ndmx))
  allocate( jl(ndmx))
  allocate( qradb2(nbrx,nbrx,lx,nsp))
  allocate( ylm(ngw, lqmax*lqmax))

  CALL ylmr2( lqmax*lqmax, ngw, gx, g, ylm )

  qradb2 = 0.0d0
     
  do is=1,nsp
     do ia=1,nas
        do jv=1,nhm
           do iv=1,nhm
              gqq(iv,jv,ia,is)=(0.,0.)
              gqqm(iv,jv,ia,is)=(0.,0.)
           enddo
        enddo
     enddo
  enddo
  
  if(ipol.eq.1) then
     gmes=a1(1)**2+a1(2)**2+a1(3)**2
     gmes=2*pi/SQRT(gmes)
  endif
  if(ipol.eq.2) then
     gmes=a2(1)**2+a2(2)**2+a2(3)**2
     gmes=2*pi/SQRT(gmes)
  endif
  if(ipol.eq.3) then
     gmes=a3(1)**2+a3(2)**2+a3(3)**2
     gmes=2*pi/SQRT(gmes)
  endif    
  ! only for Vanderbilt species 
  do is=1,nvb
     c=fpi                 !/omegab
     !
     ALLOCATE ( qrl(kkbeta(is), nbeta(is)*(nbeta(is)+1)/2, nqlc(is)) )
     !
     call fill_qrl (is, SIZE(qrl, 1), SIZE(qrl, 2), SIZE(qrl, 3), qrl)
     ! now the radial part
     do l=1,nqlc(is)                        
        xg= gmes !only orthorombic cells
        call bess(xg,l,kkbeta(is),r(1,is),jl)
        ijv = 0
        do iv= 1,nbeta(is)
           do jv=iv,nbeta(is)
              ijv = ijv + 1
!     
!     note qrl(r)=r^2*q(r)
!
              do ir=1,kkbeta(is)
                 fint(ir)=qrl(ir,ijv,l)*jl(ir)
              end do
              if (oldvan(is)) then
                 call herman_skillman_int     &
                      &  (kkbeta(is),fint,rab(1,is),qradb2(iv,jv,l,is))
              else
                 call simpson (kkbeta(is),fint,rab(1,is),qradb2(iv,jv,l,is))
              endif
              qradb2(iv,jv,l,is)=  c*qradb2(iv,jv,l,is)
              qradb2(jv,iv,l,is)=  qradb2(iv,jv,l,is)
           end do
        end do
     end do
     DEALLOCATE ( qrl )    
  enddo

  igi=-1
  do ig=1,ngw
     if(ipol.eq.1 ) then
        if(mill_l(1,ig).eq.1 .and. mill_l(2,ig).eq.0  .and. mill_l(3,ig).eq. 0) igi=ig
     endif
     if(ipol.eq.2 ) then
        if(mill_l(1,ig).eq.0 .and. mill_l(2,ig).eq.1  .and. mill_l(3,ig).eq. 0) igi=ig
     endif
     if(ipol.eq.3 ) then
        if(mill_l(1,ig).eq.0 .and. mill_l(2,ig).eq.0   .and. mill_l(3,ig).eq. 1) igi=ig
     endif
  enddo
  if( igi.ne.-1) then

!setting array beigr

      
             
     do is=1,nvb
        do iv= 1,nh(is)
           do jv=iv,nh(is)
              ivs=indv(iv,is)
              jvs=indv(jv,is)
              ivl=nhtolm(iv,is)
              jvl=nhtolm(jv,is)
!
!     lpx = max number of allowed y_lm
!     lp  = composite lm to indentify them
!
              qgbs=(0.,0.)
              do i=1,lpx(ivl,jvl)
                 lp=lpl(ivl,jvl,i)
!
!     extraction of angular momentum l from lp:  
!
                 if (lp.eq.1) then
                    l=1         
                 else if ((lp.ge.2) .and. (lp.le.4)) then
                    l=2
                 else if ((lp.ge.5) .and. (lp.le.9)) then
                    l=3
                 else if ((lp.ge.10).and.(lp.le.16)) then
                    l=4
                 else if ((lp.ge.17).and.(lp.le.25)) then
                    l=5
                 else if (lp.ge.26) then 
                    call errore(' qvanb ',' lp.ge.26 ',lp)
                 endif
!
!       sig= (-i)^l
!
                 sig=(0.,-1.)**(l-1)
                  
                 sig=sig*ap(lp,ivl,jvl)
                 qgbs=qgbs+sig*ylm(igi,lp)*qradb2(ivs,jvs,l,is)
                
              end do
   

              
              do ia=1,na(is)
                     
                 gqqm(iv,jv,ia,is)=qgbs
                 gqqm(jv,iv,ia,is)=qgbs
                 gqq(iv,jv,ia,is)=CONJG(gqqm(iv,jv,ia,is))
                 gqq(jv,iv,ia,is)=CONJG(gqqm(iv,jv,ia,is))
              end do
           end do
        enddo
     enddo
  endif

  call mp_sum(gqq(:,:,:,:))
  call mp_sum(gqqm(:,:,:,:))

  deallocate( fint)
  deallocate( jl)
  deallocate(qradb2)
  deallocate(ylm)
  
  return
end subroutine qqberry2






! this subroutine updates gqq and gqqm to the 
! (new) atomic position


subroutine qqupdate(eigr, gqqm0, gqq, gqqm, ipol)

!   gqq output: as defined above

  use cvan
  use gvecw, only: ngw
  use ions_base, only : nas => nax, nat, na, nsp
  use reciprocal_vectors, only: mill_l
  use uspp_param, only: nh, nhm
  use mp, only: mp_sum

  implicit none

 
  complex(8) eigr(ngw,nat)
  complex(8) gqq(nhm,nhm,nas,nsp)
  complex(8) gqqm(nhm,nhm,nas,nsp)
  complex(8) gqqm0(nhm,nhm,nas,nsp)

  integer ipol
  
  integer igi,ig,is,iv,jv,ia,isa


  do is=1,nsp
     do ia=1,nas
        do jv=1,nhm
           do iv=1,nhm
              gqq(iv,jv,ia,is)=(0.,0.)
              gqqm(iv,jv,ia,is)=(0.,0.)
           enddo
        enddo
     enddo
  enddo

  igi=-1
  do ig=1,ngw
     if(ipol.eq.1 ) then
        if(mill_l(1,ig).eq.1 .and. mill_l(2,ig).eq.0  .and. mill_l(3,ig).eq. 0) igi=ig
     endif
     if(ipol.eq.2 ) then
        if(mill_l(1,ig).eq.0 .and. mill_l(2,ig).eq.1  .and. mill_l(3,ig).eq. 0) igi=ig
     endif
     if(ipol.eq.3 ) then
        if(mill_l(1,ig).eq.0 .and. mill_l(2,ig).eq.0  .and. mill_l(3,ig).eq. 1) igi=ig
     endif
  enddo
  if( igi.ne.-1) then

  
     isa = 1
     do is=1,nvb
        do ia=1,na(is)
           do iv= 1,nh(is)
              do jv=iv,nh(is)
                 gqqm(iv,jv,ia,is)= gqqm0(iv,jv,ia,is)*eigr(igi,isa)
                 gqqm(jv,iv,ia,is)= gqqm0(iv,jv,ia,is)*eigr(igi,isa)
                 gqq(iv,jv,ia,is)=CONJG(gqqm(iv,jv,ia,is))
                 gqq(jv,iv,ia,is)=CONJG(gqqm(iv,jv,ia,is))
              enddo
           enddo
           isa = isa + 1
        enddo
     enddo
  endif
  call mp_sum(gqq(:,:,:,:))
  call mp_sum(gqqm(:,:,:,:))
  return
end subroutine qqupdate
      
!-----------------------------------------------------------------------
      subroutine bess(xg,l,mmax,r,jl)
!-----------------------------------------------------------------------
!     calculates spherical bessel functions j_l(qr)
!     NOTA BENE: it is assumed that r(1)=0 always
!
      implicit none
      integer l, mmax
      real(8)    xg, jl(mmax), r(mmax)
! local variables
      real(8)    eps, xrg, xrg2
      parameter(eps=1.e-8)
      integer i, ir
!
!    l=-1 (for derivative  calculations)
!
      if(l.eq.0) then
         if(xg.lt.eps) then
            do i=1,mmax
               jl(i)=0.0
            end do
         else
            jl(1)=0.
            do ir=2,mmax
               xrg=r(ir)*xg
               jl(ir)=cos(xrg)/xrg
            end do
         end if
      end if
!
!    s part
!
      if(l.eq.1) then
         if(xg.lt.eps) then
            do i=1,mmax
               jl(i)=1.0
            end do
         else
            jl(1)=1.
            do ir=2,mmax
               xrg=r(ir)*xg
               jl(ir)=sin(xrg)/xrg
            end do
         endif
      endif
!
!     p-part
!
      if(l.eq.2) then
         if(xg.lt.eps) then
            do i=1,mmax
               jl(i)=0.0
            end do
         else
            jl(1)=0.
            do ir=2,mmax
               xrg=r(ir)*xg
               jl(ir)=(sin(xrg)/xrg-cos(xrg))/xrg
            end do
         endif
      endif
!
!     d part
!
      if(l.eq.3) then
         if(xg.lt.eps) then
            do i=1,mmax
               jl(i)=0.0
            end do
         else
            jl(1)=0.
            do ir=2,mmax
               xrg=r(ir)*xg
               jl(ir)=(sin(xrg)*(3./(xrg*xrg)-1.)                       &
     &              -3.*cos(xrg)/xrg) /xrg
            end do
         endif
      endif
!
!     f part
!
      if(l.eq.4) then
         if(xg.lt.eps) then
            do i=1,mmax
               jl(i)=0.0
            end do
         else
            jl(1)=0.
            do ir=2,mmax
               xrg=r(ir)*xg
               xrg2=xrg*xrg
               jl(ir)=( sin(xrg)*(15./(xrg2*xrg)-6./xrg)                &
     &              +cos(xrg)*(1.-15./xrg2)           )/xrg
            end do
         endif
      endif
!
!     g part
!
      if(l.eq.5) then
         if(xg.lt.eps) then
            do i=1,mmax
               jl(i)=0.0
            end do
         else
            jl(1)=0.
            do ir=2,mmax
               xrg=r(ir)*xg
               xrg2=xrg*xrg
               jl(ir)=( sin(xrg)*(105./(xrg2*xrg2)-45./xrg2+1.)         &
     &              +cos(xrg)*(10./xrg-105./(xrg2*xrg)) )/xrg
            end do
         endif
      endif
!
      return
    end subroutine bess
