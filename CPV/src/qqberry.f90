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

  use kinds,              only: dp
  use uspp_param,         only: upf, lmaxq, nbetam, nh, nhm, oldvan, nvb
  use uspp,               only: indv, lpx, lpl, ap,nhtolm
  use atom,               only: rgrid
  use core
  use gvecw,              only: ngw
  use gvect,              only: mill
  use constants
  use ions_base,          only: nax, na, nsp
  use cell_base,          only: at, alat
  use gvect,              only: g, gg
  use mp,                 only: mp_sum
  use mp_global,          only: intra_bgrp_comm
  use cp_interfaces,      only: fill_qrl
  
  implicit none

  complex(dp) gqq(nhm,nhm,nax,nsp)
  complex(dp) gqqm(nhm,nhm,nax,nsp)
  real(dp) gmes
  real(dp), external :: g_mes
  integer :: ipol

! local variables

  integer :: ndm, ig, is, iv, jv, i, istart, il,l,ir, igi,ia
  real(dp), allocatable:: fint(:),jl(:)
  real(dp), allocatable:: qrl(:,:,:), qradb2(:,:,:,:) 
  real(dp) c, xg
  complex(dp) qgbs,sig
  integer :: ivs, jvs, ivl, jvl, lp, ijv
  real(dp), allocatable:: ylm(:,:)

  IF( .NOT. ALLOCATED( rgrid ) ) &
     CALL errore( ' qqberry2 ', ' rgrid not allocated ', 1 )
  IF( .NOT. ALLOCATED( upf ) ) &
     CALL errore( ' qqberry2 ', ' upf not allocated ', 1 )

  ndm = MAXVAL (upf(1:nsp)%kkbeta)
  allocate( fint( ndm), jl(ndm))
  allocate( qradb2(nbetam,nbetam,lmaxq,nsp))
  allocate( ylm(ngw, lmaxq*lmaxq))

  CALL ylmr2( lmaxq*lmaxq, ngw, g, gg, ylm )

  qradb2 = 0.0d0
     
  do is=1,nsp
     do ia=1,nax
        do jv=1,nhm
           do iv=1,nhm
              gqq(iv,jv,ia,is)=(0.d0,0.d0)
              gqqm(iv,jv,ia,is)=(0.d0,0.d0)
           enddo
        enddo
     enddo
  enddo
  
  gmes = g_mes ( ipol, at, alat )

  ! only for Vanderbilt species 
  do is=1,nvb
     c=fpi                 !/omegab
     !
     ALLOCATE ( qrl( upf(is)%kkbeta, upf(is)%nbeta*(upf(is)%nbeta+1)/2, &
                     upf(is)%nqlc ) )
     !
     call fill_qrl ( is, qrl )
     ! now the radial part
     do l=1,upf(is)%nqlc
        xg= gmes !only orthorombic cells
        !!!call bess(xg,l,upf(is)%kkbeta,rgrid(is)%r,jl)
        call sph_bes ( upf(is)%kkbeta, rgrid(is)%r, xg, l-1, jl )
        do iv= 1,upf(is)%nbeta
           do jv=iv,upf(is)%nbeta
              ijv = (jv-1)*jv/2 + iv
!     
!     note qrl(r)=r^2*q(r)
!
              do ir=1,upf(is)%kkbeta
                 fint(ir)=qrl(ir,ijv,l)*jl(ir)
              end do
              if (oldvan(is)) then
                 call herman_skillman_int ( upf(is)%kkbeta,fint,rgrid(is)%rab,&
                                            qradb2(iv,jv,l,is) )
              else
                 call simpson ( upf(is)%kkbeta,fint,rgrid(is)%rab,&
                                qradb2(iv,jv,l,is) )
              endif
              qradb2(iv,jv,l,is)=  c*qradb2(iv,jv,l,is)
              if ( iv /= jv ) qradb2(jv,iv,l,is)=  qradb2(iv,jv,l,is)
           end do
        end do
     end do
     DEALLOCATE ( qrl )    
  enddo

  igi=-1
  do ig=1,ngw
     if(ipol.eq.1 ) then
        if(mill(1,ig).eq.1 .and. mill(2,ig).eq.0  .and. mill(3,ig).eq. 0) igi=ig
     endif
     if(ipol.eq.2 ) then
        if(mill(1,ig).eq.0 .and. mill(2,ig).eq.1  .and. mill(3,ig).eq. 0) igi=ig
     endif
     if(ipol.eq.3 ) then
        if(mill(1,ig).eq.0 .and. mill(2,ig).eq.0   .and. mill(3,ig).eq. 1) igi=ig
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
              qgbs=(0.d0,0.d0)
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
                 sig=(0.d0,-1.d0)**(l-1)
                  
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

  call mp_sum(gqq(:,:,:,:),intra_bgrp_comm)
  call mp_sum(gqqm(:,:,:,:),intra_bgrp_comm)

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

  use kinds, only : dp
  use gvecw, only: ngw
  use ions_base, only : nax, nat, na, nsp
  use gvect, only: mill
  use uspp_param, only: nh, nhm, nvb, ish
  use mp, only: mp_sum
  use mp_global, only: intra_bgrp_comm

  implicit none

 
  complex(dp) eigr(ngw,nat)
  complex(dp) gqq(nhm,nhm,nax,nsp)
  complex(dp) gqqm(nhm,nhm,nax,nsp)
  complex(dp) gqqm0(nhm,nhm,nax,nsp)

  integer ipol
  
  integer igi,ig,is,iv,jv,ia,isa


  do is=1,nsp
     do ia=1,nax
        do jv=1,nhm
           do iv=1,nhm
              gqq(iv,jv,ia,is)=(0.d0,0.d0)
              gqqm(iv,jv,ia,is)=(0.d0,0.d0)
           enddo
        enddo
     enddo
  enddo

  igi=-1
  do ig=1,ngw
     if(ipol.eq.1 ) then
        if(mill(1,ig).eq.1 .and. mill(2,ig).eq.0  .and. mill(3,ig).eq. 0) igi=ig
     endif
     if(ipol.eq.2 ) then
        if(mill(1,ig).eq.0 .and. mill(2,ig).eq.1  .and. mill(3,ig).eq. 0) igi=ig
     endif
     if(ipol.eq.3 ) then
        if(mill(1,ig).eq.0 .and. mill(2,ig).eq.0  .and. mill(3,ig).eq. 1) igi=ig
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
  call mp_sum(gqq(:,:,:,:),intra_bgrp_comm)
  call mp_sum(gqqm(:,:,:,:),intra_bgrp_comm)
  return
end subroutine qqupdate
