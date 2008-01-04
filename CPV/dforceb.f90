!
! Copyright (C) 2002-2008 Quantum-ESPRESS0 groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine dforceb(c0, i, betae, ipol, bec0, ctabin, gqq, gqqm, qmat, dq2, df)

! this subroutine computes the force for electrons
! in case of Berry,s phase like perturbation
! it gives the force for the i-th state

! c0 input: Psi^0_i
! c1 input: Psi^1_i
! i  input: ot computes the force for the i-th state
! v0      input: the local zeroth order potential
! v1      input: the local first order potential
! betae  input: the functions beta_iR
! ipol   input:the polarization of nabla_k
! bec0   input: the factors <beta_iR|Psi^0_v>
! bec1   input: the factors <beta_iR|Psi^1_v>
! ctabin input: the inverse-correspondence array g'+(-)1=g
! gqq    input: the factors int dr Beta_Rj*Beta_Ri exp(iGr)
! gqqm   input: the factors int dr Beta_Rj*Beta_Ri exp(iGr)
! qmat   input: 
!   dq2  input: factors d^2hxc_ijR
!   df     output: force for the i-th state


  use  gvecs
  use gvecw, only: ngw
  use  parameters
  use electrons_base, only: nx => nbspx, n => nbsp, nspin, f
  use  constants
  use  cvan
  use ions_base, only : nat, nas => nax, na, nsp
  use cell_base, only: a1, a2, a3
  use uspp_param, only: nh, nhm
  use uspp, only : nhsa=> nkb



  implicit none
      
      
  complex(8) c0(ngw, n), betae(ngw,nhsa), df(ngw),&
       &   gqq(nhm,nhm,nas,nsp),gqqm(nhm,nhm,nas,nsp),&
       &   qmat(nx,nx)
  real(8) bec0(nhsa,n),&
       &   dq2(nat,nhm,nhm,nspin),  gmes

  integer i, ipol, ctabin(ngw,2)

! local variables

  integer j,k,ig,iv,jv,ix,jx,is,ia, isa,iss,iss1,mism
  integer ir,ism,itemp,itempa,jnl,inl
  complex(8) ci ,fi, fp, fm
  real(8) afr(nhsa), dd
  complex(8)  afrc(nhsa)
  complex(8), allocatable::  dtemp(:)

  allocate( dtemp(ngw))


  ci=(0.d0,1.d0)
 

! now the interaction term 
! first the norm-conserving part 

  do ig=1,ngw
     dtemp(ig)=(0.d0,0.d0)
  enddo

  do j=1,n
     do ig=1,ngw
        if(ctabin(ig,2) .ne. (ngw+1)) then
           if(ctabin(ig,2).ge.0) then
              dtemp(ig)=dtemp(ig)+c0(ctabin(ig,2),j)*qmat(j,i)
           else
              dtemp(ig)=dtemp(ig)+CONJG(c0(-ctabin(ig,2),j))*qmat(j,i)
           endif
        endif
     enddo
     do ig=1,ngw
        if(ctabin(ig,1) .ne. (ngw+1)) then
           if(ctabin(ig,1).ge.0) then
              dtemp(ig)=dtemp(ig)-c0(ctabin(ig,1),j)*CONJG(qmat(j,i))
           else
              dtemp(ig)=dtemp(ig)-CONJG(c0(-ctabin(ig,1),j))*conjg(qmat(j,i))
           endif
        endif
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

  fi=f(i)*ci/(2.d0*gmes)

  do ig=1,ngw
     df(ig)= fi*dtemp(ig)
  end do
      
! now the interacting Vanderbilt term
! the term (-ie/|G|)(-beta_i'R>gqq(i',j')bec0_jRj'Q^-1_ji+
! +beta_i'R>gqqm(i',j')bec0jRj'Q^-1_ij*



  if(nhsa.gt.0) then
     do inl=1,nhsa
        afrc(inl)=(0.d0,0.d0)
     end do
 
     do is=1,nvb!loop on species
        do iv=1,nh(is)      !loop on projectors           
           do jv=1,nh(is)   !loop on projectors                               
               do ia=1,na(is)
                  inl=ish(is)+(iv-1)*na(is)+ia
                  jnl=ish(is)+(jv-1)*na(is)+ia              
                  do j=1,n  !loop on states
                     afrc(inl)=afrc(inl)+gqq(iv,jv,ia,is)*bec0(jnl,j)*qmat(j,i)&
                          &     -CONJG(gqq(jv,iv,ia,is))*bec0(jnl,j)*conjg(qmat(i,j))
  

                  end do
               end do
            end do
         end do
      enddo

      do ig=1,ngw
         dtemp(ig)=(0.d0,0.d0)
      end do
      do inl=1,nhsa
         do ig=1,ngw
            dtemp(ig)=dtemp(ig)+afrc(inl)*betae(ig,inl)
         enddo
      enddo
!         call MXMA
!     &        (betae,1,2*ngw,afr,1,nhsax,dtemp,1,2*ngw,2*ngw,nhsa,1)
      do ig=1,ngw
         df(ig)=df(ig)+fi*dtemp(ig)
      end do
   endif
      
   deallocate( dtemp)
   return
 end subroutine dforceb



 subroutine enberry( detq,  ipol, enb)

   use constants
   use parameters
   use cell_base, only: a1, a2, a3

   implicit none

   complex(8) detq
   real(8) enb
   integer ipol
   real(8) gmes

     
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

   
   enb = 2.d0*AIMAG(log(detq))/gmes!take care of sign
   
   return
 end subroutine enberry



