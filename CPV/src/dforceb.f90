!
! Copyright (C) 2002-2008 Quantum ESPRESS0 groups
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


  use kinds, only : DP
  use  gvecs
  use gvecw, only: ngw
  use  parameters
  use electrons_base, only: nx => nbspx, n => nbsp, nspin, f
  use  constants
  use ions_base, only : nat, nax, na, nsp
  use cell_base, only: at, alat
  use uspp_param, only: nh, nhm, nvb, ish
  use uspp, only : nhsa=> nkb
  use efield_module, ONLY : ctabin_missing_1,ctabin_missing_2,n_g_missing_m,&
       &      ctabin_missing_rev_1,ctabin_missing_rev_2
  use mp_global, only: intra_bgrp_comm, nproc_bgrp
  use mp, only: mp_alltoall
  use parallel_include


  implicit none
      
      
  complex(DP) c0(ngw, n), betae(ngw,nhsa), df(ngw),&
       &   gqq(nhm,nhm,nax,nsp),gqqm(nhm,nhm,nax,nsp),&
       &   qmat(nx,nx)
  real(DP) bec0(nhsa,n), dq2(nat,nhm,nhm,nspin),  gmes
  real(DP), EXTERNAL :: g_mes

  integer i, ipol, ctabin(ngw,2)

! local variables

  integer j,k,ig,iv,jv,ix,jx,is,ia, isa,iss,iss1,mism
  integer ir,ism,itemp,itempa,jnl,inl
  complex(DP) ci ,fi, fp, fm
  real(DP) afr(nhsa), dd
  complex(DP)  afrc(nhsa)
  complex(DP), allocatable::  dtemp(:)
  complex(DP), allocatable :: sndbuf(:,:,:),rcvbuf(:,:,:)
  integer :: ierr, ip

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

#if defined(__MPI)

     if(ipol/=3) then
!allocate arrays
             allocate(sndbuf(n_g_missing_m(ipol),2,nproc_bgrp))
             sndbuf(:,:,:)=(0.d0,0.d0)
             allocate(rcvbuf(n_g_missing_m(ipol),2,nproc_bgrp))
!copy arrays to snd buf
             do ip=1,nproc_bgrp
                do ig=1,n_g_missing_m(ipol)
                   if(ipol==1) then
                      if(ctabin_missing_rev_1(ig,1,ip)>0) then
                         sndbuf(ig,1,ip)=-c0(ctabin_missing_rev_1(ig,1,ip),j)*CONJG(qmat(j,i))
                      elseif(ctabin_missing_rev_1(ig,1,ip)<0) then
                         sndbuf(ig,1,ip)=-conjg(c0(-ctabin_missing_rev_1(ig,1,ip),j))*CONJG(qmat(j,i))
                      endif
                   else
                     if(ctabin_missing_rev_2(ig,1,ip)>0) then
                         sndbuf(ig,1,ip)=-c0(ctabin_missing_rev_2(ig,1,ip),j)*CONJG(qmat(j,i))
                      elseif(ctabin_missing_rev_2(ig,1,ip)<0) then
                         sndbuf(ig,1,ip)=-conjg(c0(-ctabin_missing_rev_2(ig,1,ip),j))*CONJG(qmat(j,i))
                      endif
                   endif
                enddo
                do ig=1,n_g_missing_m(ipol)
                   if(ipol==1) then
                      if(ctabin_missing_rev_1(ig,2,ip)>0) then
                         sndbuf(ig,2,ip)=c0(ctabin_missing_rev_1(ig,2,ip),j)*qmat(j,i)
                      elseif(ctabin_missing_rev_1(ig,2,ip)<0) then
                         sndbuf(ig,2,ip)=conjg(c0(-ctabin_missing_rev_1(ig,2,ip),j))*qmat(j,i)
                      endif
                   else
                      if(ctabin_missing_rev_2(ig,2,ip)>0) then
                         sndbuf(ig,2,ip)=c0(ctabin_missing_rev_2(ig,2,ip),j)*qmat(j,i)
                      elseif(ctabin_missing_rev_2(ig,2,ip)<0) then
                         sndbuf(ig,2,ip)=conjg(c0(-ctabin_missing_rev_2(ig,2,ip),j))*qmat(j,i)
                      endif
                   endif
                enddo
             enddo


             CALL mp_alltoall( sndbuf, rcvbuf, intra_bgrp_comm )           

!update sca
             do ip=1,nproc_bgrp
                do ig=1,n_g_missing_m(ipol)
                   if(ipol==1) then
                      if(ctabin_missing_1(ig,1,ip)/=0) then
                         dtemp(ctabin_missing_1(ig,1,ip))=dtemp(ctabin_missing_1(ig,1,ip))+rcvbuf(ig,1,ip)
                      endif
                       if(ctabin_missing_1(ig,2,ip)/=0) then
                         dtemp(ctabin_missing_1(ig,2,ip))=dtemp(ctabin_missing_1(ig,2,ip))+rcvbuf(ig,2,ip)
                      endif
                   else
                      if(ctabin_missing_2(ig,1,ip)/=0) then
                         dtemp(ctabin_missing_2(ig,1,ip))=dtemp(ctabin_missing_2(ig,1,ip))+rcvbuf(ig,1,ip)
                      endif
                       if(ctabin_missing_2(ig,2,ip)/=0) then
                         dtemp(ctabin_missing_2(ig,2,ip))=dtemp(ctabin_missing_2(ig,2,ip))+rcvbuf(ig,2,ip)
                      endif
                   endif
                enddo
             enddo
               
            
          


             

!deallocate arrays
             deallocate(rcvbuf,sndbuf)
          endif
      
#endif
    enddo
  
  gmes = g_mes ( ipol, at, alat ) 

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



function enberry( detq,  ipol )

   use constants
   use kinds, only: dp
   use cell_base, only: alat, at
   USE electrons_base, ONLY : nspin

   implicit none

   complex(dp), intent (in) :: detq
   real(dp) :: enberry
 
   integer ipol
   real(dp) gmes
   real(dp), external :: g_mes

   gmes = g_mes ( ipol, at, alat )
   enberry = 2.d0/REAL(nspin,DP)*AIMAG(log(detq))/gmes ! take care of sign
   
   return
 end function enberry


!
! Copyright (C) 2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

FUNCTION g_mes ( ipol, at, alat )

  USE kinds, ONLY : dp
  USE constants, ONLY : pi

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ipol
  REAL(dp), INTENT(IN) :: at(3,3), alat
  REAL(dp) :: g_mes

  IF ( ipol < 1 .OR. ipol > 3) CALL errore ( 'gmes','incorrect ipol', 1)
  g_mes = 2.0_dp*pi/alat/SQRT(at(1,ipol)**2+at(2,ipol)**2+at(3,ipol)**2)
  
END FUNCTION g_mes
