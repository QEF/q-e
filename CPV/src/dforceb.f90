!
! Copyright (C) 2002-2008 Quantum ESPRESS0 groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine dforceb(c0, i, betae, ipol, bec0, ctabin, gqq, gqqm, qmat, dq2, df)

  !! This subroutine computes the force for electrons in case of Berry phase
  !! like perturbation.
  !! It gives the force for the i-th state.

! c1 input: Psi^1_i
! v0      input: the local zeroth order potential
! v1      input: the local first order potential
! bec1   input: the factors <beta_iR|Psi^1_v>

  use kinds, only : DP
  use gvecw, only: ngw
  use  parameters
  use electrons_base, only: nx => nbspx, n => nbsp, nspin, f
  use constants, only :
  use ions_base, only : nat, nax, nsp, ityp
  use cell_base, only: at, alat
  use uspp_param, only: nh, nhm, upf
  use uspp, only : nkb, ofsbeta
  use efield_module, ONLY : ctabin_missing_1,ctabin_missing_2,n_g_missing_m,&
       &      ctabin_missing_rev_1,ctabin_missing_rev_2
  use mp_global, only: intra_bgrp_comm, nproc_bgrp
  use mp, only: mp_alltoall
  use parallel_include

  implicit none
      
  complex(DP) :: c0(ngw, n)
  !! input: Psi^0_i
  complex(DP) :: betae(ngw,nkb)
  !! input: the functions beta_iR
  complex(DP) :: df(ngw)
  !! output: force for the i-th state
  complex(DP) :: gqq(nhm,nhm,nax,nsp)
  !! input: the factors int dr Beta_Rj*Beta_Ri exp(iGr)
  complex(DP) :: gqqm(nhm,nhm,nax,nsp)
  !! input: the factors int dr Beta_Rj*Beta_Ri exp(iGr)
  complex(DP) :: qmat(nx,nx)
  !! input
  real(DP) :: bec0(nkb,n)
  !! input: the factors <beta_iR|Psi^0_v>
  real(DP) :: dq2(nat,nhm,nhm,nspin)
  !! input: factors d^2hxc_ijR
  integer :: i
  !! input: ot computes the force for the i-th state
  integer :: ctabin(ngw,2)
  !! input: the inverse-correspondence array g'+(-)1=g
  integer :: ipol
  !! input:the polarization of nabla_k

  real(DP), EXTERNAL :: g_mes
  
  ! ... local variables

  integer :: j,k,ig,iv,jv,ix,jx,is,ia, iss,iss1,mism
  integer :: ir,ism,itemp,itempa,jnl,inl
  complex(DP) :: ci ,fi, fp, fm
  real(DP) :: afr(nkb), dd, gmes
  complex(DP) :: afrc(nkb)
  complex(DP), allocatable :: dtemp(:)
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



  if(nkb > 0) then
     do inl=1,nkb
        afrc(inl)=(0.d0,0.d0)
     end do
 
     do ia=1,nat
        is=ityp(ia)
        IF(upf(is)%tvanp) THEN
           do iv=1,nh(is)      !loop on projectors           
              do jv=1,nh(is)   !loop on projectors                               
                  inl = ofsbeta(ia) + iv
                  jnl = ofsbeta(ia) + jv
                  do j=1,n  !loop on states
                     afrc(inl)=afrc(inl)+gqq(iv,jv,ia,is)*bec0(jnl,j)*qmat(j,i)&
                          &     -CONJG(gqq(jv,iv,ia,is))*bec0(jnl,j)*conjg(qmat(i,j))
  

                  end do
               end do
            end do
         ENDIF
      end do

      do ig=1,ngw
         dtemp(ig)=(0.d0,0.d0)
      end do
      do inl=1,nkb
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

   use constants, only :
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
