!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

subroutine bforceion(fion,tfor,ipol,qmatinv,bec0,becdr,gqq,evalue)

! this subroutine compute the part of force for the ions due to
! electronic berry phase( see internal notes)
! it needs becdr

! fion       : input, forces on ions
! tfor       : input, if true it computes force
! a1,a2,a3   : input, direct lattice vectors
! ipol       : input, electric field polarization
! qmatinv    : input, inverse of Q matrix: Q_i,j=<Psi_i|exp(iG*r)|Psi_j>
! bec0       : input, factors <beta_iR|Psi_j>
! becdr      : input, factors d<beta_iR|Psi_j>/dR
! gqq        : input, Int_e exp(iG*r)*q_ijR(r)
! evalue     : input, scale of electric field
 
  use ions_base, only : nas => nax, na, nsp
  use cvan
  use parameters
  use constants
  use cell_base, only: a1, a2, a3
  use uspp_param, only: nh, nhm
  use uspp, only : nhsa=> nkb
  use electrons_base, only: n => nbsp, nx => nbspx


  implicit none

  real(8) evalue
  complex(8) qmatinv(nx,nx),gqq(nhm,nhm,nas,nsp)
  real(8) bec0(nhsa,n),becdr(nhsa,n,3)
  real(8) fion(3,*)
  integer ipol
  logical tfor

!local variables

  complex(8) ci, temp, temp1,temp2,temp3
  real(8) gmes
  integer iv,jv,ia,is,k,i,j,isa,ilm,jlm,inl,jnl,ism
      
  if(.not. tfor) return

  ci = (0.d0,1.d0)
     

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



  isa = 0
  do is=1,nvb
   do ia=1,na(is)
     isa = isa + 1
     do iv= 1,nh(is)
        do jv=1,nh(is)         
              inl=ish(is)+(iv-1)*na(is)+ia
              jnl=ish(is)+(jv-1)*na(is)+ia
             
              temp=(0.d0,0.d0)
              temp1=(0.d0,0.d0)
              temp2=(0.d0,0.d0)
              temp3=(0.d0,0.d0)
              do i=1,n
                 do j=1,n

                    temp = temp + ci*gmes*gqq(iv,jv,ia,is)* &!TAKECARE: sign + due to exp(+iGr) in gqq
                         &        bec0(inl,i)*bec0(jnl,j)*qmatinv(j,i)

                    temp1 = temp1 + gqq(iv,jv,ia,is)*&
     &  (  becdr(inl,i,1)*bec0(jnl,j)+bec0(inl,i)*becdr(jnl,j,1))*qmatinv(j,i)

                    temp2 = temp2 + gqq(iv,jv,ia,is)*&
     &  (  becdr(inl,i,2)*bec0(jnl,j)+bec0(inl,i)*becdr(jnl,j,2))*qmatinv(j,i)

                    temp3 = temp3 + gqq(iv,jv,ia,is)*&
     &  (  becdr(inl,i,3)*bec0(jnl,j)+bec0(inl,i)*becdr(jnl,j,3))*qmatinv(j,i)


                 enddo
              enddo

              fion(ipol,isa) = fion(ipol,isa) -   2.d0*evalue*AIMAG(temp)/gmes
              fion(1,isa) = fion(1,isa) -   2.d0*evalue*AIMAG(temp1)/gmes
              fion(2,isa) = fion(2,isa) -   2.d0*evalue*AIMAG(temp2)/gmes
              fion(3,isa) = fion(3,isa) -   2.d0*evalue*AIMAG(temp3)/gmes
           end do
        end do
     end do
  end do

  return
end subroutine bforceion
