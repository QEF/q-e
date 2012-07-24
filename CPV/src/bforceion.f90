!
! Copyright (C) 2002-2005 Quantum ESPRESSO group
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
! at         : input, direct lattice vectors, divided by alat
! ipol       : input, electric field polarization
! qmatinv    : input, inverse of Q matrix: Q_i,j=<Psi_i|exp(iG*r)|Psi_j>
! bec0       : input, factors <beta_iR|Psi_j>
! becdr      : input, factors d<beta_iR|Psi_j>/dR
! gqq        : input, Int_e exp(iG*r)*q_ijR(r)
! evalue     : input, scale of electric field
 
  use ions_base, only : nax, na, nsp
  use uspp_param, only: nvb, ish
  use kinds, only : dp
  use constants
  use cell_base, only: at, alat
  use uspp_param, only: nh, nhm
  use uspp, only : nhsa=> nkb
  use electrons_base, only: nbsp, nbspx, nspin, nbspx_bgrp
  use mp_global, only: nbgrp


  implicit none

  real(dp) evalue
  complex(dp) qmatinv(nbspx,nbspx),gqq(nhm,nhm,nax,nsp)
  real(dp) bec0(nhsa,nbspx),becdr(nhsa,nbspx,3)
  real(dp) fion(3,*)
  integer ipol
  logical tfor

!local variables

  complex(dp) ci, temp, temp1,temp2,temp3
  real(dp) :: gmes
  real(dp), external :: g_mes
  integer iv,jv,ia,is,k,i,j,isa,ilm,jlm,inl,jnl,ism
      
  if(.not. tfor) return

  if( nbgrp > 1 ) &
     call errore(' bforceion ', ' parallelization over bands not yet implemented ', 1 )

  ci = (0.d0,1.d0)
  gmes = g_mes (ipol, at, alat) 

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
              do i=1,nbsp
                 do j=1,nbsp

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
