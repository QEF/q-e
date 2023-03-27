!
! Copyright (C) 2002-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

subroutine bforceion(fion,tfor,ipol,qmatinv,bec0,becdr,gqq,evalue)

  !! This subroutine computes the part of force for the ions due to
  !! electronic Berry phase (see internal notes).  
  !! It needs \(\text{becdr}\) as input variable.
 
  use ions_base, only : nax, na, nsp, nat, ityp
  use uspp_param, only: upf, nh, nhm
  use uspp, only : nkb, ofsbeta
  use kinds, only : dp
  use constants, only :
  use cell_base, only: at, alat ! at: direct lattice vectors, divided by alat
  use electrons_base, only: nbsp, nbspx, nspin, nbspx_bgrp
  use mp_global, only: nbgrp

  implicit none

  real(dp) :: evalue
  !! input, scale of electric field
  complex(dp) :: qmatinv(nbspx,nbspx)
  !! input, inverse of Q matrix:
  !! \( Q_{i,j}=\langle\Psi_i|\exp(iG\ r)|\Psi_j\rangle\)
  complex(dp) :: gqq(nhm,nhm,nax,nsp)
  !! input, \( \text{Int}_e \exp(iG\ r)\ q_ijR(r) \)
  real(dp) :: bec0(nkb,nbspx)
  !! input, factors \(\langle \beta_{iR}|\Psi_j\rangle\)
  real(dp) :: becdr(nkb,nbspx,3)
  !! input, factors \(d\langle\beta_{iR}|\Psi_j\rangle/dR \)
  real(dp) :: fion(3,*)
  !! input, forces on ions
  integer :: ipol
  !! input, electric field polarization
  logical :: tfor
  !! input, if TRUE it computes force
  
  ! ... local variables

  complex(dp) :: ci, temp, temp1,temp2,temp3
  real(dp) :: gmes
  real(dp), external :: g_mes
  integer :: iv,jv,ia,is,k,i,j,ilm,jlm,inl,jnl,ism
      
  if(.not. tfor) return

  if( nbgrp > 1 ) &
     call errore(' bforceion ', ' parallelization over bands not yet implemented ', 1 )

  ci = (0.d0,1.d0)
  gmes = g_mes (ipol, at, alat) 

  do ia=1,nat
     is=ityp(ia)
     IF(upf(is)%tvanp) THEN
        do iv= 1,nh(is)
           do jv=1,nh(is)         
              inl = ofsbeta(ia) + iv
              jnl = ofsbeta(ia) + jv
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

              fion(ipol,ia) = fion(ipol,ia) -   2.d0*evalue*AIMAG(temp)/gmes
              fion(1,ia) = fion(1,ia) -   2.d0*evalue*AIMAG(temp1)/gmes
              fion(2,ia) = fion(2,ia) -   2.d0*evalue*AIMAG(temp2)/gmes
              fion(3,ia) = fion(3,ia) -   2.d0*evalue*AIMAG(temp3)/gmes
           end do
        end do
     END IF
  end do

  return
end subroutine bforceion
