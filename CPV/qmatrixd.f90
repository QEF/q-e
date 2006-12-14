!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

subroutine qmatrixd(c0, bec0,ctable, gqq, qmat, detq)

  ! this subroutine computes the inverse of the matrix Q
  ! Q_ij=<Psi_i^0|e^iG_ipol.r|Psi_j^0>
  ! and det Q
  ! Matrix Q is symmetric, and we make us of it

  !    c0 input: the unperturbed wavefunctions
  !    bec0 input: the coefficients <Phi_Rj|Psi_i^0>
  !    ctable input: the coorespondence array
  !    gqq input: the intqq(r) exp(iG_ipol*r) array
  !    qmat output: the inverse q matrix
  !    detq output: det Q

  use gvecs
  use gvecw, only: ngw
  use parameters
  use constants
  use cvan
  use ions_base, only : nas => nax, nsp, na
  use cell_base, only: a1, a2, a3
  use reciprocal_vectors, only: ng0 => gstart
  use uspp_param, only: nh, nhm
  use uspp, only : nhsa=> nkb
  use electrons_base, only: nx => nbspx, n => nbsp, ispin
  use mp, only: mp_sum
  use mp_global, only: intra_image_comm
  
  implicit none
  
  real(8) :: bec0(nhsa,n)
  complex(8) ::  gqq(nhm,nhm,nas,nsp)
  complex(8) :: c0(ngw,nx),  qmat(nx,nx), detq
  integer :: ctable(ngw,2)
  ! local variables
  integer ig,ix,jx, iv,jv,is,ia, inl,jnl
  complex(8) :: sca
  real(8) :: ar(nx,nx),ai(nx,nx),wr(nx),wi(nx),zr(nx,nx), &
       zi(nx,nx),fv1(nx),fv2(nx),fv3(nx)
  integer :: ipiv(nx,nx),info
  complex(8) :: work(nx)

  qmat(:,:)=(0.d0,0.d0)

  do ix=1,n
     do jx=ix,n

! first the local part

        sca=(0.d0,0.d0)
        if(ispin(ix) == ispin(jx) ) then
       
           do ig=1,ngw
              if(ctable(ig,1).ne.(ngw+1))then
                 if(ctable(ig,1).ge.0) then
                    sca=sca+CONJG(c0(ctable(ig,1),ix))*c0(ig,jx)
                 endif
              endif
           enddo

          do ig=1,ngw
             if(ctable(ig,1).ne.(ngw+1))then
                if(ctable(ig,1).lt. 0) then
                   sca=sca+c0(-ctable(ig,1),ix)*c0(ig,jx)
                endif
             endif
          enddo

          do ig=ng0,ngw
             if(ctable(ig,2).ne.(ngw+1)) then
                if(ctable(ig,2).lt.0) then
                   sca=sca+c0(-ctable(ig,2),ix)*CONJG(c0(ig,jx))
                endif
             endif
          enddo
          
          do ig=ng0,ngw
             if(ctable(ig,2).ne.(ngw+1)) then
                if(ctable(ig,2).ge.0) then
                   sca=sca+CONJG(c0(ctable(ig,2),ix))*conjg(c0(ig,jx))
                endif
             endif
          enddo
        
          call mp_sum( sca, intra_image_comm )
       endif
       qmat(ix,jx)=sca
       
       !  now the non local vanderbilt part
       
       sca =(0.d0,0.d0)
       if(ispin(ix)==ispin(jx)) then
          do is=1,nvb!loop on vanderbilt species
             do ia=1,na(is)!loop on atoms
                do iv=1,nh(is)!loop on projectors
                   do jv=1,nh(is)
                      inl=ish(is)+(iv-1)*na(is)+ia
                      jnl=ish(is)+(jv-1)*na(is)+ia                
                      sca=sca+gqq(iv,jv,ia,is)*bec0(inl,ix)*bec0(jnl,jx)
                   enddo
                enddo
             enddo
          enddo
          qmat(ix,jx)=qmat(ix,jx)+sca
       endif
       qmat(jx,ix)=qmat(ix,jx)       
     enddo
  enddo

  !LAPACK
  call zgetrf (n,n,qmat,nx,ipiv,info)
  !  write(6,*) 'info trf', info
  detq=(1.d0,0.d0)
 
  do ix=1,n
     if(ix.ne.ipiv(ix,1)) detq=-detq
  enddo
  do ix=1,n
     detq = detq*qmat(ix,ix)
  enddo

  call zgetri (n,qmat,nx,ipiv,work,nx,info)

! force qmat to be symmetric

  do ix=1,n
     do jx=ix+1,n
        qmat(jx,ix)=0.5d0*(qmat(ix,jx)+qmat(jx,ix))
        qmat(ix,jx)=qmat(jx,ix)
     enddo
  enddo

  return
end subroutine qmatrixd

            
   

