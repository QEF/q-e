!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!-------------------------------------------------------------------------
      subroutine calphiid(c0,bec,betae,phi)
!-----------------------------------------------------------------------
!     input: c0 (orthonormal with s(r(t)), bec=<c0|beta>, betae=|beta>
!     computes the matrix phi (with the old positions)
!       where  |phi> = s'|c0> = |c0> + sum q_ij |i><j|c0>
!     where s'=s(r(t))  
!     
      use ions_base, only: na, nsp
      use io_global, only: stdout
      use cvan
      use uspp_param, only: nh
      use uspp, only :nhsa=>nkb, nhsavb=>nkbus, qq
      use electrons_base, only: n => nbsp
      use gvecw, only: ngw
      use constants, only: pi, fpi
      use control_flags, only: iprint, iprsta
      use mp, only: mp_sum
!
      implicit none
      complex(8) c0(ngw,n), phi(ngw,n), betae(ngw,nhsa)
      real(8)    bec(nhsa,n), emtot
! local variables
      integer is, iv, jv, ia, inl, jnl, i, j
      real(8) qtemp(nhsavb,n) ! automatic array
!
      phi(1:ngw,1:n) = 0.0d0
!
      if (nvb.gt.0) then
         qtemp = 0.0d0
         do is=1,nvb
            do iv=1,nh(is)
               do jv=1,nh(is)
                  if(abs(qq(iv,jv,is)).gt.1.e-5) then
                     do ia=1,na(is)
                        inl=ish(is)+(iv-1)*na(is)+ia
                        jnl=ish(is)+(jv-1)*na(is)+ia
                        do i=1,n
                           qtemp(inl,i) = qtemp(inl,i) +                &
     &                                    qq(iv,jv,is)*bec(jnl,i)
                        end do
                     end do
                  endif
               end do
            end do
         end do
         !
         call MXMA                                                     &
     &       (betae,1,2*ngw,qtemp,1,nhsavb,phi,1,2*ngw,2*ngw,nhsavb,n)
      end if
!
      do j=1,n
         do i=1,ngw
            phi(i,j)=(phi(i,j)+c0(i,j))
         end do
      end do
!     =================================================================
      if(iprsta.gt.2) then
         emtot=0.
         do j=1,n
            do i=1,ngw
               emtot=emtot                                              &
     &        +2.*DBLE(phi(i,j)*CONJG(c0(i,j)))
            end do
         end do
         emtot=emtot/n
         call mp_sum(emtot)
         WRITE( stdout,*) 'in calphi sqrt(emtot)=',sqrt(emtot)
         WRITE( stdout,*)
         do is=1,nsp
            if(nsp.gt.1) then
               WRITE( stdout,'(33x,a,i4)') ' calphi: bec (is)',is
               WRITE( stdout,'(8f9.4)')                                       &
     &            ((bec(ish(is)+(iv-1)*na(is)+1,i),iv=1,nh(is)),i=1,n)
            else
               do ia=1,na(is)
                  WRITE( stdout,'(33x,a,i4)') ' calphi: bec (ia)',ia
                  WRITE( stdout,'(8f9.4)')                                    &
     &               ((bec(ish(is)+(iv-1)*na(is)+ia,i),iv=1,nh(is)),i=1,n)
               end do
            end if
         end do
      endif
!
      return
    end subroutine calphiid



!-----------------------------------------------------------------------
      subroutine calcmt(fdiag,zmat,fmat)
!-----------------------------------------------------------------------
!
!  constructs fmat=zmat^t.fdiag.zmat
!
      use electrons_base, only: nudx, nspin, nupdwn, iupdwn, nx => nbspx

      implicit none
      integer iss, nss, istart, i, j, k
      real(8) zmat(nudx,nudx,nspin), fmat(nudx,nudx,nspin),        &
   &                fdiag(nx)

      do iss=1,nspin
       nss=nupdwn(iss)
       istart=iupdwn(iss)
       do i=1,nss
        do k=1,nss
         fmat(k,i,iss)=0.0
         do j=1,nss
          fmat(k,i,iss)=fmat(k,i,iss)+                                   &
    &         zmat(j,k,iss)*fdiag(j+istart-1)*zmat(j,i,iss)
         end do
        end do
       end do
      end do

      return
      end subroutine calcmt


!-----------------------------------------------------------------------
      subroutine rotate(z0,c0,bec,c0diag,becdiag)
!-----------------------------------------------------------------------
      use kinds, only: dp
      use cvan
      use electrons_base, only: nudx, nspin, nupdwn, iupdwn, nx => nbspx, n => nbsp
      use uspp_param, only: nh
      use uspp, only :nhsa=>nkb, nhsavb=>nkbus, qq
      use gvecw, only: ngw
      use ions_base, only: nas => nax, nsp, na

      implicit none
      integer iss, nss, istart, i, j, k, ni, nj, is, jv, ia, jnl
      real(8) z0(nudx,nudx,nspin)
      real(8) bec(nhsa,n), becdiag(nhsa,n)
      complex(8) c0(ngw,nx), c0diag(ngw,nx)

      c0diag(1:ngw,1:nx)=0.d0
        do iss=1,nspin
         nss=nupdwn(iss)
         istart=iupdwn(iss)
          call MXMA(c0(1,istart),1,2*ngw,z0(1,1,iss),nudx,1,c0diag(1,istart),1,2*ngw,2*ngw,nss,nss)
        end do

        do iss=1,nspin
        nss=nupdwn(iss)
        istart=iupdwn(iss)
        do ni=1,nss
         do is=1,nsp
          do jv=1,nh(is)
           do ia=1,na(is)
            jnl=ish(is)+(jv-1)*na(is)+ia
             becdiag(jnl,ni+istart-1)=0.0
             do nj=1,nss
              becdiag(jnl,ni+istart-1)=becdiag(jnl,ni+istart-1)+        &
     &        CMPLX(z0(ni,nj,iss),0.d0)*bec(jnl,nj+istart-1)
             end do
           end do
          end do
         end do
        end do
        end do
      return
      end subroutine rotate

!-----------------------------------------------------------------------
      subroutine ddiag(nx,n,amat,dval,dvec,iflag)
!-----------------------------------------------------------------------
!
      implicit none
      integer nx,n,naux,ndim,iopt,iflag,k,i,j,info
      real(8)   dval(n)
      real(8) amat(nx,n), dvec(nx,n)
      real(8), allocatable::  ap(:)
      real(8), allocatable::  aux(:)


      ndim=(n*(n+1))/2
      naux=3*n
      allocate(ap(ndim))
      allocate(aux(naux))
      ap(:)=0.d0
      aux(:)=0.d0

      if (iflag.eq.1) then
       iopt=1
      else if (iflag.eq.0) then
       iopt=0
      else
       write(*,*) 'ERROR: diag, iflag not allowed'
       stop
      end if

      k=0
      do j=1,n
       do i=1,j
        k=k+1
!       ap(i + (j-1)*j/2)=amat(i,j)
        ap(k)=amat(i,j)
       end do
      end do

      call dspev('V','U',n,ap,dval,dvec,nx,aux,info)
      if(info.ne.0) write(6,*) 'Problems with ddiag'

      deallocate(ap)
      deallocate(aux)

      return
    end subroutine ddiag

!-----------------------------------------------------------------------
      subroutine calcm(fdiag,zmat,fmat)
!-----------------------------------------------------------------------
!
!  constructs fmat=zmat.fdiag.zmat^t
!
      use electrons_base, only: nudx, nspin, nupdwn, iupdwn, nx => nbspx, n => nbsp

      implicit none
      integer iss, nss, istart, i, j, k
      real(8) zmat(nudx,nudx,nspin), fmat(nudx,nudx,nspin),         &
    &   fdiag(nx)

      do iss=1,nspin
       nss=nupdwn(iss)
       istart=iupdwn(iss)
       do i=1,nss
        do k=1,nss
         fmat(k,i,iss)=0.0
         do j=1,nss
          fmat(k,i,iss)=fmat(k,i,iss)+                                  &
    &            zmat(k,j,iss)*fdiag(j+istart-1)*zmat(i,j,iss)
         end do
        end do
       end do
      end do

      return
      end subroutine calcm

    subroutine minparabola(ene0,dene0,ene1,passop,passo,stima)
!this subroutines finds the minimum of a quadratic real function
      
      implicit none
      real(8) ene0,dene0,ene1,passop,passo,stima
      real(8) a,b,c!a*x^2+b*x+c
      
      c=ene0
      b=dene0
      a=(ene1-b*passop-c)/(passop**2.d0)
      
      passo = -b/(2.d0*a)
      if( a.lt.0.d0) then
         if(ene1.lt.ene0) then
            passo=passop
         else 
            passo=0.5*passop
         endif
      endif


      stima=a*passo**2.d0+b*passo+c


      return
    end subroutine minparabola


subroutine pc2(a,beca,b,becb)      
               
! this function applies the operator Pc
            
!    this subroutine applies the Pc operator
!    a input :unperturbed wavefunctions
!    b input :first order wavefunctions
!    b output:b_i =b_i-a_j><a_j|S|b_i>
     
      use ions_base, only: na, nsp
      use io_global, only: stdout
      use cvan 
      use gvecw, only: ngw
      use constants, only: pi, fpi
      use control_flags, only: iprint, iprsta
      use reciprocal_vectors, only: ng0 => gstart
      use mp, only: mp_sum
      use electrons_base, only: n => nbsp, fspin
      use uspp_param, only: nh
      use uspp, only :nhsa=>nkb
      use uspp, only :qq
                           
      implicit none        
                           
      complex(8) a(ngw,n), b(ngw,n)
                     
      real(8)    beca(nhsa,n),becb(nhsa,n)
! local variables
      integer is, iv, jv, ia, inl, jnl, i, j,ig
      real(8) sca
      do i=1,n
         do j=1,n
            sca=0.
            if (ng0.eq.2) then
               b(1,i)=0.5d0*(b(1,i)+CONJG(b(1,i)))
            endif
            if(fspin(j) == fspin(i)) then
               do  ig=1,ngw           !loop on g vectors
                  sca=sca+2.d0*DBLE(CONJG(a(ig,j))*b(ig,i)) !2. for real wavefunctions
               enddo
               if (ng0.eq.2) then
                  sca=sca-DBLE(CONJG(a(1,j))*b(1,i))
               endif

               call mp_sum( sca )

               if (nvb.gt.0) then

                  do is=1,nvb
                     do iv=1,nh(is)
                        do jv=1,nh(is)
                           do ia=1,na(is)
                              inl=ish(is)+(iv-1)*na(is)+ia
                              jnl=ish(is)+(jv-1)*na(is)+ia
                              sca=sca+ qq(iv,jv,is)*beca(inl,j)*becb(jnl,i)  
                           end do
                        end do
                     end do
                  end do
               end if

               do ig=1,ngw
                  b(ig,i)=b(ig,i)-sca*a(ig,j)
               enddo
               ! this to prevent numerical errors
               if (ng0.eq.2) then
                  b(1,i)=0.5d0*(b(1,i)+CONJG(b(1,i)))
               endif
            endif
         enddo

      enddo
      return
      end subroutine pc2


      subroutine pcdaga2(a,as ,b )

! this function applies the operator Pc

!    this subroutine applies the Pc^dagerr operator
!    a input :unperturbed wavefunctions
!    b input :first order wavefunctions
!    b output:b_i =b_i - S|a_j><a_j|b_i>

      use ions_base, only: na, nsp
      use io_global, only: stdout
      use cvan
      use gvecw, only: ngw
      use constants, only: pi, fpi
      use control_flags, only: iprint, iprsta
      use reciprocal_vectors, only: ng0 => gstart
      use mp, only: mp_sum
      use electrons_base, only: n => nbsp, fspin
      use uspp_param, only: nh
      use uspp, only :nhsa=>nkb

      implicit none

      complex(8) a(ngw,n), b(ngw,n), as(ngw,n)
      real(8)    beca(nhsa,n),becb(nhsa,n)
      ! local variables
      integer is, iv, jv, ia, inl, jnl, i, j,ig
      real(8) sca
      !
      do i=1,n
         do j=1,n
            sca=0.
            if (ng0.eq.2) then
               b(1,i)=0.5d0*(b(1,i)+CONJG(b(1,i)))
            endif
            if(fspin(i) == fspin(j)) then
               do  ig=1,ngw           !loop on g vectors
                  sca=sca+2.*DBLE(CONJG(a(ig,j))*b(ig,i)) !2. for real weavefunctions
               enddo
               if (ng0.eq.2) then
                  sca=sca-DBLE(CONJG(a(1,j))*b(1,i))
               endif
               call mp_sum(sca)
               do ig=1,ngw
                  b(ig,i)=b(ig,i)-sca*as(ig,j)
               enddo
               ! this to prevent numerical errors
               if (ng0.eq.2) then
                  b(1,i)=0.5d0*(b(1,i)+CONJG(b(1,i)))
               endif
            endif
         enddo
      enddo
      return
    end subroutine pcdaga2

