      subroutine compute_rho(rho,rhoout,segni,nrxx)
!
!   This subroutine diagonalize the spin density matrix and gives as output
!   the spin up and spin down compotents of the charge
!
#include "f_defs.h"
   USE kinds, ONLY : dp
   USE io_global,  ONLY :  stdout

   implicit none
   integer :: nrxx  ! input: the dimension of the mesh
         
   real(kind=dp) ::  rho(nrxx,4),segni(nrxx),rhoout(nrxx,2)
                 ! input: the four components of the charge 
                 ! output: keep track of the spin direction
                 ! output: the spin up and spin down charge

   complex(kind=dp) :: umat(2,2),rhom(2,2),re(2,2)
                    ! the rotation matrix
                    ! the density matrix
                    ! the rotated density matrix

   real(kind=dp) :: theta,phi        ! the two angles of the rotation matrix
   real(kind=dp) :: ct,st,cp,sp      ! sinus and cosinus of theta or phi
   real(kind=dp) :: pi,eps,amag      ! pi, a small number and the atan function 
   real(kind=dp) :: mx,my,mz         ! magnetization
   real(kind=dp) :: ux,uy,uz         ! magnetization
   real(kind=dp) :: ux0,uy0,uz0      ! magnetization

   logical :: negative
   integer :: ir           ! counter on mesh points
   integer :: it,count(3)
   integer :: i,j,s1,s2    ! counter on spin

      pi=4.d0*atan(1.d0)
      eps=1.d-12
     
      ux=0.d0
      uy=0.d0
      uz=0.d0
      ux0=0.d0
      uy0=0.d0
      uz0=0.d0
      do ir=min(1,nrxx/2),nrxx
         amag=sqrt(rho(ir,2)**2+rho(ir,3)**2+rho(ir,4)**2)
         if (amag.gt.eps) then
             ux=rho(ir,2)/amag
             uy=rho(ir,3)/amag
             uz=rho(ir,4)/amag
             ux0=ux
             uy0=uy
             uz0=uz
             goto 120
         endif
      enddo
120   continue
      ux=1.d0
      uy=2.d0
      uz=3.d0
      ux0=1.d0
      uy0=2.d0
      uz0=3.d0

      count(3)=0

      do it=1,3 
         count(it)=0
         do ir=1,nrxx
            amag=sqrt(rho(ir,2)**2+rho(ir,3)**2+rho(ir,4)**2)
            if (amag.gt.eps) then
               mx=rho(ir,2)/amag
               my=rho(ir,3)/amag
               mz=rho(ir,4)/amag
               if (it.eq.2) then
                  if (abs(ux*mx+uy*my+uz*mz).lt.1.d-3) then
                     ux=ux+0.5d0*mx
                     uy=uy+0.5d0*my
                     uz=uz+0.5d0*mz
                     amag=sqrt(ux**2+uy**2+uz**2)
                     ux=ux/amag
                     uy=uy/amag
                     uz=uz/amag
                  endif
               else
                  if (abs(ux*mx+uy*my+uz*mz).lt.1.d-3) &
                                       count(it)=count(it)+1
               endif
            endif
         enddo
         if (count(1).eq.0) goto 100
      enddo
      if (count(1).lt.count(3)) then
         ux=ux0
         uy=uy0
         uz=uz0
      endif
100   WRITE( stdout,*) 'it, count', it,count(1),count(3)
      WRITE( stdout,*) ux,uy,uz

      do ir=1,nrxx
         amag=sqrt(rho(ir,2)**2+rho(ir,3)**2+rho(ir,4)**2)
         re(1,1)=0.5d0*(rho(ir,1)+amag)
         re(2,2)=0.5d0*(rho(ir,1)-amag)
         if (amag.lt.eps) then
            negative=.false.
         else
            mx=rho(ir,2)/amag
            my=rho(ir,3)/amag
            mz=rho(ir,4)/amag
            negative=(mx*ux+my*uy+mz*uz).gt.0.d0
         endif
         

         if (negative) then
            rhoout(ir,2)=DREAL(re(1,1))
            rhoout(ir,1)=DREAL(re(2,2))
            segni(ir)=-1.d0
         else
            rhoout(ir,1)=DREAL(re(1,1))
            rhoout(ir,2)=DREAL(re(2,2))
            segni(ir)=1.d0
         endif
      enddo


      return
      end
