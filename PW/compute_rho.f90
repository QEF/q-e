!
! Copyright (C) 2005 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
      subroutine compute_rho(rho,rhoout,segni,nrxx)
!
!   This subroutine diagonalize the spin density matrix and gives as output
!   the spin up and spin down compotents of the charge
!
   USE kinds, ONLY : dp
   USE cell_base,  ONLY : at, bg
   USE symme,      ONLY : symv
   USE constants, ONLY: pi
   USE io_global,  ONLY :  stdout

   implicit none
   integer :: nrxx  ! input: the dimension of the mesh
         
   real(DP) ::  rho(nrxx,4),segni(nrxx),rhoout(nrxx,2)
                 ! input: the four components of the charge 
                 ! output: keep track of the spin direction
                 ! output: the spin up and spin down charge

   complex(DP) :: re(2,2) ! the rotated density matrix

   real(DP) :: eps,amag,aux  
   real(DP) :: mx,my,mz         ! magnetization
   real(DP) :: ux(3)         ! magnetization
   real(DP) :: ux0(3)      ! magnetization

   logical :: negative
   integer :: ir           ! counter on mesh points
   integer :: it,ipol,count(3)

      eps=1.d-12
     
      ux=0.d0
      ux0=0.d0
      do ir=min(1,nrxx/2),nrxx
         amag=sqrt(rho(ir,2)**2+rho(ir,3)**2+rho(ir,4)**2)
         if (amag.gt.eps) then
             DO ipol=1,3
                ux(ipol)=rho(ir,ipol+1)/amag
             ENDDO
             ux0=ux
             goto 120
         endif
      enddo
120   continue
      ux(1)=1.d0
      ux(2)=2.d0
      ux(3)=3.d0
      ux0(1)=1.d0
      ux0(2)=2.d0
      ux0(3)=3.d0

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
                  if (abs(ux(1)*mx+ux(2)*my+ux(3)*mz).lt.1.d-3) then
                     ux(1)=ux(1)+0.5d0*mx
                     ux(2)=ux(2)+0.5d0*my
                     ux(3)=ux(3)+0.5d0*mz
                     amag=sqrt(ux(1)**2+ux(2)**2+ux(3)**2)
                     ux=ux/amag
                  endif
               else
                  if (abs(ux(1)*mx+ux(2)*my+ux(3)*mz).lt.1.d-3) &
                                       count(it)=count(it)+1
               endif
            endif
         enddo
         if (count(1).eq.0) goto 100
      enddo
      if (count(1).lt.count(3)) THEN 
         ux=ux0
         WRITE( stdout,'(5x,"it, count:",3i5)') it,count(1),count(3)
      ENDIF
100   continue
      call symv (ux)
      aux=sqrt(ux(1)**2+ux(2)**2+ux(3)**2)
      ux=ux/aux
      
      WRITE( stdout,'(/,5x, "GGA quantization axis:&
           & (",f8.4,",",f8.4,",",f8.4,"  )")') ux(1),ux(2),ux(3)

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
            negative=(mx*ux(1)+my*ux(2)+mz*ux(3)).gt.0.d0
         endif
         

         if (negative) then
            rhoout(ir,2)= DBLE(re(1,1))
            rhoout(ir,1)= DBLE(re(2,2))
            segni(ir)=-1.d0
         else
            rhoout(ir,1)= DBLE(re(1,1))
            rhoout(ir,2)= DBLE(re(2,2))
            segni(ir)=1.d0
         endif
      enddo


      return
      end subroutine compute_rho
