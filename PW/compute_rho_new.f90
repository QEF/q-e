!
! Copyright (C) 2005 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
      SUBROUTINE compute_rho_new(rho,rhoout,segni,nrxx,ux)
!
!   This subroutine diagonalizes the spin density matrix and gives as output
!   the spin up and spin down components of the charge
!
#include "f_defs.h"
   USE kinds, ONLY : dp
   USE constants, ONLY: pi
   USE io_global,  ONLY :  stdout

   implicit none

   REAL(DP) :: ux(3)     ! fixed direction to calculate signs
   INTEGER :: nrxx       ! input: the dimension of the mesh
         
   REAL(DP) ::  rho(nrxx,4),segni(nrxx),rhoout(nrxx,2)
                 ! input: the four components of the charge 
                 ! output: keep track of the spin direction
                 ! output: the spin up and spin down charge

   COMPLEX(DP) :: re(2,2) ! the rotated density matrix

   REAL(DP) :: eps,amag,uxmod ! pi, a small number and the atan function 
   REAL(DP) :: m(3)       ! magnetization
   REAL(DP) :: ux0(3)     ! magnetization

   LOGICAL :: negative
   INTEGER :: ir, i       ! counter on mesh points
   INTEGER :: it,counter(3)

   eps=1.d-12
   uxmod=sqrt(ux(1)**2+ux(2)**2+ux(3)**2) 
   if (uxmod<eps) then

      do i=1,3
         ux(i)=i
         ux0(i)=i
      enddo
      counter(3)=0

      do it=1,3 
         counter(it)=0
         do ir=1,nrxx
            amag=SQRT(rho(ir,2)**2+rho(ir,3)**2+rho(ir,4)**2)
            if (amag > eps) then
               do i=1,3
                  m(i)=rho(ir,i+1)/amag
               enddo
               if (it==2) then
                  if (ABS(ux(1)*m(1)+ux(2)*m(2)+ux(3)*m(3)).lt.1.d-3) then
                     do i=1,3
                        ux(i)=ux(i)+0.5d0*m(i)
                     enddo
                     uxmod=SQRT(ux(1)**2+ux(2)**2+ux(3)**2)
                     do i=1,3
                        ux(i)=ux(i)/uxmod
                     enddo
                  endif
               else
                  if (ABS(ux(1)*m(1)+ux(2)*m(2)+ux(3)*m(3)).lt.1.d-3) &
                                               counter(it)=counter(it)+1
               endif
            endif
         enddo
         if (counter(1)==0) goto 100
      enddo
      if (counter(1).lt.counter(3)) then
         do i=1,3
            ux(i)=ux0(i)
         enddo
      endif
!     WRITE( stdout,'(5x,"it, counter:",3i5)') it,counter(1),counter(3)
      WRITE( stdout,'(5x,"Fixed quantization axis: ", 3f12.6)') &
                          ux(1), ux(2), ux(3)
   endif
100 continue
   do ir=1,nrxx
      amag=SQRT(rho(ir,2)**2+rho(ir,3)**2+rho(ir,4)**2)
      re(1,1)=0.5_DP*(rho(ir,1)+amag)
      re(2,2)=0.5_DP*(rho(ir,1)-amag)
      if (amag.lt.eps) then
         negative=.false.
      else
         do i=1,3
            m(i)=rho(ir,i+1)/amag
         enddo
         negative=(m(1)*ux(1)+m(2)*ux(2)+m(3)*ux(3))>0.d0
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
   end subroutine compute_rho_new
