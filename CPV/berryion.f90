!#include "../include/machine.h"
! this subroutine returns the berry phase energy
! = L/2*Pi*Im(log Sum_R exp(i*(2pi/L)*R_i*rho_i))
! of the ions and the constant force on the ions
! now only for orthorombic primitive cell

subroutine berryion( tau0,fion, tfor,ipol,evalue,enbi)

!  tau0    : input, positions of ions
!  fion    : input,output, forces on ions
!  tfor    : input, flag for force calculation
!  ipol    : input, electric field polarization
!  evalue  : input, scale for electric field
!  enbi    : output, berry phase energy of the ions

     

  
  use elct
  use constants, only: pi, fpi
  use ions_base
  use parameters, only: natx
  use cell_base, only: a1, a2, a3

  implicit none

  real(kind=8) tau0(3,*)
  real(kind=8) fion(3,*)
  real(kind=8) enbi, evalue
  integer ipol, isa
  logical tfor

!local variables
  real(kind=8) gmes
  real(kind=8) pola
  integer is, ia
  complex(kind=8) temp, ci

  temp = (0.,0.)
  ci = (0.,1.)

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
  pola=0.d0
  isa = 0
  do is=1,nsp
     do ia=1,na(is)
        isa = isa + 1

!this force term is along ipol-direction
        if( tfor) then
           fion(ipol,isa)=fion(ipol,isa)+evalue*zv(is)
        endif
            
        temp = temp - ci*gmes*tau0(ipol,isa)*zv(is)
        pola=pola+evalue*zv(is)*tau0(ipol,isa)!this is just the center of ionic charge
     enddo
  enddo

  enbi=aimag(log(exp(temp)))/gmes!this sounds stupid it's just a Riemann plane
!  write(6,*) 'Pola  :', pola!ATTENZIONE
  return
end subroutine berryion

            
!-------------------------------------------------------------------------
      subroutine cofcharge(tau,cdz)
!-----------------------------------------------------------------------
!this subroutine gives the center of the ionic charge


      use ions_base, only: na, nsp, zv
!
      implicit none
      real(kind=8) tau(3,*), cdz(3)
! local variables
      real(kind=8) zmas
      integer is,i,ia,isa
!
      zmas=0.0
      do is=1,nsp
         zmas=zmas+na(is)*zv(is)
      end do
!
      isa = 0
      do i=1,3
         cdz(i)=0.0
         do is=1,nsp
            do ia=1,na(is)
               isa = isa + 1
               cdz(i)=cdz(i)+tau(i,isa)*zv(is)
            end do
         end do
         cdz(i)=cdz(i)/zmas
      end do
      write(6,*) 'Center of charge', cdz(3)!ATTENZIONE
!
      return
      end
!



!----------------------------------------------------
        subroutine noforce(fion, ipol)
!----------------------------------------------------

! this subroutine adds a electric force, in order
! to keep steady the center of mass along the electric
! field direction

          use ions_base

          implicit none

          real(kind=8) fion(3,*)
          integer ipol!el. field polarization


          integer i,ia,is,isa
          real(kind=8) fcm!force appplied on center of mass
          real(kind=8) tch!total charge

          fcm=0.d0
          tch=0.d0
          isa = 0
          do is=1,nsp
             do ia=1,na(is)
                isa = isa + 1
                fcm=fcm+fion(ipol,isa)     
                tch=tch+zv(is)
             enddo             
          enddo
          write(6,*) 'Forza su ioni in ipol:', fcm!ATTENZIONE 
          fcm=fcm/tch
          isa = 0
          do is=1,nsp
             do ia=1,na(is)
                isa = isa + 1
                fion(ipol,isa)=fion(ipol,isa)-fcm*zv(is)
             enddo
          enddo
   
          return
        end subroutine noforce

  
