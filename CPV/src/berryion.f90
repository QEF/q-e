!
! Copyright (C) 2002-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

subroutine berryion( tau0,fion, tfor,ipol,evalue,enbi)

! this subroutine returns the berry phase energy
! = L/2*Pi*Im(log Sum_R exp(i*(2pi/L)*R_i*rho_i))
! of the ions and the constant force on the ions
! now only for orthorombic primitive cell

!  tau0    : input, positions of ions
!  fion    : input,output, forces on ions
!  tfor    : input, flag for force calculation
!  ipol    : input, electric field polarization
!  evalue  : input, scale for electric field
!  enbi    : output, berry phase energy of the ions

  use kinds,      only : dp
  use constants,  only : pi
  use ions_base,  ONLY : nsp, nat, zv, ityp
  use cell_base,  only : alat, at

  implicit none

  real(dp) tau0(3,*)
  real(dp) fion(3,*)
  real(dp) enbi, evalue
  integer ipol
  logical tfor

!local variables
  real(dp) :: gmes, pola
  integer is, ia
  complex(dp) temp, ci
  real(dp), external:: g_mes

  temp = (0.0_dp,0.0_dp)
  ci = (0.0_dp,1.0_dp)

  gmes = g_mes ( ipol, at, alat)
  pola=0.0_dp
  do ia=1,nat
     is = ityp(ia)
     !this force term is along ipol-direction
     if( tfor) then
        fion(ipol,ia)=fion(ipol,ia)+evalue*zv(is)
     endif
     temp = temp - ci*gmes*tau0(ipol,ia)*zv(is)
     pola=pola+evalue*zv(is)*tau0(ipol,ia)!this is just the center of ionic charge
  enddo

  enbi=AIMAG(log(exp(temp)))/gmes!this sounds stupid it's just a Riemann plane
  return
end subroutine berryion

            
!-------------------------------------------------------------------------
      subroutine cofcharge(tau,cdz)
!-----------------------------------------------------------------------
!this subroutine gives the center of the ionic charge

      use kinds, only : dp
      use ions_base, only: na, nsp, zv, ityp, nat
!
      implicit none
      real(dp) tau(3,*), cdz(3)
! local variables
      real(dp) zmas
      integer is,i,ia
!
      zmas=0.0d0
      do is=1,nsp
         zmas=zmas+na(is)*zv(is)
      end do
!
      do i=1,3
         cdz(i)=0.0d0
         do ia=1,nat
            is=ityp(ia)
            cdz(i)=cdz(i)+tau(i,ia)*zv(is)
         end do
         cdz(i)=cdz(i)/zmas
      end do
!      write(6,*) 'Center of charge', cdz(3)!ATTENZIONE
!
      return
      end subroutine cofcharge
!



!----------------------------------------------------
        subroutine noforce(fion, ipol)
!----------------------------------------------------

! this subroutine adds an electric force, in order
! to keep steady the center of mass along the electric
! field direction

          use kinds,     only : dp
          use ions_base, ONLY : zv, nat, ityp

          implicit none

          real(dp) fion(3,*)
          integer ipol!el. field polarization


          integer i,ia
          real(dp) fcm!force appplied on center of mass
          real(dp) tch!total charge

          fcm=0.d0
          tch=0.d0
          do ia=1,nat
             fcm=fcm+fion(ipol,ia)     
             tch=tch+zv(ityp(ia))
          enddo
          fcm=fcm/tch
          do ia=1,nat
             fion(ipol,ia)=fion(ipol,ia)-fcm*zv(ityp(ia))
          enddo
   
          return
        end subroutine noforce

  
