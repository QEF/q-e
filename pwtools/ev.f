
!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
      program ev
!
!      fit of E(v) to an equation of state (EOS)
!
!      Input data file format for cubic systems:
!         a0(1)  Etot(1)
!         ...
!         a0(n)  Etot(n)
!      Input data file format for hexagonal systems:
!         V0(1)  Etot(1)
!         ...
!         V0(n)  Etot(n)
!      where V0(i)  = sqrt(3)/2 * a^2 * c   is the unit-cell volume
!            Etot(i)= min Etot(c)   for the given volume V0(i)
!      V0, a0, etot in atomic (Rydberg) units ; bulk modulus in kbar
!
!      Output data file format  for cubic systems:
!         a0(1)  Etot(1) Efit(1)  Pfit(1)  Etot(1)-Efit(1)
!         ...
!         a0(n)  Etot(n) Efit(n)  Pfit(n)  Etot(n)-Efit(n)
!      Output data file format  for hexagonal systems:
!         V0(1)  Etot(1) Efit(1)  Pfit(1)  Etot(1)-Efit(1)
!         ...
!         V0(n)  Etot(n) Efit(n)  Pfit(n)  Etot(n)-Efit(n)
!
!      where Efit is the fitted value from the EOS
!            Pfit is the corresponding pressure from the EOS
!
      implicit none
      integer nmaxpar, nmaxpt, nseek, nmin, npar,npt,istat
      parameter(nmaxpar=4, nmaxpt=100, nseek=10000, nmin=4)
      character bravais*3, filin*20
      real*8 par(nmaxpar), deltapar(nmaxpar), parmin(nmaxpar),
     &       parmax(nmaxpar), v0(nmaxpt), etot(nmaxpt), efit(nmaxpt), 
     &       fac, emin, chisq, a
      common v0, etot, emin, efit, istat, npt
!
      print '(5x,''Enter type of bravais lattice '',
     &           ''(fcc, bcc, sc, hex) > '',$)'
      read '(a)',bravais
!
      if(bravais.eq.'fcc'.or.bravais.eq.'FCC') then
         fac = 0.25
      else if(bravais.eq.'bcc'.or.bravais.eq.'BCC') then
         fac = 0.50
      else if(bravais.eq.'sc'.or.bravais.eq.'SC') then
         fac = 1.0
      else if(bravais.eq.'hex'.or.bravais.eq.'HEX') then
         fac = sqrt(3d0)/2d0
      else
         print '(5x,''ev: unexpected lattice '',a3)', bravais
         stop
      endif
!
      print '(5x,''Enter type of equation of state :''/
     &        5x,''1=birch1, 2=birch2, 3=keane, 4=murnaghan > '',$)'
      read *,istat
      if(istat.eq.1 .or. istat.eq.4) then
         npar=3
      else if(istat.eq.2 .or. istat.eq.3) then
         npar=4
      else
         print '(5x,''Unexpected eq. of state '',i2)', istat
         stop
      end if
      print '(5x,''Input file > '',$)'
      read '(a)',filin
      open(unit=2,file=filin,status='old',form='formatted')
  10  continue
      emin=1d10
      do npt=1,nmaxpt
         if (bravais.eq.'hex'.or.bravais.eq.'HEX') then
            read(2,*,err=10,end=20) v0(npt), etot(npt) 
         else
            read(2,*,err=10,end=20) a, etot(npt) 
            v0  (npt) = fac*a**3
         end if
         if(etot(npt).lt.emin) then
            par(1) = v0(npt)
            emin = etot(npt)
         end if
      enddo

      npt = nmaxpt+1
  20  npt = npt-1
!
! par(1) = V, Volume of the unit cell in (a.u.^3)
! par(2) = B, Bulk Modulus (in KBar)
! par(3) = dB/dP (adimensional)
! par(4) = d^2B/dP^2 (in KBar^(-1), used only by 2nd order formulae)
!
      par(2) = 500.0          
      par(3) = 5.0
      par(4) = -0.01
!
      parmin(1) = 0.0
      parmin(2) = 0.0
      parmin(3) = 1.0
      parmin(4) = -1.0
!
      parmax(1) = 100000.
      parmax(2) = 100000.
      parmax(3) = 15.0
      parmax(4) = 0.0
!
      deltapar(1) = 0.1
      deltapar(2) = 100.
      deltapar(3) = 1.0
      deltapar(4) = 0.01
!
      call find_minimum
     &     (npar,par,deltapar,parmin,parmax,nseek,nmin,chisq)
!
      call write_results
     &     (npt,bravais,fac,v0,etot,efit,istat,par,npar,emin,chisq)
!
      stop
      end   
!
!-----------------------------------------------------------------------
      subroutine eqstate(npar,par,chisq)
!-----------------------------------------------------------------------
!
      implicit none
      integer nmaxpt, npt, npar, i, istat
      parameter( nmaxpt=100 )
      real*8 par(npar), k0, dk0, d2k0, c0, c1, x,
     &       vol0, v0(nmaxpt), etot(nmaxpt), efit(nmaxpt),
     &       ddk, emin, conv_atomic_unit, chisq
      common v0, etot, emin, efit, istat, npt
      parameter ( conv_atomic_unit=6.79777e-6 )
!
      vol0 = par(1)
      k0   = par(2)*conv_atomic_unit ! converts k0 to atomic units...
      dk0  = par(3)
      d2k0 = par(4)/conv_atomic_unit ! and d2k0/dp2 ito (a.u.)^(-1)
!
      if(istat.eq.1.or.istat.eq.2) then
         if(istat.eq.1) then
            c0 = 0.0
         else
            c0 = ( 9.0*k0*d2k0 + 9.0*dk0**2 - 63.0*dk0 + 143.0 ) / 48.0
         endif
         c1 = 3.0*(dk0-4.0)/8.0
         do i=1,npt
            x = vol0/v0(i)
            efit(i) = 9.0*k0*vol0*( (-0.5+  c1-  c0)*x**(2d0/3d0)/2.0
     &                             +( 0.5-2*c1+3*c0)*x**(4d0/3d0)/4.0
     &                             +(       c1-3*c0)*x**(6d0/3d0)/6.0
     &                             +(            c0)*x**(8d0/3d0)/8.0
     &                             -(-1d0/8.0+c1/6.0-c0/8.0) )
         enddo
      else
         if(istat.eq.3) then
            ddk = dk0 + k0*d2k0/dk0
         else
            ddk = dk0
         endif
         do i=1,npt
            efit(i) = - k0*dk0/ddk*vol0/(ddk-1.0)
     &      + v0(i)*k0*dk0/ddk**2*( (vol0/v0(i))**ddk/(ddk-1.0) + 1.0 )
     &      - k0*(dk0-ddk)/ddk*( v0(i)*log(vol0/v0(i)) + v0(i)-vol0 )
         enddo
      endif
!
!      emin = equilibrium energy obtained by minimizing chi**2
!
      emin = 0.0
      do i = 1,npt
         emin = emin + etot(i)-efit(i)
      enddo
      emin = emin/npt
!
      chisq = 0.0
      do i = 1,npt
          efit(i) = efit(i)+emin
          chisq   = chisq + (etot(i)-efit(i))**2
      enddo
      chisq = chisq/npt
!
      return
      end
!
!-----------------------------------------------------------------------
      subroutine write_results
     &      (npt,bravais,fac,v0,etot,efit,istat,par,npar,emin,chisq)
!-----------------------------------------------------------------------
!
      implicit none
      integer npt, istat, npar, i, iun, ios
      character filout*20, bravais*3
      real*8 v0(npt), etot(npt), efit(npt), par(npar), emin, chisq, fac
      real*8 p(npt), birch, keane


 10   continue
      print '(5x,''Output file > '',$)'
      read '(a)',filout
      if(filout.ne.' ') then
         iun=8
         open(unit=iun,file=filout,form='formatted',status='new',
     $        iostat=ios)
         if (ios.ne.0) then
            print '(5x,''Cannot open, try again! '')'
            go to 10
         end if
      else
         iun=6
      end if

      if(istat.eq.1) then
         write(iun,'(''# equation of state: '',
     &        ''birch 1st order.  chisq = '',d10.4)') chisq
      else if(istat.eq.2) then
         write(iun,'(''# equation of state: '',
     &        ''birch 2nd order.  chisq = '',d10.4)') chisq
      else if(istat.eq.3) then
         write(iun,'(''# equation of state: '',
     &        ''keane.            chisq = '',d10.4)') chisq
      else if(istat.eq.4) then
         write(iun,'(''# equation of state: '',
     &        ''murnaghan.        chisq = '',d10.4)') chisq
      end if

      if(istat.eq.1 .or. istat.eq.4) par(4) = 0.0

      if(istat.eq.1 .or. istat.eq.2) then
         do i=1,npt
            p(i)=birch(v0(i)/par(1),par(2),par(3),par(4))
         end do
      else
         do i=1,npt
            p(i)=keane(v0(i)/par(1),par(2),par(3),par(4))
         end do
      end if

      if(bravais.ne.'hex'.and.bravais.ne.'HEX') then
         write(iun,'(''# a0 ='',f6.2,''  k0 ='',i5,'' kbar,  dk0 ='',
     &                f6.2,''  d2k0 ='',f7.3,''  emin ='',f11.5/)')
     &      (par(1)/fac)**(1d0/3d0), int(par(2)), par(3), par(4), emin
         write(iun,'(f7.3,2f12.5,3x,f8.2,3x,f12.5)') 
     &        ( (v0(i)/fac)**(1d0/3d0), etot(i), efit(i), p(i), 
     &        etot(i)-efit(i), i=1,npt ) 
      else
         write(iun,'(''# V0 ='',f8.2,''  k0 ='',i5,'' kbar,  dk0 ='',
     &                f6.2,''  d2k0 ='',f7.3,''  emin ='',f11.5/)')
     &               par(1), int(par(2)), par(3), par(4), emin
         write(iun,'(f8.2,2f12.5,3x,f8.2,3x,f12.5)') ( v0(i),
     &        etot(i), efit(i), p(i), etot(i)-efit(i), i=1,npt )
      end if

      if(filout.ne.' ') close(unit=iun)
      return
      end
!
!-----------------------------------------------------------------------
      subroutine find_minimum
     &   (npar,par,deltapar,parmin,parmax,nseek,nmin,chisq)
!-----------------------------------------------------------------------
!
!     Very Stupid Minimization
!
      implicit none
      integer maxpar, nseek, npar, nmin, n,j,i
      parameter (maxpar=4)
      real*8 par(npar), deltapar(npar), parmin(npar), parmax(npar), 
     &       parnew(maxpar), chisq, chinew, bidon
      real*8 rndm
      external rndm
!
!      various initializations
!
      chisq = 1.0e30
      chinew= 1.0e30
      call eqstate(npar,par,chisq)
      do j = 1,nmin
         do i = 1,nseek
            do n = 1,npar
  10           parnew(n) = par(n) + (0.5 - rndm())*deltapar(n)
               if(parnew(n).gt.parmax(n) .or. parnew(n).lt.parmin(n)) 
     &         go to 10
            enddo
!
            call eqstate(npar,parnew,chinew)
!
            if(chinew.lt.chisq) then
               do n = 1,npar
                  par(n) = parnew(n)
               enddo
               chisq = chinew
            endif
         enddo
         do n = 1,npar
            deltapar(n) = deltapar(n)/10.0
         enddo
      enddo
!
      call eqstate(npar,par,chisq)
!
      return
      end
!
      function birch(x,k0,dk0,d2k0)
!
      implicit none
      real*8 birch, x, k0,dk0, d2k0
      real*8 c0, c1

      if(d2k0.ne.0.0) then
         c0 = ( 9.0*k0*d2k0 + 9.0*dk0**2 - 63.0*dk0 + 143.0 ) / 48.0
      else
         c0 = 0.0
      endif
      c1 = 3.0*(dk0-4.0)/8.0
      birch = 3.0*k0*( (-0.5+  c1-  c0)*x**( -5.0/3d0)
     &     +( 0.5-2*c1+3*c0)*x**( -7.0/3d0)
     &     +(       c1-3*c0)*x**( -9.0/3d0)
     &     +(            c0)*x**(-11.0/3d0) )
      return
      end
!     
      function keane(x,k0,dk0,d2k0)
!     
      implicit none
      real*8 keane, x, k0, dk0, d2k0, ddk
      
      ddk = dk0 + k0*d2k0/dk0
      keane = k0*dk0/ddk**2*( x**(-ddk) - 1d0 ) + (dk0-ddk)/ddk*log(x)
      
      return
      end
