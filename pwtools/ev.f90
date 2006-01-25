
!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
      PROGRAM ev
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
      IMPLICIT NONE
      INTEGER nmaxpar, nmaxpt, nseek, nmin, npar,npt,istat
      PARAMETER(nmaxpar=4, nmaxpt=100, nseek=10000, nmin=4)
      CHARACTER bravais*3, filin*20
      REAL*8 par(nmaxpar), deltapar(nmaxpar), parmin(nmaxpar), &
             parmax(nmaxpar), v0(nmaxpt), etot(nmaxpt), efit(nmaxpt), &
             fac, emin, chisq, a
      COMMON v0, etot, emin, efit, istat, npt
!
      PRINT '(5x,''Enter type of bravais lattice '', &
             &   ''(fcc, bcc, sc, hex) > '',$)'
      READ '(a)',bravais
!
      IF(bravais.EQ.'fcc'.OR.bravais.EQ.'FCC') THEN
         fac = 0.25
      ELSE IF(bravais.EQ.'bcc'.OR.bravais.EQ.'BCC') THEN
         fac = 0.50
      ELSE IF(bravais.EQ.'sc'.OR.bravais.EQ.'SC') THEN
         fac = 1.0
      ELSE IF(bravais.EQ.'hex'.OR.bravais.EQ.'HEX') THEN
         fac = SQRT(3d0)/2d0
      ELSE
         PRINT '(5x,''ev: unexpected lattice '',a3)', bravais
         STOP
      ENDIF
!
      PRINT '(5x,''Enter type of equation of state :''/&
             &5x,''1=birch1, 2=birch2, 3=keane, 4=murnaghan > '',$)'
      READ *,istat
      IF(istat.EQ.1 .OR. istat.EQ.4) THEN
         npar=3
      ELSE IF(istat.EQ.2 .OR. istat.EQ.3) THEN
         npar=4
      ELSE
         PRINT '(5x,''Unexpected eq. of state '',i2)', istat
         STOP
      END IF
      PRINT '(5x,''Input file > '',$)'
      READ '(a)',filin
      OPEN(unit=2,file=filin,status='old',form='formatted')
  10  CONTINUE
      emin=1d10
      DO npt=1,nmaxpt
         IF (bravais.EQ.'hex'.OR.bravais.EQ.'HEX') THEN
            READ(2,*,err=10,END=20) v0(npt), etot(npt) 
         ELSE
            READ(2,*,err=10,END=20) a, etot(npt) 
            v0  (npt) = fac*a**3
         END IF
         IF(etot(npt).LT.emin) THEN
            par(1) = v0(npt)
            emin = etot(npt)
         END IF
      ENDDO

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
      CALL find_minimum &
           (npar,par,deltapar,parmin,parmax,nseek,nmin,chisq)
!
      CALL write_results &
           (npt,bravais,fac,v0,etot,efit,istat,par,npar,emin,chisq)
!
      STOP
      END   
!
!-----------------------------------------------------------------------
      SUBROUTINE eqstate(npar,par,chisq)
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
      INTEGER nmaxpt, npt, npar, i, istat
      PARAMETER( nmaxpt=100 )
      REAL*8 par(npar), k0, dk0, d2k0, c0, c1, x, &
             vol0, v0(nmaxpt), etot(nmaxpt), efit(nmaxpt), &
             ddk, emin, conv_atomic_unit, chisq
      COMMON v0, etot, emin, efit, istat, npt
      PARAMETER ( conv_atomic_unit=6.79777e-6 )
!
      vol0 = par(1)
      k0   = par(2)*conv_atomic_unit ! converts k0 to atomic units...
      dk0  = par(3)
      d2k0 = par(4)/conv_atomic_unit ! and d2k0/dp2 ito (a.u.)^(-1)
!
      IF(istat.EQ.1.OR.istat.EQ.2) THEN
         IF(istat.EQ.1) THEN
            c0 = 0.0
         ELSE
            c0 = ( 9.0*k0*d2k0 + 9.0*dk0**2 - 63.0*dk0 + 143.0 ) / 48.0
         ENDIF
         c1 = 3.0*(dk0-4.0)/8.0
         DO i=1,npt
            x = vol0/v0(i)
            efit(i) = 9.0*k0*vol0*( (-0.5+  c1-  c0)*x**(2d0/3d0)/2.0 &
                                   +( 0.5-2*c1+3*c0)*x**(4d0/3d0)/4.0 &
                                   +(       c1-3*c0)*x**(6d0/3d0)/6.0 &
                                   +(            c0)*x**(8d0/3d0)/8.0 &
                                   -(-1d0/8.0+c1/6.0-c0/8.0) )
         ENDDO
      ELSE
         IF(istat.EQ.3) THEN
            ddk = dk0 + k0*d2k0/dk0
         ELSE
            ddk = dk0
         ENDIF
         DO i=1,npt
            efit(i) = - k0*dk0/ddk*vol0/(ddk-1.0) &
            + v0(i)*k0*dk0/ddk**2*( (vol0/v0(i))**ddk/(ddk-1.0) + 1.0 ) &
            - k0*(dk0-ddk)/ddk*( v0(i)*LOG(vol0/v0(i)) + v0(i)-vol0 )
         ENDDO
      ENDIF
!
!      emin = equilibrium energy obtained by minimizing chi**2
!
      emin = 0.0
      DO i = 1,npt
         emin = emin + etot(i)-efit(i)
      ENDDO
      emin = emin/npt
!
      chisq = 0.0
      DO i = 1,npt
          efit(i) = efit(i)+emin
          chisq   = chisq + (etot(i)-efit(i))**2
      ENDDO
      chisq = chisq/npt
!
      RETURN
      END
!
!-----------------------------------------------------------------------
      SUBROUTINE write_results &
            (npt,bravais,fac,v0,etot,efit,istat,par,npar,emin,chisq)
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
      INTEGER npt, istat, npar, i, iun, ios
      CHARACTER filout*20, bravais*3
      REAL*8 v0(npt), etot(npt), efit(npt), par(npar), emin, chisq, fac
      REAL*8 p(npt), birch, keane


 10   CONTINUE
      PRINT '(5x,''Output file > '',$)'
      READ '(a)',filout
      IF(filout.NE.' ') THEN
         iun=8
         OPEN(unit=iun,file=filout,form='formatted',status='new', &
              iostat=ios)
         IF (ios.NE.0) THEN
            PRINT '(5x,''Cannot open, try again! '')'
            go to 10
         END IF
      ELSE
         iun=6
      END IF

      IF(istat.EQ.1) THEN
         WRITE(iun,'(''# equation of state: '', &
                    &''birch 1st order.  chisq = '',d10.4)') chisq
      ELSE IF(istat.EQ.2) THEN
         WRITE(iun,'(''# equation of state: '', &
                    &''birch 2nd order.  chisq = '',d10.4)') chisq
      ELSE IF(istat.EQ.3) THEN
         WRITE(iun,'(''# equation of state: '', &
                    &''keane.            chisq = '',d10.4)') chisq
      ELSE IF(istat.EQ.4) THEN
         WRITE(iun,'(''# equation of state: '', &
                    &''murnaghan.        chisq = '',d10.4)') chisq
      END IF

      IF(istat.EQ.1 .OR. istat.EQ.4) par(4) = 0.0

      IF(istat.EQ.1 .OR. istat.EQ.2) THEN
         DO i=1,npt
            p(i)=birch(v0(i)/par(1),par(2),par(3),par(4))
         END DO
      ELSE
         DO i=1,npt
            p(i)=keane(v0(i)/par(1),par(2),par(3),par(4))
         END DO
      END IF

      IF(bravais.NE.'hex'.AND.bravais.NE.'HEX') THEN
         WRITE(iun,'(''# a0 ='',f6.2,''  k0 ='',i5,'' kbar,  dk0 ='', &
                    &f6.2,''  d2k0 ='',f7.3,''  emin ='',f11.5/)') &
            (par(1)/fac)**(1d0/3d0), INT(par(2)), par(3), par(4), emin
         WRITE(iun,'(f7.3,2f12.5,3x,f8.2,3x,f12.5)') &
              ( (v0(i)/fac)**(1d0/3d0), etot(i), efit(i), p(i), &
              etot(i)-efit(i), i=1,npt ) 
      ELSE
         WRITE(iun,'(''# V0 ='',f8.2,''  k0 ='',i5,'' kbar,  dk0 ='', &
                    & f6.2,''  d2k0 ='',f7.3,''  emin ='',f11.5/)') &
                     par(1), INT(par(2)), par(3), par(4), emin
         WRITE(iun,'(f8.2,2f12.5,3x,f8.2,3x,f12.5)') &
             ( v0(i), etot(i), efit(i), p(i), etot(i)-efit(i), i=1,npt )
      END IF

      IF(filout.NE.' ') CLOSE(unit=iun)
      RETURN
      END
!
!-----------------------------------------------------------------------
      SUBROUTINE find_minimum &
         (npar,par,deltapar,parmin,parmax,nseek,nmin,chisq)
!-----------------------------------------------------------------------
!
!     Very Stupid Minimization
!
      USE random_numbers, ONLY : rndm
      IMPLICIT NONE
      INTEGER maxpar, nseek, npar, nmin, n,j,i
      PARAMETER (maxpar=4)
      REAL*8 par(npar), deltapar(npar), parmin(npar), parmax(npar), &
             parnew(maxpar), chisq, chinew, bidon
!
!      various initializations
!
      chisq = 1.0e30
      chinew= 1.0e30
      CALL eqstate(npar,par,chisq)
      DO j = 1,nmin
         DO i = 1,nseek
            DO n = 1,npar
  10           parnew(n) = par(n) + (0.5 - rndm())*deltapar(n)
               IF(parnew(n).GT.parmax(n) .OR. parnew(n).LT.parmin(n)) &
               go to 10
            ENDDO
!
            CALL eqstate(npar,parnew,chinew)
!
            IF(chinew.LT.chisq) THEN
               DO n = 1,npar
                  par(n) = parnew(n)
               ENDDO
               chisq = chinew
            ENDIF
         ENDDO
         DO n = 1,npar
            deltapar(n) = deltapar(n)/10.0
         ENDDO
      ENDDO
!
      CALL eqstate(npar,par,chisq)
!
      RETURN
      END
!
      FUNCTION birch(x,k0,dk0,d2k0)
!
      IMPLICIT NONE
      REAL*8 birch, x, k0,dk0, d2k0
      REAL*8 c0, c1

      IF(d2k0.NE.0.0) THEN
         c0 = ( 9.0*k0*d2k0 + 9.0*dk0**2 - 63.0*dk0 + 143.0 ) / 48.0
      ELSE
         c0 = 0.0
      ENDIF
      c1 = 3.0*(dk0-4.0)/8.0
      birch = 3.0*k0*( (-0.5+  c1-  c0)*x**( -5.0/3d0) &
           +( 0.5-2*c1+3*c0)*x**( -7.0/3d0) &
           +(       c1-3*c0)*x**( -9.0/3d0) &
           +(            c0)*x**(-11.0/3d0) )
      RETURN
      END
!     
      FUNCTION keane(x,k0,dk0,d2k0)
!     
      IMPLICIT NONE
      REAL*8 keane, x, k0, dk0, d2k0, ddk
      
      ddk = dk0 + k0*d2k0/dk0
      keane = k0*dk0/ddk**2*( x**(-ddk) - 1d0 ) + (dk0-ddk)/ddk*LOG(x)
      
      RETURN
      END
