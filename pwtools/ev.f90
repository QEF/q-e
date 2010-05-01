!
! Copyright (C) 2003-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Last edition:
! Date: April 30, 2010
! Author: Eyvaz Isaev 
! Department of Physics, Chemistry and Biology (IFM), Linkoping University, Sweden
! 
! a) Input: Add lattice parameters units: au or Ang
! b) Output: More info printed out
! c) Output: Additional output file with E+PV
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
      USE kinds, only: DP
      USE constants, only: bohr_radius_angs, ry_kbar
      IMPLICIT NONE
      INTEGER nmaxpar, nmaxpt, nseek, nmin, npar,npt,istat
      PARAMETER(nmaxpar=4, nmaxpt=100, nseek=10000, nmin=4)
      CHARACTER bravais*3, au_unit*3, filin*256
      REAL(DP) par(nmaxpar), deltapar(nmaxpar), parmin(nmaxpar), &
             parmax(nmaxpar), v0(nmaxpt), etot(nmaxpt), efit(nmaxpt), &
             fac, emin, chisq, a
      LOGICAL :: convert_to_au
!
!
      PRINT '(5x,''Lattice parameter or Volume are in '', &
             &   ''(au, Ang) > '',$)' 
      READ '(a)', au_unit
      convert_to_au = au_unit.eq.'Ang' .or. au_unit.eq.'ANG' .or. &
                        au_unit.eq.'ang'
      IF (convert_to_au) PRINT '(5x,"Assuming Angstrom")'
 
      PRINT '(5x,''Enter type of bravais lattice '', &
             &   ''(fcc, bcc, sc, hex) > '',$)'
      READ '(a)',bravais
!
      IF(bravais.EQ.'fcc'.OR.bravais.EQ.'FCC') THEN
         fac = 0.25d0
      ELSE IF(bravais.EQ.'bcc'.OR.bravais.EQ.'BCC') THEN
         fac = 0.50d0
      ELSE IF(bravais.EQ.'sc'.OR.bravais.EQ.'SC') THEN
         fac = 1.0d0
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
	    IF (convert_to_au) v0(npt)=v0(npt)/bohr_radius_angs**3
         ELSE
            READ(2,*,err=10,END=20) a, etot(npt)
	    IF (convert_to_au) a = a/bohr_radius_angs
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
      par(2) = 500.0d0          
      par(3) = 5.0d0
      par(4) = -0.01d0
!
      parmin(1) = 0.0d0
      parmin(2) = 0.0d0
      parmin(3) = 1.0d0
      parmin(4) = -1.0d0
!
      parmax(1) = 100000.d0
      parmax(2) = 100000.d0
      parmax(3) = 15.0d0
      parmax(4) = 0.0d0
!
      deltapar(1) = 0.1d0
      deltapar(2) = 100.d0
      deltapar(3) = 1.0d0
      deltapar(4) = 0.01d0
!
      CALL find_minimum &
           (npar,par,deltapar,parmin,parmax,nseek,nmin,chisq)
!
      CALL write_results &
           (npt,bravais,fac,v0,etot,efit,istat,par,npar,emin,chisq)
!
      STOP
    CONTAINS  
!
!-----------------------------------------------------------------------
      SUBROUTINE eqstate(npar,par,chisq)
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
      INTEGER nmaxpt, npar, i
      PARAMETER( nmaxpt=100 )
      REAL(DP) par(npar), k0, dk0, d2k0, c0, c1, x, vol0, ddk, chisq
!
      vol0 = par(1)
      k0   = par(2)/ry_kbar ! converts k0 to Ry atomic units...
      dk0  = par(3)
      d2k0 = par(4)*ry_kbar ! and d2k0/dp2 to (Ry a.u.)^(-1)
!
      IF(istat.EQ.1.OR.istat.EQ.2) THEN
         IF(istat.EQ.1) THEN
            c0 = 0.0d0
         ELSE
            c0 = ( 9.d0*k0*d2k0 + 9.d0*dk0**2-63.d0*dk0+143.d0 )/48.d0
         ENDIF
         c1 = 3.d0*(dk0-4.d0)/8.d0
         DO i=1,npt
            x = vol0/v0(i)
            efit(i) = 9.d0*k0*vol0*( (-0.5d0+c1-c0)*x**(2.d0/3.d0)/2.d0 &
                         +( 0.50-2.d0*c1+3.d0*c0)*x**(4.d0/3.d0)/4.d0 &
                         +(       c1-3.d0*c0)*x**(6.d0/3.d0)/6.d0 &
                         +(            c0)*x**(8.d0/3.d0)/8.d0 &
                         -(-1.d0/8.d0+c1/6.d0-c0/8.d0) )
         ENDDO
      ELSE
         IF(istat.EQ.3) THEN
            ddk = dk0 + k0*d2k0/dk0
         ELSE
            ddk = dk0
         ENDIF
         DO i=1,npt
            efit(i) = - k0*dk0/ddk*vol0/(ddk-1.d0) &
            + v0(i)*k0*dk0/ddk**2*( (vol0/v0(i))**ddk/(ddk-1.d0)+1.d0) &
            - k0*(dk0-ddk)/ddk*( v0(i)*LOG(vol0/v0(i)) + v0(i)-vol0 )
         ENDDO
      ENDIF
!
!      emin = equilibrium energy obtained by minimizing chi**2
!
      emin = 0.0d0
      DO i = 1,npt
         emin = emin + etot(i)-efit(i)
      ENDDO
      emin = emin/npt
!
      chisq = 0.0d0
      DO i = 1,npt
          efit(i) = efit(i)+emin
          chisq   = chisq + (etot(i)-efit(i))**2
      ENDDO
      chisq = chisq/npt
!
      RETURN
    END SUBROUTINE eqstate
!
!-----------------------------------------------------------------------
      SUBROUTINE write_results &
            (npt,bravais,fac,v0,etot,efit,istat,par,npar,emin,chisq)
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
      INTEGER npt, istat, npar, i, iun, ios
      CHARACTER filout*256, bravais*3
      REAL(DP) v0(npt), etot(npt), efit(npt), par(npar), emin, chisq, fac
      REAL(DP) p(npt), epv(npt)

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

      IF(istat.EQ.1 .OR. istat.EQ.4) par(4) = 0.0d0

      IF(istat.EQ.1 .OR. istat.EQ.2) THEN
         DO i=1,npt
            p(i)=birch(v0(i)/par(1),par(2),par(3),par(4))
         END DO
      ELSE
         DO i=1,npt
            p(i)=keane(v0(i)/par(1),par(2),par(3),par(4))
         END DO
      END IF
      
      do i=1,npt
         epv(i) = etot(i) + p(i)*v0(i) / ry_kbar
      enddo

      IF(bravais.NE.'hex'.AND.bravais.NE.'HEX') THEN
         WRITE(iun,'(''# a0 ='',f6.2,'' a.u., k0 ='',i5,'' kbar, dk0 ='', &
                    &f6.2,'' d2k0 ='',f7.3,'' emin ='',f11.5)') &
            (par(1)/fac)**(1d0/3d0), INT(par(2)), par(3), par(4), emin
         WRITE(iun,'(''# a0 ='',f7.3,'' Ang, k0 ='', f6.1,'' GPa,  V0 = '',f7.3,'' (a.u.)^3,  V0 ='', &
                    &f7.3,'' A^3 '',/)') &
              (par(1)/fac)**(1d0/3d0)*bohr_radius_angs, par(2)/10, par(1), &
           &  par(1)*bohr_radius_angs**3


	write(iun,'(73("#"))') 
	write(iun,'("# Lat.Par", 5x, "E_calc", 8x, "E_fit", 7x, &
	     & "E_diff", 4x, "Pressure", 6x, "Enthalpy")')
	write(iun,'("# Ang", 11x, "Ry", 11x, "Ry", 12x, &
	     & "Ry", 8x, "GPa", 11x, "Ry")')
	write(iun,'(73("#"))') 


         WRITE(iun,'(f7.3,2x,f12.5, 2x,f12.5, f12.5, 3x, f8.2, 3x,f12.5)') &
              & ( (v0(i)/fac)**(1d0/3d0)*bohr_radius_angs, etot(i), efit(i),  &
              & etot(i)-efit(i), p(i), epv(i), i=1,npt ) 
      ELSE
         WRITE(iun,'(''# V0 ='',f8.2,'' a.u.^3,  k0 ='',i5,'' kbar,  dk0 ='', &
                    & f6.2,''  d2k0 ='',f7.3,''  emin ='',f11.5)') &
                    & par(1), INT(par(2)), par(3), par(4), emin


         WRITE(iun,'("# V0 =",f8.2,"  Ang^3,  k0 =",f6.1," GPa"/)') &
                    & par(1)*bohr_radius_angs**3, par(2)/10 


	write(iun,'(74("#"))') 
	write(iun,'("# Vol.", 8x, "E_calc", 8x, "E_fit", 7x, &
	     & "E_diff", 4x, "Pressure", 6x, "Enthalpy")')
	write(iun,'("# Ang^3", 9x, "Ry", 11x, "Ry", 12x, &
	     & "Ry", 8x, "GPa", 11x, "Ry")')
	write(iun,'(74("#"))') 


         WRITE(iun,'(f8.2,2x,f12.5, 2x,f12.5, f12.5, 3x, f8.2, 3x,f12.5)') &
              ( v0(i)*bohr_radius_angs**3, etot(i), efit(i),  &
               etot(i)-efit(i), p(i), epv(i), i=1,npt ) 


      END IF

      IF(filout.NE.' ') CLOSE(unit=iun)
      RETURN
    END SUBROUTINE write_results
!
!-----------------------------------------------------------------------
      SUBROUTINE find_minimum &
         (npar,par,deltapar,parmin,parmax,nseek,nmin,chisq)
!-----------------------------------------------------------------------
!
!     Very Stupid Minimization
!
      USE random_numbers, ONLY : randy
      IMPLICIT NONE
      INTEGER maxpar, nseek, npar, nmin, n,j,i
      PARAMETER (maxpar=4)
      REAL(DP) par(npar), deltapar(npar), parmin(npar), parmax(npar), &
             parnew(maxpar), chisq, chinew, bidon
!
!      various initializations
!
      chisq = 1.0d30
      chinew= 1.0d30
      CALL eqstate(npar,par,chisq)
      DO j = 1,nmin
         DO i = 1,nseek
            DO n = 1,npar
  10           parnew(n) = par(n) + (0.5d0 - randy())*deltapar(n)
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
            deltapar(n) = deltapar(n)/10.d0
         ENDDO
      ENDDO
!
      CALL eqstate(npar,par,chisq)
!
      RETURN
    END SUBROUTINE find_minimum
!
      FUNCTION birch(x,k0,dk0,d2k0)
!
      IMPLICIT NONE
      REAL(DP) birch, x, k0,dk0, d2k0
      REAL(DP) c0, c1

      IF(d2k0.NE.0.d0) THEN
         c0 = (9.d0*k0*d2k0 + 9.d0*dk0**2 - 63.d0*dk0 + 143.d0 )/48.d0
      ELSE
         c0 = 0.0d0
      ENDIF
      c1 = 3.d0*(dk0-4.d0)/8.d0
      birch = 3.d0*k0*( (-0.5d0+  c1-  c0)*x**( -5.d0/3.d0) &
           +( 0.5d0-2.d0*c1+3.0d0*c0)*x**( -7.d0/3.d0) &
           +(       c1-3*c0)*x**( -9.0d0/3d0) &
           +(            c0)*x**(-11.0d0/3d0) )
      RETURN
    END FUNCTION birch
!     
      FUNCTION keane(x,k0,dk0,d2k0)
!     
      IMPLICIT NONE
      REAL(DP) keane, x, k0, dk0, d2k0, ddk
      
      ddk = dk0 + k0*d2k0/dk0
      keane = k0*dk0/ddk**2*( x**(-ddk) - 1d0 ) + (dk0-ddk)/ddk*LOG(x)
      
      RETURN
    END FUNCTION keane
  END PROGRAM ev
