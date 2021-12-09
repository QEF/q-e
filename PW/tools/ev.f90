!
! Copyright (C) 2003-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributions by Eyvaz Isaev
! Dept of Physics, Chemistry and Biology (IFM), Linkoping University, Sweden
!
! a) Input: Add lattice parameters units: au or Ang
! b) Output: More info printed out
! c) Output: Additional output file with E+PV
!
PROGRAM ev
!
!      fit of E(v) to an equation of state (EOS)
!
!      Interactive input:
!         au or Ang
!         structure
!         equation of state
!         input data file
!         output data file
!
!      Input data file format for cubic systems:
!         a0(1)  Etot(1)
!         ...
!         a0(n)  Etot(n)
!      where a0 is the lattice parameter (a.u. or Ang)
!      Input data file format for noncubic (e.g. hexagonal) systems:
!         V0(1)  Etot(1)
!         ...
!         V0(n)  Etot(n)
!      where V0 is the unit-cell volume (a.u.^3 or Ang^3)
!      e.g. for an hexagonal cell,
!         V0(i)  = sqrt(3)/2 * a^2 * c    unit-cell volume
!         Etot(i)= min Etot(c)   for the given volume V0(i)
!      Etot in atomic (Rydberg) units
!
!      Output data file format  for cubic systems:
!      # a0=... a.u., K0=... kbar, dk0=..., d2k0=... kbar^-1, Emin=... Ry
!      # a0=... Ang,  K0=... GPa , V0=... (a.u.)^3, V0 = Ang^3
!         a0(1)  Etot(1) Efit(1)  Etot(1)-Efit(1)  Pfit(1)  Enth(1)
!         ...
!         a0(n)  Etot(n) Efit(n)  Etot(n)-Efit(n)  Pfit(n)  Enth(n)
!      Output data file format  for noncubic systems:
!      # V0=...(a.u.)^3, K0=... kbar, dk0=..., d2k0=... kbar^-1, Emin=... Ry
!      # V0=...Ang^3,  K0=... GPa
!         V0(1)  Etot(1) Efit(1)  Etot(1)-Efit(1)  Pfit(1)  Enth(1)
!         ...
!         V0(n)  Etot(n) Efit(n)  Etot(n)-Efit(n)  Pfit(n)  Enth(n)
!      where
!            a0(i), V0(i), Etot(i) as in input
!            Efit(i) is the fitted value from the EOS
!            Pfit(i) is the corresponding pressure from the EOS (GPa)
!            Enth(i)=Efit(i)+Pfit(i)*V0(i) is the enthalpy (Ry)
!!
      USE kinds, ONLY: DP
      USE constants, ONLY: bohr_radius_angs, ry_kbar
      USE mp_global, ONLY : mp_startup, mp_global_end
      USE mp_world,  ONLY : world_comm
      USE mp,        ONLY : mp_bcast
      USE io_global, ONLY : ionode, ionode_id
      IMPLICIT NONE
      INTEGER, PARAMETER:: nmaxpar=4, nmaxpt=100 
      INTEGER :: npar,npt,istat, ierr
      CHARACTER :: bravais*8, au_unit*3, filin*256
      REAL(DP) :: par(nmaxpar), v0(nmaxpt), etot(nmaxpt), efit(nmaxpt), &
             fac, emin, chisq, a
      REAL(DP), PARAMETER :: gpa_kbar = 10.0_dp
      LOGICAL :: in_angstrom
      CHARACTER(LEN=256) :: fileout
  !
  CALL mp_startup ( )
  !
  IF ( ionode ) THEN

      WRITE(*,'(5x,"Lattice parameter or Volume are in (au, Ang) > ")', advance="NO")
      READ(5,'(a)') au_unit
      in_angstrom = au_unit=='Ang' .or. au_unit=='ANG' .or. &
                    au_unit=='ang'
      IF (in_angstrom) WRITE(*,'(5x,"Assuming Angstrom")')
      WRITE(*,'(5x,"Enter type of bravais lattice (fcc, bcc, sc, noncubic) > ")', advance="NO")
      READ(5,'(a)') bravais
!
      IF(trim(bravais)=='fcc'.or.trim(bravais)=='FCC') THEN
         fac = 0.25d0
      ELSEIF(trim(bravais)=='bcc'.or.trim(bravais)=='BCC') THEN
         fac = 0.50d0
      ELSEIF(trim(bravais)=='sc'.or.trim(bravais)=='SC') THEN
         fac = 1.0d0
      ELSEIF(bravais=='noncubic'.or.bravais=='NONCUBIC' .or.  &
             trim(bravais)=='hex'.or.trim(bravais)=='HEX' ) THEN
!         fac = sqrt(3d0)/2d0 ! not used
         fac = 0.0_DP ! not used
      ELSE
         WRITE(*,'(5x,"ev: unexpected lattice <",a,">")') trim(bravais)
         STOP
      ENDIF
!
      WRITE(*,'(5x,"Enter type of equation of state :"/&
            &  5x,"1=birch1, 2=birch2, 3=keane, 4=murnaghan > ")', advance="NO")
      READ(5,*) istat
      IF(istat==1 .or. istat==4) THEN
         npar=3
      ELSEIF(istat==2 .or. istat==3) THEN
         npar=4
      ELSE
         WRITE(*,'(5x,"Unexpected eq. of state ",i2)') istat
         STOP
      ENDIF
      WRITE(*,'(5x,"Input file > ")', advance="NO")
      READ(5,'(a)') filin
      OPEN(unit=2,file=filin,status='old',form='formatted',iostat=ierr)
      IF (ierr/=0) THEN
         ierr= 1 
         GO TO 99
      END IF
  10  CONTINUE
      emin=1d10
      DO npt=1,nmaxpt
         IF (bravais=='noncubic'.or.bravais=='NONCUBIC' .or. &
             bravais=='hex'.or.bravais=='HEX' ) THEN
            READ(2,*,err=10,END=20) v0(npt), etot(npt)
            IF (in_angstrom) v0(npt)=v0(npt)/bohr_radius_angs**3
         ELSE
            READ(2,*,err=10,END=20) a, etot(npt)
            IF (in_angstrom) a = a/bohr_radius_angs
            v0  (npt) = fac*a**3
         ENDIF
         IF(etot(npt)<emin) THEN
            par(1) = v0(npt)
            emin = etot(npt)
         ENDIF
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

      CALL find_minimum &
           (npar,par,chisq)
!
      CALL write_results &
           (npt,in_angstrom,fac,v0,etot,efit,istat,par,npar,emin,chisq, &
            fileout)
!
      CALL write_evdata_xml  &
           (npt,fac,v0,etot,efit,istat,par,npar,emin,chisq,fileout, ierr)

      IF (ierr /= 0) GO TO 99
    ENDIF
99  CALL mp_bcast ( ierr, ionode_id, world_comm )
    IF ( ierr == 1) THEN
       CALL errore( 'ev', 'file '//trim(filin)//' cannot be opened', ierr )
    ELSE IF ( ierr == 2 ) THEN
       CALL errore( 'ev', 'file '//trim(fileout)//' cannot be opened', ierr )
    ELSE IF ( ierr == 11 ) THEN
       CALL errore( 'write_evdata_xml', 'no free units to write ', ierr )
    ELSE IF ( ierr == 12 ) THEN
       CALL errore( 'write_evdata_xml', 'error opening the xml file ', ierr )
    ENDIF

  CALL mp_global_end()

      STOP
    CONTAINS
!
!-----------------------------------------------------------------------
      SUBROUTINE eqstate(npar,par,chisq,ediff)
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
      INTEGER, INTENT(in) :: npar
      REAL(DP), INTENT(in) :: par(npar)
      REAL(DP), INTENT(out):: chisq
      REAL(DP), OPTIONAL, INTENT(out):: ediff(npt)
      INTEGER :: i
      REAL(DP) :: k0, dk0, d2k0, c0, c1, x, vol0, ddk
!
      vol0 = par(1)
      k0   = par(2)/ry_kbar ! converts k0 to Ry atomic units...
      dk0  = par(3)
      d2k0 = par(4)*ry_kbar ! and d2k0/dp2 to (Ry a.u.)^(-1)
!
      IF(istat==1.or.istat==2) THEN
         IF(istat==1) THEN
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
         IF(istat==3) THEN
            ddk = dk0 + k0*d2k0/dk0
         ELSE
            ddk = dk0
         ENDIF
         DO i=1,npt
            efit(i) = - k0*dk0/ddk*vol0/(ddk-1.d0) &
            + v0(i)*k0*dk0/ddk**2*( (vol0/v0(i))**ddk/(ddk-1.d0)+1.d0) &
            - k0*(dk0-ddk)/ddk*( v0(i)*log(vol0/v0(i)) + v0(i)-vol0 )
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
          IF(present(ediff)) ediff(i) = efit(i)-etot(i)
       ENDDO
      chisq = chisq/npt
!
      RETURN
    END SUBROUTINE eqstate
!
!-----------------------------------------------------------------------
      SUBROUTINE write_results &
            (npt,in_angstrom,fac,v0,etot,efit,istat,par,npar,emin,chisq, &
             filout)
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
      INTEGER, INTENT(in) :: npt, istat, npar
      REAL(DP), INTENT(in):: v0(npt), etot(npt), efit(npt), emin, chisq, fac
      REAL(DP), INTENT(inout):: par(npar)
      REAL(DP), EXTERNAL :: keane, birch
      LOGICAL, INTENT(in) :: in_angstrom
      CHARACTER(len=256), intent(inout) :: filout
      !
      REAL(DP) :: p(npt), epv(npt)
      INTEGER :: i, iun
      LOGICAL :: exst

      WRITE(*,'(5x,"Output file > ")', advance="NO")
      READ (5,'(a)') filout
      IF(filout/=' ') THEN
         iun=8
         INQUIRE(file=filout,exist=exst)
         IF (exst) WRITE(*,'(5x,"Beware: file ",A," will be overwritten")')&
                  trim(filout)
         OPEN(unit=iun,file=filout,form='formatted',status='unknown', &
              iostat=ierr)
         IF (ierr/=0) THEN
            ierr= 2 
            GO TO 99
         END IF
      ELSE
         iun=6
      ENDIF

      IF(istat==1) THEN
         WRITE(iun,'("# equation of state: birch 1st order.  chisq = ", &
                   & d12.4)') chisq
      ELSEIF(istat==2) THEN
         WRITE(iun,'("# equation of state: birch 3rd order.  chisq = ", &
                   & d12.4)') chisq
      ELSEIF(istat==3) THEN
         WRITE(iun,'("# equation of state: keane.            chisq = ", &
                   & d12.4)') chisq
      ELSEIF(istat==4) THEN
         WRITE(iun,'("# equation of state: murnaghan.        chisq = ", &
                   & d12.4)') chisq
      ENDIF

      IF(istat==1 .or. istat==4) par(4) = 0.0d0

      IF(istat==1 .or. istat==2) THEN
         DO i=1,npt
            p(i)=birch(v0(i)/par(1),par(2),par(3),par(4))
         ENDDO
      ELSE
         DO i=1,npt
            p(i)=keane(v0(i)/par(1),par(2),par(3),par(4))
         ENDDO
      ENDIF

      DO i=1,npt
         epv(i) = etot(i) + p(i)*v0(i) / ry_kbar
      ENDDO

      IF ( fac /= 0.0_dp ) THEN
! cubic case
         WRITE(iun,'("# a0 =",f8.4," a.u., k0 =",i5," kbar, dk0 =", &
                    &f6.2," d2k0 =",f7.3," emin =",f11.5)') &
            (par(1)/fac)**(1d0/3d0), int(par(2)), par(3), par(4), emin
         WRITE(iun,'("# a0 =",f9.5," Ang, k0 =", f6.1," GPa,  V0 = ", &
                  & f8.2," (a.u.)^3,  V0 =", f8.2," A^3 ",/)') &
           & (par(1)/fac)**(1d0/3d0)*bohr_radius_angs, par(2)/gpa_kbar, &
             par(1), par(1)*bohr_radius_angs**3

        WRITE(iun,'(73("#"))')
        WRITE(iun,'("# Lat.Par", 7x, "E_calc", 8x, "E_fit", 7x, &
             & "E_diff", 4x, "Pressure", 6x, "Enthalpy")')
        IF (in_angstrom) THEN
           WRITE(iun,'("# Ang", 13x, "Ry", 11x, "Ry", 12x, &
             & "Ry", 8x, "GPa", 11x, "Ry")')
           WRITE(iun,'(73("#"))')
           WRITE(iun,'(f9.5,2x,f12.5, 2x,f12.5, f12.5, 3x, f8.2, 3x,f12.5)') &
              & ( (v0(i)/fac)**(1d0/3d0)*bohr_radius_angs, etot(i), efit(i),  &
              & etot(i)-efit(i), p(i)/gpa_kbar, epv(i), i=1,npt )
        ELSE
           WRITE(iun,'("# a.u.",12x, "Ry", 11x, "Ry", 12x, &
             & "Ry", 8x, "GPa", 11x, "Ry")')
           WRITE(iun,'(73("#"))')
           WRITE(iun,'(f9.5,2x,f12.5, 2x,f12.5, f12.5, 3x, f8.2, 3x,f12.5)') &
              & ( (v0(i)/fac)**(1d0/3d0), etot(i), efit(i),  &
              & etot(i)-efit(i), p(i)/gpa_kbar, epv(i), i=1,npt )
        ENDIF

      ELSE
! noncubic case
         WRITE(iun,'("# V0 =",f8.2," a.u.^3,  k0 =",i5," kbar,  dk0 =", &
                    & f6.2,"  d2k0 =",f7.3,"  emin =",f11.5)') &
                    & par(1), int(par(2)), par(3), par(4), emin

         WRITE(iun,'("# V0 =",f8.2,"  Ang^3,  k0 =",f6.1," GPa"/)') &
                    & par(1)*bohr_radius_angs**3, par(2)/gpa_kbar

        WRITE(iun,'(74("#"))')
        WRITE(iun,'("# Vol.", 8x, "E_calc", 8x, "E_fit", 7x, &
             & "E_diff", 4x, "Pressure", 6x, "Enthalpy")')
        IF (in_angstrom) THEN
          WRITE(iun,'("# Ang^3", 9x, "Ry", 11x, "Ry", 12x, &
             & "Ry", 8x, "GPa", 11x, "Ry")')
          WRITE(iun,'(74("#"))')
           WRITE(iun,'(f8.2,2x,f12.5, 2x,f12.5, f12.5, 3x, f8.2, 3x,f12.5)') &
              ( v0(i)*bohr_radius_angs**3, etot(i), efit(i),  &
               etot(i)-efit(i), p(i)/gpa_kbar, epv(i), i=1,npt )
         else
          WRITE(iun,'("# a.u.^3",8x, "Ry", 11x, "Ry", 12x, &
             & "Ry", 8x, "GPa", 11x, "Ry")')
          WRITE(iun,'(74("#"))')
           WRITE(iun,'(f8.2,2x,f12.5, 2x,f12.5, f12.5, 3x, f8.2, 3x,f12.5)') &
              ( v0(i), etot(i), efit(i),  &
               etot(i)-efit(i), p(i)/gpa_kbar, epv(i), i=1,npt )
         end if

      ENDIF
      IF(filout/=' ') CLOSE(unit=iun)
 99   RETURN
    END SUBROUTINE write_results
!
      
      ! This subroutine is passed to LMDIF to be minimized
      ! LMDIF takes as input the difference between f_fit and f_real
      !       and computes the chi^2 internally.
      SUBROUTINE EOSDIFF(m_, n_, par_, f_, i_)
       IMPLICIT NONE
         INTEGER,INTENT(in)  :: m_, n_
         INTEGER,INTENT(inout)   :: i_
         REAL(DP),INTENT(in)    :: par_(n_)
         REAL(DP),INTENT(out)   :: f_(m_)
         REAL(DP) :: chisq_
         !
         CALL eqstate(n_,par_,chisq_, f_)
      END SUBROUTINE

!-----------------------------------------------------------------------
      SUBROUTINE find_minimum(npar,par,chisq)
!-----------------------------------------------------------------------
!
      USE lmdif_module, ONLY : lmdif0
      IMPLICIT NONE
      INTEGER ,INTENT(in)  :: npar
      REAL(DP),INTENT(out) :: par(nmaxpar)
      REAL(DP),INTENT(out) :: chisq
      !
      REAL(DP) :: vchisq(npar)
      REAL(DP) :: ediff(npt)
      INTEGER :: i
      !
      par(1) = v0(npt/2)
      par(2) = 500.0d0
      par(3) = 5.0d0
      par(4) = -0.01d0 ! unused for some eos
      !      
      CALL lmdif0(EOSDIFF, npt, npar, par, ediff, 1.d-12, i)
      !
      IF(i>0 .and. i<5) THEN
         PRINT*, "Minimization succeeded"
      ELSEIF(i>=5) THEN
         PRINT*, "Minimization stopped before convergence"
      ELSEIF(i<=0) THEN 
        PRINT*, "Minimization error"
        STOP
      ENDIF
      !
      CALL eqstate(npar,par,chisq)

      END SUBROUTINE
!-----------------------------------------------------------------------

  END PROGRAM ev

      FUNCTION birch(x,k0,dk0,d2k0)
!
      USE kinds, ONLY : DP
      IMPLICIT NONE
      REAL(DP) birch, x, k0,dk0, d2k0
      REAL(DP) c0, c1

      IF(d2k0/=0.d0) THEN
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
      USE kinds, ONLY : DP
      IMPLICIT NONE
      REAL(DP) keane, x, k0, dk0, d2k0, ddk

      ddk = dk0 + k0*d2k0/dk0
      keane = k0*dk0/ddk**2*( x**(-ddk) - 1d0 ) + (dk0-ddk)/ddk*log(x)

      RETURN
    END FUNCTION keane
!-----------------------------------------------------------------------
    SUBROUTINE write_evdata_xml &
        (npt,fac,v0,etot,efit,istat,par,npar,emin,chisq,filout, ierr)
!-----------------------------------------------------------------------
!
  USE kinds,     ONLY : dp
  USE constants, ONLY : ry_kbar, bohr_radius_angs
  IMPLICIT NONE
  INTEGER, INTENT(in) :: npt, istat, npar
  REAL(DP), INTENT(in):: v0(npt), etot(npt), efit(npt), emin, chisq, fac
  REAL(DP), INTENT(in):: par(npar)
  CHARACTER(len=256), INTENT(IN) :: filout
  INTEGER, INTENT(out) :: ierr
  !
  INTEGER :: iunout = 11
  REAL(DP) :: p(npt), volume(2), a0(2), alldata(6,npt)
  INTEGER :: i, iun
  CHARACTER(len=256) :: filename
  REAL(DP), EXTERNAL :: birch, keane

  IF (filout/=' ') THEN
     filename = TRIM(filout) // '.xml'
  ELSE
     filename = 'ev.xml'
  ENDIF
  !
  ! ... open XML descriptor
  !
  OPEN ( UNIT=iunout, FILE = TRIM( filename ), FORM='formatted', IOSTAT = ierr )
  IF ( ierr /= 0 ) THEN
     WRITE (6,*) 'Failed opening file ' // TRIM(filename)
     RETURN
  END IF

  WRITE (iunout,'("<xml>")')
  WRITE (iunout,'("<EQUATIONS_OF_STATE>")')
  WRITE (iunout,'("<EQUATION_TYPE>")')
  IF (istat==1) THEN
     WRITE (iunout,'("Birch 1st order")')
  ELSEIF (istat==2) THEN
     WRITE (iunout,'("Birch 2nd order")')
  ELSEIF (istat==3) THEN
     WRITE (iunout,'("Keane")')
  ELSEIF (istat==4) THEN
     WRITE (iunout,'("Murnaghan")')
  ENDIF
  WRITE (iunout,'("</EQUATION_TYPE>")')
  WRITE (iunout,'("<CHI_SQUARE>")')
  WRITE (iunout,'(1pe25.12)') chisq
  WRITE (iunout,'("</CHI_SQUARE>")')
  WRITE (iunout,'("</EQUATIONS_OF_STATE>")')

  IF (istat==1 .or. istat==2) THEN
     DO i=1,npt
        p(i)=birch(v0(i)/par(1),par(2),par(3),par(4))
     ENDDO
  ELSE
     DO i=1,npt
        p(i)=keane(v0(i)/par(1),par(2),par(3),par(4))
     ENDDO
  ENDIF

  DO i=1,npt
     alldata (1,i) = v0(i)
     alldata (2,i) = etot(i) 
     alldata (3,i) = efit(i)
     alldata (4,i) = etot(i) - efit(i)
     alldata (5,i) = p(i) 
     alldata (6,i) = etot(i) + p(i) * v0(i) / ry_kbar
  ENDDO

  WRITE (iunout,'("<EQUATIONS_PARAMETERS>")')

  volume(1)=par(1)
  volume(2)=par(1)*bohr_radius_angs**3
  WRITE (iunout, '("<EQUILIBRIUM_VOLUME_AU_A>")')
  WRITE (iunout, '(2(1pe25.15))') volume(:)
  WRITE (iunout, '("</EQUILIBRIUM_VOLUME_AU_A>")')
  WRITE (iunout, '("<BULK_MODULUS_KBAR>")')
  WRITE (iunout, '(1pe25.15)') par(2)
  WRITE (iunout, '("</BULK_MODULUS_KBAR>")')
  WRITE (iunout, '("<DERIVATIVE_BULK_MODULUS>")')
  WRITE (iunout, '(1pe25.15)') par(3)
  WRITE (iunout, '("</DERIVATIVE_BULK_MODULUS>")')
  WRITE (iunout, '("<SECOND_DERIVATIVE_BULK_MODULUS>")')
  WRITE (iunout, '(1pe25.15)') par(4)
  WRITE (iunout, '("</SECOND_DERIVATIVE_BULK_MODULUS>")')
  WRITE (iunout, '("<MINIMUM_ENERGY_RY>")')
  WRITE (iunout, '(1pe25.15)') emin
  WRITE (iunout, '("</MINIMUM_ENERGY_RY>")')
  WRITE (iunout, '("<CELL_FACTOR>")')
  WRITE (iunout, '(1pe25.15)') fac
  WRITE (iunout, '("</CELL_FACTOR>")')
  IF (fac /= 0.0_DP) THEN
     a0(1) = (par(1)/fac)**(1d0/3d0)
     a0(2) = (par(1)/fac)**(1d0/3d0) * bohr_radius_angs
     WRITE (iunout, '("<CELL_PARAMETER_AU_A>")')
     WRITE (iunout, '(2(1pe25.15))') a0
     WRITE (iunout, '("</CELL_PARAMETER_AU_A>")')
  ENDIF
  WRITE (iunout,'("</EQUATIONS_PARAMETERS>")')

  WRITE (iunout,'("<FIT_CHECK>")')
  WRITE (iunout,'("<NUMBER_OF_DATA>")')
  WRITE (iunout,'(i8)') npt
  WRITE (iunout,'("</NUMBER_OF_DATA>")')
  WRITE (iunout,'("<VOL_ENE_EFIT_DELTA_P_GIBBS>")')
  WRITE (iunout,'(6(1pe25.15))') alldata(:,:)
  WRITE (iunout,'("</VOL_ENE_EFIT_DELTA_P_GIBBS>")')

  WRITE (iunout,'("</FIT_CHECK>")')
  CLOSE (unit=iunout, status='keep')
  
  RETURN
END SUBROUTINE write_evdata_xml

