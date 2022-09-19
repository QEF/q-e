!
! Copyright (C) 2002-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE EFERMI( NEL, NBANDS, DEL, NKPTS, OCC, EF, EIGVAL, &
                   entropy, ismear, nspin )
  !-----------------------------------------------------------------------
  !! Fermi energy & smearing package written by A. De Vita in July 1992
  !! from R.J. Needs original version.
  !
  !! This version is by N. Marzari. Cold smearing added in Oct. 1995.
  !
  !! Given a set of weights and the eigenvalues associated to their
  !! k-points for Bz sampling, this subroutine performs two tasks: 
  !
  !! * determines the fermi level and the occupancy of the states
  !!   according to the chosen (ismear) thermal broadening;
  !! * calculates -TS (for schemes 1, 2 and 4, one-half of -TS
  !!   is the entropy correction that should be added to the
  !!   total energy to recover the zero-broadening energy
  !!   (i.e. the true ground-state energy). This is not really
  !!   necessary anymore using schemes 3, 5, and 6: the free
  !!   energy E-TS is automatically independent of the temperature (!)
  !!   up to the fourth (3) or third order (5, 6) in T. Note that
  !!   (5, 6) do not have negative occupation numbers, at variance
  !!   with (3).
  !
  !! The suggested smearing scheme is ismear=6 (cold smearing II,
  !! Marzari et al., Phys. Rev. Lett. 82, 3296 (1999) )
  !
  !! The six smearing schemes (choose one with parameter \(\text{ismear}\))
  !! are:
  !
  !! (1) GAUSSIAN  
  !!  See: c-l fu and k-m ho, Phys. Rev. B 28, 5480 (1983).
  !!  Their implementation was variational but *not* corrected for
  !!  second order deviation in sigma, as also was the similar scheme
  !!  (with opposite sign deviation) in: R.J.Needs, R.M.Martin and O.H.
  !!  Nielsen, Phys. Rev. B 33, 3778 (1986).  
  !!  Using the correction calculated hereafter everything should be ok.
  !!  The smearing function is a gaussian normalised to 2.  
  !!  The occupation function is the associated complementary error
  !!  function. 
  !
  !! (2) FERMI-DIRAC  
  !!  See: M.J.Gillan, J. Phys. Condens. Matter 1, 689 (1989), following 
  !!  the scheme outlined in j.callaway and n.h.march, solid state phys. 38,
  !!  136 (1984), after d.n.mermin, phys. rev 137, a1441 (1965).
  !!  the occupation function is twice the single electron 
  !!  fermi-dirac distribution.
  !
  !! (3) HERMITE-DELTA EXPANSION (1ST ORD. - right delcor now, nmar)  
  !!  See: Methfessel and Paxton, Phys. Rev. B 40, 3616 (1989).
  !!  the smearing function is a truncated expansion of Dirac's delta
  !!  in Hermite polinomials. 
  !!  For the smearing function implemented here the truncation is 
  !!  at the first non-trivial expansion term \(\text{d1}(x)\).
  !!  The occupation function is the associated primitive (note: 
  !!  the occupation function is neither monotonic nor limited 
  !!  between 0.0 and 2.0: please check the compatibility of this with 
  !!  your code version and verify in a test calculation that the 
  !!  Fermi level is uniquely determined).  
  !!  The entropy correction holds up to the third order in delta at least,
  !!  and it is not necessary (put = 0.) for the Hermite delta expansion,
  !!  since the linear entropy term in sigma is zero by construction
  !!  in that case (well, we still need the correct free energy, hence
  !!  \(\text{delcor}\) is set to its true value, \(\text{nmar}\)).
  !
  !! (4) GAUSSIAN SPLINES (nmar)  
  !!  Similar to a gaussian smearing, but does not require the
  !!  function inversion to calculate the gradients on the occupancies. 
  !!  It is thus to be preferred in a scheme in which the occ. are 
  !!  independent variables (N. Marzari).
  !
  !! (5) COLD SMEARING I (nmar)  
  !!  Similar to Methfessel-Paxton (zeroes the linear order in the entropy), 
  !!  but now with positive-definite occupations (note they can be greater 
  !!  than 1). This version has a=-0.5634 (minimization of the bump), not 
  !!  a=-0.8165 (monotonic function in the tail)- N. Marzari.
  !
  !! (6) COLD SMEARING II (nmar)  
  !!  The one to use. (5) and (6) are practically identical; this is more elegant. 
  !!  For a discussion, see Marzari et al., Phys. Rev. Lett. 82, 3296 (1999), 
  !!  or Marzari's PhD thesis (Univ. of Cambridge, 1996), at 
  !!  http://quasiamore.mit.edu/phd
  !
  !-----------------------------------------------------------------------
  !     PLEASE INQUIRE WITH ADV/NMAR FOR REFERENCE & SUGGESTIONS IF
  !     YOU PLAN TO USE THE PRESENT CORRECTED BZ SAMPLING SCHEME 
  !-----------------------------------------------------------------------
  !
  !     ANOTHER NOTE:
  !     Thanks to the possible > 2 or < 0 
  !     orbital occupancies in the general case of smearing function,
  !     (e.g. in the M-P case) the algorithm to find EF has been
  !     chosen to be the robust bisection method (from Numerical
  !     Recipes) to allow for non monotonic relation between total 
  !     NEL (see above) and EF. One value for EF which solves
  !     NEL(EF) - Z = 0   is always found.
  !-----------------------------------------------------------------------

  USE kinds, ONLY : DP

  implicit none
 
  INTEGER, INTENT(IN) :: nel
  !! number of electrons for each unit cell
  INTEGER, INTENT(IN) :: nbands
  !! number of bands for each k-point
  INTEGER, INTENT(IN) :: nkpts
  !! number of k-points
  REAL(kind=DP), INTENT(OUT) :: occ(nbands,nkpts)
  !! the occupancy of each state
  REAL(kind=DP), INTENT(OUT) :: ef
  !! the Fermi energy
  REAL(kind=DP), INTENT(IN) :: eigval(nbands,nkpts)
  !! eigenvalues
  REAL(kind=DP), INTENT(OUT) :: entropy
  !! -TS (such that the variational functional, i.e.
  !! the free energy, is E-TS)
  INTEGER, INTENT(IN) :: ismear
  !! smearing scheme (see routine comments)
  INTEGER, INTENT(IN) :: nspin
  !! 1:spin restricted, 2:spin unrestricted
  REAL(kind=DP), INTENT(IN)  :: del
  !! width of gaussian smearing function

  ! ... local variables
  
  REAL(kind=DP) :: weight(nkpts)       ! The weight of each k-point.
  REAL(kind=DP) :: sort(nbands*nkpts)  ! The eigenvalues are written into sort which is
                                       ! then sorted into ascending numerical value, from
                                       ! which bounds on ef can easily be obtained.
  REAL(kind=DP), EXTERNAL :: FERMID,DELTHM,POSHM,POSHM2, EFERMI_SPLINE
  
  INTEGER, PARAMETER :: JMAX =300             !The max number of bisections to get Ef 
  REAL(kind=DP), PARAMETER :: XACC=1.0D-17    !The desired accuracy on Ef

  INTEGER :: isppt,j,nkp,neig,nn,n, inel, nel2, j2
  REAL(kind=DP) :: fspin, entrofac,entrospin
  REAL(kind=DP) :: pi,ee,eesh,sq2i,piesqq,z,en
  REAL(kind=DP) :: eigmin, eigmax, xe1,xe2,z1
  REAL(kind=DP) :: x,fmid, f, rtbis,dx,xmid,fi,a
  REAL(kind=DP) :: zeta,elow, test
  REAL(kind=DP) :: delcor              ! The correction -0.5*T*S (the correction is needed 
                                       ! for ismear 1,2 and 4 only)

  if ((nspin == 1).or.(nspin == 2)) then
     continue
  else
     write(*,*) 'ERROR: EFERMI with nspin different from 1 or 2'
     stop
  end if
  
  fspin=DBLE(nspin)
  entrofac=3.d0-fspin
  entrospin=2.d0/fspin

!  if ((nspin == 2).and.(ismear == 2)) then
!     write(*,*) 'ERROR: EFERMI with nspin.eq.2 and ismear.ne.2'
!     stop
!  end if
  
  if (nspin == 1) then
     if (2*nbands == nel) then
        DO ISPPT = 1, NKPTS
           DO J = 1,NBANDS
              OCC(J,ISPPT) = 2.0d0 
           end do
        end do
        return
     end if
  else
     if (nbands == nel) then
        DO ISPPT = 1, NKPTS
           DO J = 1,NBANDS
              OCC(J,ISPPT) = 1.0d0 
           end do
        end do
        return
     end if
  end if
  
  pi=acos(0.d0)*2.d0
  ee=exp(1.d0)
  eesh=sqrt(ee)*0.5d0
  sq2i=sqrt(2.0d0)*0.5d0
  piesqq=sqrt(ee*pi)*0.25d0
 
! note that this has to be changed if k-points are introduced !

  do nkp=1,nkpts
     weight(nkp)=1.d0/DBLE(nkpts)
  end do

  Z    = DBLE (NEL)
 
! COPY EIGVAL INTO SORT ARRAY.
 
  NEIG = 0
  DO  ISPPT = 1,NKPTS
     DO   J = 1, NBANDS
        NEIG = NEIG + 1
        SORT(NEIG) = EIGVAL(J,ISPPT)
     enddo
  enddo
 
!-----------------------------------------------------------------------
!  THE ARRAY IS ORDERED INTO ASCENDING ORDER OF EIGENVALUE
!-----------------------------------------------------------------------
 
  DO N=2,NKPTS*NBANDS
     EN=SORT(N)
     DO NN=N-1,1,-1
        IF (SORT(NN).LE.EN) THEN
           EXIT
        ENDIF
        SORT(NN+1)=SORT(NN)
     enddo
     SORT(NN+1)=EN
  end do
  eigmin=sort(1)
  eigmax=sort(NKPTS*nbands) 
 
!-----------------------------------------------------------------------
!  if the temperature is 0 (well, le.1d-9) then set manually the
!  Fermi energy between the HOMO and LUMO
!-----------------------------------------------------------------------

 
  if ((abs(del).le.1.d-9).and.(nspin.eq.1)) then
     if ((2*(nel/2)).ne.nel) then
        write(*,*) 'EFERMI: etemp=0.0 but nel is odd !'
        stop
     end if
     nel2=nel/2
     entropy=0.d0
     ef=0.5d0*(sort(NKPTS*nel2)+sort(NKPTS*nel2+1))
     DO ISPPT = 1,NKPTS
        DO J = 1, NBANDS
           if (eigval(J,ISPPT).le.ef) then
              occ(j,isppt)=2.d0
           else
              occ(j,isppt)=0.d0
           end if
        end do
     end do
     TEST = 0.d0
!     write(*,'(a8,f12.6)') 'Efermi: ',ef
     DO ISPPT = 1,NKPTS
        DO J = 1,NBANDS
!     write(*,'(a8,f12.6,f10.6)') 'Eigs,f: ',&
!    &        eigval(J,ISPPT),OCC(J,ISPPT)
           TEST = TEST + WEIGHT(ISPPT)*OCC(J,ISPPT)
        end do
     end do
!      this is commented since occ is normalized to 2
!      test=test*2.0
     IF ( ABS(TEST-Z) .GT. 1.0D-5) THEN
        WRITE(*,*) '*** WARNING *** OCCUPANCIES MANUALLY SET'
        DO ISPPT = 1,NKPTS
           DO J = 1, NBANDS
              if (j.le.nel2) then
                 occ(j,isppt)=2.d0
              else
                 occ(j,isppt)=0.d0
              end if
!         write(*,'(a8,f12.6,f10.6)') 'Eigs,f: ',
!    &         eigval(J,ISPPT),OCC(J,ISPPT)
           end do
        end do
     end if
     return
  else if ((abs(del).le.1.d-9).and.(nspin.ne.1)) then
     if ((2*(nel/2)).ne.nel) then
        write(*,*) 'EFERMI: etemp=0.0 but nel is odd !'
        stop
     end if
     entropy=0.d0
     ef=0.5d0*(sort(NKPTS*nel)+sort(NKPTS*nel+1))
     DO ISPPT = 1,NKPTS
        DO J = 1, NBANDS
           if (eigval(J,ISPPT).le.ef) then
              occ(j,isppt)=1.d0
           else
              occ(j,isppt)=0.d0
           end if

        end do
     end do
     TEST = 0.d0
!     write(*,'(a8,f12.6)') 'Efermi: ',ef

     DO ISPPT = 1,NKPTS
        DO J = 1,NBANDS
!     write(*,'(a8,f12.6,f10.6)') 'Eigs,f: ', &
!    &        eigval(J,ISPPT),OCC(J,ISPPT)
           TEST = TEST + WEIGHT(ISPPT)*OCC(J,ISPPT)
        end do
     end do
!      this is commented since occ is normalized to 2
!      test=test*2.0
     IF ( ABS(TEST-Z) .GT. 1.0D-5) THEN
        WRITE(*,*) '*** WARNING *** OCCUPANCIES MANUALLY SET'
        DO ISPPT = 1,NKPTS
           DO J = 1, NBANDS
              if (j.le.nel) then
                 occ(j,isppt)=1.d0
              else
                 occ(j,isppt)=0.d0
              end if
!         write(*,'(a8,f12.6,f10.6)') 'Eigs,f: ',
!    &         eigval(J,ISPPT),OCC(J,ISPPT)
           end do
        end do
     end if
     return
  end if



!-----------------------------------------------------------------------
!     THE UPPER BOUND XE2 AND THE LOWER BOUND XE1 
!     ARE PUT TO FIRST AND LAST EIGENVALUE, THEN
!     THE ACTUAL FERMI ENERGY IS FOUND BY BISECTION
!     UPPER BOUND IS ACTUALLY UPPED A BIT, JUST IN CASE
!-----------------------------------------------------------------------
 
  XE1=SORT(1)
  XE2=SORT(NKPTS*NBANDS)+del*5.d0
!     write(*,*) NEL,NBANDS,DEL,NKPTS,ismear
!     write(*,*) xe1,xe2
!
!           WRITE(*,*) ' '
  IF(ISMEAR.EQ.1) THEN
!           WRITE(*,*) 'GAUSSIAN BROADENING'
  ELSEIF(ISMEAR.EQ.2) THEN
!           WRITE(*,*) 'FERMI-DIRAC BROADENING'
  ELSEIF(ISMEAR.EQ.3) THEN
!           WRITE(*,*) 'HERMITE-DIRAC BROADENING'
  ELSEIF(ISMEAR.EQ.4) THEN
!           WRITE(*,*) 'GAUSSIAN SPLINES BROADENING'
  ELSEIF(ISMEAR.EQ.5) THEN
!           WRITE(*,*) 'COLD SMEARING I'
  ELSEIF(ISMEAR.EQ.6) THEN
!           WRITE(*,*) 'COLD SMEARING II'
  ENDIF
!
! FMID = FUNC(X2) in Numerical Recipes. 
!
  Z1=0.D0
  DO  ISPPT = 1,NKPTS
     DO  J = 1,NBANDS
        X = (XE2 - EIGVAL(J,ISPPT))/DEL
        IF(ISMEAR.EQ.1) THEN
           Z1 = Z1 + WEIGHT(ISPPT)*( 2.d0 - erfc(X) )/fspin
        ELSEIF(ISMEAR.EQ.2) THEN
           Z1 = Z1 + WEIGHT(ISPPT)*FERMID(-X)/fspin
        ELSEIF(ISMEAR.EQ.3) THEN
           Z1 = Z1 + WEIGHT(ISPPT)*DELTHM(X)/fspin
        ELSEIF(ISMEAR.EQ.4) THEN
           Z1 = Z1 + WEIGHT(ISPPT)*EFERMI_SPLINE(-X)/fspin
        ELSEIF(ISMEAR.EQ.5) THEN
           Z1 = Z1 + WEIGHT(ISPPT)*POSHM(X)/fspin
        ELSEIF(ISMEAR.EQ.6) THEN
           Z1 = Z1 + WEIGHT(ISPPT)*POSHM2(X)/fspin
        ENDIF
     END DO
  END DO 

 
  FMID= Z1-Z
!     write(*,*) fmid,z1,z
 
! F = FUNC(X1)
 
  Z1=0.D0
  DO ISPPT = 1,NKPTS
     DO J = 1,NBANDS
        X = (XE1 - EIGVAL(J,ISPPT))/DEL
        IF(ISMEAR.EQ.1) THEN
           Z1 = Z1 + WEIGHT(ISPPT)*( 2.d0 - erfc(X) )/fspin
        ELSEIF(ISMEAR.EQ.2) THEN
           Z1 = Z1 + WEIGHT(ISPPT)*FERMID(-X)/fspin
        ELSEIF(ISMEAR.EQ.3) THEN
           Z1 = Z1 + WEIGHT(ISPPT)*DELTHM(X)/fspin
        ELSEIF(ISMEAR.EQ.4) THEN
           Z1 = Z1 + WEIGHT(ISPPT)*EFERMI_SPLINE(-X)/fspin
        ELSEIF(ISMEAR.EQ.5) THEN
           Z1 = Z1 + WEIGHT(ISPPT)*POSHM(X)/fspin
        ELSEIF(ISMEAR.EQ.6) THEN
           Z1 = Z1 + WEIGHT(ISPPT)*POSHM2(X)/fspin
        ENDIF
     END DO
  END DO

 
  F= Z1-Z
!     write(*,*) f,z1,z
 
  IF(F*FMID .GE. 0.D0) THEN
     WRITE(*,*) 'WARNING: NO FERMI ENERGY INSIDE EIGENVALUES ?'
  ENDIF
  IF(F .LT. 0.D0) THEN
     RTBIS = XE1
     DX = XE2 - XE1
  ELSE
     RTBIS = XE2
     DX = XE1 - XE2
  ENDIF

  DO  J = 1, JMAX
     DX = DX * 0.5D0
     XMID = RTBIS + DX
 
! FMID=FUNC(XMID)
 
     Z1=0.D0
     DO ISPPT = 1,NKPTS
        DO  J2 = 1,NBANDS
           X = (XMID - EIGVAL(J2,ISPPT))/DEL
           IF(ISMEAR.EQ.1) THEN
              Z1 = Z1 + WEIGHT(ISPPT)*( 2.d0 - erfc(X) )/fspin
           ELSEIF(ISMEAR.EQ.2) THEN
              Z1 = Z1 + WEIGHT(ISPPT)*FERMID(-X)/fspin
           ELSEIF(ISMEAR.EQ.3) THEN
              Z1 = Z1 + WEIGHT(ISPPT)*DELTHM(X)/fspin
           ELSEIF(ISMEAR.EQ.4) THEN
              Z1 = Z1 + WEIGHT(ISPPT)*EFERMI_SPLINE(-X)/fspin
           ELSEIF(ISMEAR.EQ.5) THEN
              Z1 = Z1 + WEIGHT(ISPPT)*POSHM(X)/fspin
           ELSEIF(ISMEAR.EQ.6) THEN
              Z1 = Z1 + WEIGHT(ISPPT)*POSHM2(X)/fspin
           ENDIF
        END DO
     END DO

 
     FMID= Z1-Z
 
     IF(FMID .LE. 0.D0) RTBIS=XMID
     IF(ABS(DX) .LT. XACC .OR. FMID .EQ. 0) THEN
        EXIT
     ENDIF
  ENDDO
  IF(J >= JMAX) THEN
     WRITE(*,*) 'CANNOT BISECT FOREVER, CAN I ?'
     stop
  ENDIF
  EF = RTBIS
 

  DO ISPPT = 1, NKPTS
     DO J = 1,NBANDS
        X = ( EF-EIGVAL(J,ISPPT))/DEL
        IF(ISMEAR.EQ.1) THEN
           OCC(J,ISPPT) = 2.d0 - erfc(X)  
        ELSEIF(ISMEAR.EQ.2) THEN
           OCC(J,ISPPT) = FERMID(-X)
        ELSEIF(ISMEAR.EQ.3) THEN
           OCC(J,ISPPT) = DELTHM(X)
        ELSEIF(ISMEAR.EQ.4) THEN
           OCC(J,ISPPT) = EFERMI_SPLINE(-X)
        ELSEIF(ISMEAR.EQ.5) THEN
           OCC(J,ISPPT) = POSHM(X)
        ELSEIF(ISMEAR.EQ.6) THEN
           OCC(J,ISPPT) = POSHM2(X)
        ENDIF
! occupations are normalized to two or one depending on nspin
        OCC(J,ISPPT) = OCC(J,ISPPT)/fspin
     ENDDO
  ENDDO


!-------------------------------------------------------------
! CALCULATES THE CORRECTION TERM TO GET "0 TEMPERATURE" ENERGY
!-------------------------------------------------------------

  DELCOR=0.0D0
  DO ISPPT = 1, NKPTS
     DO J = 1,NBANDS
        X = ( EF-EIGVAL(J,ISPPT))/DEL
        IF(ISMEAR.EQ.1) THEN
           DELCOR=DELCOR &
     &    -DEL*WEIGHT(ISPPT)*EXP(-X*X)/(2.D0*SQRT(pi))
        ELSEIF(ISMEAR.EQ.2) THEN
           FI=FERMID(-X)/entrospin
           IF(ABS(FI) .GT. 1.d-12) THEN
              IF(ABS(FI-1.D0) .GT. 1.d-12) THEN
                 DELCOR=DELCOR+DEL*WEIGHT(ISPPT)* &
                      &     (FI*LOG(FI)+(1.D0-FI)*LOG(1.D0-FI))
              ENDIF
           ENDIF
        ELSEIF(ISMEAR.EQ.3) THEN
           DELCOR=DELCOR+DEL/2.0d0*WEIGHT(ISPPT)  &
     &    *(2.0d0*x*x-1.d0)*exp(-x*x)/(2.0d0*sqrt(pi))
        ELSEIF(ISMEAR.EQ.4) THEN
           x=abs(x)
           zeta=eesh*abs(x)*exp(-(x+sq2i)**2)+piesqq*erfc(x+sq2i)
           delcor=delcor-del*WEIGHT(ISPPT)*zeta
        ELSEIF(ISMEAR.EQ.5) THEN
           a=-0.5634d0
!        a=-0.8165
           DELCOR=DELCOR-DEL/2.d0*WEIGHT(ISPPT) &
! NOTE g's are all intended to be normalized to 1 !
! this following line is -2*int_minf^x [t*g(t)]dt
     &   *(2.d0*a*x**3-2.d0*x*x+1 )*exp(-x*x)/(2.d0*sqrt(pi))
        ELSEIF(ISMEAR.EQ.6) THEN
         DELCOR=DELCOR-DEL/2.d0*WEIGHT(ISPPT) &
! NOTE g's are all intended to be normalized to 1 !
! this following line is -2*int_minf^x [t*g(t)]dt
     &   *(1.d0-sqrt(2.d0)*x)*exp(-(x-1.d0/sqrt(2.d0))**2)/sqrt(pi)
        ENDIF
     END DO
  END DO


!--------------------------------------------------------
!  the correction is also stored in sort, for compatibility,
!  and -TS is stored in entropy
!--------------------------------------------------------

  sort(1)=delcor
  entropy=entrospin*delcor
 
!--------------------------------------------------------
!     TEST WHETHER OCCUPANCY ADDS UP TO Z
!--------------------------------------------------------

  TEST = 0.d0
!     write(*,'(a8,f12.6)') 'Efermi: ',ef
  DO ISPPT = 1,NKPTS
     DO J = 1,NBANDS
!       write(*,'(a8,f12.6,f10.6)') 'Eigs,f: ',
!    &       eigval(J,ISPPT),OCC(J,ISPPT)
        TEST = TEST + WEIGHT(ISPPT)*OCC(J,ISPPT)
     end do
  end do

  IF ( ABS(TEST-Z) .GT. 1.0D-5) THEN
     WRITE(*,*) '*** WARNING ***'
     WRITE(*,220) TEST,NEL
220     FORMAT(' SUM OF OCCUPANCIES =',F30.20 ,' BUT NEL =',I5)
!      ELSE
!
!230     FORMAT(' TOTAL CHARGE = ',F15.8)
  ENDIF
!
!     TEST WHETHER THE MATERIAL IS A SEMICONDUCTOR
!
  IF ( MOD( NEL, 2) .EQ. 1) RETURN
  INEL = NEL/2
  ELOW = EIGVAL(INEL+1,1)
  DO  ISPPT = 2,NKPTS
     ELOW =MIN( ELOW, EIGVAL(INEL+1,ISPPT))
  ENDDO

  DO  ISPPT = 1,NKPTS
     IF (ELOW .LT. EIGVAL(INEL,ISPPT)) RETURN
  END DO

  if (NKPTS.gt.1) then
     WRITE (*,*) 'MATERIAL MAY BE A SEMICONDUCTOR'
  end if
!
  RETURN
END SUBROUTINE efermi
!-----------------------------------------------------------------------
FUNCTION fermid(xx)
  
  USE kinds, ONLY : DP

  implicit none

  REAL(kind=DP), INTENT(in) :: xx
  
  REAL(kind=DP) :: fermid

  IF(XX .GT. 30.D0) THEN
     FERMID=0.D0
  ELSEIF(XX .LT. -30.D0) THEN
     FERMID=2.D0
  ELSE
     FERMID=2.D0/(1.D0+EXP(XX))
  ENDIF
!
  RETURN
END FUNCTION fermid
!-----------------------------------------------------------------------
FUNCTION delthm(xx)
!
  USE kinds, ONLY : DP

  implicit none

  REAL(kind=DP) :: delthm
  REAL(kind=DP), INTENT(in) :: xx

  REAL(kind=DP) :: pi

  pi=3.14159265358979d0
  IF(XX .GT. 10.D0) THEN
     DELTHM=2.D0
  ELSEIF(XX .LT. -10.D0) THEN
     DELTHM=0.D0
  ELSE
     DELTHM=(2.D0-erfc(XX))+XX*EXP(-XX*XX)/SQRT(PI)
  ENDIF
!
  RETURN
END FUNCTION delthm
!-----------------------------------------------------------------------
FUNCTION efermi_spline(x)

  USE kinds,  ONLY : DP
  
  implicit none

  REAL(kind=DP) :: efermi_spline
  REAL(kind=DP), INTENT(in) :: x

  REAL(kind=DP) :: eesqh,sq2i,fx

  eesqh=sqrt(exp(1.d0))*0.5d0
  sq2i=sqrt(2.d0)*0.5d0
  if (x.ge.0.d0) then
     fx=eesqh*exp(-(x+sq2i)**2)
  else
     fx=1.d0-eesqh*exp(-(x-sq2i)**2)
  endif
  efermi_spline=2.d0*fx
!
  return
END FUNCTION efermi_spline
!-----------------------------------------------------------------------
FUNCTION poshm(x)
!
! NOTE g's are all intended to be normalized to 1 !
! function = 2 * int_minf^x [g(t)] dt
!

  USE kinds,    ONLY : DP

  implicit none

  REAL(kind=DP) :: poshm
  REAL(kind=DP), INTENT(in) :: x

  REAL(kind=DP) :: pi,a

  pi=3.141592653589793238d0
  a=-0.5634d0
  !     a=-0.8165
  IF(X .GT. 10.D0) THEN
     POSHM=2.D0
  ELSEIF(X .LT. -10.D0) THEN
     POSHM=0.D0
  ELSE
     POSHM=(2.D0-erfc(X))+(-2.d0*a*x*x+2*x+a)*EXP(-X*X)/SQRT(PI)/2.d0
  ENDIF
!
  RETURN
END FUNCTION poshm
!-----------------------------------------------------------------------
FUNCTION poshm2(x)
!
! NOTE g's are all intended to be normalized to 1 !
! function = 2 * int_minf^x [g(t)] dt
!

  USE kinds,   ONLY : DP

  implicit none

  REAL(kind=DP) :: poshm2
  REAL(kind=DP), INTENT(in) :: x

  REAL(kind=DP) :: pi

  pi=3.141592653589793238d0
  IF(X .GT. 10.D0) THEN
     POSHM2=2.D0
  ELSEIF(X .LT. -10.D0) THEN
     POSHM2=0.D0
  ELSE
     POSHM2=(2.D0-erfc(X-1.d0/sqrt(2.d0)))+ &
          &  sqrt(2.d0)*exp(-x*x+sqrt(2.d0)*x-0.5d0)/sqrt(pi)
  ENDIF
!
  RETURN
END FUNCTION poshm2
!-----------------------------------------------------------------------
