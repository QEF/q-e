!
! Copyright (C) 2002-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!  AB INITIO COSTANT PRESSURE MOLECULAR DYNAMICS
!  ----------------------------------------------
!  Car-Parrinello Parallel Program
!  Carlo Cavazzoni - Gerardo Ballabio
!  SISSA, Trieste, Italy - 1997-99
!  Last modified: Fri Dec  3 10:42:00 MET 1999
!  ----------------------------------------------
!  routines in this file:
!  SUBROUTINE adjef(nk,e,wk,wke,fke,ef,qtot,ne,temp,sume,nspin)
!  REAL(DP) FUNCTION stepf(x)
!  ----------------------------------------------
!  BEGIN manual

      SUBROUTINE adjef(nk,e,wk,wke,fke,ef,qtot,ne,temp,sume,nspin)

!  this routine computes Fermi energy and weights of occupied states
!  using an improved Gaussian-smearing method
!  refs: C.L.Fu and K.M.Ho, Phys.Rev. B28, 5480 (1983)
!        M.Methfessel and A.T.Paxton Phys.Rev. B40 (15 aug. 89).
!
!  taken from APW code by J. Soler and A. Williams (jk+ss)
!  added computation of occupation numbers without k-point weight
!  ----------------------------------------------
!  END manual

      USE kinds
      USE io_global, ONLY: stdout

      IMPLICIT NONE

! ... declare subroutine arguments
      INTEGER ne,nk,nspin
      REAL(DP) ef,qtot,temp,sume
      REAL(DP) e(ne,nk,nspin),wke(ne,nk,nspin)
      REAL(DP) wk(nk),fke(ne,nk,nspin)
      REAL(DP), PARAMETER  :: tol = 1.d-10
      INTEGER, PARAMETER :: nitmax = 100

! ... declare functions
      REAL(DP) stepf

! ... declare other variables
      REAL(DP) sumq,emin,emax,fac,t,drange
      INTEGER ik,ispin,ie,iter

!  end of declarations
!  ----------------------------------------------

!     qtot=DBLE(nel)
      sumq=0.d0
      sume=0.d0
      emin=e(1,1,1)
      emax=e(1,1,1)
      fac=2.d0
      IF (nspin.EQ.2) fac=1.d0

      DO ik=1,nk
        DO ispin=1,nspin
          DO ie=1,ne
            wke(ie,ik,ispin)=wk(ik)*fac
            fke(ie,ik,ispin)=fac
            sumq=sumq+wke(ie,ik,ispin)
            sume=sume+wke(ie,ik,ispin)*e(ie,ik,ispin)
            emin= MIN (emin,e(ie,ik,ispin))
            emax= MAX (emax,e(ie,ik,ispin))
          END DO
        END DO
      END DO
      ef=emax
      IF (dabs(sumq-qtot).LT.tol) RETURN
      IF (sumq.LT.qtot) THEN
        WRITE( stdout,*) 'FERMIE: NOT ENOUGH STATES'
        WRITE( stdout,*) 'FERMIE: QTOT,SUMQ=',qtot,sumq
        STOP
      END IF

      t=DMAX1(temp,1.d-6)
      drange=t*dsqrt(-dlog(tol*.01d0))
      emin=emin-drange
      emax=emax+drange
      DO iter=1,nitmax
        ef=0.5d0*(emin+emax)
        sumq=0.d0
        sume=0.d0
        DO ik=1,nk
          DO ispin=1,nspin 
            DO ie=1,ne
              wke(ie,ik,ispin)=fac/2.d0*                &
                  wk(ik)*stepf((e(ie,ik,ispin)-ef)/t)
              fke(ie,ik,ispin)=fac/2.d0*                &
                  stepf((e(ie,ik,ispin)-ef)/t)
              sumq=sumq+wke(ie,ik,ispin)
              sume=sume+wke(ie,ik,ispin)*e(ie,ik,ispin)
            END DO
          END DO
        END DO
        IF (dabs(sumq-qtot).LT.tol) RETURN
        IF (sumq.LE.qtot) emin=ef
        IF (sumq.GE.qtot) emax=ef
      END DO

      WRITE( stdout,*) 'FERMIE: ITERATION HAS NOT CONVERGED.'
      WRITE( stdout,*) 'FERMIE: QTOT,SUMQ=',qtot,sumq
      STOP

      END SUBROUTINE adjef

!  ----------------------------------------------

      DOUBLE PRECISION FUNCTION stepf(x)
      USE kinds
      IMPLICIT NONE
      REAL(DP) :: x
      REAL(DP), EXTERNAL :: erfc
      REAL(DP), PARAMETER :: c=0.5641895835D0
!     stepf=erfc(x)
      stepf=1.d0/(exp(min(x,100.d0))+1.d0)
      END FUNCTION stepf


      SUBROUTINE adjef_s(e,fke,ef,nel,nx,temp,sume)

! e(nstati)
! fke(nstati) = f(nstati)       (output)
! ef = fermi energy             (output)
! nel = n. electrons
! ne = nstati
! temp = broadening (au)
! sume = sum e(nstati)          (output)


!  CALCULATES FERMI ENERGY AND WEIGHTS OF OCCUPIED STATES USING
!  AN IMPROVED GAUSSIAN-SMEARING METHOD
!  REFS: C.L.FU AND K.M.HO, PHYS.REV. B28, 5480 (1983)
!        M.METHFESSEL AND A.T.PAXTON PHYS.REV. B40 (15 AUG. 89).
!
!  Taken from APW code by J. Soler and A. Williams (jk+ss)
!  Added computation of occupation numbers without k-point weight

      use kinds
      USE io_global, ONLY: stdout
      IMPLICIT NONE

      integer nx,nel
      real(DP) E(nx),FKE(nx),temp,sume,ef,tol
      integer nitmax
      PARAMETER (TOL=1.D-10,NITMAX=100)
      integer iter,ie
      real(DP) t,emin,emax,stepf
      real(DP) sumq,fac,qtot,drange
      QTOT=DBLE(NEL)
      SUMQ=0.D0
      SUME=0.D0
      EMIN=E(1)
      EMAX=E(1)
      fac=2.d0
      do ie=1,nx
        FKE(IE)=fac
        SUMQ=SUMQ+FKE(IE)
        SUME=SUME+E(IE)
        EMIN=MIN(EMIN,E(IE))
        EMAX=MAX(EMAX,E(IE))
      end do
      EF=EMAX
      IF (DABS(SUMQ-QTOT).LT.TOL) RETURN
      IF (SUMQ.LT.QTOT) THEN
         WRITE( stdout,*) 'FERMIE: NOT ENOUGH STATES'
         WRITE( stdout,*) 'FERMIE: QTOT,SUMQ=',QTOT,SUMQ
         STOP
      ENDIF
      T=MAX(TEMP,1.D-6)
      DRANGE=T*DSQRT(-DLOG(TOL*.01D0))
      EMIN=EMIN-DRANGE
      EMAX=EMAX+DRANGE
      DO ITER=1,NITMAX
         EF=0.5D0*(EMIN+EMAX)
         SUMQ=0.D0
         SUME=0.D0
         do ie=1,nx
           FKE(IE)=fac*STEPF((E(IE)-EF)/T)
           SUMQ=SUMQ+FKE(IE)
           SUME=SUME+FKE(IE)*E(IE)
         enddo
         IF (DABS(SUMQ-QTOT).LT.TOL) RETURN
         IF (SUMQ.LE.QTOT) EMIN=EF
         IF (SUMQ.GE.QTOT) EMAX=EF
      ENDDO
      WRITE( stdout,*) 'FERMIE: ITERATION HAS NOT CONVERGED.'
      WRITE( stdout,*) 'FERMIE: QTOT,SUMQ=',QTOT,SUMQ
      STOP
      END SUBROUTINE adjef_s
