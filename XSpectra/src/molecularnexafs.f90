!
! Copyright (C) 2013-2016 Guido Fratesi
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! MOLECULARNEXAFS
!
! This code combines the output of different codes (pw.x, projwfc.x,
! and mainly xspectra.x) to evaluate the NEXAFS spectrum of a molecule
! containing inequivalent atoms belonging to the same species.
!
! The theoretical framework adopted is explained in the following
! publication, which should be cited in papers using this tool
! (additionally to relevant papers for the NEXAFS calculation):
!
! * Azimuthal Dichroism in Near-Edge X‑ray Absorption Fine Structure
!   Spectra of Planar Molecules
!   Guido Fratesi, Valeria Lanzilotto, Luca Floreano, and Gian Paolo Brivio
!   J. Phys. Chem. C 2013, 117, 6632−6638
!   http://dx.doi.org/10.1021/jp312569q
!
! Further examples are also described in:
!   Phys. Chem. Chem. Phys., 2014, 16, 14834
!   http://dx.doi.org/10.1039/c4cp01625d
!
! The spectrum for each inequivalent atom of the same species sums up
! into a total spectrum, where the contributions from each atom (as
! evaluated by xspectra.x) have to be properly aligned before the
! summation. For example, a molecule like pyridine contains C atoms in
! three different sites, C1, C2, and C3. This example is adopted below.
!
! For each inequivalent atom in the molecule, the code requires:
!
! 1) A full-core-hole (FCH) calculation of the total energy by pw.x
!    (by a FCH pseudopotential at C_i). The code will return the
!    core-level shifts (CLS) as a by-product.
!
! 2) Half-core-hole (HCH) calculations of the NEXAFS spectrum by
!    xspectra.x for the three directions X-Y-Z.  Those will be read as
!    prefix(iat)//suffix(idir)
!
! 3) Since xspectra.x energy scales may not refer to the vacuum level
!    (more commonly the Fermi level or the middle of the HOMO-LUMO gap
!    may have been used), one also need to input the reference energy
!    level of the HCH xspectra.x calculation, and the corresponding
!    vacuum level.
!
!    The vacuum level might be estimated by the Makov-Payne correction
!    directly in the HCH pw.x SCF calculation if the unit cell is
!    cubic, or by running additional calculations accessing the vacuum
!    level, then by evaluating the vacuum-to-HOMO energy difference in
!    this case, summed to the HOMO of the original calculation.
!
!    For adsorbed molecules, one should refer to the Fermi level both
!    for NEXAFS and the XPS, as is generally done in the experiments.
!    So we should set both evacuum=0 and efermi=0, and the new zero
!    will be the average core level binding energy with respect to the
!    Fermi level.
!
! 4) Eventually, energies from xspectra.x output files will be shifted
!    (E->E+eshift) so that the new zero of energy is the weighted
!    average of core-level binding energies:
!
!       eshift(iat)= -(evacuum(iat)-efermi(iat)) + dcorebe(iat)
!
!    where the first term sets the vacuum level at zero and the second
!    one adds to the transition energy the core level shift.
!
!    Recall that by the pseudopotential approach binding energies are
!    defined only up to a constant.
!
! 5) Optionally, the PDOS of p symmetry can be taken as an approximation
!    to the spectrum for a K edge: the comparison is useful for checking
!    the results.
!
!=======================================================================
! NAMELIST / CONTROL /	
!-----------------------------------------------------------------------
!
!    syslabel     character(len=80)     DEFAULT='MOLECULE'
!    System label for output files
!
!    doxps     logical     DEFAULT=.TRUE.
!    Plot XPS spectrum
!
!    donexafs     logical     DEFAULT=.TRUE.
!    Compute and plot NEXAFS spectrum
!
!    nat     integer     DEFAULT=0
!    The number of inequivalent atoms in the molecule
!
!    atweight(1:nat)     real(DP)     DEFAULT=1, 1, . . .
!    The multiplicities of the inequivalent atoms in the molecule
!
!-----------------------------------------------------------------------
!NAMELIST / XPS /	
!-----------------------------------------------------------------------
!
!    erangexps(1:2)     real(DP)     DEFAULT=(-5:5)
!    Energy range for plotting
!
!    nptxps     integer     DEFAULT=501
!    Number of points for plotting
!
!    delorentz, degauss, lorentzratio     real(DP)     DEFAULT=0.2, 0.2, 0.5
!    Pseudovoigt parameters.
!    The spectrum is plotted as:
!    lorentzratio * lorentzian + (1-lorentzratio)*gaussian
!    delorentz=HWHM of the Lorentzian
!    degauss=STDDEV of the Gaussian
!
!    etotfch(1:nat)     real(DP)     DEFAULT=0
!    Total energy (Ry) with a FCH in that given atom
!
!-----------------------------------------------------------------------
!NAMELIST / NEXAFS /	
!-----------------------------------------------------------------------
!
!    dosingleatoms     logical     DEFAULT=.TRUE.
!    Set to true to also output the results of each atom in a separate file
!
!    dopdosp     logical     DEFAULT=.FALSE.
!    Set to true to also read and sum the PDOS (p) of each atom, with
!    analogous procedure as for the NEXAFS spectrum.
!
!    erangenexafs(1:2)     real(DP)     DEFAULT=(-10:20)
!    Energy range for plotting
!
!    nptnexafs     integer     DEFAULT=3001
!    Number of points for plotting
!
!    efermi(1:nat)     real(DP)     DEFAULT=0
!    Reference energy level of the HCH xspectra.x calculation. Should
!    be left to zero in case of molecules adsorbed on metals (assuming
!    the Fermi level was taken as a reference for xspectra.x).
!
!    evacuum(1:nat)     real(DP)     DEFAULT=0
!    Vacuum level of the HCH calculation. Should be left to zero in
!    case of molecules adsorbed on metals.
!
!    prefix(1:nat)     character(len=80)     DEFAULT='xanes.dat.001', 'xanes.dat.002', . . .
!    The files containing the spectrum computed by xspectra.x should
!    follow this naming convention:
!    prefix(iat)//suffix(nu), where
!    Iat=1,nat identifies the inequivalent atom
!    nu=1,2,3 the three directions x,y,z
!
!    suffix(1:3)     character(len=80)     DEFAULT='.xspectraX.dat', '.xspectraY.dat', '.xspectraZ.dat'
!    Direction-dependent suffix of xspectra.x spectrum file
!
!    atlabel(0:nat)     character(len=80)     DEFAULT='SUM', '001', '002', . . .
!    Atomic labels (the 0-th element is for the summed spectrum).
!
!    pdospfile(1:nat)     character(len=80)     DEFAULT='pdosp.dat.001', 'pdosp.dat.002', . . .
!    The files containing the p-pdos from projwfc.x.
!
!=======================================================================
!
PROGRAM nexafsanalysis
  IMPLICIT NONE
  !
  INTEGER, PARAMETER  :: DP=SELECTED_REAL_KIND(14,200)
  REAL(DP), PARAMETER :: PI=3.14159265358979323846_DP 
  REAL(DP), PARAMETER :: RY2EV=13.605698066_DP
  !
  INTEGER, PARAMETER  :: MAXAT=20
  CHARACTER(LEN=6), PARAMETER :: FMTSINGLEAT="(i3.3)"
  !
  ! CONTROL namelist
  CHARACTER(LEN=80) :: syslabel  ! system label for output files
  LOGICAL  :: doxps  ! plot XPS-CLS?
  LOGICAL  :: donexafs  ! compute and plot NEXAFS?
  INTEGER  :: nat  ! number of inequivalent atoms
  REAL(DP) :: atweight(MAXAT)  ! atom weights (multiplicities)
  !
  ! Variables for XPS (CLS)
  ! First the CLS namelist
  REAL(DP) :: erangexps(2)  ! emin emax for plotting
  INTEGER  :: nptxps  ! number of points for plotting
  REAL(DP) :: delorentz, degauss, lorentzratio  ! pseudo-voigt parameters
  REAL(DP) :: etotfch(MAXAT) ! total energy in Ry
  ! Others
  REAL(DP) :: dcorebe(0:MAXAT) ! CLS in eV
  !
  ! Variables for NEXAFS
  ! First the NEXAFS namelist
  LOGICAL  :: dosingleatoms  ! plot each atom individually?
  LOGICAL  :: dopdosp  ! align (p) PDOS as if it were the NEXAFS spectrum
  REAL(DP) :: erangenexafs(2)  ! emin emax for plotting
  INTEGER  :: nptnexafs  ! number of points for plotting
  REAL(DP) :: efermi(MAXAT)  ! Fermi level used in xspectra.x
  REAL(DP) :: evacuum(MAXAT)  ! Vacuum level
  CHARACTER (LEN=80) :: prefix(MAXAT)  ! atom-dependent prefix of xspectra.x file
  CHARACTER (LEN=80) :: suffix(3)  ! direction-dependent suffix of xspectra.x file
  CHARACTER (LEN=80) :: atlabel(0:MAXAT)  ! atom-dependent label
  CHARACTER (LEN=80) :: pdospfile(MAXAT)  ! atom-dependent file with p-pdos
  !
  ! I/O units
  INTEGER, PARAMETER :: STDOUT=6, STDIN=5, STDERR=0
  INTEGER, PARAMETER :: INPUNIT=100
  INTEGER, PARAMETER :: OUTUNIT=200
  !
  ! Local variables
  INTEGER :: iat
  !
  NAMELIST /CONTROL/ syslabel, doxps, donexafs, nat, atweight
  NAMELIST /XPS/ erangexps, nptxps, delorentz, degauss, lorentzratio, etotfch
  NAMELIST /NEXAFS/ dosingleatoms, erangenexafs, nptnexafs, efermi, evacuum, &
       &            prefix, suffix, atlabel, dopdosp, pdospfile
  !
  ! Some defaults
  ! CONTROL
  syslabel='MOLECULE'
  doxps=.true.
  donexafs=.true.
  nat=0
  atweight(:)=1._DP
  ! XPS
  erangexps=(/ -5 , 5 /)
  nptxps=501
  delorentz=0.2_DP
  degauss=0.2_DP
  lorentzratio=0.5_DP
  etotfch(:)=0._DP
  ! NEXAFS
  dosingleatoms=.true.
  dopdosp=.false.
  erangenexafs=(/ -10 , 20 /)
  nptnexafs=3001
  efermi(:)=0._DP
  evacuum(:)=0._DP
  suffix=(/ '.xspectraX.dat', '.xspectraY.dat', '.xspectraZ.dat' /)
  atlabel(0)="SUM"
  DO iat=1,MAXAT
     WRITE (atlabel(iat),FMTSINGLEAT) iat
     prefix(iat)='xanes.dat.'//atlabel(iat)
     pdospfile(iat)='pdosp.dat.'//atlabel(iat)
  END DO
  !
  ! Read input namelists
  READ (STDIN,NML=CONTROL)
  READ (STDIN,NML=XPS)
  READ (STDIN,NML=NEXAFS)
  !
  IF (nat>MAXAT) STOP "Too many atoms: change MAXAT and recompile."
  !
  ! Write the actually used variables to file
  WRITE (STDERR,'(A)') "# NEXAFS analysis running with the following input:"
  WRITE (STDERR,NML=CONTROL)
  WRITE (STDERR,NML=XPS)
  WRITE (STDERR,NML=NEXAFS)
  !
  ! XPS to be done in any case, also if plotting not required
  CALL getxps ()
  !
  ! NEXAFS is optional
  IF (donexafs) CALL getnexafs ()
  !
  CALL byebye ()
  !
!!$============================================================================
  !
CONTAINS
  !
  SUBROUTINE getxps()
    IMPLICIT NONE
    CHARACTER(LEN=80) :: xpsout = ' '  ! output file for XPS spectrum
    REAL(DP), ALLOCATABLE :: xps(:)
    REAL(DP) :: etotav, xpstot, w
    INTEGER :: ie
    INTEGER :: iat
    !
    etotav=SUM(atweight(1:nat)*etotfch(1:nat))/SUM(atweight(1:nat))
    dcorebe(1:nat)=(etotfch(1:nat)-etotav)*RY2EV
    !dcorebe(0)=SUM(atweight(1:nat)*dcorebe(1:nat))/SUM(atweight(1:nat)) ! =0
    dcorebe(0)=0
    !
    ! Output multiplicities and energies
    WRITE (STDOUT,'(A)') "## CLSs follow"
    DO iat=1,nat
       WRITE (STDOUT,'(A,i3,f10.3,f15.7,1x,A)') "#cls ", iat, &
          REAL(atweight(iat),8)/MAXVAL(atweight(1:nat)), dcorebe(iat), &
          TRIM(atlabel(iat))
    END DO
    !
    IF (doxps) THEN
       xpsout=TRIM(syslabel)//".XPS.dat"
       OPEN (UNIT=OUTUNIT, FILE=TRIM(xpsout))
       ! Output multiplicities and energies in the output file
       DO iat=1,nat
          WRITE (OUTUNIT,'(A,i3,f10.3,f15.7,1x,A)') "#cls ", iat, &
             REAL(atweight(iat),8)/MAXVAL(atweight(1:nat)), dcorebe(iat), &
             TRIM(atlabel(iat))
       END DO
       !
       ! output spectrum
       ALLOCATE ( xps(nat) )
       DO ie=1,nptxps
          w=erangexps(1)+(erangexps(2)-erangexps(1))*(ie-1)/(nptxps-1)
          DO iat=1,nat
             xps(iat)=atweight(iat)* &
                  pseudovoigt(w-dcorebe(iat),delorentz,degauss,lorentzratio)
          END DO
          xpstot=SUM(xps(1:nat))
          WRITE (OUTUNIT,'(102f15.7)') w, xpstot, xps(1:nat)
       END DO
       DEALLOCATE (xps)
       !
       CLOSE (OUTUNIT)
    END IF
    !
  END SUBROUTINE getxps
  !
!!$============================================================================
  !
  SUBROUTINE getnexafs()
    IMPLICIT NONE
    REAL(DP) :: eshift(MAXAT) ! in eV
    INTEGER  :: iat
    !
    IF (donexafs) THEN
       !
       ! Compute the shifts
       WRITE (STDOUT,'(A)') '## eshifts follow'
       DO iat=1, nat
          eshift(iat)= -(evacuum(iat)-efermi(iat)) + dcorebe(iat)
          WRITE (STDOUT,'(A,i3,3f20.10,1x,A)') "#es ", iat, eshift(iat), evacuum(iat), efermi(iat),TRIM(atlabel(iat))
       END DO
       !
       CALL interpolate_sum_output ('NEXAFS',eshift)
       !
    END IF
    !
    IF (dopdosp) THEN
       !
       ! Compute the shifts
       !WRITE (STDOUT,'(A)') '## eshifts follow'
       DO iat=1, nat
          eshift(iat)= -(evacuum(iat)            ) + dcorebe(iat)
          !WRITE (STDOUT,'(A,i3,3f20.10,1x,A)') "#es ", iat, eshift(iat), evacuum(iat), efermi(iat),TRIM(atlabel(iat))
       END DO
       !
       CALL interpolate_sum_output ('PDOS-p',eshift)
       !
    END IF
    !
  END SUBROUTINE getnexafs
  !
  SUBROUTINE interpolate_sum_output (mode,eshift)
    IMPLICIT NONE
    CHARACTER(LEN=6) :: mode
    REAL(DP) :: eshift(MAXAT) ! in eV
    !
    CHARACTER(LEN=80) :: nexafsout = ' '  ! output file for NEXAFS spectrum
    INTEGER, PARAMETER :: MAXDATA=10000
    INTEGER, PARAMETER :: whichcol(3)=(/4,5,3/)
    REAL(DP), ALLOCATABLE :: nexafsspec(:,:,:)
    !
    REAL(DP) :: w, xdum(MAXDATA), ydum(MAXDATA), pdos(5)
    INTEGER :: ie, iat, idir, nl, jje, j
    CHARACTER (LEN=200) :: line
    CHARACTER (LEN=200) :: fileinp
    !
    IF ((TRIM(mode)/='NEXAFS').AND.(TRIM(mode)/='PDOS-p')) THEN
       STOP 'error in mode'
    END IF
    !
    ALLOCATE (nexafsspec(nptnexafs,0:3,0:nat))
    !
    ! Open and reinterpolate the files
    DO iat=1, nat
       DO idir=1, 3
          IF (TRIM(suffix(idir))=='') THEN
             nexafsspec(:,idir,iat)=0
          ELSE
             !
             ! Read data file
             IF (TRIM(mode)=='NEXAFS') THEN
                fileinp=TRIM(prefix(iat))//TRIM(suffix(idir))
             ELSE IF (TRIM(mode)=='PDOS-p') THEN
                fileinp=TRIM(pdospfile(iat))
             END IF
             OPEN(UNIT=INPUNIT,FILE=fileinp,ACTION='READ')
             !
             nl=0
             DO j=1, MAXDATA
                READ (INPUNIT, '(A)', IOSTAT=jje) line
                IF (jje/=0) EXIT
                IF (INDEX(line,"#")==0) THEN
                   nl=nl+1
                   IF (TRIM(mode)=='NEXAFS') THEN
                      READ (line,*) xdum(nl), ydum(nl)
                   ELSE IF (TRIM(mode)=='PDOS-p') THEN
                      READ (line,*) pdos(1:5)
                      xdum(nl)=pdos(1)
                      ydum(nl)=pdos(whichcol(idir))
                   END IF
                END IF
             END DO
             CLOSE (INPUNIT)
             !
             ! Apply shift to energy
             xdum(:)=xdum(:)+eshift(iat)
             !
             ! Reinterpolate
             DO ie=1,nptnexafs
                w=erangenexafs(1)+(erangenexafs(2)-erangenexafs(1))*(ie-1)/(nptnexafs-1)
                IF ((w<xdum(1)).OR.(w>xdum(nl))) THEN
                   nexafsspec(ie,idir,iat)=0
                ELSE
                   DO j=2,nl
                      IF (w<=xdum(j)) THEN
                         nexafsspec(ie,idir,iat)=ydum(j-1)+ &
                              (w-xdum(j-1)) * (ydum(j)-ydum(j-1))/(xdum(j)-xdum(j-1))
                         EXIT
                      END IF
                   END DO
                END IF
             END DO
             !
          END IF
       END DO
    END DO
    !
    ! Sum over the directions
    DO iat=1,nat
       DO ie=1,nptnexafs
          nexafsspec(ie,0,iat)=SUM(nexafsspec(ie,1:3,iat))
       END DO
    END DO
    !
    ! Sum over the atoms
    DO ie=1,nptnexafs
       DO idir=0,3
          nexafsspec(ie,idir,0)=DOT_PRODUCT(nexafsspec(ie,idir,1:nat),atweight(1:nat))
       END DO
    END DO
    !
    ! Output the spectrum
    nexafsout=TRIM(syslabel)//"."//TRIM(mode)//".dat"
    OPEN (UNIT=OUTUNIT, FILE=nexafsout)
    DO iat=0,nat
       WRITE (OUTUNIT,'(A15,4A20)',ADVANCE='NO') '#'//TRIM(atlabel(iat))//'#E-dBE(eV)', 'TOT', 'X', 'Y', 'Z'
    END DO
    WRITE (OUTUNIT,*)
    DO ie=1,nptnexafs
       w=erangenexafs(1)+(erangenexafs(2)-erangenexafs(1))*(ie-1)/(nptnexafs-1)
       WRITE (OUTUNIT,'(1f15.6,4g20.9)',ADVANCE='NO') w, nexafsspec(ie,0:3,0)
       DO iat=1,nat
          WRITE (OUTUNIT,'(1f15.6,4g20.9)',ADVANCE='NO') w-dcorebe(iat), nexafsspec(ie,0:3,iat)
       END DO
       WRITE (OUTUNIT,*)
    END DO
    CLOSE (OUTUNIT)
    !
    ! Output the spectrum for single atoms
    IF (dosingleatoms) THEN
       DO iat=0,nat
          nexafsout=TRIM(syslabel)//"."//TRIM(mode)//"."//TRIM(atlabel(iat))//".dat"
          OPEN (UNIT=OUTUNIT, FILE=nexafsout)
          WRITE (OUTUNIT,'(A15,4A20,A15)') '#'//TRIM(atlabel(iat))//'#E(eV)', 'TOT', 'X', 'Y', 'Z', 'E-dBE'
          DO ie=1,nptnexafs
             w=erangenexafs(1)+(erangenexafs(2)-erangenexafs(1))*(ie-1)/(nptnexafs-1)
             WRITE (OUTUNIT,'(1f15.6,4g20.9,1f15.6)') &
                  w, nexafsspec(ie,0:3,iat), w-dcorebe(iat)
          END DO
          CLOSE (OUTUNIT)
       END DO
    END IF
    !
    DEALLOCATE (nexafsspec)
    !
  END SUBROUTINE interpolate_sum_output
  !
  !
  !
!!$============================================================================
!!$   PLOTTING FUNCTIONS  
!!$============================================================================
  !
  FUNCTION pseudovoigt(w,l,s,x)
    REAL(DP) pseudovoigt
    REAL(DP), INTENT(in) :: w,l,s,x
    pseudovoigt=x*lorentzian(w,l)+(1-x)*gaussian(w,s)
  END FUNCTION pseudovoigt
  !
  FUNCTION lorentzian(w,l)
    REAL(DP) lorentzian
    REAL(DP), INTENT(in) :: w   ! distance from the center of the distribution
    REAL(DP), INTENT(in) :: l   ! half width at half maximum
    lorentzian=l/PI/(w**2+l**2)
  END FUNCTION lorentzian
  !
  FUNCTION gaussian(w,s)
    REAL(DP) gaussian
    REAL(DP), INTENT(in) :: w   ! distance from the center of the distribution
    REAL(DP), INTENT(in) :: s   ! standard deviation
    gaussian=1/(SQRT(2*PI)*s)*EXP(-w**2/2/s**2)
  END FUNCTION gaussian
  !
!!$============================================================================
  !
  SUBROUTINE byebye()
    !
    WRITE (STDERR,*)
    WRITE (STDERR,'(A)')" N   N EEEEE X   X   A   FFFFF  SSSS     DDDD   OOO  N   N EEEEE "
    WRITE (STDERR,'(A)')" NN  N E      X X   A A  F     S         D   D O   O NN  N E     "
    WRITE (STDERR,'(A)')" N N N EEE     X    A A  FFF    SSS      D   D O   O N N N EEE   "
    WRITE (STDERR,'(A)')" N  NN E      X X  AAAAA F         S     D   D O   O N  NN E     "
    WRITE (STDERR,'(A)')" N   N EEEEE X   X A   A F     SSSS      DDDD   OOO  N   N EEEEE "
    WRITE (STDERR,*)
    !
  END SUBROUTINE byebye
  !
END PROGRAM nexafsanalysis
