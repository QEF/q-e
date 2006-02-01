!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!--------------------------------------------------------------------
PROGRAM dos
  !--------------------------------------------------------------------
  !
  ! Calculates the Density of States (DOS),
  ! separated into up and down components for LSDA
  !
  ! Input (namelist &inputpp ... &end):                   Default value
  !
  !    prefix        prefix of input file produced by pw.x    'pwscf'
  !                    (wavefunctions are not needed)
  !    outdir        directory containing the input file       ./
  !    ngauss        Type of gaussian broadening (optional)    0
  !            =  0  Simple Gaussian (default)
  !            =  1  Methfessel-Paxton of order 1
  !            = -1  Marzari-Vanderbilt "cold smearing"
  !            =-99  Fermi-Dirac function
  !    degauss       gaussian broadening, Ry (not eV!)          see below
  !    Emin, Emax    min, max energy (eV) for DOS plot          band extrema
  !    DeltaE        energy grid step (eV)                      none
  !    fildos        output file containing DOS(E)              "prefix".dos
  !
  ! Output:
  !
  !   The total DOS (states/eV plotted vs E in eV) is written to file "fildos"
  !
  ! Important notice:
  !    The tetrahedron method is used if
  !    - the input data file has been produced by pw.x using the option
  !      occupations='tetrahedra', AND
  !    - a value for degauss is not given as input to namelist &inputpp
  !    Gaussian broadening is used in all other cases:
  !    - if degauss is set to some value in namelist &inputpp, that value
  !      (and the optional value for ngauss) is used
  !    - if degauss is NOT set to any value in namelist &inputpp, the 
  !      value of degauss and of ngauss are read from the input data
  !      file (they will be the same used in the pw.x calculations)
  !    - if degauss is NOT set to any value in namelist &inputpp, AND
  !      there is no value of degauss and of ngauss in the input data
  !      file, degauss=DeltaE (in Ry) and ngauss=0 will be used
  !
  USE io_global,  ONLY : stdout, ionode, ionode_id
  USE io_files,   ONLY : nd_nmbr, prefix, tmp_dir, trimcheck
  USE constants,  ONLY : rytoev
  USE kinds,      ONLY : DP
  USE klist,      ONLY : xk, wk, degauss, ngauss, lgauss, nks, nkstot
  USE ktetra,     ONLY : ntetra, tetra, ltetra
  USE wvfct,      ONLY : nbnd, et
  USE lsda_mod,   ONLY : nspin
  USE noncollin_module, ONLY: noncolin
  USE mp,         ONLY : mp_bcast
  !
  IMPLICIT NONE
  CHARACTER(len=256) :: fildos, outdir
  REAL(DP) :: E, DOSofE (2), DOSint, Elw, Eup, DeltaE, Emin, Emax, &
                   degauss1
  INTEGER :: nrot, ik, n, ndos, ngauss1, ios
  NAMELIST /inputpp/ outdir, prefix, fildos, degauss, ngauss, &
       Emin, Emax, DeltaE
  LOGICAL :: minus_q
                                                    
  !
  CALL start_postproc (nd_nmbr)
  !
  IF ( ionode ) THEN
     !
     !   set default values for variables in namelist
     !
     outdir='./'
     prefix ='pwscf'
     fildos =' '
     Emin   =-1000000.
     Emax   = 1000000.
     DeltaE = 0.01
     ngauss = 0
     degauss= 0.d0
     !
     CALL input_from_file ( )
     !
     READ (5, inputpp, err=200, iostat=ios )
200  CALL errore('dos','reading inputpp namelist',ABS(ios))
     !
     tmp_dir = trimcheck (outdir)
     ! save the value of degauss and ngauss: they are read from file
     degauss1 = degauss
     ngauss1  = ngauss
     !
  END IF
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( tmp_dir, ionode_id )
  CALL mp_bcast( prefix, ionode_id )
  !
  CALL read_file( )
  !
  IF ( ionode ) THEN
     !
     IF (nks.NE.nkstot) CALL errore ('dos', 'not implemented', 1)
     !
     IF (degauss1.NE.0.d0) THEN
        degauss=degauss1
        ngauss =ngauss1
        WRITE( stdout,'(/5x,"Gaussian broadening (read from input): ",&
             &        "ngauss,degauss=",i4,f12.6/)') ngauss,degauss
        ltetra=.FALSE.
        lgauss=.TRUE.
     ELSE IF (ltetra) THEN
        WRITE( stdout,'(/5x,"Tetrahedra used"/)')
     ELSE IF (lgauss) THEN
        WRITE( stdout,'(/5x,"Gaussian broadening (read from file): ",&
             &        "ngauss,degauss=",i4,f12.6/)') ngauss,degauss
     ELSE
        degauss=DeltaE/rytoev
        ngauss =0
        WRITE( stdout,'(/5x,"Gaussian broadening (default values): ",&
             &        "ngauss,degauss=",i4,f12.6/)') ngauss,degauss
        ltetra=.FALSE.
        lgauss=.TRUE.
     END IF
     !
     ! find band extrema
     !
     Elw = et (1, 1)
     Eup = et (nbnd, 1)
     DO ik = 2, nks
        Elw = MIN (Elw, et (1, ik) )
        Eup = MAX (Eup, et (nbnd, ik) )
     ENDDO
     IF (degauss.NE.0.d0) THEN
        Eup = Eup + 3.d0 * degauss
        Elw = Elw - 3.d0 * degauss
     ENDIF
     !
     Emin=MAX(Emin/rytoev,Elw)
     Emax=MIN(Emax/rytoev,Eup)
     DeltaE = DeltaE / rytoev
     ndos = NINT ( (Emax - Emin) / DeltaE+0.500001)
     DOSint = 0.d0
     !
     IF ( fildos == ' ' ) fildos = TRIM(prefix)//'.dos'
     OPEN (unit = 4, file = fildos, status = 'unknown', form = 'formatted')
     IF (nspin.EQ.1.OR.nspin.EQ.4) THEN
        WRITE(4,'("#  E (eV)   dos(E)     Int dos(E)")')
     ELSE
        WRITE(4,'("#  E (eV)   dosup(E)     dosdw(E)   Int dos(E)")')
     END IF
     DO n= 1, ndos
        E = Emin + (n - 1) * DeltaE
        IF (ltetra) THEN
           CALL dos_t(et,nspin,nbnd, nks,ntetra,tetra, E, DOSofE)
        ELSE
           CALL dos_g(et,nspin,nbnd, nks,wk,degauss,ngauss, E, DOSofE)
        ENDIF
        IF (nspin.EQ.1.OR.nspin.EQ.4) THEN
           DOSint = DOSint + DOSofE (1) * DeltaE
           WRITE (4, '(f7.3,2e12.4)') E * rytoev, DOSofE(1)/rytoev, DOSint
        ELSE
           DOSint = DOSint + (DOSofE (1) + DOSofE (2) ) * DeltaE
           WRITE (4, '(f7.3,3e12.4)') E * rytoev, DOSofE/rytoev, DOSint
        ENDIF
     ENDDO

     CLOSE (unit = 4)
     !
  END IF
  !
  CALL stop_pp
  !
END PROGRAM dos

