!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
program dos
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
  !    ngauss        Type pf gaussian broadening (optional)    0
  !            =  0  Simple Gaussian (default)
  !            =  1  Methfessel-Paxton of order 1
  !            = -1  Marzari-Vanderbilt "cold smearing"
  !            = 99  Fermi-Dirac function
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
#include "f_defs.h"
  USE io_global,  ONLY : stdout, ionode_id
  USE io_files,   ONLY : nd_nmbr, prefix, tmp_dir
  USE constants,  ONLY : rytoev
  USE kinds, ONLY : DP
  USE klist,      ONLY : xk, wk, degauss, ngauss, lgauss, nks, nkstot
  USE ktetra,     ONLY : ntetra, tetra, ltetra
  USE wvfct,      ONLY : nbnd, et
  USE lsda_mod,   ONLY : nspin
#ifdef __PARA
  use para,       ONLY : me
  use mp,         ONLY : mp_bcast
#endif
  implicit none
  character(len=256) :: fildos
  character(len=256) :: outdir
  real(kind=DP) :: E, DOSofE (2), DOSint, Elw, Eup, DeltaE, Emin, Emax, &
       degauss1
  integer :: nrot, ik, n, ndos, ngauss1, ios
  namelist /inputpp/ outdir, prefix, fildos, degauss, ngauss, &
       Emin, Emax, DeltaE
  logical :: minus_q
                                                    
  CHARACTER (LEN=256) :: input_file
  INTEGER             :: nargs, iiarg, ierr, ilen
  INTEGER, EXTERNAL   :: iargc
  !
  call start_postproc (nd_nmbr)
#ifdef __PARA
  if (me == 1) then
#endif
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
  !
  ! ... Input from file ?
  !
  nargs = iargc()
  !
  DO iiarg = 1, ( nargs - 1 )
     !
     CALL getarg( iiarg, input_file )
     IF ( TRIM( input_file ) == '-input' .OR. &
          TRIM( input_file ) == '-inp'   .OR. &
          TRIM( input_file ) == '-in' ) THEN
        !
        CALL getarg( ( iiarg + 1 ) , input_file )
        OPEN ( UNIT = 5, FILE = input_file, FORM = 'FORMATTED', &
               STATUS = 'OLD', IOSTAT = ierr )
        CALL errore( 'iosys', 'input file ' // TRIM( input_file ) // &
                   & ' not found' , ierr )
        !
     END IF
     !
  END DO

  read (5, inputpp, err=200, iostat=ios )
200 call errore('dos','reading inputpp namelist',abs(ios))
  !
  tmp_dir = trim(outdir)
  ! save the value of degauss and ngauss: they are read from file
  degauss1 = degauss
  ngauss1  = ngauss
#ifdef __PARA
  end if
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( tmp_dir, ionode_id )
  CALL mp_bcast( prefix, ionode_id )
#endif
  !
  call read_file( )
  !
#ifdef __PARA
  if (me == 1) then
#endif
  if (nks.ne.nkstot) call errore ('dos', 'not implemented', 1)
  !
  if (degauss1.ne.0.d0) then
     degauss=degauss1
     ngauss =ngauss1
     WRITE( stdout,'(/5x,"Gaussian broadening (read from input): ",&
          &        "ngauss,degauss=",i4,f12.6/)') ngauss,degauss
     ltetra=.false.
     lgauss=.true.
  else if (ltetra) then
     WRITE( stdout,'(/5x,"Tetrahedra used"/)')
  else if (lgauss) then
     WRITE( stdout,'(/5x,"Gaussian broadening (read from file): ",&
          &        "ngauss,degauss=",i4,f12.6/)') ngauss,degauss
  else
     degauss=DeltaE/rytoev
     ngauss =0
     WRITE( stdout,'(/5x,"Gaussian broadening (default values): ",&
          &        "ngauss,degauss=",i4,f12.6/)') ngauss,degauss
     ltetra=.false.
     lgauss=.true.
  end if
  !
  ! find band extrema
  !
  Elw = et (1, 1)
  Eup = et (nbnd, 1)
  do ik = 2, nks
     Elw = min (Elw, et (1, ik) )
     Eup = max (Eup, et (nbnd, ik) )
  enddo
  if (degauss.ne.0.d0) then
     Eup = Eup + 3d0 * degauss
     Elw = Elw - 3d0 * degauss
  endif
  !
  Emin=max(Emin/rytoev,Elw)
  Emax=min(Emax/rytoev,Eup)
  DeltaE = DeltaE / rytoev
  ndos = nint ( (Emax - Emin) / DeltaE+0.500001)
  DOSint = 0.0
  !
  if ( fildos == ' ' ) fildos = trim(prefix)//'.dos'
  open (unit = 4, file = fildos, status = 'unknown', form = 'formatted')
  if (nspin.eq.1) then
     write(4,'("#  E (eV)   dos(E)     Int dos(E)")')
  else
     write(4,'("#  E (eV)   dosup(E)     dosdw(E)   Int dos(E)")')
  end if
  do n= 1, ndos
     E = Emin + (n - 1) * DeltaE
     if (ltetra) then
        call dos_t(et,nspin,nbnd, nks,ntetra,tetra, E, DOSofE)
     else
        call dos_g(et,nspin,nbnd, nks,wk,degauss,ngauss, E, DOSofE)
     endif
     if (nspin.eq.1) then
        DOSint = DOSint + DOSofE (1) * DeltaE
        write (4, '(f7.3,2e12.4)') E * rytoev, DOSofE(1)/rytoev, DOSint
     else
        DOSint = DOSint + (DOSofE (1) + DOSofE (2) ) * DeltaE
        write (4, '(f7.3,3e12.4)') E * rytoev, DOSofE/rytoev, DOSint
     endif
  enddo

  close (unit = 4)
#ifdef __PARA
  end if
#endif
  call stop_pp
end program dos

