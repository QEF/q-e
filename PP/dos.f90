!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
program dos_e

  character :: nodenumber * 3
  call start_postproc (nodenumber)
  call dos (nodenumber)

  call stop_pp
  stop
end program dos_e
!
!--------------------------------------------------------------------
subroutine dos (nodenumber)
  !--------------------------------------------------------------------
  !
  ! Input (namelist &inputpp ... &end):                Default value
  !
  !    prefix        prefix of input file produced by pw.x  'pwscf'
  !                    (wavefunctions are not needed)
  !    tmp_dir       directory where filpun is searched   ./
  !    ngauss1,      gaussian broadening parameters       0
  !     degauss1       if absent, read from file          0.d0
  !    Emin, Emax    min, max energy (eV) for DOS plot    band extrema
  !    DeltaE        energy grid step (eV)                none
  !    fildos        output file containing DOS(E)        dos.out
  !
  ! In order to use tetrahedron method, the input file must
  ! contain information on tetrahedra (option ltetra=.true.)
  ! if degauss1, ngauss1 are specified they override what is
  ! specified in the input file wrt summation method
  !
  ! Parallel case implemented only for one processor!
  !
  use pwcom
  use io

  character(len=3) :: nodenumber
  character(len=14) :: fildos
  real(kind=DP) :: E, DOSofE (2), DOSint, Elw, Eup, DeltaE, Emin, Emax, &
       degauss1
  integer :: nrot, ik, n, ndos, ngauss1
  namelist /inputpp/ tmp_dir, prefix, fildos, degauss1,ngauss1,&
       Emin, Emax, DeltaE
  logical :: minus_q
  !
  nd_nmbr = nodenumber
  !
  !   set default values for variables in namelist
  !
  tmp_dir='./'
  prefix ='pwscf'
  fildos =' '
  Emin   =-1000000.
  Emax   = 1000000.
  DeltaE = 0.01
  ngauss1= 0
  degauss1=0.d0
  !
  read (5, inputpp, err=200, iostat=ios )
200 call errore('dos','reading inputpp namelist',abs(ios))
  !
  call read_file
  !
  if (nks.ne.nkstot) call errore ('dos', 'not implemented', 1)
  !
  if (degauss1.ne.0.d0) then
     degauss=degauss1
     ngauss =ngauss1
     write(6,'(/5x,"Gaussian broadening (read from input): ",&
          &        "ngauss,degauss=",i4,f12.6/)') ngauss,degauss
     ltetra=.false.
     lgauss=.true.
  else if (ltetra) then
     write(6,'(/5x,"Tetrahedra used"/)')
  else if (lgauss) then
     write(6,'(/5x,"Gaussian broadening (read from file): ",&
          &        "ngauss,degauss=",i4,f12.6/)') ngauss,degauss
  else
     if (degauss1.eq.0.d0) call errore('dos','I need a gaussian broadening!',1)
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
  Emin=max(Emin/13.6058,Elw)
  Emax=min(Emax/13.6058,Eup)
  DeltaE = DeltaE / 13.6058
  ndos = nint ( (Emax - Emin) / DeltaE+0.500001)
  DOSint = 0.0
  !
  open (unit = 4, file = fildos, status = 'unknown', form = 'formatted')
  if (nspin.eq.1) then
     write(4,'("#  E (eV)   dos(E)     Int dos(E)")')
  else
     write(4,'("#  E (eV)   dosup(E)     dosdw(E)   Int dos(E)")')
  end if
  do n= 1, ndos
     E = Emin + (n - 1) * DeltaE
     if (ltetra) then
        call dos_t(et,nspin,nbndx,nbnd, nks,ntetra,tetra, E, DOSofE)
     else
        call dos_g(et,nspin,nbndx,nbnd, nks,wk,degauss,ngauss, E, DOSofE)
     endif
     if (nspin.eq.1) then
        DOSint = DOSint + DOSofE (1) * DeltaE
        write (4, '(f7.3,2e12.4)') E * 13.6058, DOSofE(1)/13.6058, DOSint
     else
        DOSint = DOSint + (DOSofE (1) + DOSofE (2) ) * DeltaE
        write (4, '(f7.3,3e12.4)') E * 13.6058, DOSofE/13.6058, DOSint
     endif
  enddo

  close (unit = 4)
  return
end subroutine dos
