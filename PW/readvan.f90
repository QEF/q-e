!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
subroutine readvan (is, iunps)
  !---------------------------------------------------------------------
  !
  !     This routine reads the quantities which defines the Vanderbilt
  !     pseudopotential from the file produced by the atomic program.
  !     It is compatible only with the last versions of the atomic code
  !     6.0.0 or later.
  !     In particular the Herman Skillman mesh is no more supported.
  !     It assume multiple rinner values.
  !
  USE kinds, only: dp
  USE parameters, ONLY: nchix, lmaxx, nbrx, ndmx, npsx, nqfx, lqmax
  use atom,  only: zmesh, mesh,dx, r, rab, chi, oc, nchi, lchi, &
       rho_at, rho_atc, nlcc
  use pseud, only: zp, lmax, lloc
  use uspp_param, only: vloc_at, dion, betar, qqq, qfcoef, qfunc, nqf, nqlc, &
       rinner, nbeta, kkbeta, lll, iver, ifqopt, psd, tvanp
  use funct, only: dft, which_dft, ismeta, ishybrid
  implicit none

  !
  !    First the arguments passed to the subroutine
  !
  integer :: is, iunps
  ! The number of the pseudopotential
  ! The unit of the pseudo file
  !
  !   The local variables which are used to read the Vanderbilt file.
  !   They are used only for the pseudopotential report. They are no
  !   longer used in the rest of the code.
  !

  real(kind=DP) :: exfact, etotpseu, wwnl (nchix), ee (nchix), eloc, &
       dummy, rc (lmaxx + 1), eee (nbrx), ddd (nbrx, nbrx), rcloc, ru(ndmx)
  ! index of the exchange and correlation use
  ! total pseudopotential energy
  ! the occupation of the valence states
  ! the energy of the valence states
  ! energy of the local potential
  ! dummy real variable
  ! the cut-off radii of the pseudopotential
  ! energies of the beta function
  ! the screened D_{\mu,\nu} parameters
  ! the cut-off radius of the local potential
  ! the total charge density

  integer :: ios, idmy (3), i, nnlz (nchix), keyps, ifpcor, irel, &
       iptype (nbrx), npf, lp, l, mb, nb, ir
  ! integer variable for I/O control
  ! Contains the date of creation of the pseu
  ! Dummy counter
  ! The nlm values of the valence states
  ! the type of pseudopotential. Only US allo
  ! if = 1 the pseudopotential has nlcc
  ! it says if the pseudopotential is relativ
  ! more recent parameters
  ! as above
  ! counter on Q angular momenta
  ! counter on angular momenta
  ! beta function counter
  ! beta function counter
  ! mesh points counter
  character (len=20) :: line, xctit
  ! The title of the pseudopotential
  ! Name of the xctype
  !
  !     We first check the input variables
  !
  if (is.lt.0.or.is.gt.npsx) call errore ('readvan', 'Wrong is number &', 1)

  read (iunps, '(6i5)', err = 100, iostat = ios) (iver (i, is) , i = 1, 3) ,&
       (idmy (i) , i = 1, 3)
  if (iver (1, is) .lt.6) then
     call errore ('readvan', 'This version is too old', iver (1, is) )
  elseif (iver (1, is) .gt.7) then
     call errore ('readvan', 'This version is too new', iver (1, is) )

  endif

  read (iunps, '(a20,3f15.9)', err = 100, iostat = ios) line, zmesh &
       (is) , zp (is) , exfact
  if (exfact.eq.0) then
     dft = 'PZ'
  elseif (exfact.eq.1) then
     dft = 'BLYP'
  elseif (exfact.eq. - 5.or.exfact.eq.3) then
     dft = 'BP'
  elseif (exfact.eq. - 6.or.exfact.eq.4) then
     dft = 'PW91'
  elseif (exfact.eq.5) then
     dft = 'PBE'
  else
     call errore ('readvan', 'Wrong xc in pseudopotential', 1)
  endif
  !
  call which_dft( dft )
  !
  IF ( ismeta ) &
    CALL errore( 'readvan ', 'META-GGA not implemented in PWscf', 1 )
#if defined (EXX)
#else
  IF ( ishybrid ) &
    CALL errore( 'readvan ', 'HYBRID XC not implemented in PWscf', 1 )
#endif

  !
  read (iunps, '(2i5,1pe19.11)', err = 100, iostat = ios) nchi (is) &
       , mesh (is) , etotpseu
  !
  !    Test the input values
  !
  if (nchi (is) .gt.nchix) call errore ('readvan', 'nchi> nchix', &
       nchi (is) )
  if (nchi (is) .lt.0) call errore ('readvan', 'wrong nchi ', is)
  if (mesh (is) > ndmx .or. mesh (is) <= 0) call errore ('readvan', &
       'wrong mesh', is)
  if (zp (is) .le.0.d0) call errore ('readvan', 'wrong zp', is)
  !
  !   Set the pseudopotential name
  !
  psd (is) = line (1:2)
  !
  !    The next 3 parameters are not passed to the rest.
  !    No control on them
  !
  read (iunps, '(i5,2f15.9)', err = 100, iostat = ios) (nnlz (nb) , &
       wwnl (nb) , ee (nb) , nb = 1, nchi (is) )
  read (iunps, '(2i5,f15.9)', err = 100, iostat = ios) keyps, &
       ifpcor, dummy
  if (keyps.ne.3) call errore ('readvan', 'keyps .ne. 3', keyps)
  !
  !   If this atom has nlcc set the appropriate variables
  !
  if (ifpcor.eq.1) then
     nlcc (is) = .true.
  else
     nlcc (is) = .false.
  endif
  !
  !     Read information on the angular momenta, and on Q pseudization
  !
  read (iunps, '(2i5,f9.5,2i5,f9.5)', err = 100, iostat = ios) lmax &
       (is) , lloc (is) , eloc, ifqopt (is) , nqf (is) , dummy
  !
  !    NB: In the Vanderbilt atomic code the angular momenta goes
  !        from 1 to lmax+1. In our code from 0 to lmax. In this routine
  !        we use the Vanderbilt convention up to the end. Then we change
  !        to interface with the rest of the code.
  !
  if (lmax (is) .gt.lmaxx + 1.or.lmax (is) .le.0) call errore (' read &
       &van', 'Wrong lmax', lmax (is) )
  if (lloc (is) .eq. - 1) lloc (is) = lmax (is) + 1
  if (lloc (is) .gt.lmax (is) + 1.or.lloc (is) .lt.0) call errore ( &
       'readvan', 'wrong lloc', is)
  if (nqf (is) > nqfx .or. nqf (is) < 0) call errore (' readvan', &
       'Wrong nqf', nqf (is) )
  if (ifqopt (is) .lt.0) call errore ('readvan', 'wrong ifqopt', is)
  !
  !     Reads and test the values of rinner
  !
  read (iunps, *, err = 100, iostat = ios) (rinner (lp, is), &
       lp = 1, lmax (is) * 2 - 1)
  do lp = 1, lmax (is) * 2 - 1
     if (rinner (lp, is) .lt.0.d0) call errore ('readvan', 'Wrong rinner', is)
  enddo
  read (iunps, '(i5)', err = 100, iostat = ios) irel
  !
  !       set the number of angular momentum terms in q_ij to read in
  !

  nqlc (is) = 2 * lmax (is) - 1

  if (nqlc (is) .gt.lqmax.or.nqlc (is) .lt.0) call errore (' readvan', &
       &'Wrong  nqlc', nqlc (is) )
  read (iunps, '(1p4e19.11)', err = 100, iostat = ios) (rc (l) , l = &
       1, lmax (is) )
  !
  !     reads the number of beta functions
  !

  read (iunps, '(2i5)', err = 100, iostat = ios) nbeta (is) , &
       kkbeta (is)
  if (nbeta (is) .gt.nbrx.or.nbeta (is) .lt.0) call errore ( &
       'readvan', 'wrong nbeta', is)
  if (kkbeta (is) .gt.mesh (is) .or.kkbeta (is) .lt.0) call errore ( &
       'readvan', 'wrong kkbeta', is)
  !
  !    Now reads the main Vanderbilt parameters
  !
  do nb = 1, nbeta (is)
     read (iunps, '(i5)', err = 100, iostat = ios) lll (nb, is)
     read (iunps, '(1p4e19.11)', err = 100, iostat = ios) eee (nb) , &
          (betar (ir, nb, is) , ir = 1, kkbeta (is) )
     if (lll (nb, is) .gt.lmaxx.or.lll (nb, is) .lt.0)&
          call errore ('readvan', ' wrong lll ', is)
     do mb = nb, nbeta (is)
        read (iunps, '(1p4e19.11)', err = 100, iostat = ios) &
             dion (nb, mb, is) , ddd (nb, mb) , qqq (nb, mb, is), &
             (qfunc (ir, nb, mb, is), ir = 1, kkbeta (is) ) , &
             ( (qfcoef (i, lp, nb, mb, is) , i = 1, nqf(is) ) , lp = 1, nqlc(is) )
        !
        !     Use the symmetry of the coefficients
        !
        dion (mb, nb, is) = dion (nb, mb, is)

        qqq (mb, nb, is) = qqq (nb, mb, is)
        do ir = 1, kkbeta (is)
           qfunc (ir, mb, nb, is) = qfunc (ir, nb, mb, is)

        enddo
        do i = 1, nqf (is)
           do lp = 1, nqlc (is)
              qfcoef (i, lp, mb, nb, is) = qfcoef (i, lp, nb, mb, is)
           enddo
        enddo
     enddo
  enddo
  !
  !    for versions later than 7.2
  !
  if (10 * iver (1, is) + iver (2, is) .ge.72) then
     read (iunps, '(6i5)', err = 100, iostat = ios) (iptype (nb) , &
          nb = 1, nbeta (is) )
     read (iunps, '(i5,f15.9)', err = 100, iostat = ios) npf, dummy
  endif
  !
  !   reads the local potential of the vanderbilt
  !
  read (iunps, '(1p4e19.11)', err = 100, iostat = ios) rcloc, &
       (vloc_at (ir, is) , ir = 1, mesh (is) )
  !
  !   If present reads the core charge
  !
  if (nlcc (is) ) then
     if (iver (1, is) .ge.7) read (iunps, '(1p4e19.11)', err = 100, &
          iostat = ios) dummy
     read (iunps, '(1p4e19.11)', err = 100, iostat = ios) &
          (rho_atc ( ir, is) , ir = 1, mesh (is) )
  endif
  !
  !     Reads the total atomic charge and forgets
  !
  read (iunps, '(1p4e19.11)', err = 100, iostat = ios) (ru (ir) , &
       ir = 1, mesh (is) )
  !
  !     Reads the valence atomic charge
  !
  read (iunps, '(1p4e19.11)', err = 100, iostat = ios) (rho_at (ir,is) , &
       ir = 1, mesh (is) )
  !
  !     Reads the mesh
  !
  read (iunps, '(1p4e19.11)', err = 100, iostat = ios) (r (ir, is) , &
       ir = 1, mesh (is) )
  read (iunps, '(1p4e19.11)', err = 100, iostat = ios) (rab (ir, is), &
       ir = 1, mesh (is) )
  !
  !    For compatibility with rho_atc in the non-US case
  !
  if (nlcc (is) ) then
     do ir = 2, mesh (is)
        rho_atc (ir, is) = rho_atc (ir, is) / 4.0 / 3.14159265 / r(ir,is)**2
     enddo
     rho_atc (1, is) = 0.d0
  endif
  !
  !     Put the local potential in the variables of the code, with the same
  !     units
  !
  lloc (is) = 0
  do ir = 2, mesh (is)
     vloc_at (ir, is) = vloc_at (ir, is) / r (ir, is)
  enddo
  vloc_at (1, is) = vloc_at (2, is)
  !
  !    Set lmax in the range 0-lmax
  !
  lmax (is) = lmax (is) - 1
  !
  !    Reads the wavefunctions of the atom
  !
  if (iver (1, is) .ge.7) then
     read (iunps, *, err = 100, iostat = ios) i
     if (i.ne.nchi (is) ) call errore ('readvan', &
          & 'unexpected or unimplemented case', 1)

  endif
  read (iunps, *, err = 100, iostat = ios) ( (chi (ir, nb, is), &
       ir = 1, mesh (is) ), nb = 1, nchi (is) )
  do nb = 1, nchi (is)
     i = nnlz (nb) / 100
     lchi (nb, is) = nnlz (nb) / 10 - i * 10
  enddo
100 call errore ('readvan', 'error reading pseudo file', abs (ios) )
  !
  !    Here we write on output information on the pseudopotential
  !
  if (exfact.eq.0.) xctit = '      ceperley-alder'
  if (exfact.eq. - 1.) xctit = '              wigner'
  if (exfact.eq. - 2.) xctit = '     hedin-lundqvist'
  if (exfact.eq. - 3.) xctit = ' gunnarson-lundqvist'

  if (exfact.gt.0.) xctit = '      slater x-alpha'
  !      WRITE( stdout,200) is
200 format (/4x,60('=')/4x,'|  pseudopotential report', &
       &        ' for atomic species:',i3,11x,'|')
  !      WRITE( stdout,300) 'pseudo potential version', iver(1,is),
  !     +     iver(2,is), iver(3,is)
300 format (4x,'|  ',1a30,3i4,13x,' |' /4x,60('-'))
  !      WRITE( stdout,400) line, xctit
400 format (4x,'|  ',2a20,' exchange-corr  |')
  !      WRITE( stdout,500) zmesh(is), is, zp(is), exfact
500 format (4x,'|  z =',f5.0,4x,'zv(',i2,') =',f5.0,4x,'exfact =', &
       &     f10.5, 9x,'|')
  !      WRITE( stdout,600) ifpcor, etotpseu
600 format (4x,'|  ifpcor = ',i2,10x,' atomic energy =',f10.5, &
       &       ' Ry',6x,'|')
  !      WRITE( stdout,700)
700 format(4x,'|  index    orbital      occupation    energy',14x,'|')
  !      WRITE( stdout,800) ( nb, nnlz(nb), wwnl(nb), ee(nb), nb=1,nchi(is) )
800 format(4x,'|',i5,i11,5x,f10.2,f12.2,15x,'|')
900 format('(4x,"|  rinner =",',i1,'f8.4,',i2,'x,"|")')
  !      WRITE( stdout,1000)
1000 format(4x,'|    new generation scheme:',32x,'|')
  !      WRITE( stdout,1100) nbeta(is),kkbeta(is),rcloc
1100 format(4x,'|    nbeta = ',i2,5x,'kkbeta =',i5,5x, &
       &     'rcloc =',f10.4,4x,'|'/ &
       &     4x,'|    ibeta    l     epsilon   rcut',25x,'|')
  do nb = 1, nbeta (is)
     lp = lll (nb, is) + 1
     !         WRITE( stdout,1200) nb,lll(nb,is),eee(nb),rc(lp)
1200 format    (4x,'|',5x,i2,6x,i2,4x,2f7.2,25x,'|')
  enddo
  !      WRITE( stdout,1300)

1300 format (4x,60('='))
  return
end subroutine readvan

