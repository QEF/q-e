!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine d3_summary
  !-----------------------------------------------------------------------
  !
  !    This routine writes on output the quantities which have been read
  !    from the punch file, and the quantities computed in the d3_setup
  !    file.
  !
  !    if iverbosity = 0 only a partial summary is done.
  !
  USE kinds, only : DP
  USE constants,  ONLY : amu_ry
  USE ions_base,  ONLY : nat, ityp, ntyp => nsp, atm, tau, amass
  USE run_info, ONLY : title
  USE io_global,  ONLY : stdout
  USE symm_base,   ONLY : s, sr, sname, ftau
  USE control_flags, ONLY : iverbosity
  USE fft_base, ONLY : dffts, dfftp
  USE gvecw, ONLY : ecutwfc
  use pwcom
  use qpoint, ONLY: xq
  use phcom
  use d3com
!
  implicit none
  integer :: i, l, nt, mu, nu, ipol, apol, na, isymq, isym, nsymtot, &
       ik, ib, irr, imode0, iaux
  ! generic counter
  ! counter on angular momenta
  ! counter on atomic types
  ! counter on modes
  ! counter on modes
  ! counter on polarizations
  ! counter on polarizations
  ! counter on atoms
  ! counter on symmetries
  ! counter on symmetries
  ! counter on symmetries
  ! counter on k points
  ! counter on beta functions
  ! counter on irreducible representation
  ! the first mode

  real (DP) :: ft1, ft2, ft3, xkg (3)
  ! fractionary translation
  ! k point in crystal coordinates

  WRITE( stdout, 100) title, ibrav, alat, omega, nat, ntyp, &
       ecutwfc, ecutwfc * dual
100 format (/,5x,a75,/,/, &
       &     'bravais-lattice index     = ',i12,/,5x, &
       &     'lattice parameter (a_0)   = ',f12.4,'  a.u.',/,5x, &
       &     'unit-cell volume          = ',f12.4,' (a.u.)^3',/,5x, &
       &     'number of atoms/cell      = ',i12,/,5x, &
       &     'number of atomic types    = ',i12,/,5x, &
       &     'kinetic-energy cut-off    = ',f12.4,'  Ry',/,5x, &
       &     'charge density cut-off    = ',f12.4,'  Ry',/,5x,/)
  !
  !    and here more detailed information. Description of the unit cell
  !
  WRITE( stdout, '(2(3x,3(2x,"celldm(",i1,")=",f11.5),/))') (i, &
       celldm (i) , i = 1, 6)
  WRITE( stdout, '(5x, &
       &  "crystal axes: (cart. coord. in units of a_0)",/, &
       &         3(15x,"a(",i1,") = (",3f8.4," )  ",/ ) )')  (apol, &
       & (at (ipol, apol) , ipol = 1, 3) , apol = 1, 3)
  WRITE( stdout, '(5x, &
       &"reciprocal axes: (cart. coord. in units 2 pi/a_0)",/, &
       &         3(15x,"b(",i1,") = (",3f8.4," )  ",/ ) )')  (apol, &
       & (bg (ipol, apol) , ipol = 1, 3) , apol = 1, 3)
  !
  !    description of the atoms inside the unit cell
  !
  WRITE( stdout, '(/, 5x,"Atoms inside the unit cell: ")')
  WRITE( stdout, '(/,3x,"Cartesian axes")')
  WRITE( stdout, '(/,5x,"site n.  atom      mass ", &
       &                "          positions (a_0 units)")')

  WRITE( stdout, '(7x,i2,5x,a6,f8.4,"   tau(",i2, ") = (",3f11.5,"  )")') &
       (na, atm (ityp (na) ) , amass (ityp (na) ) / amu_ry, na, &
       (tau (ipol, na ) , ipol = 1, 3) , na = 1, nat)
  WRITE( stdout, '(/,5x,"Computing dynamical matrix for ")')
  WRITE( stdout, '(20x,"q = (",3f11.5," )")') (xq (ipol) , ipol = 1, 3)
  if (q0mode_todo (1) .le.0) then
     WRITE( stdout, '(/,5x,"Computing all the modes ")')
  else
     WRITE( stdout, '(/,5x,"Computing only selected modes: ")')
     do i = 1, 3 * nat
        if (q0mode (i) ) WRITE( stdout, '(5x,"Mode to be computed: ",i5)') i
     enddo
  endif
  !
  !   description of symmetries
  !
  WRITE( stdout, * )
  if (nsymg0.le.1) then
     WRITE( stdout, '(5x,"No symmetry! for q=0 ")')
  else
     WRITE( stdout, '(5x,i2," + 1 = ",i2," q=0 Sym.Ops. ",/)') &
          nsymg0, nsymg0 + 1
  endif
  if (.not.lgamma) then
     WRITE( stdout, * )
     if (nsymq.le.1.and..not.minus_q) then
        WRITE( stdout, '(5x,"No symmetry!")')
     else
        if (minus_q) then
           WRITE( stdout, '(5x,i2," Sym.Ops. (with q -> -q+G )",/)') &
                nsymq + 1
        else
           WRITE( stdout, '(5x,i2," Sym.Ops. (no q -> -q+G )",/)') &
                nsymq
        endif
     endif
  endif
  if (iverbosity.eq.1) then

     WRITE( stdout, '(36x,"s",24x,"frac. trans.")')
     if (minus_q) then
        iaux = 0
     else
        iaux = 1
     endif
     do isymq = iaux, nsymg0
        if (isymq.eq.0) then
           isym = irotmq
           WRITE( stdout, '(/,5x,"This transformation sends q -> -q+G")')
        else
           !
           ! It seems to me variable irgq is useless !
           !               isym = irgq(isymq)
           isym = isymq
        endif
        if (isymq.eq.nsymq + 1) then
           WRITE( stdout, '(//,5x,&
                &"In the following are listed symmetries of the crystal")')
           WRITE( stdout, '(5x,"not belonging to the small group of q")')
        endif
        WRITE( stdout, '(/6x,"isym = ",i2,5x,a45/)') isymq, sname (isym)
        if (ftau (1, isym) .ne.0.or.ftau (2, isym) .ne.0.or.ftau (3, &
             isym) .ne.0) then
           ft1 = at (1, 1) * ftau (1, isym) / dfftp%nr1 + at (1, 2) * ftau ( &
                2, isym) / dfftp%nr2 + at (1, 3) * ftau (3, isym) / dfftp%nr3
           ft2 = at (2, 1) * ftau (1, isym) / dfftp%nr1 + at (2, 2) * ftau ( &
                2, isym) / dfftp%nr2 + at (2, 3) * ftau (3, isym) / dfftp%nr3
           ft3 = at (3, 1) * ftau (1, isym) / dfftp%nr1 + at (3, 2) * ftau ( &
                2, isym) / dfftp%nr2 + at (3, 3) * ftau (3, isym) / dfftp%nr3
           WRITE( stdout, '(1x,"cryst.",3x,"s(",i2,") = (",3(i6,5x), &
                &" )    f =( ",f10.7," )")') isymq,  (s (1, ipol, isym),&
                ipol = 1, 3) , DBLE (ftau (1, isym) )  / DBLE (dfftp%nr1)
           WRITE( stdout, '(17x," (",3(i6,5x), &
                &                    " )       ( ",f10.7," )")') &
                (s (2, ipol, &
                &isym) , ipol = 1, 3) , DBLE (ftau (2, isym) )  / DBLE (dfftp%nr2)
           WRITE( stdout, '(17x," (",3(i6,5x), &
                &                    " )       ( ",f10.7," )"/)')  (s (3, ipol, &
                & isym) , ipol = 1, 3) , DBLE (ftau (3, isym) )  / DBLE (dfftp%nr3)
           WRITE( stdout,'(1x,"cart.",3x,"s(",i2,") = (",3f11.7, &
                &         " )    f =( ", f10.7," )")') &
                          isymq, (sr (1,ipol,isym), ipol=1,3), ft1
           WRITE( stdout, '(17x," (",3f11.7, " )       ( ",f10.7," )")') &
                           (sr (2,ipol,isym) , ipol = 1, 3) , ft2
           WRITE( stdout, '(17x," (",3f11.7, " )       ( ",f10.7," )"/)')&
                           (sr (3,ipol,isym) , ipol = 1, 3) , ft3
        else
           WRITE( stdout, '(1x,"cryst.",3x,"s(",i2,") = (",3(i6,5x), " )")') &
                isymq,  (s (1, ipol, isym) , ipol = 1, 3)
           WRITE( stdout, '(17x," (",3(i6,5x)," )")') &
                        (s (2, ipol, isym), ipol = 1, 3)
           WRITE( stdout, '(17x," (",3(i6,5x)," )"/)') &
                        (s (3, ipol, isym) , ipol = 1, 3)
           WRITE( stdout, '(1x,"cart.",3x,"s(",i2,") = (",3f11.7, " )")') &
                isymq,  (sr (1, ipol, isym) , ipol = 1, 3)
           WRITE( stdout, '(17x," (",3f11.7," )")') &
                        (sr (2, ipol, isym) , ipol = 1, 3)
           WRITE( stdout, '(17x," (",3f11.7," )"/)') &
                        (sr (3, ipol, isym) , ipol = 1, 3)
        endif
     enddo
  endif
  !
  !     Description of the reciprocal lattice vectors
  !
  WRITE( stdout, '(/5x,"G cutoff =",f10.4,"  (", &
       &       i7," G-vectors)","     FFT grid: (",i3, &
       &       ",",i3,",",i3,")")') gcutm, ngm, dfftp%nr1, dfftp%nr2, dfftp%nr3

  if (doublegrid) WRITE( stdout, '(5x,"G cutoff =",f10.4,"  (", &
       &                      i7," G-vectors)","  smooth grid: (",i3, &
       &                      ",",i3,",",i3,")")') gcutms, ngms,  &
       &dffts%nr1, dffts%nr2, dffts%nr3
  if (degauss.eq.0.d0) then
     WRITE( stdout, '(5x,"number of k points=",i5)') nkstot
  else
     WRITE( stdout, '(5x,"number of k points=",i5, &
          &               "  gaussian broad. (Ry)=",f8.4,5x, &
          &               "ngauss = ",i3)') nkstot, degauss, ngauss
  endif
  WRITE( stdout, '(23x,"cart. coord. in units 2pi/a_0")')
  do ik = 1, nkstot
     WRITE( stdout, '(8x,"k(",i5,") = (",3f12.7,"), wk =",f12.7)') ik, &
          (xk (ipol, ik) , ipol = 1, 3) , wk (ik)
  enddo
  if (iverbosity.eq.1) then
     WRITE( stdout, '(/23x,"cryst. coord.")')
     do ik = 1, nkstot
        do ipol = 1, 3
           ! xkg are the compone
           xkg (ipol) = at (1, ipol) * xk (1, ik) + at (2, ipol) * xk (2, &
                ik) + at (3, ipol) * xk (3, ik)
           ! of xk in the crysta
           ! rec. lattice basis
        enddo
        WRITE( stdout, '(8x,"k(",i5,") = (",3f12.7,"), wk =",f12.7)') &
             ik, (xkg (ipol) , ipol = 1, 3) , wk (ik)
     enddo

  endif
  !
  CALL print_ps_info ( )
  !
  ! Representation for q=0
  !
  if (.not.lgamma) then
     WRITE( stdout, '(//5x,"Atomic displacements (q=0 Repr):")')
     WRITE( stdout, '(5x,"There are ",i5, &
          &    " irreducible representations")') nirrg0
     imode0 = 0
     do irr = 1, nirrg0
        WRITE( stdout, '(/, 5x,"Representation ",i5,i7, &
             &             " modes - To be done")') irr, npertg0 (irr)
        if (iverbosity.eq.1) then

           WRITE( stdout, '(5x,"Phonon polarizations are as follows:",/)')
           if (npertg0 (irr) .eq.1) then
              WRITE( stdout, '(20x," mode # ",i3)') imode0 + 1
              WRITE( stdout, '(20x," (",2f10.5,"   ) ")')  ( (ug0 (mu, nu) , nu = &
                   & imode0 + 1, imode0 + npertg0 (irr) ) , mu = 1, 3 * nat)
           elseif (npertg0 (irr) .eq.2) then
              WRITE( stdout, '(2(10x," mode # ",i3,16x))') imode0 + 1, &
                   imode0 + 2
              WRITE( stdout, '(2(10x," (",2f10.5,"   ) "))')  ( (ug0 (mu, nu),&
                   nu = imode0 + 1, imode0 + npertg0 (irr) ) , mu = 1, 3 * nat)
           else
              WRITE( stdout, '(4x,3(" mode # ",i3,13x))') imode0 + 1, &
                   imode0 + 2, imode0 + 3
              WRITE( stdout, '((5x,3("(",2f10.5," ) ")))') ( (ug0 (mu, &
                   nu) , nu = imode0 + 1, imode0 + npertg0 (irr) ) , mu = 1, &
                   3 * nat)
           endif
           imode0 = imode0 + npertg0 (irr)
        endif
     enddo
  endif
  !
  ! Representation for a generic q
  !
  WRITE( stdout, '(//5x,"Atomic displacements:")')
  WRITE( stdout, '(5x,"There are ",i5," irreducible representations") &
       &') nirr
  imode0 = 0
  do irr = 1, nirr
     WRITE( stdout, '(/, 5x,"Representation ",i5,i7, &
          &                " modes - To be done")') irr, npert (irr)
     if (iverbosity.eq.1) then

        WRITE( stdout, '(5x,"Phonon polarizations are as follows:",/)')
        if (npert (irr) .eq.1) then
           WRITE( stdout, '(20x," mode # ",i3)') imode0 + 1
           WRITE( stdout, '(20x," (",2f10.5,"   ) ")')  ( (u (mu, nu) , nu = &
                imode0 + 1, imode0 + npert (irr) ) , mu = 1, 3 * nat)
        elseif (npert (irr) .eq.2) then
           WRITE( stdout, '(2(10x," mode # ",i3,16x))') imode0 + 1, &
                imode0 + 2
           WRITE( stdout, '(2(10x," (",2f10.5,"   ) "))')  ( (u (mu, nu) , &
                nu = imode0 + 1, imode0 + npert (irr) ) , mu = 1, 3 * nat)
        else
           WRITE( stdout, '(4x,3(" mode # ",i3,13x))') imode0 + 1, imode0 &
                + 2, imode0 + 3
           WRITE( stdout, '((5x,3("(",2f10.5," ) ")))') ( (u (mu, nu) , &
                nu = imode0 + 1, imode0 + npert (irr) ) , mu = 1, 3 * nat)
        endif
        imode0 = imode0 + npert (irr)
     endif

  enddo
  WRITE( stdout, '(/20x,"**   Complex  Version    **")')
  !
  FLUSH( stdout )
  !
  return
end subroutine d3_summary
