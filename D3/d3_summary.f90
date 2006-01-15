!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
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
  USE ions_base,  ONLY : nat, ityp, ntyp => nsp, atm, tau, amass
  USE io_global,  ONLY : stdout
  USE kinds, only : DP
  use pwcom
  USE uspp_param, ONLY : rinner, nqlc, nqf, lll, nbeta, psd, iver, tvanp
  USE atom, ONLY: mesh, numeric, xmin, dx, nlcc
  USE control_flags, ONLY : iverbosity
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

  real (DP) :: ft1, ft2, ft3, sr (3, 3), xkg (3)
  ! fractionary translation
  ! fractionary translation
  ! fractionary translation
  ! the symmetry matrix in cartesian coord
  ! k point in crystal coordinates

  character :: ps * 5
  ! the name of the pseudo
  WRITE( stdout, 100) title_ph, crystal, ibrav, alat, omega, nat, ntyp, &
       ecutwfc, ecutwfc * dual
100 format (/,5x,a75,/,/,5x, 'crystal is ',a20,/,/,5x, &
       &     'bravais-lattice index     = ',i12,/,5x, &
       &     'lattice parameter (a_0)   = ',f12.4,'  a.u.',/,5x, &
       &     'unit-cell volume          = ',f12.4,' (a.u.)^3',/,5x, &
       &     'number of atoms/cell      = ',i12,/,5x, &
       &     'number of atomic types    = ',i12,/,5x, &
       &     'kinetic-energy cut-off    = ',f12.4,'  Ry',/,5x, &
       &     'charge densisty cut-off   = ',f12.4,'  Ry',/,5x,/)
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
       (na, atm (ityp (na) ) , amass (ityp (na) ) / amconv, na, &
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
        call s_axis_to_cart (s (1, 1, isym), sr, at, bg)
        if (ftau (1, isym) .ne.0.or.ftau (2, isym) .ne.0.or.ftau (3, &
             isym) .ne.0) then
           ft1 = at (1, 1) * ftau (1, isym) / nr1 + at (1, 2) * ftau ( &
                2, isym) / nr2 + at (1, 3) * ftau (3, isym) / nr3
           ft2 = at (2, 1) * ftau (1, isym) / nr1 + at (2, 2) * ftau ( &
                2, isym) / nr2 + at (2, 3) * ftau (3, isym) / nr3
           ft3 = at (3, 1) * ftau (1, isym) / nr1 + at (3, 2) * ftau ( &
                2, isym) / nr2 + at (3, 3) * ftau (3, isym) / nr3
           WRITE( stdout, '(1x,"cryst.",3x,"s(",i2,") = (",3(i6,5x), &
                &" )    f =( ",f10.7," )")') isymq,  (s (1, ipol, isym),&
                ipol = 1, 3) , DBLE (ftau (1, isym) )  / DBLE (nr1)
           WRITE( stdout, '(17x," (",3(i6,5x), &
                &                    " )       ( ",f10.7," )")') &
                (s (2, ipol, &
                &isym) , ipol = 1, 3) , DBLE (ftau (2, isym) )  / DBLE (nr2)
           WRITE( stdout, '(17x," (",3(i6,5x), &
                &                    " )       ( ",f10.7," )"/)')  (s (3, ipol, &
                & isym) , ipol = 1, 3) , DBLE (ftau (3, isym) )  / DBLE (nr3)
           WRITE( stdout, '(1x,"cart.",3x,"s(",i2,") = (",3f11.7, &
                &                    " )    f =( ",f10.7," )")') isymq,  (sr (1 &
                &, ipol) , ipol = 1, 3) , ft1
           WRITE( stdout, '(17x," (",3f11.7, &
                &                    " )       ( ",f10.7," )")')  (sr (2, ipol) &
                & , ipol = 1, 3) , ft2
           WRITE( stdout, '(17x," (",3f11.7, &
                &                    " )       ( ",f10.7," )"/)')  (sr (3, ipol &
                &) , ipol = 1, 3) , ft3
        else
           WRITE( stdout, '(1x,"cryst.",3x,"s(",i2,") = (",3(i6,5x), &
                &                    " )")') isymq,  (s (1, ipol, isym) , ipol = &
                &1, 3)
           WRITE( stdout, '(17x," (",3(i6,5x)," )")') (s (2, ipol, isym) &
                , ipol = 1, 3)
           WRITE( stdout, '(17x," (",3(i6,5x)," )"/)') (s (3, ipol, &
                isym) , ipol = 1, 3)
           WRITE( stdout, '(1x,"cart.",3x,"s(",i2,") = (",3f11.7, &
                &                    " )")') isymq,  (sr (1, ipol) , ipol = 1, 3)
           WRITE( stdout, '(17x," (",3f11.7," )")') (sr (2, ipol) , &
                ipol = 1, 3)
           WRITE( stdout, '(17x," (",3f11.7," )"/)') (sr (3, ipol) , &
                ipol = 1, 3)
        endif
     enddo
  endif
  !
  !     Description of the reciprocal lattice vectors
  !
  WRITE( stdout, '(/5x,"G cutoff =",f10.4,"  (", &
       &       i7," G-vectors)","     FFT grid: (",i3, &
       &       ",",i3,",",i3,")")') gcutm, ngm, nr1, nr2, nr3

  if (doublegrid) WRITE( stdout, '(5x,"G cutoff =",f10.4,"  (", &
       &                      i7," G-vectors)","  smooth grid: (",i3, &
       &                      ",",i3,",",i3,")")') gcutms, ngms,  &
       &nr1s, nr2s, nr3s
  if (degauss.eq.0.d0) then
     WRITE( stdout, '(5x,"number of k points=",i5)') nkstot
  else
     WRITE( stdout, '(5x,"number of k points=",i5, &
          &               "  gaussian broad. (ryd)=",f8.4,5x, &
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
  do nt = 1, ntyp
     if (tvanp (nt) ) then
        ps = '(US)'
        WRITE( stdout, '(/5x,"pseudo",i2," is ",a2, &
             &        1x,a5,"   zval =",f5.1,"   lmax=",i2, &
             &        "   lloc=",i2)') nt, psd (nt) , ps, zp (nt) , lmax (nt) &
             &, lloc (nt)
        WRITE( stdout, '(5x,"Version ", 3i3, " of US pseudo code")') &
             (iver (i, nt) , i = 1, 3)
        WRITE( stdout, '(/,5x,"Using log mesh of ", i3, " points",/)') &
             mesh (nt)
        WRITE( stdout, '(5x,"The pseudopotential has ",i2, &
             &       " beta functions with: ",/)') nbeta (nt)
        do ib = 1, nbeta (nt)
           WRITE( stdout, '(15x," l(",i1,") = ",i3)') ib, lll (ib, nt)

        enddo
        WRITE( stdout, '(/,5x,"Q(r) pseudized with ", &
             &          i2," coefficients,  rinner = ",3f8.3, /, &
             &          58x,2f8.3)') nqf (nt) ,  (rinner (i, nt) , i = 1, nqlc ( &
             &nt) )
     else
        if (nlc (nt) .eq.1.and.nnl (nt) .eq.1) then
           ps = '(vbc)'
        elseif (nlc (nt) .eq.2.and.nnl (nt) .eq.3) then
           ps = '(bhs)'
        elseif (nlc (nt) .eq.1.and.nnl (nt) .eq.3) then
           ps = '(our)'
        else
           ps = '     '

        endif

        WRITE( stdout, '(/5x,"pseudo",i2," is ",a2, &
             &        1x,a5,"   zval =",f5.1,"   lmax=",i2, &
             &        "   lloc=",i2)') nt, psd (nt) , ps, zp (nt) , lmax (nt) &
             &, lloc (nt)
        if (numeric (nt) ) then
           WRITE( stdout, '(5x,"(in numerical form: ",i3, &
                &" grid points",", xmin = ",f5.2,", dx = ", &
                &f6.4,")")') mesh (nt) , xmin (nt) , dx (nt)
        else
           WRITE( stdout, '(/14x,"i=",7x,"1",13x,"2",10x,"3")')
           WRITE( stdout, '(/5x,"core")')
           WRITE( stdout, '(5x,"alpha =",4x,3g13.5)') (alpc (i, nt) , i = &
                1, 2)
           WRITE( stdout, '(5x,"a(i)  =",4x,3g13.5)')  (cc (i, nt) , i = 1, 2)
           do l = 0, lmax (nt)
              WRITE( stdout, '(/5x,"l = ",i2)') l
              WRITE( stdout, '(5x,"alpha =",4x,3g13.5)') (alps (i, l, nt) , &
                   i = 1, 3)
              WRITE( stdout, '(5x,"a(i)  =",4x,3g13.5)')  (aps (i, l, nt) , i = 1, &
                   &3)
              WRITE( stdout, '(5x,"a(i+3)=",4x,3g13.5)') (aps (i, l, nt) , i &
                   = 4, 6)
           enddo
           if (nlcc (nt) ) WRITE( stdout, 200) a_nlcc (nt), b_nlcc (nt), &
                alpha_nlcc (nt)
200        format(/5x,'nonlinear core correction: ', &
                &          'rho(r) = ( a + b r^2) exp(-alpha r^2)', &
                &          /,5x,'a    =',4x,g11.5, &
                &          /,5x,'b    =',4x,g11.5, &
                &          /,5x,'alpha=',4x,g11.5)
        endif
     endif
  enddo
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
  CALL flush_unit( stdout )
  !
  return
end subroutine d3_summary
