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
subroutine phq_summary
  !-----------------------------------------------------------------------
  !
  !    This routine writes on output the quantities which have been read
  !    from the punch file, and the quantities computed in the phq_setup
  !    file.
  !
  !    if iverbosity = 0 only a partial summary is done.
  !
  !
  USE ions_base,     ONLY : nat, ityp, atm, tau, ntyp => nsp, amass
  USE io_global,     ONLY : stdout
  USE char,          ONLY : crystal, sname
  USE cell_base,     ONLY : at, bg, ibrav, alat, omega, celldm
  USE klist,         ONLY : lgauss, degauss, ngauss, nkstot, xk, wk
  USE gvect,         ONLY : ecutwfc, dual, nr1, nr2, nr3, gcutm, ngm
  USE gsmooth,       ONLY : doublegrid, nr1s, nr2s, nr3s, gcutms, ngms
  USE symme,         ONLY : s, ftau
  USE constants,     ONLY : amconv
  USE noncollin_module, ONLY : noncolin
  USE spin_orb,      ONLY : lspinorb, domag
  USE funct,         ONLY : write_dft_name
  USE printout_base, ONLY : title
  use phcom
  USE control_flags, ONLY : iverbosity
  
  implicit none

  integer :: i, l, nt, mu, nu, ipol, apol, na, isymq, isym, nsymtot, &
       ik, ib, irr, imode0
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

  real(DP) :: ft1, ft2, ft3, sr (3, 3), xkg (3)
  ! fractionary translation
  ! fractionary translation
  ! fractionary translation
  ! the symmetry matrix in cartesian coord
  ! k point in crystal coordinates

  !
  WRITE( stdout, 100) title, crystal, ibrav, alat, omega, nat, ntyp, &
       ecutwfc, ecutwfc * dual, tr2_ph, alpha_mix (1), &
       nmix_ph
100 format (/,5x,a75,/,/,5x, 'crystal is ',a20,/,/,5x, &
       &     'bravais-lattice index     = ',i12,/,5x, &
       &     'lattice parameter (a_0)   = ',f12.4,'  a.u.',/,5x, &
       &     'unit-cell volume          = ',f12.4,' (a.u.)^3',/,5x, &
       &     'number of atoms/cell      = ',i12,/,5x, &
       &     'number of atomic types    = ',i12,/,5x, &
       &     'kinetic-energy cut-off    = ',f12.4,'  Ry',/,5x, &
       &     'charge density cut-off    = ',f12.4,'  Ry',/,5x, &
       &     'convergence threshold     = ',1pe12.1,/,5x, &
       &     'beta                      = ',0pf12.4,/,5x, &
       &     'number of iterations used = ',i12)

  CALL write_dft_name ( )

  !
  !  Here add a message if this is a noncollinear or a spin_orbit calculation
  !
  IF (noncolin) THEN
     IF (lspinorb) THEN
        IF (domag) THEN 
           WRITE( stdout, '(5x, "Noncollinear calculation with spin-orbit",/)')
        ELSE
           WRITE( stdout, '(5x, "Non magnetic calculation with spin-orbit",/)')
        ENDIF
     ELSE
        WRITE( stdout, '(5x, "Noncollinear calculation without spin-orbit",/)')
     END IF
  END IF
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

  WRITE( stdout, '(7x,i2,5x,a6,f8.4,"   tau(",i2, &
       &                              ") = (",3f11.5,"  )")')  &
       &(na, atm (ityp (na) ) , amass (ityp (na) )  / amconv, na,  &
       &(tau (ipol, na) , ipol = 1, 3) , na = 1, nat)
  WRITE( stdout, '(/,5x,"Computing dynamical matrix for ")')
  WRITE( stdout, '(20x,"q = (",3f11.5," )")') (xq (ipol) , ipol = 1, 3)
  !
  !   description of symmetries
  !
  WRITE( stdout, * )
  if (nsymq.le.1.and..not.minus_q) then
     WRITE( stdout, '(5x,"No symmetry!")')
  else
     if (minus_q) then
        WRITE( stdout, '(5x,i2," Sym.Ops. (with q -> -q+G )",/)') &
             nsymq + 1
     else
        WRITE( stdout, '(5x,i2," Sym.Ops. (no q -> -q+G )",/)') nsymq
     endif

  endif
  if (iverbosity.eq.1) then

     WRITE( stdout, '(36x,"s",24x,"frac. trans.")')
     if (minus_q) then
        nsymtot = nsymq + 1
     else
        nsymtot = nsymq
     endif
     do isymq = 1, nsymtot
        if (isymq.gt.nsymq) then
           isym = irotmq
           WRITE( stdout, '(/,5x,"This transformation sends q -> -q+G")')
        else
           isym = irgq (isymq)
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
                &                    " )    f =( ",f10.7," )")') isymq,  (s (1, &
                & ipol, isym) , ipol = 1, 3) , DBLE (ftau (1, isym) )  / DBLE (nr1)
           WRITE( stdout, '(17x," (",3(i6,5x), &
                &                    " )       ( ",f10.7," )")')  (s (2, ipol, &
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
       &                      ",",i3,",",i3,")")') gcutms, ngms, nr1s, nr2s, nr3s
  if (.NOT.lgauss) then
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
           xkg (ipol) = at (1, ipol) * xk (1, ik) + at (2, ipol) * xk (2, &
                ik) + at (3, ipol) * xk (3, ik)
           ! xkg are the components  of xk in the reciprocal lattice basis
        enddo
        WRITE( stdout, '(8x,"k(",i5,") = (",3f12.7,"), wk =",f12.7)') &
             ik, (xkg (ipol) , ipol = 1, 3) , wk (ik)
     enddo

  endif

  CALL print_ps_info ( )

  WRITE( stdout, '(//5x,"Atomic displacements:")')
  WRITE( stdout, '(5x,"There are ",i3," irreducible representations")') nirr
  imode0 = 0
  do irr = 1, nirr
     if (done_irr (irr) .eq.1) then
        WRITE( stdout, '(/, 5x,"Representation ",i5,i7, &
             &                  " modes -  Done")') irr, npert (irr)
     elseif (comp_irr (irr) .eq.1) then
        WRITE( stdout, '(/, 5x,"Representation ",i5,i7, &
             &             " modes - To be done")') irr, npert (irr)
     elseif (comp_irr (irr) .eq.0) then
        if (lgamma_gamma) then
            if ((irr-1)/3+1==nasr) then
               WRITE( stdout, '(/, 5x,"Representation ",i5,i7, &
                 &     " modes - Calculated using asr")') irr, npert (irr)
            else
               WRITE( stdout, '(/, 5x,"Representation ",i5,i7, &
                 &     " modes - Calculated using symmetry")') irr, npert (irr)
            endif
        else
            WRITE( stdout, '(/, 5x,"Representation ",i5,i7, &
             &     " modes - Not done in this run")') irr, npert (irr)
        endif
     endif
     if (iverbosity.eq.1) then

        WRITE( stdout, '(5x,"Phonon polarizations are as follows:",/)')
        if (npert (irr) .eq.1) then
           WRITE( stdout, '(20x," mode # ",i3)') imode0 + 1
           WRITE( stdout, '(20x," (",2f10.5,"   ) ")')  ( (u (mu, nu) ,&
                &nu = imode0 + 1, imode0 + npert (irr) ) , mu = 1, 3 * nat)
        elseif (npert (irr) .eq.2) then
           WRITE( stdout, '(2(10x," mode # ",i3,16x))') imode0 + 1, &
                imode0 + 2
           WRITE( stdout, '(2(10x," (",2f10.5,"   ) "))')  ( (u (mu, nu) , nu &
                &= imode0 + 1, imode0 + npert (irr) ) , mu = 1, 3 * nat)
        else
           WRITE( stdout, '(4x,3(" mode # ",i3,13x))') imode0 + 1, imode0 &
                + 2, imode0 + 3
           WRITE( stdout, '((5x,3("(",2f10.5," ) ")))') ( (u (mu, nu) , &
                nu = imode0 + 1, imode0 + npert (irr) ) , mu = 1, 3 * nat)
        endif
        imode0 = imode0 + npert (irr)
     endif
  enddo
  if (.not.all_comp) then
     WRITE( stdout, '(/,5x,"Compute atoms: ",8(i5,","))') (atomo (na) &
          , na = 1, nat_todo)
  endif
  !
  CALL flush_unit( stdout )
  !
  return
end subroutine phq_summary
