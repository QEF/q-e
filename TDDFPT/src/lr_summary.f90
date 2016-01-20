!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------------
SUBROUTINE lr_summary
  !------------------------------------------------------------------------------
  !
  ! This subroutine writes on output the quantities which have been read
  ! from the punch file, and the quantities computed in the lr_init_nfo and
  ! lr_nscf_setup files.
  ! If lr_verbosity = 0 only a partial summary is done.
  ! Inspired by PH/phq_summary.f90
  !
  ! Created by Iurii Timrov (2013)
  !
  USE kinds,                  ONLY : DP
  USE ions_base,              ONLY : nat, ityp, atm, tau, ntyp => nsp, amass
  USE io_global,              ONLY : stdout
  USE cell_base,              ONLY : at, bg, ibrav, alat, omega, celldm
  USE klist,                  ONLY : lgauss, smearing, degauss, ngauss, nkstot, xk, wk
  USE fft_base,               ONLY : dfftp, dffts
  USE gvect,                  ONLY : gcutm, ngm
  USE gvecs,                  ONLY : doublegrid, dual, gcutms, ngms
  USE symm_base,              ONLY : s, sr, ftau, sname, t_rev
  USE noncollin_module,       ONLY : noncolin
  USE spin_orb,               ONLY : lspinorb, domag
  USE funct,                  ONLY : write_dft_name
  USE run_info,               ONLY : title
  USE qpoint,                 ONLY : xq
  USE gvecw,                  ONLY : ecutwfc
  USE lr_variables,           ONLY : lr_verbosity

  USE lr_symm_base,           ONLY : irotmq, irgq, minus_q, nsymq
  !
  IMPLICIT NONE
  !
  INTEGER :: i, ipol, apol, na, isymq, isym, nsymtot, ik
  ! generic counter
  ! counter on polarizations
  ! counter on polarizations
  ! counter on atoms
  ! counter on symmetries
  ! counter on symmetries
  ! counter on symmetries
  ! counter on k points
  ! counter on beta functions
  !
  REAL(DP) :: ft1, ft2, ft3, xkg (3)
  ! fractionary translations
  ! k point in crystal coordinates
  !
  CALL start_clock ('lr_summary')
  !
  WRITE( stdout, 100) title, ibrav, alat, omega, nat, ntyp, &
                    & ecutwfc, ecutwfc * dual
100 FORMAT (/,5x,a75,/,/,5x, &
       &     'bravais-lattice index     = ',i12,/,5x, &
       &     'lattice parameter (a_0)   = ',f12.4,'  a.u.',/,5x, &
       &     'unit-cell volume          = ',f12.4,' (a.u.)^3',/,5x, &
       &     'number of atoms/cell      = ',i12,/,5x, &
       &     'number of atomic types    = ',i12,/,5x, &
       &     'kinetic-energy cut-off    = ',f12.4,'  Ry',/,5x, &
       &     'charge density cut-off    = ',f12.4,'  Ry')
  !
  CALL write_dft_name ( )
  !
  !  Here add a message if this is a noncollinear or a spin_orbit calculation.
  !
  IF (noncolin) THEN
     IF (lspinorb) THEN
        IF (domag) THEN
           WRITE( stdout, '(5x, "Magnetic calculation with spin-orbit",/)')
        ELSE
           WRITE( stdout, '(5x, "Non magnetic calculation with spin-orbit",/)')
        ENDIF
     ELSE
        WRITE( stdout, '(5x, "Noncollinear calculation without spin-orbit",/)')
     END IF
  ELSE
     WRITE(stdout,'(/)')
  END IF
  !
  ! And here more detailed information. Description of the unit cell.
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
  ! Description of the atoms inside the unit cell.
  !
  WRITE( stdout, '(/, 5x,"Atoms inside the unit cell: ")')
  WRITE( stdout, '(/,3x,"Cartesian axes")')
  WRITE( stdout, '(/,5x,"site n.  atom      mass ", &
       &                "          positions (a_0 units)")')

  WRITE( stdout, '(7x,i2,5x,a6,f8.4,"   tau(",i2, &
       &                              ") = (",3f11.5,"  )")')  &
       &(na, atm (ityp (na) ) , amass (ityp (na) ), na,  &
       &(tau (ipol, na) , ipol = 1, 3) , na = 1, nat)
  WRITE( stdout, '(/,5x,"Linear response calculation for ")')
  WRITE( stdout, '(20x,"q = (",3f12.7," )")') (xq (ipol) , ipol = 1, 3)
  !
  ! Description of symmetries
  !
  WRITE( stdout, * )
  !
  IF (nsymq.le.1.and..not.minus_q) THEN
     WRITE( stdout, '(5x,"No symmetry!")')
  ELSE
     IF (minus_q) THEN
        WRITE( stdout, '(5x,i2," Sym.Ops. (with q -> -q+G )",/)') nsymq + 1
     ELSE
        WRITE( stdout, '(5x,i2," Sym.Ops. (no q -> -q+G )",/)') nsymq
     ENDIF
  ENDIF
  !
  !IF (lr_verbosity > 1) THEN
     !
     WRITE( stdout, '(36x,"s",24x,"frac. trans.")')
     !
     IF (minus_q) THEN
        nsymtot = nsymq + 1
     ELSE
        nsymtot = nsymq
     ENDIF
     !
     DO isymq = 1, nsymtot
        !
        IF (isymq.gt.nsymq) THEN
           isym = irotmq
           WRITE( stdout, '(/,5x,"This transformation sends q -> -q+G")')
        ELSE
           !isym = irgq (isymq)
           isym = isymq
        ENDIF
        !
        WRITE( stdout, '(/6x,"isym = ",i2,5x,a45/)') isymq, sname (isym)
        !
        IF (noncolin.and.domag) &
            WRITE(stdout,'(1x, "Time Reversal",i3)') t_rev(isym)
        !
        IF (ftau(1,isym).ne.0 .or. ftau(2,isym).ne.0 .or. ftau(3,isym).ne.0) THEN
           !
           ft1 = at (1, 1) * ftau (1, isym) / dfftp%nr1 + at (1, 2) * ftau ( &
                2, isym) / dfftp%nr2 + at (1, 3) * ftau (3, isym) / dfftp%nr3
           ft2 = at (2, 1) * ftau (1, isym) / dfftp%nr1 + at (2, 2) * ftau ( &
                2, isym) / dfftp%nr2 + at (2, 3) * ftau (3, isym) / dfftp%nr3
           ft3 = at (3, 1) * ftau (1, isym) / dfftp%nr1 + at (3, 2) * ftau ( &
                2, isym) / dfftp%nr2 + at (3, 3) * ftau (3, isym) / dfftp%nr3
           !
           WRITE( stdout, '(1x,"cryst.",3x,"s(",i2,") = (",3(i6,5x)," )  f =( ",f10.7," )")') &
                 isymq,  (s (1, ipol, isym) , ipol = 1, 3) , DBLE (ftau (1, isym) )  / DBLE (dfftp%nr1)
           WRITE( stdout, '(17x," (",3(i6,5x), " )     ( ",f10.7," )")')  &
                         (s (2, ipol, isym) , ipol = 1, 3) , DBLE (ftau (2, isym) )  / DBLE (dfftp%nr2)
           WRITE( stdout, '(17x," (",3(i6,5x)," )       ( ",f10.7," )"/)') &
                         (s (3, ipol, isym) , ipol = 1, 3) , DBLE (ftau (3, isym) )  / DBLE (dfftp%nr3)
           !
           WRITE( stdout, '(1x,"cart.",4x,"s(",i2,") = (",3f11.7, " )  f =( ",f10.7," )")') &
                 isymq, (sr (1, ipol,isym) , ipol = 1, 3) , ft1
           WRITE( stdout, '(17x," (",3f11.7, " )       ( ",f10.7," )")') &
                        (sr (2, ipol,isym) , ipol = 1, 3) , ft2
           WRITE( stdout, '(17x," (",3f11.7, " )       ( ",f10.7," )"/)') &
                        (sr (3, ipol,isym) , ipol = 1, 3) , ft3
        ELSE
           WRITE( stdout, '(1x,"cryst.",3x,"s(",i2,") = (",3(i6,5x)," )")') &
                        isymq,  (s (1, ipol, isym) , ipol = 1, 3)
           WRITE( stdout, '(17x," (",3(i6,5x)," )")') &
                        (s (2, ipol, isym), ipol = 1, 3)
           WRITE( stdout, '(17x," (",3(i6,5x)," )"/)') &
                        (s (3, ipol, isym) , ipol = 1, 3)
           !
           WRITE( stdout, '(1x,"cart.",4x,"s(",i2,") = (",3f11.7, " )")') &
                        isymq,  (sr (1, ipol,isym) , ipol = 1, 3)
           WRITE( stdout, '(17x," (",3f11.7," )")') &
                        (sr (2, ipol,isym) , ipol = 1, 3)
           WRITE( stdout, '(17x," (",3f11.7," )"/)') &
                        (sr (3, ipol,isym) , ipol = 1, 3)
        ENDIF
     ENDDO
  !ENDIF
  !
  ! Description of the reciprocal lattice vectors
  !
  WRITE( stdout, '(/5x,"G cutoff =",f10.4,"  (", &
       &       i7," G-vectors)","     FFT grid: (",i3, &
       &       ",",i3,",",i3,")")') gcutm, ngm, dfftp%nr1, dfftp%nr2, dfftp%nr3
  !
  IF (doublegrid) WRITE( stdout, '(5x,"G cutoff =",f10.4,"  (", &
       &                      i7," G-vectors)","  smooth grid: (",i3, &
       &                      ",",i3,",",i3,")")') gcutms, ngms, dffts%nr1, dffts%nr2, dffts%nr3
  IF (.NOT.lgauss) THEN
     WRITE( stdout, '(5x,"number of k points=",i6)') nkstot
  ELSE
     WRITE( stdout, '(/5x,"number of k points=", i6, 2x, &
          &             a," smearing, width (Ry)=",f8.4)') &
          &             nkstot, TRIM(smearing), degauss
  ENDIF
  !
  IF ( nkstot < 1000 ) THEN
     WRITE( stdout, '(23x,"cart. coord. in units 2pi/a_0")')
     DO ik = 1, nkstot
        WRITE( stdout, '(8x,"k(",i5,") = (",3f12.7,"), wk =",f12.7)') ik, &
             (xk (ipol, ik) , ipol = 1, 3) , wk (ik)
     ENDDO
  ENDIF
  !
  IF ( lr_verbosity > 1 ) THEN
     WRITE( stdout, '(/23x,"cryst. coord.")')
     DO ik = 1, nkstot
        DO ipol = 1, 3
           xkg (ipol) = at (1, ipol) * xk (1, ik) + at (2, ipol) * xk (2, &
                ik) + at (3, ipol) * xk (3, ik)
           ! xkg are the components  of xk in the reciprocal lattice basis
        ENDDO
        WRITE( stdout, '(8x,"k(",i5,") = (",3f12.7,"), wk =",f12.7)') &
             ik, (xkg (ipol) , ipol = 1, 3) , wk (ik)
     ENDDO
  ENDIF
  !
  ! Print information about the pseudopotential
  !
  CALL print_ps_info ( )
  !
  WRITE(stdout,'(/)')
  !
  FLUSH( stdout )
  !
  CALL stop_clock ('lr_summary')
  !
  RETURN
  !
END SUBROUTINE lr_summary
