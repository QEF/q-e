  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! Adapted from the code PH/phq_summary - Quantum-ESPRESSO group     
  !
  !-----------------------------------------------------------------------
  SUBROUTINE epw_summary
  !-----------------------------------------------------------------------
  !!
  !!    Output symmetry informations
  !!
  USE ions_base,     ONLY : nat, ityp, atm, tau, ntyp => nsp, amass
  USE io_global,     ONLY : stdout
  USE cell_base,     ONLY : at, bg, ibrav, alat, omega, celldm
  USE klist,         ONLY : lgauss, degauss, ngauss, nkstot, xk, wk
  USE gvect,         ONLY : gcutm, ngm
  USE gvecs,         ONLY : dual, doublegrid, gcutms, ngms
  USE gvecw,         ONLY : ecutwfc
  USE symm_base,     ONLY : s, sname, ftau, s_axis_to_cart,sr
  USE funct,         ONLY : write_dft_name
  USE epwcom,        ONLY : title
  USE phcom,         ONLY : DP, tr2_ph, nmix_ph, alpha_mix
  USE lr_symm_base,  ONLY : irotmq, minus_q, nsymq, irgq
  USE control_flags, ONLY : iverbosity
  USE fft_base,      ONLY : dfftp, dffts
#if defined(__NAG)
  USE f90_unix_io,    ONLY : flush
#endif
  implicit none

  integer :: i, ipol, apol, na, isymq, isym, nsymtot, ik
  ! generic counter
  ! counter on angular momenta
  ! counter on atomic types
  ! counter on modes
  ! counter on symmetries
  ! counter on symmetries
  ! counter on symmetries
  ! counter on k points

  real(DP) :: ft1, ft2, ft3, xkg (3)
  ! fractionary translation
  ! fractionary translation
  ! fractionary translation
  ! the symmetry matrix in cartesian coord
  ! k point in crystal coordinates

  !
  WRITE( stdout, 100) title, ibrav, alat, omega, nat, ntyp, &
       ecutwfc, ecutwfc * dual, tr2_ph, alpha_mix (1), &
       nmix_ph
100 format (/,5x,a75,/,/,5x, &
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
     WRITE(stdout,'(/)') 
  !
  !    Description of the unit cell
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
       &(na, atm (ityp (na) ) , amass (ityp (na) )/ (1.66042d-24/9.1095d-28*0.5d0), na,  &
       &(tau (ipol, na) , ipol = 1, 3) , na = 1, nat)
  !
  !
  !   description of symmetries
  !
  WRITE( stdout, * )
  IF (nsymq.le.1.and..not.minus_q) THEN
     WRITE( stdout, '(5x,"No symmetry!")')
  ELSE
     IF (minus_q) THEN
        WRITE( stdout, '(5x,i2," Sym.Ops. (with q -> -q+G )",/)') &
             nsymq + 1
     ELSE
        WRITE( stdout, '(5x,i2," Sym.Ops. (no q -> -q+G )",/)') nsymq
     ENDIF

  ENDIF
  IF (iverbosity.eq.1) THEN

     WRITE( stdout, '(36x,"s",24x,"frac. trans.")')
     IF (minus_q) THEN
        nsymtot = nsymq + 1
     ELSE
        nsymtot = nsymq
     ENDIF
     DO isymq = 1, nsymtot
        IF (isymq.gt.nsymq) THEN
           isym = irotmq
           WRITE( stdout, '(/,5x,"This transformation sends q -> -q+G")')
        ELSE
           isym = irgq (isymq)
        ENDIF
        !WRITE( stdout, '(6x,"isym = ",i2,5x,a)') isymq, sname (isym)
        WRITE( stdout, *) 'isym = ', isymq, sname (isym)
      !  CALL s_axis_to_cart (s (1, 1, isym), sr, at, bg)
        CALL s_axis_to_cart ()
        IF (ftau (1, isym) .ne.0.or.ftau (2, isym) .ne.0.or.ftau (3, &
             isym) .ne.0) THEN
           ft1 = at (1, 1) * ftau (1, isym) / dfftp%nr1 + at (1, 2) * ftau ( &
                2, isym) / dfftp%nr2 + at (1, 3) * ftau (3, isym) / dfftp%nr3
           ft2 = at (2, 1) * ftau (1, isym) / dfftp%nr1 + at (2, 2) * ftau ( &
                2, isym) / dfftp%nr2 + at (2, 3) * ftau (3, isym) / dfftp%nr3
           ft3 = at (3, 1) * ftau (1, isym) / dfftp%nr1 + at (3, 2) * ftau ( &
                2, isym) / dfftp%nr2 + at (3, 3) * ftau (3, isym) / dfftp%nr3
           WRITE( stdout, '(1x,"cryst.",3x,"s(",i2,") = (",3i6, &
                &                    " )    f =( ",f10.7," )")') isymq,  (s (1 &
                &, ipol,isym) , ipol = 1, 3) , DBLE(ftau (1, isym) )/DBLE (dfftp%nr1)
           WRITE( stdout, '(17x," (",3(i6,5x), &
                &                    " )       ( ",f10.7," )")')  (s (2, ipol, &
                &isym) , ipol = 1, 3) , DBLE(ftau (2, isym) )  / DBLE(dfftp%nr2)
           WRITE( stdout, '(17x," (",3(i6,5x), &
                &                    " )       ( ",f10.7," )"/)')  (s (3, ipol, &
                & isym) , ipol = 1, 3) , DBLE (ftau (3, isym) )  / DBLE(dfftp%nr3)
           WRITE( stdout, '(1x,"cart.",3x,"s(",i2,") = (",3f11.7, &
                &                    " )    f =( ",f10.7," )")') isymq,  (sr (1 &
                &, ipol,isym) , ipol = 1, 3) , ft1
           WRITE( stdout, '(17x," (",3f11.7, &
                &                    " )       ( ",f10.7," )")')  (sr (2, ipol,isym) &
                & , ipol = 1, 3) , ft2
           WRITE( stdout, '(17x," (",3f11.7, &
                &                    " )       ( ",f10.7," )"/)')  (sr (3, ipol,isym &
                &) , ipol = 1, 3) , ft3
        ELSE
           WRITE( stdout, '(1x,"cryst.",3x,"s(",i2,") = (",3i6, &
                &                    " )")') isymq,  (s (1, ipol, isym) , ipol = &
                &1, 3)
           WRITE( stdout, '(17x," (",3(i6,5x)," )")') (s (2, ipol, isym) &
                , ipol = 1, 3)
           WRITE( stdout, '(17x," (",3(i6,5x)," )"/)') (s (3, ipol, &
                isym) , ipol = 1, 3)
           WRITE( stdout, '(1x,"cart.",3x,"s(",i2,") = (",3f11.7, &
                &                    " )")') isymq,  (sr (1, ipol,isym) , ipol = 1, 3)
           WRITE( stdout, '(17x," (",3f11.7," )")') (sr (2, ipol,isym) , &
                ipol = 1, 3)
           WRITE( stdout, '(17x," (",3f11.7," )"/)') (sr (3, ipol,isym) , &
                ipol = 1, 3)
        ENDIF
     ENDDO
  ENDIF
  !
  !     Description of the reciprocal lattice vectors
  !
  WRITE( stdout, '(/5x,"G cutoff =",f10.4,"  (", &
       &       i7," G-vectors)","     FFT grid: (",i3, &
       &       ",",i3,",",i3,")")') gcutm, ngm, dfftp%nr1, dfftp%nr2, dfftp%nr3

  IF (doublegrid) WRITE( stdout, '(5x,"G cutoff =",f10.4,"  (", &
         &                   i7," G-vectors)","  smooth grid: (",i3, &
         &                   ",",i3,",",i3,")")') gcutms, ngms, dffts%nr1, dffts%nr1, dffts%nr1
  IF (.NOT.lgauss) THEN
     WRITE( stdout, '(5x,"number of k points=",i5)') nkstot
  ELSE
     WRITE( stdout, '(5x,"number of k points=",i5, &
          &               "  gaussian broad. (Ry)=",f8.4,5x, &
          &               "ngauss = ",i3)') nkstot, degauss, ngauss
  ENDIF
  IF (iverbosity==1 .or. nkstot < 10000 ) THEN
     WRITE( stdout, '(23x,"cart. coord. in units 2pi/a_0")')
     DO ik = 1, nkstot
        WRITE( stdout, '(8x,"k(",i5,") = (",3f12.7,"), wk =",f12.7)') ik, &
             (xk (ipol, ik) , ipol = 1, 3) , wk (ik)
     ENDDO
  ENDIF
  IF (iverbosity.eq.1) THEN
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
  CALL print_ps_info ( )
  !
  CALL flush( stdout )
  !
  END SUBROUTINE epw_summary
