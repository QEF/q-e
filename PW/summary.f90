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
SUBROUTINE summary()
  !-----------------------------------------------------------------------
  !
  !    This routine writes on output all the information obtained from
  !    the input file and from the setup routine, before starting the
  !    self-consistent calculation.
  !
  !    if iverbosity = 0 only a partial summary is done.
  !
  USE io_global,       ONLY : stdout
  USE kinds,           ONLY : DP
  USE constants,       ONLY : amconv
  USE atom,            ONLY : numeric, xmin, dx, nlcc, mesh
  USE cell_base,       ONLY : alat, ibrav, omega, at, bg, celldm
  USE ions_base,       ONLY : nat, atm, zv, tau, ntyp => nsp, ityp
  USE char,            ONLY : title, sname
  USE cellmd,          ONLY : calc, cmass
  USE ions_base,       ONLY : amass
  USE gvect,           ONLY : nr1, nr2, nr3, dual, ecutwfc, ecfixed, q2sigma, &
                              ngm, gcutm, qcutz
  USE gsmooth,         ONLY : nr1s, nr2s, nr3s, doublegrid, ngms, gcutms
  USE lsda_mod,        ONLY : lsda, starting_magnetization
  USE klist,           ONLY : degauss, ngauss, lgauss, nkstot, xk, wk
  USE ktetra,          ONLY : ltetra
  USE pseud,           ONLY : zp, alps, alpc, cc, aps, nlc, nnl, lmax, lloc, &
                              a_nlcc, b_nlcc, alpha_nlcc
  USE symme,           ONLY : nsym, invsym, s, t_rev, ftau
  USE control_flags,   ONLY : imix, nmix, mixing_beta, nstep, diis_ndim, &
                              tr2, isolve, lmd, lbfgs, lpath, iverbosity
  USE uspp_param,      ONLY : nqf, rinner, nqlc, nbeta, iver, lll, psd, tvanp
  USE noncollin_module,ONLY : noncolin
  USE spin_orb,        ONLY : domag, lspinorb
  USE funct,           ONLY : write_dft_name
  USE bp,              ONLY : lelfield, gdir, nppstr, efield, nberrycyc
  USE fixed_occ,       ONLY : f_inp, tfixed_occ
  USE wvfct,           ONLY : nbnd
  USE lsda_mod,        ONLY : nspin
  !
  IMPLICIT NONE
  !
  ! ... declaration of the local variables
  !
  INTEGER :: i, ipol, apol, na, isym, ik, ib, nt, l, ngmtot, ibnd
    ! counter on the celldm elements
    ! counter on polarizations
    ! counter on direct or reciprocal lattice vect
    ! counter on atoms
    ! counter on symmetries
    ! counter on k points
    ! counter on beta functions
    ! counter on types
    ! counter on angular momenta
    ! total number of G-vectors (parallel executio
  REAL(DP) :: sr(3,3), ft1, ft2, ft3
    ! symmetry matrix in real axes
    ! fractionary translation
  REAL(DP), ALLOCATABLE :: xau(:,:)
    ! atomic coordinate referred to the crystal axes
  REAL(DP) :: xkg(3)
    ! coordinates of the k point in crystal axes
  CHARACTER :: mixing_style * 9
  CHARACTER :: ps * 5
    ! name of pseudo type
  REAL(DP) :: xp
    ! fraction contributing to a given atom type (obsolescent)
  !
  ! ... we start with a general description of the run
  !
  IF ( imix ==  0 ) mixing_style = 'plain'
  IF ( imix ==  1 ) mixing_style = 'TF'
  IF ( imix ==  2 ) mixing_style = 'local-TF'
  !
  IF ( title /= ' ') WRITE( stdout, "(/,5X,'Title: ',/,5X,A75)" ) title
  !
  WRITE( stdout, 100) ibrav, alat, omega, nat, ntyp, &
                      ecutwfc, dual * ecutwfc, tr2,  &
                      mixing_beta, nmix, mixing_style
  !
100 FORMAT( /,/,5X, &
       &     'bravais-lattice index     = ',I12,/,5X, &
       &     'lattice parameter (a_0)   = ',F12.4,'  a.u.',/,5X, &
       &     'unit-cell volume          = ',F12.4,' (a.u.)^3',/,5X, &
       &     'number of atoms/cell      = ',I12,/,5X, &
       &     'number of atomic types    = ',I12,/,5X, &
       &     'kinetic-energy cutoff     = ',F12.4,'  Ry',/,5X, &
       &     'charge density cutoff     = ',F12.4,'  Ry',/,5X, &
       &     'convergence threshold     = ',1PE12.1,/,5X, &
       &     'beta                      = ',0PF12.4,/,5X, &
       &     'number of iterations used = ',I12,2X,A,' mixing')
  !
  call write_dft_name
  !
  IF ( lmd .OR. lbfgs .OR. lpath ) &
     WRITE( stdout, '(5X,"nstep                     = ",I12,/)' ) nstep
  !
  IF ( lspinorb ) &
     WRITE( stdout, '(5X,"Noncollinear calculation with spin-orbit",/)' )
  !
  IF ( qcutz > 0.D0 ) THEN
     !
     WRITE( stdout, 110 ) ecfixed, qcutz, q2sigma
     !
110  FORMAT( 5X,'A smooth kinetic-energy cutoff is imposed at ', &
          &  F12.4,' Ry',/5X,'height of the smooth ', &
          &  'step-function =',F21.4,' Ry',/5X, &
          &  'width of the smooth step-function  =',F21.4,' Ry',/ )
     !
  END IF


  IF ( lelfield ) THEN !here informations for berry's phase el. fields calculations
     WRITE(stdout, *)
     WRITE(stdout, '(''     Using Berry phase electric field'')')
     WRITE(stdout, '(''     Direction :'', i4)') gdir
     WRITE(stdout, '(''     Intensity (a.u.) :'', f13.10)') efield
     WRITE(stdout, '(''     Strings composed by:'', i5,'' k-points'')') nppstr
     WRITE(stdout, '(''     Number of iterative cycles:'', i4)') nberrycyc
     WRITE(stdout, *)
  ENDIF
  !
  ! ... and here more detailed information. Description of the unit cell
  !
  WRITE( stdout, '(2(3X,3(2X,"celldm(",I1,")=",F11.6),/))' ) &
       ( i, celldm(i), i = 1, 6 )
  !
  WRITE( stdout, '(5X, &
       &     "crystal axes: (cart. coord. in units of a_0)",/, &
       &       3(15x,"a(",i1,") = (",3f10.6," )  ",/ ) )')  (apol,  &
       (at (ipol, apol) , ipol = 1, 3) , apol = 1, 3)
  !
  WRITE( stdout, '(5x, &
       &   "reciprocal axes: (cart. coord. in units 2 pi/a_0)",/, &
       &            3(15x,"b(",i1,") = (",3f10.6," )  ",/ ) )')  (apol,&
       &  (bg (ipol, apol) , ipol = 1, 3) , apol = 1, 3)
  !
  DO nt = 1, ntyp
     !
     IF ( tvanp(nt) ) THEN
        !
        ps = '(US)'
        WRITE( stdout, '(/5x,"PSEUDO",i2," is ",a2, &
             &        1x,a5,"   zval =",f5.1,"   lmax=",i2, &
             &        "   lloc=",i2)') nt, psd (nt) , ps, zp (nt) , lmax (nt) &
             &, lloc (nt)
        WRITE( stdout, '(5x,"Version ", 3i3, " of US pseudo code")') &
             (iver (i, nt) , i = 1, 3)
        WRITE( stdout, '(5x,"Using log mesh of ", i5, " points")') mesh (nt)
        WRITE( stdout, '(5x,"The pseudopotential has ",i2, &
             &       " beta functions with: ")') nbeta (nt)
        DO ib = 1, nbeta (nt)
           WRITE( stdout, '(15x," l(",i1,") = ",i3)') ib, lll (ib, nt)
        END DO
        WRITE( stdout, '(5x,"Q(r) pseudized with ", &
             &          i2," coefficients,  rinner = ",3f8.3,/ &
             &          52x,3f8.3,/ &
             &          52x,3f8.3)') nqf(nt), (rinner(i,nt), i=1,nqlc(nt) )
        !
     ELSE
        !
        IF ( nlc(nt) == 1 .AND. nnl(nt) == 1 ) THEN
           !
           ps = '(vbc)'
           !
        ELSE IF ( nlc(nt) == 2 .AND. nnl(nt) == 3 ) THEN
           !
           ps = '(bhs)'
           !
        ELSE IF ( nlc(nt) == 1 .AND. nnl(nt) == 3 ) THEN
           !
           ps = '(our)'
           !
        ELSE
           !
           ps = '     '
           !
        END IF
        !
        WRITE( stdout, '(/5x,"PSEUDO",i2," is ",a2, 1x,a5,"   zval =",f5.1,&
             &      "   lmax=",i2,"   lloc=",i2)') &
                        nt, psd(nt), ps, zp(nt), lmax(nt), lloc(nt)
        IF (numeric (nt) ) THEN
           WRITE( stdout, '(5x,"(in numerical form: ",i5,&
                &" grid points",", xmin = ",f5.2,", dx = ",f6.4,")")')&
                & mesh (nt) , xmin (nt) , dx (nt)
        ELSE
           WRITE( stdout, '(/14x,"i=",7x,"1",13x,"2",10x,"3")')
           WRITE( stdout, '(/5x,"core")')
           WRITE( stdout, '(5x,"alpha =",4x,3g13.5)') (alpc (i, nt) , i = 1, 2)
           WRITE( stdout, '(5x,"a(i)  =",4x,3g13.5)') (cc (i, nt) , i = 1, 2)
           DO l = 0, lmax(nt)
              WRITE( stdout, '(/5x,"l = ",i2)') l
              WRITE( stdout, '(5x,"alpha =",4x,3g13.5)') (alps (i, l, nt) , &
                   i = 1, 3)
              WRITE( stdout, '(5x,"a(i)  =",4x,3g13.5)') (aps (i, l, nt) , i = 1,3)
              WRITE( stdout, '(5x,"a(i+3)=",4x,3g13.5)') (aps (i, l, nt) , i= 4, 6)
           ENDDO
           IF ( nlcc(nt) ) WRITE( stdout, 200) a_nlcc(nt), b_nlcc(nt), alpha_nlcc(nt)
200        FORMAT(/5x,'nonlinear core correction: ', &
                &     'rho(r) = ( a + b r^2) exp(-alpha r^2)', &
                & /,5x,'a    =',4x,g11.5, &
                & /,5x,'b    =',4x,g11.5, &
                & /,5x,'alpha=',4x,g11.5)
        ENDIF
     ENDIF

  ENDDO
  WRITE( stdout, '(/5x, "atomic species   valence    mass     pseudopotential")')
  xp = 1.d0
  DO nt = 1, ntyp
     IF (calc.EQ.' ') THEN
        WRITE( stdout, '(5x,a6,6x,f10.2,2x,f10.5,5x,5 (a2,"(",f5.2,")"))') &
                   atm(nt), zv(nt), amass(nt), psd(nt), xp
     ELSE
        WRITE( stdout, '(5x,a6,6x,f10.2,2x,f10.5,5x,5 (a2,"(",f5.2,")"))') &
                   atm(nt), zv(nt), amass(nt)/amconv, psd(nt), xp
     END IF
  ENDDO

  IF (calc.EQ.'cd' .OR. calc.EQ.'cm' ) &
     WRITE( stdout, '(/5x," cell mass =", f10.5, " AMU ")') cmass/amconv
  IF (calc.EQ.'nd' .OR. calc.EQ.'nm' ) &
     WRITE( stdout, '(/5x," cell mass =", f10.5, " AMU/(a.u.)^2 ")') cmass/amconv

  IF (lsda) THEN
     WRITE( stdout, '(/5x,"Starting magnetic structure ", &
          &      /5x,"atomic species   magnetization")')
     DO nt = 1, ntyp
        WRITE( stdout, '(5x,a6,9x,f6.3)') atm(nt), starting_magnetization(nt)
     ENDDO
  ENDIF
  !
  !   description of symmetries
  !
  IF (nsym.LE.1) THEN
     WRITE( stdout, '(/5x,"No symmetry!")')
  ELSE
     IF (invsym) THEN
        WRITE( stdout, '(/5x,i2," Sym.Ops. (with inversion)",/)') nsym
     ELSE
        WRITE( stdout, '(/5x,i2," Sym.Ops. (no inversion)",/)') nsym
     ENDIF
  ENDIF
  IF (iverbosity.EQ.1) THEN
     WRITE( stdout, '(36x,"s",24x,"frac. trans.")')
     DO isym = 1, nsym
        WRITE( stdout, '(/6x,"isym = ",i2,5x,a45/)') isym, sname(isym)
        IF(noncolin.and.domag) WRITE(stdout,*) 'Time Reversal ', t_rev(isym)
        CALL s_axis_to_cart (s(1,1,isym), sr, at, bg)
        IF (ftau(1,isym).NE.0.OR.ftau(2,isym).NE.0.OR.ftau(3,isym).NE.0) THEN
           ft1 = at(1,1)*ftau(1,isym)/nr1 + at(1,2)*ftau(2,isym)/nr2 + &
                 at(1,3)*ftau(3,isym)/nr3
           ft2 = at(2,1)*ftau(1,isym)/nr1 + at(2,2)*ftau(2,isym)/nr2 + &
                 at(2,3)*ftau(3,isym)/nr3
           ft3 = at(3,1)*ftau(1,isym)/nr1 + at(3,2)*ftau(2,isym)/nr2 + &
                 at(3,3)*ftau(3,isym)/nr3
           WRITE( stdout, '(1x,"cryst.",3x,"s(",i2,") = (",3(i6,5x), &
                 &        " )    f =( ",f10.7," )")') &
                 isym, (s(1,ipol,isym),ipol=1,3), DBLE(ftau(1,isym))/DBLE(nr1)
           WRITE( stdout, '(17x," (",3(i6,5x), " )       ( ",f10.7," )")') &
                       (s(2,ipol,isym),ipol=1,3), DBLE(ftau(2,isym))/DBLE(nr2)
           WRITE( stdout, '(17x," (",3(i6,5x), " )       ( ",f10.7," )"/)') &
                       (s(3,ipol,isym),ipol=1,3), DBLE(ftau(3,isym))/DBLE(nr3)
           WRITE( stdout, '(1x,"cart. ",3x,"s(",i2,") = (",3f11.7, &
                 &        " )    f =( ",f10.7," )")') &
                 isym, (sr(1,ipol),ipol=1,3), ft1
           WRITE( stdout, '(17x," (",3f11.7, " )       ( ",f10.7," )")') &
                       (sr(2,ipol),ipol=1,3), ft2
           WRITE( stdout, '(17x," (",3f11.7, " )       ( ",f10.7," )"/)') &
                       (sr(3,ipol),ipol=1,3), ft3
        ELSE
           WRITE( stdout, '(1x,"cryst.",3x,"s(",i2,") = (",3(i6,5x), " )")') &
                                     isym,  (s (1, ipol, isym) , ipol = 1,3)
           WRITE( stdout, '(17x," (",3(i6,5x)," )")')  (s(2,ipol,isym), ipol=1,3)
           WRITE( stdout, '(17x," (",3(i6,5x)," )"/)') (s(3,ipol,isym), ipol=1,3)
           WRITE( stdout, '(1x,"cart. ",3x,"s(",i2,") = (",3f11.7," )")') &
                                         isym,  (sr (1, ipol) , ipol = 1, 3)
           WRITE( stdout, '(17x," (",3f11.7," )")')  (sr (2, ipol) , ipol = 1, 3)
           WRITE( stdout, '(17x," (",3f11.7," )"/)') (sr (3, ipol) , ipol = 1, 3)
        ENDIF
     ENDDO

  ENDIF
  !
  !    description of the atoms inside the unit cell
  !
  WRITE( stdout, '(/,3x,"Cartesian axes")')
  WRITE( stdout, '(/,5x,"site n.     atom                  positions (a_0 units)")')

  WRITE( stdout, '(7x,i3,8x,a6," tau(",i3,") = (",3f12.7,"  )")') &
             (na, atm(ityp(na)), na, (tau(ipol,na), ipol=1,3), na=1,nat)
  !
  !  output of starting magnetization
  !
  IF (iverbosity.EQ.1) THEN
     !
     !   allocate work space
     !
     ALLOCATE (xau(3,nat))
     !
     !     Compute the coordinates of each atom in the basis of the direct la
     !     vectors
     !
     DO na = 1, nat
        DO ipol = 1, 3
           xau(ipol,na) = bg(1,ipol)*tau(1,na) + bg(2,ipol)*tau(2,na) + &
                          bg(3,ipol)*tau(3,na)
        ENDDO
     ENDDO
     !
     !   description of the atoms inside the unit cell
     !   (in crystallographic coordinates)
     !
     WRITE( stdout, '(/,3x,"Crystallographic axes")')
     WRITE( stdout, '(/,5x,"site n.     atom        ", &
          &             "          positions (cryst. coord.)")')

     WRITE( stdout, '(7x,i2,8x,a6," tau(",i3,") = (",3f11.7,"  )")') &
           (na, atm(ityp(na)), na,  (xau(ipol,na), ipol=1,3), na=1,nat)
     !
     !   deallocate work space
     !
     DEALLOCATE(xau)
  ENDIF

  IF (lgauss) THEN
     WRITE( stdout, '(/5x,"number of k points=",i5, &
          &               "  gaussian broad. (ryd)=",f8.4,5x, &
          &               "ngauss = ",i3)') nkstot, degauss, ngauss
  ELSE IF (ltetra) THEN
     WRITE( stdout,'(/5x,"number of k points=",i5, &
          &        " (tetrahedron method)")') nkstot
  ELSE
     WRITE( stdout, '(/5x,"number of k points=",i5)') nkstot

  ENDIF
  WRITE( stdout, '(23x,"cart. coord. in units 2pi/a_0")')
  DO ik = 1, nkstot
     WRITE( stdout, '(8x,"k(",i5,") = (",3f12.7,"), wk =",f12.7)') ik, &
          (xk (ipol, ik) , ipol = 1, 3) , wk (ik)
  ENDDO
  IF (iverbosity.EQ.1) THEN
     WRITE( stdout, '(/23x,"cryst. coord.")')
     DO ik = 1, nkstot
        DO ipol = 1, 3
           xkg(ipol) = at(1,ipol)*xk(1,ik) + at(2,ipol)*xk(2,ik) + &
                       at(3,ipol)*xk(3,ik)
           ! xkg are the component in the crystal RL basis
        ENDDO
        WRITE( stdout, '(8x,"k(",i5,") = (",3f12.7,"), wk =",f12.7)') &
             ik, (xkg (ipol) , ipol = 1, 3) , wk (ik)
     ENDDO
  ENDIF
  ngmtot = ngm
#ifdef __PARA
  CALL ireduce (1, ngmtot)
#endif
  WRITE( stdout, '(/5x,"G cutoff =",f10.4,"  (", &
       &       i7," G-vectors)","     FFT grid: (",i3, &
       &       ",",i3,",",i3,")")') gcutm, ngmtot, nr1, nr2, nr3
  IF (doublegrid) THEN
     ngmtot = ngms
     !
     CALL ireduce (1, ngmtot)
     !
     WRITE( stdout, '(5x,"G cutoff =",f10.4,"  (", &
          &    i7," G-vectors)","  smooth grid: (",i3, &
          &    ",",i3,",",i3,")")') gcutms, ngmtot, nr1s, nr2s, nr3s
  ENDIF

  IF ( isolve == 2 ) &
     WRITE( stdout, '(/,5X,"reduced basis size: ",I5)' ) diis_ndim

  IF (tfixed_occ) THEN
     WRITE( stdout, '(/,5X,"Occupations read from input ")' ) 
     IF (nspin==2) THEN
        WRITE(stdout, '(/,5X," Spin-up")' ) 
        WRITE(stdout, '(/,(5X,8f9.4))') (f_inp(ibnd,1),ibnd=1,nbnd)
        WRITE(stdout, '(/,5X," Spin-down")' ) 
        WRITE(stdout, '(/,(5X,8f9.4))') (f_inp(ibnd,2),ibnd=1,nbnd)
     ELSE
        WRITE(stdout, '(/,(5X,8f9.4))') (f_inp(ibnd,1), ibnd=1,nbnd)
     END IF
  END IF

  !
  CALL flush_unit( stdout )
  !
  RETURN
  !
END SUBROUTINE summary
