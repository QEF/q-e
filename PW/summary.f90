!
! Copyright (C) 2001-2007 PWSCF group
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
  USE cell_base,       ONLY : alat, ibrav, omega, at, bg, celldm
  USE ions_base,       ONLY : nat, atm, zv, tau, ntyp => nsp, ityp
  USE char,            ONLY : title, sname
  USE cellmd,          ONLY : calc, cmass
  USE ions_base,       ONLY : amass
  USE gvect,           ONLY : nr1, nr2, nr3, dual, ecutwfc, ecfixed, q2sigma, &
                              ngm, gcutm, qcutz
  USE gsmooth,         ONLY : nr1s, nr2s, nr3s, doublegrid, ngms, gcutms
  USE lsda_mod,        ONLY : lsda, starting_magnetization
  USE ldaU,            ONLY : lda_plus_U, Hubbard_u, Hubbard_alpha, &
                              Hubbard_l, Hubbard_lmax
  USE klist,           ONLY : degauss, ngauss, lgauss, nkstot, xk, wk, &
                              nelec, nelup, neldw, two_fermi_energies
  USE ktetra,          ONLY : ltetra
  USE symme,           ONLY : nsym, invsym, s, t_rev, ftau
  USE rap_point_group, ONLY : code_group, nclass, nelem, elem, which_irr, &
                              char_mat, name_rap, name_class, gname, ir_ram
  USE rap_point_group_so, ONLY : nrap, nelem_so, elem_so, has_e, which_irr_so, &
                              char_mat_so, name_rap_so, name_class_so, d_spin, &
                              name_class_so1
  USE rap_point_group_is, ONLY : nsym_is, sr_is, ftau_is, d_spin_is, gname_is, &
                                 sname_is, code_group_is
  USE control_flags,   ONLY : imix, nmix, mixing_beta, nstep, &
                              tr2, isolve, lmd, lbfgs, lpath, iverbosity
  USE noncollin_module,ONLY : noncolin
  USE spin_orb,        ONLY : domag, lspinorb
  USE funct,           ONLY : write_dft_name
  USE bp,              ONLY : lelfield, gdir, nppstr, efield, nberrycyc
  USE fixed_occ,       ONLY : f_inp, tfixed_occ
  USE uspp_param,      ONLY : upf
  USE wvfct,           ONLY : nbnd
  USE lsda_mod,        ONLY : nspin
  USE mp_global,       ONLY : intra_pool_comm
  USE mp,              ONLY : mp_sum
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
    ! total number of G-vectors (parallel execution)
  INTEGER :: &
          nclass_ref   ! The number of classes of the point group
  LOGICAL :: is_complex, is_complex_so
  INTEGER :: irot, iclass
    !
  REAL(DP) :: sr(3,3,48), ft1, ft2, ft3
    ! symmetry matrix in real axes
    ! fractionary translation
  REAL(DP), ALLOCATABLE :: xau(:,:)
    ! atomic coordinate referred to the crystal axes
  REAL(DP) :: xkg(3)
    ! coordinates of the k point in crystal axes
  CHARACTER :: mixing_style * 9
  CHARACTER :: group_name*11
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
  WRITE( stdout, 100) ibrav, alat, omega, nat, ntyp
  IF ( two_fermi_energies ) THEN
     WRITE( stdout, 101) nelec, nelup, neldw
  ELSE
     WRITE( stdout, 102) nelec
  END IF
  WRITE( stdout, 103) nbnd, ecutwfc, dual*ecutwfc,  &
                     tr2, mixing_beta, nmix, mixing_style
  !
100 FORMAT( /,/,5X, &
       &     'bravais-lattice index     = ',I12,/,5X, &
       &     'lattice parameter (a_0)   = ',F12.4,'  a.u.',/,5X, &
       &     'unit-cell volume          = ',F12.4,' (a.u.)^3',/,5X, &
       &     'number of atoms/cell      = ',I12,/,5X, &
       &     'number of atomic types    = ',I12)
101 FORMAT(5X, &
       &     'number of electrons       = ',F12.2,' (up:',f7.2,', down:',f7.2,')')
102 FORMAT(5X, &
       &     'number of electrons       = ',f12.2)
103 FORMAT(5X, &
       &     'number of Kohn-Sham states= ',I12,/,5X, &
       &     'kinetic-energy cutoff     = ',F12.4,'  Ry',/,5X, &
       &     'charge density cutoff     = ',F12.4,'  Ry',/,5X, &
       &     'convergence threshold     = ',1PE12.1,/,5X, &
       &     'mixing beta               = ',0PF12.4,/,5X, &
       &     'number of iterations used = ',I12,2X,A,' mixing')
  !
  call write_dft_name ( ) 
  !
  IF ( lmd .OR. lbfgs .OR. lpath ) &
     WRITE( stdout, '(5X,"nstep                     = ",I12,/)' ) nstep
  !
  IF (noncolin) THEN
     IF (lspinorb) THEN
        IF (domag) THEN
           WRITE( stdout, '(5X, "Noncollinear calculation with spin-orbit",/)')
        ELSE
           WRITE( stdout, '(5X, "Non magnetic calculation with spin-orbit",/)')
        ENDIF
     ELSE
        WRITE( stdout, '(5X, "Noncollinear calculation without spin-orbit",/)')
     END IF
  END IF
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
  WRITE( stdout, '(/2(3X,3(2X,"celldm(",I1,")=",F11.6),/))' ) &
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
  CALL print_ps_info ( )
  !
  WRITE( stdout, '(/5x, "atomic species   valence    mass     pseudopotential")')
  xp = 1.d0
  DO nt = 1, ntyp
     IF (calc.EQ.' ') THEN
        WRITE( stdout, '(5x,a6,6x,f10.2,2x,f10.5,5x,5 (a2,"(",f5.2,")"))') &
                   atm(nt), zv(nt), amass(nt), upf(nt)%psd, xp
     ELSE
        WRITE( stdout, '(5x,a6,6x,f10.2,2x,f10.5,5x,5 (a2,"(",f5.2,")"))') &
                   atm(nt), zv(nt), amass(nt)/amconv, upf(nt)%psd, xp
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
  IF ( lda_plus_U ) THEN
     WRITE( stdout, '(/5x,"LDA+U calculation, Hubbard_lmax = ",i1)') &
                    Hubbard_lmax
     WRITE( stdout, '(5x,"atomic species  L   Hubbard U  Hubbard alpha")' ) 
     DO nt = 1, ntyp
        IF ( Hubbard_U(nt) /= 0.D0 .OR. Hubbard_alpha(nt) /= 0.D0 ) THEN
           WRITE( stdout,'(5x,a6,5x,i6,2f12.6)') &
             atm(nt), Hubbard_L(nt), Hubbard_U(nt), Hubbard_alpha(nt)
        END IF
     END DO
  END IF
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
     nsym_is=0
     DO isym = 1, nsym
        WRITE( stdout, '(/6x,"isym = ",i2,5x,a45/)') isym, sname(isym)
        CALL s_axis_to_cart (s(1,1,isym), sr(1,1,isym), at, bg)
        IF (noncolin) THEN
           IF (domag) THEN
              WRITE(stdout,*) 'Time Reversal ', t_rev(isym)
              IF (t_rev(isym)==0) THEN
                 nsym_is=nsym_is+1
                 CALL s_axis_to_cart (s(1,1,isym), sr_is(1,1,nsym_is), at, bg)
                 CALL find_u(sr_is(1,1,nsym_is), d_spin_is(1,1,nsym_is))
                 ftau_is(:,nsym_is)=ftau(:,isym)
                 sname_is(nsym_is)=sname(isym)
              ENDIF
           ELSE
              CALL find_u(sr(1,1,isym),d_spin(1,1,isym))
           END IF
        END IF
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
                 isym, (sr(1,ipol,isym),ipol=1,3), ft1
           WRITE( stdout, '(17x," (",3f11.7, " )       ( ",f10.7," )")') &
                       (sr(2,ipol,isym),ipol=1,3), ft2
           WRITE( stdout, '(17x," (",3f11.7, " )       ( ",f10.7," )"/)') &
                       (sr(3,ipol,isym),ipol=1,3), ft3
        ELSE
           WRITE( stdout, '(1x,"cryst.",3x,"s(",i2,") = (",3(i6,5x), " )")') &
                                     isym,  (s (1, ipol, isym) , ipol = 1,3)
           WRITE( stdout, '(17x," (",3(i6,5x)," )")')  (s(2,ipol,isym), ipol=1,3)
           WRITE( stdout, '(17x," (",3(i6,5x)," )"/)') (s(3,ipol,isym), ipol=1,3)
           WRITE( stdout, '(1x,"cart. ",3x,"s(",i2,") = (",3f11.7," )")') &
                                         isym,  (sr (1, ipol,isym) , ipol = 1, 3)
           WRITE( stdout, '(17x," (",3f11.7," )")')  (sr (2, ipol,isym) , ipol = 1, 3)
           WRITE( stdout, '(17x," (",3f11.7," )"/)') (sr (3, ipol,isym) , ipol = 1, 3)
        END IF
     END DO
     CALL find_group(nsym,sr,gname,code_group)
     IF (noncolin.AND.domag) THEN
        CALL find_group(nsym_is,sr_is,gname_is,code_group_is)
        CALL set_irr_rap_so(code_group_is,nclass_ref,nrap,char_mat_so, &
                            name_rap_so,name_class_so,name_class_so1)
        CALL divide_class_so(code_group_is,nsym_is,sr_is,d_spin_is, &
                             has_e,nclass,nelem_so,elem_so,which_irr_so)
        IF (nclass.ne.nclass_ref) CALL errore('summary', &
                                                 'point double group ?',1)
     ELSE
        IF (noncolin) THEN
           CALL set_irr_rap_so(code_group,nclass_ref,nrap,char_mat_so, &
                               name_rap_so,name_class_so,name_class_so1)
           CALL divide_class_so(code_group,nsym,sr,d_spin,has_e,nclass,  &
                                nelem_so, elem_so,which_irr_so)
           IF (nclass.ne.nclass_ref) CALL errore('summary', &
                                                'point double group ?',1)
        ELSE
           CALL set_irr_rap(code_group,nclass_ref,char_mat,name_rap, &
                            name_class,ir_ram)
           CALL divide_class(code_group,nsym,sr,nclass,nelem,elem,which_irr)
           IF (nclass.ne.nclass_ref) CALL errore('summary','point group ?',1)
        ENDIF
     ENDIF
     CALL write_group_info(.true.)
  END IF
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

     WRITE( stdout, '(7x,i3,8x,a6," tau(",i3,") = (",3f11.7,"  )")') &
           (na, atm(ityp(na)), na,  (xau(ipol,na), ipol=1,3), na=1,nat)
     !
     !   deallocate work space
     !
     DEALLOCATE(xau)
  ENDIF

  IF (lgauss) THEN
     WRITE( stdout, '(/5x,"number of k points=",i5, &
          &               "  gaussian broad. (Ry)=",f8.4,5x, &
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
  CALL mp_sum (ngmtot, intra_pool_comm)
  WRITE( stdout, '(/5x,"G cutoff =",f10.4,"  (", &
       &       i7," G-vectors)","     FFT grid: (",i3, &
       &       ",",i3,",",i3,")")') gcutm, ngmtot, nr1, nr2, nr3
  IF (doublegrid) THEN
     !
     ngmtot = ngms
     CALL mp_sum (ngmtot, intra_pool_comm)
     !
     WRITE( stdout, '(5x,"G cutoff =",f10.4,"  (", &
          &    i7," G-vectors)","  smooth grid: (",i3, &
          &    ",",i3,",",i3,")")') gcutms, ngmtot, nr1s, nr2s, nr3s
  ENDIF

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
!
!-----------------------------------------------------------------------
SUBROUTINE print_ps_info
  !-----------------------------------------------------------------------
  !
  USE io_global,       ONLY : stdout
  USE io_files,        ONLY : psfile
  USE ions_base,       ONLY : ntyp => nsp
  USE atom,            ONLY : rgrid
  USE uspp_param,      ONLY : upf
  USE paw_variables,   ONLY : rad, xlm
  USE funct,           ONLY : dft_is_gradient

  !
  INTEGER :: nt, lmax
  CHARACTER :: ps*35
  !
  DO nt = 1, ntyp
     !
     IF ( upf(nt)%tpawp ) THEN
        ! Note: for PAW pseudo also tvanp is .true.
        ps="Projector augmented-wave"
     ELSE IF ( upf(nt)%tvanp ) THEN
        ps='Ultrasoft'
     ELSE
        ps='Norm-conserving'
     END IF
     !
     IF ( upf(nt)%nlcc ) ps = TRIM(ps) // ' + core correction'
     !
     WRITE( stdout, '(/5x,"PseudoPot. #",i2," for ",a2," read from file ",a)')&
             nt, upf(nt)%psd, TRIM (psfile(nt))
     !
     WRITE( stdout, '( 5x,"Pseudo is ",a,", Zval =",f5.1)') &
            TRIM (ps), upf(nt)%zp
     !
     WRITE( stdout, '(5x,A)') TRIM(upf(nt)%generated)
     !
     IF(upf(nt)%tpawp) THEN
        lmax = 2*upf(nt)%paw%lmax_rho
        IF ( dft_is_gradient() ) lmax = lmax + xlm
        WRITE( stdout, '(5x, a, i4, a, i3)') &
               "Setup to integrate on", ((lmax+1)*(lmax+2))/2, &
               " directions: integral exact up to l =", lmax
        WRITE( stdout, '(5x,a,a)') &
               "Shape of augmentation charge: ", TRIM(upf(nt)%paw%augshape)
     ENDIF
     WRITE( stdout, '(5x,"Using radial grid of ", i4, " points, ", &
         &i2," beta functions with: ")') rgrid(nt)%mesh, upf(nt)%nbeta
     DO ib = 1, upf(nt)%nbeta
        IF (ib<10) THEN
           WRITE( stdout, '(15x," l(",i1,") = ",i3)') ib, upf(nt)%lll(ib)
        ELSE
           WRITE( stdout, '(14x," l(",i2,") = ",i3)') ib, upf(nt)%lll(ib)
        ENDIF
     END DO

     IF ( upf(nt)%tvanp ) THEN
        IF (upf(nt)%nqf==0) THEN
           WRITE( stdout, '(5x,"Q(r) pseudized with 0 coefficients ",/)') 
        ELSE
           WRITE( stdout, '(5x,"Q(r) pseudized with ", &
           &          i2," coefficients,  rinner = ",3f8.3,/ &
           &          52x,3f8.3,/ 52x,3f8.3)') &
           &          upf(nt)%nqf, (upf(nt)%rinner(i), i=1,upf(nt)%nqlc)
        END IF
     ENDIF

  ENDDO
END SUBROUTINE print_ps_info
