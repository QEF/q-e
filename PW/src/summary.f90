!
! Copyright (C) 2001-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE summary()
  !-----------------------------------------------------------------------
  !
  !    This routine writes on output all the information obtained from
  !    the input file and from the setup routine, before starting the
  !    self-consistent calculation.
  !
  !    if iverbosity < 1 only a partial summary is done.
  !
  USE io_global,       ONLY : stdout
  USE kinds,           ONLY : DP
  USE run_info,        ONLY: title
  USE constants,       ONLY : amu_ry, rytoev
  USE cell_base,       ONLY : alat, ibrav, omega, at, bg, celldm
  USE ions_base,       ONLY : nat, atm, zv, tau, ntyp => nsp, ityp
  USE cellmd,          ONLY : calc, cmass
  USE ions_base,       ONLY : amass
  USE gvect,           ONLY : ecutrho, ngm, ngm_g, gcutm
  USE gvecs,           ONLY : doublegrid, ngms, ngms_g, gcutms
  USE fft_base,        ONLY : dfftp
  USE fft_base,        ONLY : dffts
  USE lsda_mod,        ONLY : lsda, starting_magnetization
  USE ldaU,            ONLY : lda_plus_U, Hubbard_u, Hubbard_j, Hubbard_alpha, &
                              Hubbard_l, lda_plus_u_kind, Hubbard_lmax,&
                              Hubbard_J0, Hubbard_beta
  USE klist,           ONLY : degauss, smearing, lgauss, nkstot, xk, wk, &
                              nelec, nelup, neldw, two_fermi_energies
  USE ktetra,          ONLY : ltetra
  USE control_flags,   ONLY : imix, nmix, mixing_beta, nstep, lscf, &
                              tr2, isolve, lmd, lbfgs, iverbosity, tqr, tq_smoothing, tbeta_smoothing
  USE noncollin_module,ONLY : noncolin
  USE spin_orb,        ONLY : domag, lspinorb
  USE funct,           ONLY : write_dft_name, dft_is_hybrid
  USE bp,              ONLY : lelfield, gdir, nppstr_3d, efield, nberrycyc, &
                              l3dstring,efield_cart,efield_cry
  USE fixed_occ,       ONLY : f_inp, tfixed_occ
  USE uspp_param,      ONLY : upf
  USE wvfct,           ONLY : nbnd
  USE gvecw,           ONLY : qcutz, ecfixed, q2sigma, ecutwfc
  USE mp_bands,        ONLY : intra_bgrp_comm
  USE mp,              ONLY : mp_sum
  USE esm,             ONLY : do_comp_esm, esm_summary
  USE martyna_tuckerman,ONLY: do_comp_mt
  USE realus,          ONLY : real_space
  USE exx,             ONLY : ecutfock
  USE fcp_variables,   ONLY : lfcpopt, lfcpdyn
  USE fcp,             ONLY : fcp_summary
  !
  IMPLICIT NONE
  !
  ! ... declaration of the local variables
  !
  INTEGER :: i, ipol, apol, na, ik, nt, ibnd, nk_
    ! counter on the celldm elements
    ! counter on polarizations
    ! counter on direct or reciprocal lattice vect
    ! counter on atoms
    ! counter on symmetries
    ! counter on k points
    ! counter on beta functions
    ! counter on types
    ! counter on bands
    ! actual number of k-points
    !
  REAL(DP), ALLOCATABLE :: xau(:,:)
    ! atomic coordinate referred to the crystal axes
  REAL(DP) :: xkg(3)
    ! coordinates of the k point in crystal axes
  CHARACTER :: mixing_style * 9
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
  WRITE( stdout, 103) nbnd, ecutwfc, ecutrho
  IF ( dft_is_hybrid () ) WRITE( stdout, 104) ecutfock
  IF ( lscf) WRITE( stdout, 105) tr2, mixing_beta, nmix, mixing_style
  !
100 FORMAT( /,/,5X, &
       &     'bravais-lattice index     = ',I12,/,5X, &
       &     'lattice parameter (alat)  = ',F12.4,'  a.u.',/,5X, &
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
       &     'charge density cutoff     = ',F12.4,'  Ry')
104 FORMAT(5X, &
       &     'cutoff for Fock operator  = ',F12.4,'  Ry')
105 FORMAT(5X, &
       &     'convergence threshold     = ',1PE12.1,/,5X, &
       &     'mixing beta               = ',0PF12.4,/,5X, &
       &     'number of iterations used = ',I12,2X,A,' mixing')
  !
  call write_dft_name ( ) 
  !
  IF ( lmd .OR. lbfgs ) &
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
  !
  CALL plugin_summary()
  !
  !
  ! ... ESM (Effective screening medium)
  !
  IF ( do_comp_esm )  CALL esm_summary()
  !
  ! ... FCP (Ficticious charge particle)
  !
  IF ( lfcpopt .or. lfcpdyn )  CALL fcp_summary()
  !
  IF ( do_comp_mt )  WRITE( stdout, &
            '(5X, "Assuming isolated system, Martyna-Tuckerman method",/)')

  IF ( lelfield ) THEN !here information for berry's phase el. fields calculations
     WRITE(stdout, *)
     WRITE(stdout, '("     Using Berry phase electric field")')
     if(.not.l3dstring) then
        WRITE(stdout, '("     Direction :", i4)') gdir
        WRITE(stdout, '("     Intensity (Ry a.u.) :", f13.10)') efield
        WRITE(stdout, '("     Strings composed by:", i5," k-points")') nppstr_3d(gdir)
     else
        write(stdout,'("     In a.u.(Ry)  cartesian system of reference" )')
        do i=1,3
           write(stdout,'(7x,f13.10)') efield_cart(i)
        enddo
        write(stdout,'("     In a.u.(Ry)  crystal system of reference" )')
        do i=1,3
           write(stdout,'(7x,f13.10)') efield_cry(i)
        enddo
     endif
     WRITE(stdout, '("     Number of iterative cycles:", i4)') nberrycyc
     WRITE(stdout, *)
  ENDIF
  !
  ! ... and here more detailed information. Description of the unit cell
  !
  WRITE( stdout, '(/2(3X,3(2X,"celldm(",I1,")=",F11.6),/))' ) &
       ( i, celldm(i), i = 1, 6 )
  !
  WRITE( stdout, '(5X, &
       &     "crystal axes: (cart. coord. in units of alat)",/, &
       &       3(15x,"a(",i1,") = (",3f11.6," )  ",/ ) )')  (apol,  &
       (at (ipol, apol) , ipol = 1, 3) , apol = 1, 3)
  !
  WRITE( stdout, '(5x, &
       &   "reciprocal axes: (cart. coord. in units 2 pi/alat)",/, &
       &            3(15x,"b(",i1,") = (",3f10.6," )  ",/ ) )')  (apol,&
       &  (bg (ipol, apol) , ipol = 1, 3) , apol = 1, 3)
  !
  CALL print_ps_info ( )
  !
  !
  ! ... print the vdw table information if needed
  CALL print_vdw_info ()
  !
  WRITE( stdout, '(/5x, "atomic species   valence    mass     pseudopotential")')
  xp = 1.d0
  DO nt = 1, ntyp
     WRITE( stdout, '(5x,a6,6x,f10.2,2x,f10.5,5x,5 (a2,"(",f5.2,")"))') &
                   atm(nt), zv(nt), amass(nt), upf(nt)%psd, xp
  ENDDO

  IF (calc.EQ.'cd' .OR. calc.EQ.'cm' ) &
     WRITE( stdout, '(/5x," cell mass =", f10.5, " AMU ")') cmass/amu_ry
  IF (calc.EQ.'nd' .OR. calc.EQ.'nm' ) &
     WRITE( stdout, '(/5x," cell mass =", f10.5, " AMU/(a.u.)^2 ")') cmass/amu_ry

  IF (lsda) THEN
     WRITE( stdout, '(/5x,"Starting magnetic structure ", &
          &      /5x,"atomic species   magnetization")')
     DO nt = 1, ntyp
        WRITE( stdout, '(5x,a6,9x,f6.3)') atm(nt), starting_magnetization(nt)
     ENDDO
  ENDIF
  !
  ! Some output for LDA+U
  !
  IF ( lda_plus_U ) THEN
     IF (lda_plus_u_kind == 0) THEN
        !
        WRITE( stdout, '(/,/,5x,"Simplified LDA+U calculation (l_max = ",i1, &
           &") with parameters (eV):")') Hubbard_lmax
        WRITE( stdout, '(5x,A)') &
           &"atomic species    L          U    alpha       J0     beta"
        DO nt = 1, ntyp
           IF ( Hubbard_U(nt) /= 0.D0 .OR. Hubbard_alpha(nt) /= 0.D0 .OR. &
                Hubbard_J0(nt) /= 0.D0 .OR. Hubbard_beta(nt) /= 0.D0 ) THEN
              WRITE( stdout,'(5x,a6,12x,i1,2x,4f9.4)') atm(nt), Hubbard_L(nt), &
                 Hubbard_U(nt)*rytoev, Hubbard_alpha(nt)*rytoev, &
                 Hubbard_J0(nt)*rytoev, Hubbard_beta(nt)*rytoev
           END IF
        END DO
        !
     ELSEIF(lda_plus_u_kind == 1) THEN
        !
        WRITE( stdout, '(/,/,5x,"Full LDA+U calculation (l_max = ",i1, &
           &") with parameters (eV):")') Hubbard_lmax
        DO nt = 1, ntyp
           IF (Hubbard_U(nt) /= 0.d0) THEN
              IF (Hubbard_l(nt) == 0) THEN
                 WRITE (stdout,'(5x,a,i2,a,f12.8)') &
                    'U(',nt,') =', Hubbard_U(nt) * rytoev
              ELSEIF (Hubbard_l(nt) == 1) THEN
                 WRITE (stdout,'(5x,2(a,i3,a,f9.4,3x))') &
                    'U(',nt,') =', Hubbard_U(nt)*rytoev, &
                    'J(',nt,') =', Hubbard_J(1,nt)*rytoev
              ELSEIF (Hubbard_l(nt) == 2) THEN
                 WRITE (stdout,'(5x,3(a,i3,a,f9.4,3x))') &
                    'U(',nt,') =', Hubbard_U(nt)*rytoev, &
                    'J(',nt,') =', Hubbard_J(1,nt)*rytoev, &
                    'B(',nt,') =', Hubbard_J(2,nt)*rytoev
              ELSEIF (Hubbard_l(nt) == 3) THEN
                 WRITE (stdout,'(5x,4(a,i3,a,f9.4,3x))') &
                    'U (',nt,') =', Hubbard_U(nt)*rytoev,   &
                    'J (',nt,') =', Hubbard_J(1,nt)*rytoev, &
                    'E2(',nt,') =', Hubbard_J(2,nt)*rytoev, &
                    'E3(',nt,') =', Hubbard_J(3,nt)*rytoev
              END IF
           END IF
        ENDDO
        IF (lspinorb) THEN
           WRITE(stdout, '(5x,"LDA+U on averaged j=l+1/2,l-1/2 radial WFs")')
        END IF
        !
      END IF
      !
      WRITE( stdout,'(/)')
  END IF
  !
  !   description of symmetries
  !
  CALL  print_symmetries ( iverbosity, noncolin, domag )
  !
  !    description of the atoms inside the unit cell
  !
  WRITE( stdout, '(/,3x,"Cartesian axes")')
  WRITE( stdout, '(/,5x,"site n.     atom                  positions (alat units)")')

  WRITE( stdout, '(6x,i4,8x,a6," tau(",i4,") = (",3f12.7,"  )")') &
             (na, atm(ityp(na)), na, (tau(ipol,na), ipol=1,3), na=1,nat)
  !
  !  output of starting magnetization
  !
  IF (iverbosity > 0) THEN
     !
     !   allocate work space
     !
     ALLOCATE (xau(3,nat))
     !
     !     Compute the coordinates of each atom in the basis of the
     !     direct lattice vectors
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

     WRITE( stdout, '(6x,i4,8x,a6," tau(",i4,") = (",3f11.7,"  )")') &
           (na, atm(ityp(na)), na,  (xau(ipol,na), ipol=1,3), na=1,nat)
     !
     !   deallocate work space
     !
     DEALLOCATE(xau)
  ENDIF

  IF ( lsda ) THEN
     !
     ! ... LSDA case: do not print replicated k-points
     !
     nk_ = nkstot/2
  ELSE
     nk_ = nkstot
  END IF

  IF (lgauss) THEN
     WRITE( stdout, '(/5x,"number of k points=", i6, 2x, &
          &             a," smearing, width (Ry)=",f8.4)') &
          &             nk_, TRIM(smearing), degauss
  ELSE IF (ltetra) THEN
     WRITE( stdout,'(/5x,"number of k points=",i6, &
          &        " (tetrahedron method)")') nk_
  ELSE
     WRITE( stdout, '(/5x,"number of k points=",i6)') nk_

  ENDIF
  IF ( iverbosity > 0 .OR. nk_ < 100 ) THEN
     WRITE( stdout, '(23x,"cart. coord. in units 2pi/alat")')
     DO ik = 1, nk_
        WRITE( stdout, '(8x,"k(",i5,") = (",3f12.7,"), wk =",f12.7)') ik, &
             (xk (ipol, ik) , ipol = 1, 3) , wk (ik)
     ENDDO
  ELSE
     WRITE( stdout, '(/5x,a)') &
     "Number of k-points >= 100: set verbosity='high' to print them."
  ENDIF
  IF ( iverbosity > 0 ) THEN
     WRITE( stdout, '(/23x,"cryst. coord.")')
     DO ik = 1, nk_
        DO ipol = 1, 3
           xkg(ipol) = at(1,ipol)*xk(1,ik) + at(2,ipol)*xk(2,ik) + &
                       at(3,ipol)*xk(3,ik)
           ! xkg are the component in the crystal RL basis
        ENDDO
        WRITE( stdout, '(8x,"k(",i5,") = (",3f12.7,"), wk =",f12.7)') &
             ik, (xkg (ipol) , ipol = 1, 3) , wk (ik)
     ENDDO
  ENDIF
  WRITE( stdout, '(/5x,"Dense  grid: ",i8," G-vectors", 5x, &
       &               "FFT dimensions: (",i4,",",i4,",",i4,")")') &
       &         ngm_g, dfftp%nr1, dfftp%nr2, dfftp%nr3
  IF (doublegrid) THEN
     WRITE( stdout, '(/5x,"Smooth grid: ",i8," G-vectors", 5x, &
       &               "FFT dimensions: (",i4,",",i4,",",i4,")")') &
       &         ngms_g, dffts%nr1, dffts%nr2, dffts%nr3
  ENDIF

  IF ( real_space ) WRITE( stdout, &
       & '(5x,"Real space treatment of Beta functions,", &
       &      " V.1 (BE SURE TO CHECK MANUAL!)")' )
  IF ( tbeta_smoothing ) WRITE( stdout, '(5x,"Beta functions are smoothed ")' )
  IF ( tqr ) WRITE( stdout, '(5x,"Real space treatment of Q(r)")' )
  IF ( tq_smoothing ) WRITE( stdout, '(5x,"Augmentation charges are smoothed ")' )

  IF (tfixed_occ) THEN
     WRITE( stdout, '(/,5X,"Occupations read from input ")' ) 
     IF ( lsda ) THEN
        WRITE(stdout, '(/,5X," Spin-up")' ) 
        WRITE(stdout, '(/,(5X,8f9.4))') (f_inp(ibnd,1),ibnd=1,nbnd)
        WRITE(stdout, '(/,5X," Spin-down")' ) 
        WRITE(stdout, '(/,(5X,8f9.4))') (f_inp(ibnd,2),ibnd=1,nbnd)
     ELSE
        WRITE(stdout, '(/,(5X,8f9.4))') (f_inp(ibnd,1), ibnd=1,nbnd)
     END IF
  END IF
  !
  FLUSH( stdout )
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
  USE io_files,        ONLY : pseudo_dir, psfile
  USE ions_base,       ONLY : ntyp => nsp
  USE atom,            ONLY : rgrid
  USE uspp_param,      ONLY : upf
  USE funct,           ONLY : dft_is_gradient
  IMPLICIT NONE
  !
  INTEGER :: nt, ib, i
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
     WRITE( stdout, '(/5x,"PseudoPot. #",i2," for ",a2," read from file:", &
           & /5x,a)') nt, upf(nt)%psd, TRIM(pseudo_dir)//TRIM (psfile(nt))
     WRITE( stdout, '(5x,"MD5 check sum: ", a )') upf(nt)%md5_cksum
     !
     WRITE( stdout, '( 5x,"Pseudo is ",a,", Zval =",f5.1)') &
            TRIM (ps), upf(nt)%zp
     !
     WRITE( stdout, '(5x,A)') TRIM(upf(nt)%generated)
     !
     IF(upf(nt)%tpawp) &
        WRITE( stdout, '(5x,a,a)') &
               "Shape of augmentation charge: ", TRIM(upf(nt)%paw%augshape)
     !
     ! info added for 1/r pseudos (AF)
     IF(upf(nt)%tcoulombp ) &
        WRITE( stdout, '(5x,a,a)') "1/r Coulomb pseudo"
     !
     WRITE( stdout, '(5x,"Using radial grid of ", i4, " points, ", &
         &i2," beta functions with: ")') rgrid(nt)%mesh, upf(nt)%nbeta
     DO ib = 1, upf(nt)%nbeta
        IF (ib < 10 ) THEN
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
!
!-----------------------------------------------------------------------
SUBROUTINE print_vdw_info
  !-----------------------------------------------------------------------
  !
  USE io_global,       ONLY : stdout
  USE io_files,        ONLY : psfile
  USE funct,           ONLY : get_inlc 
  USE kernel_table,    ONLY : vdw_table_name, vdw_kernel_md5_cksum

  integer :: inlc

  inlc = get_inlc()
  if ( inlc > 0 ) then
     WRITE( stdout, '(/5x,"vdW kernel table read from file ",a)') TRIM (vdw_table_name)
     WRITE( stdout, '(5x,"MD5 check sum: ", a )') vdw_kernel_md5_cksum
  endif 

END SUBROUTINE print_vdw_info
!
SUBROUTINE print_symmetries ( iverbosity, noncolin, domag )
  !-----------------------------------------------------------------------
  !
  USE kinds,           ONLY : dp
  USE io_global,       ONLY : stdout 
  USE symm_base,       ONLY : nsym, nsym_ns, nsym_na, invsym, s, sr, &
                              t_rev, ftau, sname
  USE rap_point_group, ONLY : code_group, nclass, nelem, elem, &
       which_irr, char_mat, name_rap, name_class, gname, ir_ram, elem_name
  USE rap_point_group_so, ONLY : nrap, nelem_so, elem_so, has_e, &
       which_irr_so, char_mat_so, name_rap_so, name_class_so, d_spin, &
       name_class_so1, elem_name_so
  USE rap_point_group_is, ONLY : nsym_is, sr_is, ftau_is, d_spin_is, &
       gname_is, sname_is, code_group_is
  USE cell_base,       ONLY : at, ibrav
  USE fft_base, ONLY : dfftp
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iverbosity
  LOGICAL, INTENT(IN) :: noncolin, domag
  !
  INTEGER :: nclass_ref   ! The number of classes of the point group
  INTEGER :: isym, ipol
  REAL (dp) :: ft1, ft2, ft3
  !
  !
  IF (nsym <= 1) THEN
     WRITE( stdout, '(/5x,"No symmetry found")')
  ELSE
     IF (invsym) THEN
        IF ( nsym_ns > 0 ) THEN
           WRITE( stdout, '(/5x,i2," Sym. Ops., with inversion, found ", &
                    &  "(",i2," have fractional translation)")' ) nsym, nsym_ns
        ELSE 
           WRITE( stdout, '(/5x,i2," Sym. Ops., with inversion, found")' )&
                         nsym
        END IF
     ELSE
        IF ( nsym_ns > 0 ) THEN
           WRITE( stdout, '(/5x,i2," Sym. Ops. (no inversion) found ",&
                    &  "(",i2," have fractional translation)")' ) nsym, nsym_ns
        ELSE
           WRITE( stdout,'(/5x,i2," Sym. Ops. (no inversion) found")' ) nsym
        END IF
     ENDIF
  ENDIF
  IF ( nsym_na > 0 ) THEN 
      WRITE( stdout, '(10x,"(note: ",i2," additional sym.ops. were found ", &
                   &   "but ignored",/,10x," their fractional translations ",&
                   &   "are incommensurate with FFT grid)",/)') nsym_na
  ELSE
      WRITE( stdout, '(/)' )
  END IF
  IF ( iverbosity > 0 ) THEN
     WRITE( stdout, '(36x,"s",24x,"frac. trans.")')
     nsym_is=0
     DO isym = 1, nsym
        WRITE( stdout, '(/6x,"isym = ",i2,5x,a45/)') isym, sname(isym)
        IF (noncolin) THEN
           IF (domag) THEN
              WRITE(stdout,*) 'Time Reversal ', t_rev(isym)
              IF (t_rev(isym)==0) THEN
                 nsym_is=nsym_is+1
                 sr_is(:,:,nsym_is) = sr(:,:,isym)
                 CALL find_u(sr_is(1,1,nsym_is), d_spin_is(1,1,nsym_is))
                 ftau_is(:,nsym_is)=ftau(:,isym)
                 sname_is(nsym_is)=sname(isym)
              ENDIF
           ELSE
              CALL find_u(sr(1,1,isym),d_spin(1,1,isym))
           END IF
        END IF
        IF ( ftau(1,isym).NE.0 .OR. ftau(2,isym).NE.0 .OR. &
             ftau(3,isym).NE.0) THEN
           ft1 = at(1,1)*ftau(1,isym)/dfftp%nr1 + at(1,2)*ftau(2,isym)/dfftp%nr2 + &
                at(1,3)*ftau(3,isym)/dfftp%nr3
           ft2 = at(2,1)*ftau(1,isym)/dfftp%nr1 + at(2,2)*ftau(2,isym)/dfftp%nr2 + &
                at(2,3)*ftau(3,isym)/dfftp%nr3
           ft3 = at(3,1)*ftau(1,isym)/dfftp%nr1 + at(3,2)*ftau(2,isym)/dfftp%nr2 + &
                at(3,3)*ftau(3,isym)/dfftp%nr3
           WRITE( stdout, '(1x,"cryst.",3x,"s(",i2,") = (",3(i6,5x), &
                &        " )    f =( ",f10.7," )")') &
                isym, (s(1,ipol,isym),ipol=1,3), DBLE(ftau(1,isym))/DBLE(dfftp%nr1)
           WRITE( stdout, '(17x," (",3(i6,5x), " )       ( ",f10.7," )")') &
                (s(2,ipol,isym),ipol=1,3), DBLE(ftau(2,isym))/DBLE(dfftp%nr2)
           WRITE( stdout, '(17x," (",3(i6,5x), " )       ( ",f10.7," )"/)') &
                (s(3,ipol,isym),ipol=1,3), DBLE(ftau(3,isym))/DBLE(dfftp%nr3)
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
     !
     CALL find_group(nsym,sr,gname,code_group)
     !
     ! ... Do not attempt calculation of classes if the lattice is provided
     ! ... in input as primitive vectors and not computed from parameters:
     ! ... the resulting vectors may not be accurate enough for the algorithm 
     !
     IF ( ibrav == 0 ) RETURN
     !
     IF (noncolin.AND.domag) THEN
        CALL find_group(nsym_is,sr_is,gname_is,code_group_is)
        CALL set_irr_rap_so(code_group_is,nclass_ref,nrap,char_mat_so, &
             name_rap_so,name_class_so,name_class_so1)
        CALL divide_class_so(code_group_is,nsym_is,sr_is,d_spin_is, &
             has_e,nclass,nelem_so,elem_so,which_irr_so)
        IF (nclass.ne.nclass_ref) CALL errore('summary', &
             'point double group ?',1)
        CALL set_class_el_name_so(nsym_is,sname_is,has_e,nclass,nelem_so, &
                                  elem_so,elem_name_so)
     ELSE
        IF (noncolin) THEN
           CALL set_irr_rap_so(code_group,nclass_ref,nrap,char_mat_so, &
                name_rap_so,name_class_so,name_class_so1)
           CALL divide_class_so(code_group,nsym,sr,d_spin,has_e,nclass,  &
                nelem_so, elem_so,which_irr_so)
           IF (nclass.ne.nclass_ref) CALL errore('summary', &
                'point double group ?',1)
           CALL set_class_el_name_so(nsym,sname,has_e,nclass,nelem_so, &
                                     elem_so,elem_name_so)
        ELSE
           CALL set_irr_rap(code_group,nclass_ref,char_mat,name_rap, &
                name_class,ir_ram)
           CALL divide_class(code_group,nsym,sr,nclass,nelem,elem,which_irr)
           IF (nclass.ne.nclass_ref) CALL errore('summary','point group ?',1)
           CALL set_class_el_name(nsym,sname,nclass,nelem,elem,elem_name)
        ENDIF
     ENDIF
     CALL write_group_info(.true.)
     !
  END IF
  !
END SUBROUTINE print_symmetries
