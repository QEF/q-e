!
! Copyright (C) 2001-2022 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE forces()
  !----------------------------------------------------------------------------
  !! This routine is a driver routine which computes the forces
  !! acting on the atoms. The complete expression of the forces
  !! contains many parts which are computed by different routines:
  !
  !! - force_lc: local potential contribution 
  !! - force_us: non-local potential contribution
  !! - (esm_)force_ew: (ESM) electrostatic ewald term
  !! - force_cc: nonlinear core correction contribution
  !! - force_corr: correction term for incomplete self-consistency
  !! - force_hub: contribution due to the Hubbard term;
  !! - force_london: Grimme DFT+D dispersion forces
  !! - force_d3: Grimme-D3 (DFT-D3) dispersion forces
  !! - force_xdm: XDM dispersion forces
  !! - more terms from external electric fields, Martyna-Tuckerman, etc.
  !! - force_sol: contribution due to 3D-RISM
  !
  USE kinds,             ONLY : DP
  USE io_global,         ONLY : stdout
  USE cell_base,         ONLY : at, bg, alat, omega  
  USE ions_base,         ONLY : nat, ntyp => nsp,nsp, ityp, tau, zv, amass, extfor, atm
  USE gvect,             ONLY : ngm, gstart, ngl, igtongl, g, gg, gcutm
  USE lsda_mod,          ONLY : nspin
  USE symme,             ONLY : symvector
  USE vlocal,            ONLY : strf, vloc
  USE force_mod,         ONLY : force, sumfor
  USE scf,               ONLY : rho
  USE ions_base,         ONLY : if_pos
  USE ldaU,              ONLY : lda_plus_u, Hubbard_projectors
  USE extfield,          ONLY : tefield, forcefield, gate, forcegate, relaxz
  USE control_flags,     ONLY : gamma_only, remove_rigid_rot, textfor, &
                                iverbosity, llondon, ldftd3, lxdm, ts_vdw, &
                                mbd_vdw, lforce => tprnfor, istep
  USE bp,                ONLY : lelfield, gdir, l3dstring, efield_cart, &
                                efield_cry,efield
  USE uspp,              ONLY : okvan
  USE martyna_tuckerman, ONLY : do_comp_mt, wg_corr_force
  USE london_module,     ONLY : force_london
  USE dftd3_api,         ONLY : get_atomic_number
  USE dftd3_qe,          ONLY : dftd3_pbc_gdisp, dftd3

  USE xdm_module,        ONLY : force_xdm
  USE tsvdw_module,      ONLY : FtsvdW
  USE libmbd_interface,  ONLY : FmbdvdW
  USE esm,               ONLY : do_comp_esm, esm_bc, esm_force_ew
  USE qmmm,              ONLY : qmmm_mode
  USE rism_module,       ONLY : lrism, force_rism
  USE extffield,         ONLY : apply_extffield_PW
  USE input_parameters,  ONLY : nextffield
  !
#if defined(__CUDA)
  USE device_fbuff_m,          ONLY : dev_buf
#endif
  !
#if defined (__ENVIRON)
  USE plugin_flags,        ONLY : use_environ
  USE environ_base_module, ONLY : calc_environ_force
  USE environ_pw_module,   ONLY : is_ms_gcs, run_ms_gcs
#endif
#if defined (__OSCDFT)
  USE plugin_flags,        ONLY : use_oscdft
  USE oscdft_base,         ONLY : oscdft_ctx
  USE oscdft_forces_subs,  ONLY : oscdft_apply_forces, oscdft_print_forces
#endif
  !
  IMPLICIT NONE
  !
  REAL(DP), ALLOCATABLE :: forcenl(:,:),         &
                           forcelc(:,:),         &
                           forcecc(:,:),         &
                           forceion(:,:),        &
                           force_disp(:,:),      &
                           force_d3(:,:),        &
                           force_disp_xdm(:,:),  &
                           force_mt(:,:),        &
                           forcescc(:,:),        &
                           forces_bp_efield(:,:),&
                           forceh(:,:), &
                           force_sol(:,:)
  ! nonlocal, local, core-correction, ewald, scf correction terms, and hubbard
  !
  ! aux is used to store a possible additional density
  ! now defined in real space
  !
  COMPLEX(DP), ALLOCATABLE :: auxg(:), auxr(:)
  !
  REAL(DP) :: sumscf, sum_mm
  REAL(DP), PARAMETER :: eps = 1.e-12_dp
  INTEGER  :: ipol, na
  ! counter on polarization
  ! counter on atoms
  !
  REAL(DP), ALLOCATABLE :: taupbc(:,:)
  INTEGER :: atnum(1:nat)
  REAL(DP) :: stress_dftd3(3,3)
  !
  INTEGER :: ierr
  !
  force(:,:)    = 0.D0
  !
  ! Early return if all forces to be set to zero
  !
  IF ( ALL( if_pos == 0 ) ) RETURN
  !
  CALL start_clock( 'forces' )
#if defined(__CUDA)
  ! Cleanup scratch space used in previous SCF iterations.
  ! This will reduce memory footprint.
  CALL dev_buf%reinit(ierr)
  IF (ierr .ne. 0) CALL infomsg('forces', 'Cannot reset GPU buffers! Some buffers still locked.')
#endif
  !
  !
  ALLOCATE( forcenl(3,nat), forcelc(3,nat), forcecc(3,nat), &
            forceh(3,nat), forceion(3,nat), forcescc(3,nat) )
  !    
  forcescc(:,:) = 0.D0
  forceh(:,:)   = 0.D0
  !
  ! ... The nonlocal contribution is computed here
  !
  call start_clock('frc_us')
  CALL force_us( forcenl )
  call stop_clock('frc_us')
  !
  ! ... The local contribution
  !
  call start_clock('frc_lc')
  CALL force_lc( nat, tau, ityp, ntyp, alat, omega, ngm, ngl, igtongl, &
                 g, rho%of_r(:,1), gstart, gamma_only, vloc, forcelc )
  call stop_clock('frc_lc')
  !
  ! ... The NLCC contribution
  !
  call start_clock('frc_cc')
  CALL force_cc( forcecc )
  call stop_clock('frc_cc')

  ! ... The Hubbard contribution
  !     (included by force_us if using beta as local projectors)
  !
  IF ( lda_plus_u .AND. Hubbard_projectors.NE.'pseudo' ) CALL force_hub( forceh )
  !
  ! ... The ionic contribution is computed here
  !
  IF( do_comp_esm ) THEN
     CALL esm_force_ew( forceion )
  ELSE
     CALL force_ew( alat, nat, ntyp, ityp, zv, at, bg, tau, omega, g, &
                    gg, ngm, gstart, gamma_only, gcutm, strf, forceion )
  ENDIF
  !
  ! ... the semi-empirical dispersion correction
  !
  IF ( llondon ) THEN
    !
    ALLOCATE( force_disp(3,nat) )
    force_disp(:,:) = 0.0_DP
    force_disp = force_london( alat , nat , ityp , at , bg , tau )
    !
  ENDIF
  !
  ! ... The Grimme-D3 dispersion correction
  !
  IF ( ldftd3 ) THEN
    !
    CALL start_clock('force_dftd3')
    ALLOCATE( force_d3(3, nat) )
    force_d3(:,:) = 0.0_DP
    ! taupbc are atomic positions in alat units, centered around r=0
    ALLOCATE ( taupbc(3,nat) )
    taupbc(:,:) = tau(:,:)
    CALL cryst_to_cart( nat, taupbc, bg, -1 ) 
    taupbc(:,:) = taupbc(:,:) - NINT(taupbc(:,:))
    CALL cryst_to_cart( nat, taupbc, at,  1 ) 
    atnum(:) = get_atomic_number(atm(ityp(:)))
    CALL dftd3_pbc_gdisp( dftd3, alat*taupbc, atnum, alat*at, &
                          force_d3, stress_dftd3 )
    force_d3 = -2.d0*force_d3
    DEALLOCATE( taupbc)
    CALL stop_clock('force_dftd3')
  ENDIF
  !
  !
  IF (lxdm) THEN
     ALLOCATE( force_disp_xdm(3,nat) )
     force_disp_xdm = 0._dp
     force_disp_xdm = force_xdm(nat)
  ENDIF
  !
  ! ... The SCF contribution
  !
  call start_clock('frc_scc')
#if defined(__CUDA)
  ! Cleanup scratch space again, next subroutines uses a lot of memory.
  ! In an ideal world this should be done only if really needed (TODO).
  CALL dev_buf%reinit(ierr)
  IF (ierr .ne. 0) CALL errore('forces', 'Cannot reset GPU buffers! Buffers still locked: ', abs(ierr))
#endif
  !
  CALL force_corr( forcescc )
  call stop_clock('frc_scc') 
  !
  IF (do_comp_mt) THEN
    !
    ALLOCATE( force_mt(3,nat) )
    CALL wg_corr_force( .TRUE., omega, nat, ntyp, ityp, ngm, g, tau, zv, strf, &
                        rho%of_g(:,1), force_mt )
  ENDIF
  !
  ! ... The solvation contribution (3D-RISM)
  !
  IF (lrism) THEN
     ALLOCATE ( force_sol ( 3 , nat ) )
     CALL force_rism( force_sol )
  END IF
  !
  ! ... call void routine for user define/ plugin patches on internal forces
  !
#if defined(__LEGACY_PLUGINS)
  CALL plugin_int_forces() 
#endif 
#if defined (__ENVIRON)
  IF (use_environ) CALL calc_environ_force(force)
#endif
#if defined (__OSCDFT)
  IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==1)) CALL oscdft_apply_forces(oscdft_ctx)
#endif
  !
  ! ... Berry's phase electric field terms
  !
  IF (lelfield) THEN
     ALLOCATE( forces_bp_efield(3,nat) )
     forces_bp_efield(:,:) = 0.d0
     IF (.NOT.l3dstring) THEN
        IF (okvan) CALL forces_us_efield( forces_bp_efield, gdir, efield )
        CALL forces_ion_efield( forces_bp_efield, gdir, efield )
     ELSE
        IF (okvan) THEN
           DO ipol = 1, 3
              CALL forces_us_efield( forces_bp_efield, ipol, efield_cry(ipol) )
           ENDDO
        ENDIF
        DO ipol = 1, 3
           CALL forces_ion_efield( forces_bp_efield, ipol, efield_cart(ipol) )
        ENDDO
     ENDIF
  ENDIF
  !
  ! ... here we sum all the contributions and compute the total force acting
  ! ... on the crystal
  !
  DO ipol = 1, 3
     !
     sumfor = 0.D0
     !
     DO na = 1, nat
        !
        force(ipol,na) = force(ipol,na)    + &
                         forcenl(ipol,na)  + &
                         forceion(ipol,na) + &
                         forcelc(ipol,na)  + &
                         forcecc(ipol,na)  + &
                         forceh(ipol,na)   + &
                         forcescc(ipol,na)
        !
        IF ( llondon )  force(ipol,na) = force(ipol,na) + force_disp(ipol,na)
        IF ( ldftd3 )   force(ipol,na) = force(ipol,na) + force_d3(ipol,na)
        IF ( lxdm )     force(ipol,na) = force(ipol,na) + force_disp_xdm(ipol,na)
        ! factor 2 converts from Ha to Ry a.u.
        ! the IF condition is to avoid double counting
        IF ( mbd_vdw ) THEN
          force(ipol, na) = force(ipol, na) + 2.0_dp*FmbdvdW(ipol, na)
        ELSE IF ( ts_vdw ) THEN
          force(ipol, na) = force(ipol, na) + 2.0_dp*FtsvdW(ipol, na)
        ENDIF
        IF ( tefield )  force(ipol,na) = force(ipol,na) + forcefield(ipol,na)
        IF ( gate )     force(ipol,na) = force(ipol,na) + forcegate(ipol,na) ! TB
        IF (lelfield)   force(ipol,na) = force(ipol,na) + forces_bp_efield(ipol,na)
        IF (do_comp_mt) force(ipol,na) = force(ipol,na) + force_mt(ipol,na) 
        IF ( lrism )    force(ipol,na) = force(ipol,na) + force_sol(ipol,na)
        !
        sumfor = sumfor + force(ipol,na)
        !
     ENDDO
     !
     !TB
     IF ((gate.AND.relaxz).AND.(ipol==3)) WRITE( stdout, '("Total force in z direction = 0 disabled")')
     !
     IF ( (do_comp_esm .AND. ( esm_bc /= 'pbc' )).OR.(gate.AND.relaxz) ) THEN
        !
        ! ... impose total force along xy = 0
        !
        DO na = 1, nat
           IF ( ipol /= 3) force(ipol,na) = force(ipol,na)  &
                                            - sumfor / DBLE( nat )
        ENDDO
        !
     ELSEIF ( qmmm_mode < 0 ) THEN
        !
        ! ... impose total force = 0 except in a QM-MM calculation
        !
        DO na = 1, nat
           force(ipol,na) = force(ipol,na) - sumfor / DBLE( nat ) 
        ENDDO
        !
     ENDIF
     !
  ENDDO
  !
  ! ... call run_extffield to apply external force fields on ions
  ! 
  IF ( nextffield > 0 ) THEN 
     tau(:,:) = tau(:,:)*alat
     CALL apply_extffield_PW(istep,nextffield,tau,force)
     tau(:,:) = tau(:,:)/alat
  END IF
  !
  ! ... resymmetrize (should not be needed, but ...)
  !
  CALL symvector( nat, force )
  !
  IF ( remove_rigid_rot ) &
     CALL remove_tot_torque( nat, tau, amass(ityp(:)), force  )
  !
  IF( textfor ) force(:,:) = force(:,:) + extfor(:,:)
  !
  ! ... call void routine for user define/ plugin patches on external forces
  !
#if defined(__LEGACY_PLUGINS)
  CALL plugin_ext_forces() 
#endif 
#if defined (__ENVIRON)
  IF (use_environ) THEN
     IF (is_ms_gcs()) CALL run_ms_gcs()
  END IF
#endif
  !
  ! ... write on output the forces
  !
  WRITE( stdout, '(/,5x,"Forces acting on atoms (cartesian axes, Ry/au):", / )')
  DO na = 1, nat
     WRITE( stdout, 9035) na, ityp(na), force(:,na)
  ENDDO
  !
  ! ... forces on fixed coordinates are set to zero ( C.S. 15/10/2003 )
  !
  force(:,:)    = force(:,:)    * DBLE( if_pos )
  forcescc(:,:) = forcescc(:,:) * DBLE( if_pos )
  !
  IF ( iverbosity > 0 ) THEN
     !
     IF ( do_comp_mt ) THEN
        WRITE( stdout, '(5x,"The Martyna-Tuckerman correction term to forces")')
        DO na = 1, nat
           WRITE( stdout, 9035) na, ityp(na), ( force_mt(ipol,na), ipol = 1, 3 )
        ENDDO
     END IF
     !
     WRITE( stdout, '(5x,"The non-local contrib.  to forces")')
     DO na = 1, nat
        WRITE( stdout, 9035) na, ityp(na), ( forcenl(ipol,na), ipol = 1, 3 )
     ENDDO
     WRITE( stdout, '(5x,"The ionic contribution  to forces")')
     DO na = 1, nat
        WRITE( stdout, 9035) na, ityp(na), ( forceion(ipol,na), ipol = 1, 3 )
     ENDDO
     WRITE( stdout, '(5x,"The local contribution  to forces")')
     DO na = 1, nat
        WRITE( stdout, 9035) na, ityp(na), ( forcelc(ipol,na), ipol = 1, 3 )
     ENDDO
     WRITE( stdout, '(5x,"The core correction contribution to forces")')
     DO na = 1, nat
        WRITE( stdout, 9035) na, ityp(na), ( forcecc(ipol,na), ipol = 1, 3 )
     ENDDO
     WRITE( stdout, '(5x,"The Hubbard contrib.    to forces")')
     DO na = 1, nat
        WRITE( stdout, 9035) na, ityp(na), ( forceh(ipol,na), ipol = 1, 3 )
     ENDDO
     WRITE( stdout, '(5x,"The SCF correction term to forces")')
     DO na = 1, nat
        WRITE( stdout, 9035) na, ityp(na), ( forcescc(ipol,na), ipol = 1, 3 )
     ENDDO
     !
     IF ( llondon) THEN
        WRITE( stdout, '(/,5x,"Dispersion contribution to forces:")')
        DO na = 1, nat
           WRITE( stdout, 9035) na, ityp(na), (force_disp(ipol,na), ipol = 1, 3)
        ENDDO
     END IF
     !
     IF ( ldftd3 ) THEN
        WRITE( stdout, '(/,5x,"DFT-D3 dispersion contribution to forces:")')
        DO na = 1, nat
           WRITE( stdout, 9035) na, ityp(na), (force_d3(ipol,na), ipol = 1, 3)
        ENDDO
     END IF
     !
     IF (lxdm) THEN
        WRITE( stdout, '(/,5x,"XDM contribution to forces:")')
        DO na = 1, nat
           WRITE( stdout, 9035) na, ityp(na), (force_disp_xdm(ipol,na), ipol = 1, 3)
        ENDDO
     END IF
     !
     ! again, as above, if condition is to avoid redundant printing
     IF ( mbd_vdw ) THEN
        WRITE( stdout, '(/,5x, "MBD contribution to forces")')
        DO na = 1, nat
           WRITE( stdout, 9035) na, ityp(na), (2.0d0*FmbdvdW(ipol, na), ipol = 1, 3)
        ENDDO
     ELSE IF ( ts_vdw ) THEN
        WRITE( stdout, '(/,5x, "TS-VDW contribution to forces")')
        DO na = 1, nat
           WRITE( stdout, 9035) na, ityp(na), (2.0d0*FtsvdW(ipol, na), ipol = 1, 3)
        ENDDO
     ENDIF

     !
     ! TB gate forces
     IF ( gate ) THEN
        WRITE( stdout, '(/,5x,"Gate contribution to forces:")')
        DO na = 1, nat
           WRITE( stdout, 9035) na, ityp(na), (forcegate(ipol,na), ipol = 1, 3)
        ENDDO
     END IF
     !
     IF ( lrism ) THEN
        WRITE( stdout, '(/,5x,"3D-RISM Solvation contribution to forces:")')
        DO na = 1, nat
           WRITE( stdout, 9035) na, ityp(na), (force_sol(ipol,na), ipol = 1, 3)
        END DO
     END IF
     !
  END IF
#if defined (__OSCDFT)
  IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==1)) CALL oscdft_print_forces(oscdft_ctx)
#endif
  !
  sumfor = 0.D0
  sumscf = 0.D0
  !
  DO na = 1, nat
     !
     sumfor = sumfor + force(1,na)**2 + force(2,na)**2 + force(3,na)**2
     sumscf = sumscf + forcescc(1,na)**2 + forcescc(2,na)**2+ forcescc(3,na)**2
     !
  ENDDO
  !
  sumfor = SQRT( sumfor )
  sumscf = SQRT( sumscf )
  !
  WRITE( stdout, '(/5x,"Total force = ",F12.6,5X, &
              &  "Total SCF correction = ",F12.6)') sumfor, sumscf
  !
  IF ( llondon .AND. iverbosity > 0 ) THEN
     !
     sum_mm = 0.D0
     DO na = 1, nat
        sum_mm = sum_mm + &
                 force_disp(1,na)**2 + force_disp(2,na)**2 + force_disp(3,na)**2
     ENDDO
     sum_mm = SQRT( sum_mm )
     WRITE ( stdout, '(/,5x, "Total Dispersion Force = ",F12.6)') sum_mm
     !
  END IF
  !
  IF ( ldftd3 .AND. iverbosity > 0 ) THEN
     !
     sum_mm = 0.D0
     DO na = 1, nat
        sum_mm = sum_mm + &
                 force_d3(1,na)**2 + force_d3(2,na)**2 + force_d3(3,na)**2
     ENDDO
     sum_mm = SQRT( sum_mm )
     WRITE ( stdout, '(/,5x, "DFT-D3 dispersion Force = ",F12.6)') sum_mm
     !
  END IF
  !
  IF ( lxdm .AND. iverbosity > 0 ) THEN
     !
     sum_mm = 0.D0
     DO na = 1, nat
        sum_mm = sum_mm + &
                 force_disp_xdm(1,na)**2 + force_disp_xdm(2,na)**2 + force_disp_xdm(3,na)**2
     ENDDO
     sum_mm = SQRT( sum_mm )
     WRITE ( stdout, '(/,5x, "Total XDM Force = ",F12.6)') sum_mm
     !
  END IF
  !
  IF ( lrism .AND. iverbosity > 0 ) THEN
     !
     sum_mm = 0.D0
     DO na = 1, nat
        sum_mm = sum_mm + &
                 force_sol(1,na)**2 + force_sol(2,na)**2 + force_sol(3,na)**2
     END DO
     sum_mm = SQRT( sum_mm )
     WRITE ( stdout, '(/,5x, "Total 3D-RISM Solvation Force = ",F12.6)') sum_mm
     !
  END IF
  !
  DEALLOCATE( forcenl, forcelc, forcecc, forceh, forceion, forcescc )
  IF ( llondon  ) DEALLOCATE( force_disp       )
  IF ( ldftd3   ) DEALLOCATE( force_d3         )
  IF ( lxdm     ) DEALLOCATE( force_disp_xdm   ) 
  IF ( lelfield ) DEALLOCATE( forces_bp_efield )
  IF ( lrism    ) DEALLOCATE( force_sol        )
  IF(ALLOCATED(force_mt))   DEALLOCATE( force_mt )
  !
  ! FIXME: what is the following line good for?
  !
  lforce = .TRUE.
  !
  CALL stop_clock( 'forces' )
  !
  IF ( ( sumfor < 10.D0*sumscf ) .AND. ( sumfor > nat*eps ) ) &
  WRITE( stdout,'(5x,"SCF correction compared to forces is large: ", &
                   &  "reduce conv_thr to get better values")')
  RETURN
  !
9035 FORMAT(5X,'atom ',I4,' type ',I2,'   force = ',3F14.8)
  !
END SUBROUTINE forces
