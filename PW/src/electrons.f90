!
! Copyright (C) 2001-2025 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE electrons()
  !----------------------------------------------------------------------------
  !! General self-consistency loop, also for hybrid functionals.  
  !! For non-hybrid functionals it just calls \(\texttt{electron_scf}\).
  !
  USE kinds,                ONLY : DP
  USE check_stop,           ONLY : check_stop_now, stopped_by_user
  USE io_global,            ONLY : stdout, ionode
  USE fft_base,             ONLY : dfftp
  USE gvecs,                ONLY : doublegrid
  USE gvect,                ONLY : ecutrho
  USE lsda_mod,             ONLY : nspin, absmag
  USE ener,                 ONLY : etot, hwf_energy, eband, deband, ehart, &
                                   vtxc, etxc, etxcc, ewld, demet, epaw, &
                                   elondon, edftd3, vsol, esol, ef_up, ef_dw
  USE tsvdw_module,         ONLY : EtsvdW
  USE scf,                  ONLY : rho, rho_core, rhog_core, v, vltot, vrs, &
                                   kedtau, vnew
  USE control_flags,        ONLY : tr2, nexxiter, conv_elec, restart, lmd, &
                                   do_makov_payne, sic
  USE sic_mod,              ONLY : sic_energy, occ_f2fn, occ_fn2f, save_rhon, sic_first
  USE io_files,             ONLY : iunres, seqopn
  USE ldaU,                 ONLY : eth
  USE extfield,             ONLY : tefield, etotefield
  USE wvfct,                ONLY : nbnd, wg, et
  USE klist,                ONLY : nks
  USE uspp,                 ONLY : okvan
  USE exx,                  ONLY : aceinit,exxinit, exxenergy2, exxbuff, &
                                   fock0, fock1, fock2, fock3, dexx, use_ace, local_thr, &
                                   domat
  USE xc_lib,               ONLY : xclib_dft_is, exx_is_active, stop_exx
  USE control_flags,        ONLY : adapt_thr, tr2_init, tr2_multi, gamma_only
  !
  USE paw_variables,        ONLY : okpaw, ddd_paw, total_core_energy, only_paw
  USE paw_onecenter,        ONLY : PAW_potential
  USE paw_symmetry,         ONLY : PAW_symmetrize_ddd
  USE ions_base,            ONLY : nat
  USE loc_scdm,             ONLY : use_scdm, localize_orbitals
  USE loc_scdm_k,           ONLY : localize_orbitals_k
  !
  USE add_dmft_occ,         ONLY : dmft
  USE rism_module,          ONLY : lrism, rism_calc3d
  USE makovpayne,           ONLY : makov_payne
  USE vlocal,               ONLY : strf, vloc
  USE ions_base,            ONLY : tau
  !
  IMPLICIT NONE
  !
  ! ... a few local variables
  !
  REAL(DP) :: charge
  !! the total charge
  REAL(DP) :: exxen
  !! used to compute exchange energy
  REAL(DP), EXTERNAL :: exxenergyace
  INTEGER :: idum
  !! dummy counter on iterations
  INTEGER :: iter
  !! counter on iterations
  INTEGER :: printout, ik, ios
  !
  REAL(DP) :: tr2_min
  !! estimated error on energy coming from diagonalization
  REAL(DP) :: tr2_final
  !! final threshold for exx minimization 
  !! when using adaptive thresholds.
  LOGICAL :: first, exst
  REAL(DP) :: etot_cmp_paw(nat,2,2)
  LOGICAL :: DoLoc
  !
  !
  DoLoc = local_thr.gt.0.0d0
  exxen = 0.0d0
  iter = 0
  first = .true.
  tr2_final = tr2
  IF ( xclib_dft_is('hybrid') ) THEN
     !printout = 0  ! do not print etot and energy components at each scf step
     printout = 1  ! print etot, not energy components at each scf step
  ELSE IF ( lmd ) THEN
     printout = 1  ! print etot, not energy components at each scf step
  ELSE
     printout = 2  ! print etot and energy components at each scf step
  ENDIF
  IF (xclib_dft_is('hybrid') .AND. adapt_thr ) tr2= tr2_init
  fock0 = 0.D0
  fock1 = 0.D0
  fock3 = 0.D0
  !force higher verbosity for the printout in DMFT case
  IF (.NOT. exx_is_active () ) fock2 = 0.D0
  !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%  Iterate hybrid functional  %%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
  IF ( restart ) THEN
     CALL seqopn (iunres, 'restart_e', 'formatted', exst)
     IF ( exst ) THEN
        ios = 0
        READ (iunres, *, iostat=ios) iter, tr2, dexx
        IF ( ios /= 0 ) THEN
           iter = 0
        ELSE IF ( iter < 0 .OR. iter > nexxiter ) THEN
           iter = 0
        ELSE 
           READ (iunres, *) exxen, fock0, fock1, fock2
           ! FIXME: et and wg should be read from xml file
           READ (iunres, *) (wg(1:nbnd,ik),ik=1,nks)
           READ (iunres, *) (et(1:nbnd,ik),ik=1,nks)
           !$acc update device(et)
           CLOSE ( unit=iunres, status='delete')
           ! ... if restarting here, exx was already active
           ! ... initialize stuff for exx
           first = .false.
!civn  see non-scf comment about this
           Call stop_exx()
!
           CALL exxinit(DoLoc)
           IF( DoLoc.and.gamma_only) THEN
             CALL localize_orbitals( )
           ELSE IF (DoLoc) THEN
             CALL localize_orbitals_k( )
           ENDIF 
           ! FIXME: ugly hack, overwrites exxbuffer from exxinit
           CALL seqopn (iunres, 'restart_exx', 'unformatted', exst)
           IF (exst) READ (iunres, iostat=ios) exxbuff
           IF (ios /= 0) WRITE(stdout,'(5x,"Error in EXX restart!")')
!civn 
           !IF (use_ace) CALL aceinit ( DoLoc )
           domat = .false.
! 
           !
           CALL v_of_rho( rho, rho_core, rhog_core, &
               ehart, etxc, vtxc, eth, etotefield, charge, v)
           IF (lrism) CALL rism_calc3d(rho%of_g(:, 1), esol, vsol, v%of_r, tr2)
           IF (okpaw) CALL PAW_potential(rho%bec, ddd_paw, epaw,etot_cmp_paw)
           CALL set_vrs( vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, &
                         nspin, doublegrid )
           !
           WRITE(stdout,'(5x,"Calculation (EXX) restarted from iteration #", &
                        & i6)') iter 
        ENDIF
     ENDIF
     CLOSE ( unit=iunres, status='delete')
  ENDIF
  !
  !  ... energy calculation of neutral case
  !  ... FIXME: these lines should be called before electrons, not inside it
  !
  IF (sic .and. sic_energy) THEN
     WRITE(stdout,'(5x,"Energy calculation for the neutral polaron")')
     IF( .not. sic_first) THEN
        CALL occ_f2fn()
        CALL potinit()
        CALL wfcinit()
     END IF
     CALL electrons_scf(printout,exxen)
     CALL save_rhon(rho)
     CALL close_files(.true.)
     CALL occ_fn2f()
     CALL potinit()
     CALL wfcinit()
     sic_first = .false.
  END IF
  !
  DO idum=1,nexxiter
     !
     iter = iter + 1
     !
     ! ... Self-consistency loop. For hybrid functionals the exchange potential
     ! ... is calculated with the orbitals at previous step (none at first step)
     !
     CALL electrons_scf ( printout, exxen )
     !
     IF ( .NOT. xclib_dft_is('hybrid') ) RETURN
     !
     ! ... From now on: hybrid DFT only
     !
     IF ( stopped_by_user .OR. .NOT. conv_elec ) THEN
        conv_elec=.FALSE.
        IF ( .NOT. first) THEN
           WRITE(stdout,'(5x,"Calculation (EXX) stopped during iteration #", &
                        & i6)') iter
           CALL seqopn (iunres, 'restart_e', 'formatted', exst)
           WRITE (iunres, *) iter-1, tr2, dexx
           WRITE (iunres, *) exxen, fock0, fock1, fock2
           WRITE (iunres, *) (wg(1:nbnd,ik),ik=1,nks)
           WRITE (iunres, *) (et(1:nbnd,ik),ik=1,nks)
           CLOSE (unit=iunres, status='keep')
           CALL seqopn (iunres, 'restart_exx', 'unformatted', exst)
           WRITE (iunres) exxbuff
           CLOSE (unit=iunres, status='keep')
        ENDIF
        RETURN
     ENDIF
     !
     first =  first .AND. .NOT. exx_is_active ( )
     !
     ! "first" is true if the scf step was performed without exact exchange
     !
     IF ( first ) THEN
        !
        first = .false.
        !
        ! Activate exact exchange, set orbitals used in its calculation,
        ! then calculate exchange energy (will be useful at next step)
        !
        CALL exxinit(DoLoc)
        IF( DoLoc.and.gamma_only) THEN
          CALL localize_orbitals( )
        ELSE IF (DoLoc) THEN
          CALL localize_orbitals_k( )
        ENDIF 
        IF (use_ace) THEN
           CALL aceinit ( DoLoc ) 
           fock2 = exxenergyace()
        ELSE
           fock2 = exxenergy2()
        ENDIF
        exxen = 0.50d0*fock2 
        etot = etot - etxc 
        !
        ! Recalculate potential because XC functional has changed,
        ! start self-consistency loop on exchange
        !
        CALL v_of_rho( rho, rho_core, rhog_core, &
             ehart, etxc, vtxc, eth, etotefield, charge, v)
        etot = etot + etxc + exxen
        !
        IF (lrism) CALL rism_calc3d(rho%of_g(:, 1), esol, vsol, v%of_r, tr2)
        !
        IF (okpaw) CALL PAW_potential(rho%bec, ddd_paw, epaw,etot_cmp_paw)
        CALL set_vrs( vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, &
             nspin, doublegrid )
        !
     ELSE
        !
        ! fock1 is the exchange energy calculated for orbitals at step n,
        !       using orbitals at step n-1 in the expression of exchange
        !
        IF (use_ace) THEN
           fock1 = exxenergyace()
        ELSE
           fock1 = exxenergy2()
        ENDIF
        !
        ! Set new orbitals for the calculation of the exchange term
        !
        CALL exxinit(DoLoc)
        IF( DoLoc.and.gamma_only) THEN
          CALL localize_orbitals( )
        ELSE IF (DoLoc) THEN
          CALL localize_orbitals_k( )
        ENDIF 
        IF (use_ace) CALL aceinit ( DoLoc, fock3 )
        !
        ! fock2 is the exchange energy calculated for orbitals at step n,
        !       using orbitals at step n in the expression of exchange 
        ! fock0 is fock2 at previous step
        !
        fock0 = fock2
        IF (use_ace) THEN
           fock2 = exxenergyace()
        ELSE
           fock2 = exxenergy2()
        ENDIF
        !
        ! check for convergence. dexx is positive definite: if it isn't,
        ! there is some numerical problem. One such cause could be that
        ! the treatment of the divergence in exact exchange has failed. 
        ! FIXME: to be properly implemented for all cases
        !
!civn 
!       IF (use_ace .AND. (nspin == 1) .AND. gamma_only) THEN
        IF ( DoLoc ) THEN
          dexx =  0.5D0 * ((fock1-fock0)+(fock3-fock2)) 
        ELSE
          dexx = fock1 - 0.5D0*(fock0+fock2)
        ENDIF 
        !
        IF ( dexx < 0.0_dp ) THEN
           !IF( Doloc ) THEN
              WRITE(stdout,'(5x,a,1e12.3)') "BEWARE: negative dexx:", dexx
              dexx = ABS ( dexx )
           !ELSE
           !   CALL errore( 'electrons', 'dexx is negative! &
           !        & Check that exxdiv_treatment is appropriate for the system,&
           !        & or ecutfock may be too low', 1 )
           !ENDIF
        ENDIF
        !
        !   remove the estimate exchange energy exxen used in the inner SCF
        !
        etot = etot + exxen + 0.5D0*fock2 - fock1
        hwf_energy = hwf_energy + exxen + 0.5D0*fock2 - fock1 ! [LP]
        exxen = 0.5D0*fock2 
        !
        IF ( dexx < tr2_final ) THEN
           WRITE( stdout, 9066 ) '!!', etot, hwf_energy
        ELSE
           WRITE( stdout, 9066 ) '  ', etot, hwf_energy
        ENDIF
        IF ( dexx>1.d-8 ) THEN
          WRITE( stdout, 9067 ) dexx
        ELSE
          WRITE( stdout, 9068 ) dexx
        ENDIF
        
        WRITE( stdout, 9062 ) - fock1
        IF (use_ace) THEN
           WRITE( stdout, 9063 ) 0.5D0*fock2
        ELSE
           WRITE( stdout, 9064 ) 0.5D0*fock2
        ENDIF
        !
        IF ( dexx < tr2_final ) THEN
           IF ( do_makov_payne ) CALL makov_payne( etot, tau, &
                rho%of_r(:,1), rho%of_g(:,1), strf, vloc, gamma_only,&
                etot_in_hartree = .false., output_in_hartree = .false., &
                vacuum_level = .true. )
           WRITE( stdout, 9101 )
           RETURN
        ENDIF
        !
        IF ( adapt_thr ) THEN
           tr2 = MAX(tr2_multi * dexx, tr2_final)
           WRITE( stdout, 9121 ) tr2
        ENDIF
     ENDIF
     !
     WRITE( stdout,'(/5x,"EXX: now go back to refine exchange calculation")')
     !
     IF ( check_stop_now() .or. (iter.ge.nexxiter) ) THEN
        WRITE(stdout,'(5x,"Calculation (EXX) stopped after iteration #", &
                        & i6)') iter
        conv_elec=.FALSE.
        IF(iter.ge.nexxiter) conv_elec=.TRUE. ! it will print ace and wfc files for restart
        CALL seqopn (iunres, 'restart_e', 'formatted', exst)
        WRITE (iunres, *) iter, tr2, dexx
        WRITE (iunres, *) exxen, fock0, fock1, fock2
        ! FIXME: et and wg are written to xml file
        WRITE (iunres, *) (wg(1:nbnd,ik),ik=1,nks)
        WRITE (iunres, *) (et(1:nbnd,ik),ik=1,nks)
        CLOSE (unit=iunres, status='keep')
        RETURN
     ENDIF
     !
  ENDDO
  !
  WRITE( stdout, 9120 ) iter
  FLUSH( stdout )
  !
  RETURN
  !
  ! ... formats
  !
9062 FORMAT( '     - averaged Fock potential =',0PF17.8,' Ry' )
9063 FORMAT( '     + Fock energy (ACE)       =',0PF17.8,' Ry' )
9064 FORMAT( '     + Fock energy (full)      =',0PF17.8,' Ry' )
9066 FORMAT(/,A2,'   total energy              =',0PF17.8,' Ry' &
            /'     Harris-Foulkes estimate   =',0PF17.8,' Ry' )
9067 FORMAT('     est. exchange err (dexx)  =',0PF17.8,' Ry' )
9068 FORMAT('     est. exchange err (dexx)  =',1PE17.1,' Ry' )
9101 FORMAT(/'     EXX self-consistency reached' )
9120 FORMAT(/'     EXX convergence NOT achieved after ',i3,' iterations: stopping' )
9121 FORMAT(/'     scf convergence threshold =',1PE17.1,' Ry' )
  !
END SUBROUTINE electrons
!
!----------------------------------------------------------------------------
SUBROUTINE electrons_scf ( printout, exxen )
  !----------------------------------------------------------------------------
  !! This routine is a driver of the self-consistent cycle.  
  !! It uses the routine \(\texttt{c_bands}\) for computing the bands at fixed
  !! Hamiltonian, the routine \(\texttt{sum_band}\) to compute the charge density,
  !! the routine \(\texttt{v_of_rho}\) to compute the new potential and the routine
  !! \(\text{mix_rho}\) to mix input and output charge densities.
  !
  USE kinds,                ONLY : DP
  USE check_stop,           ONLY : check_stop_now, stopped_by_user
  USE io_global,            ONLY : stdout, ionode
  USE cell_base,            ONLY : at, bg, alat, omega, tpiba2
  USE ions_base,            ONLY : zv, nat, nsp, ityp, tau, compute_eextfor, atm, &
                                   ntyp => nsp
  USE starting_scf,         ONLY : starting_pot
  USE bp,                   ONLY : lelfield
  USE fft_base,             ONLY : dfftp
  USE gvect,                ONLY : ngm, gstart, g, gg, gcutm
  USE gvecs,                ONLY : doublegrid, ngms
  USE klist,                ONLY : xk, wk, nelec, ngk, nks, nkstot, lgauss, &
                                   two_fermi_energies, tot_charge
  USE fixed_occ,            ONLY : one_atom_occupations
  USE lsda_mod,             ONLY : lsda, nspin, magtot, absmag, isk
  USE vlocal,               ONLY : strf, vloc
  USE wvfct,                ONLY : nbnd, et
  USE gvecw,                ONLY : ecutwfc
  USE ener,                 ONLY : etot, hwf_energy, eband, deband, ehart, &
                                   vtxc, etxc, etxcc, ewld, demet, epaw, &
                                   elondon, edftd3, ef_up, ef_dw, exdm, ef, &
                                   egrand, vsol, esol, esic, esci
  USE scf,                  ONLY : scf_type, scf_type_COPY, bcast_scf_type,&
                                   create_scf_type, destroy_scf_type, &
                                   open_mix_file, close_mix_file, &
                                   rho, rho_core, rhog_core, v, vltot, vrs, &
                                   kedtau, vnew
  USE control_flags,        ONLY : mixing_beta, tr2, ethr, niter, nmix, &
                                   conv_elec, sic, &
                                   restart, io_level, do_makov_payne,  &
                                   gamma_only, iverbosity, textfor,     &
                                   llondon, ldftd3, scf_must_converge, lxdm, ts_vdw, &
                                   mbd_vdw, use_gpu
  USE control_flags,        ONLY : n_scf_steps, scf_error, scissor, gamma_only
  USE sci_mod,              ONLY : sci_iter

  USE io_files,             ONLY : iunmix, output_drho
  USE ldaU,                 ONLY : eth, lda_plus_u, lda_plus_u_kind, &
                                   niter_with_fixed_ns, hub_pot_fix, &
                                   nsg, nsgnew, v_nsg, at_sc, neighood, &
                                   ldim_u, is_hubbard_back, apply_U, orbital_resolved
  USE extfield,             ONLY : tefield, etotefield, gate, etotgatefield !TB
  USE noncollin_module,     ONLY : noncolin, magtot_nc, i_cons,  bfield, &
                                   lambda, report, domag, nspin_mag, npol
  USE io_rho_xml,           ONLY : write_scf
  USE uspp,                 ONLY : okvan
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp_pools,             ONLY : root_pool, me_pool, my_pool_id, &
                                   inter_pool_comm, intra_pool_comm
  USE mp,                   ONLY : mp_sum, mp_bcast
  !
  USE london_module,        ONLY : energy_london
  USE dftd3_api,            ONLY : dftd3_pbc_dispersion, get_atomic_number
  USE dftd3_qe,             ONLY : dftd3
  USE xdm_module,           ONLY : energy_xdm
  USE tsvdw_module,         ONLY : EtsvdW
  !
  USE paw_variables,        ONLY : okpaw, ddd_paw, total_core_energy, only_paw
  USE paw_onecenter,        ONLY : PAW_potential
  USE paw_symmetry,         ONLY : PAW_symmetrize_ddd
  USE dfunct,               ONLY : newd
  USE esm,                  ONLY : do_comp_esm, esm_printpot, esm_ewald
  USE gcscf_module,         ONLY : lgcscf, gcscf_mu, gcscf_ignore_mun, gcscf_set_nelec
  USE clib_wrappers,        ONLY : memstat
  USE fcp_module,           ONLY : lfcp, fcp_mu
  USE rism_module,          ONLY : lrism, rism_calc3d, rism_printpot
  USE iso_c_binding,        ONLY : c_int
  !
  USE plugin_variables,     ONLY : plugin_etot
  USE libmbd_interface,     ONLY : EmbdvdW
  USE add_dmft_occ,         ONLY : dmft, dmft_update, v_dmft, dmft_updated
  !
  USE device_fbuff_m,       ONLY : dev_buf, pin_buf
  USE pwcom,                ONLY : report_mag 
  USE makovpayne,           ONLY : makov_payne
  !
#if defined (__ENVIRON)
  USE plugin_flags,         ONLY : use_environ
  USE environ_base_module,  ONLY : calc_environ_energy, print_environ_energies
  USE environ_pw_module,    ONLY : calc_environ_potential
#endif
#if defined (__OSCDFT)
  USE plugin_flags,         ONLY : use_oscdft
  USE oscdft_base,          ONLY : oscdft_ctx
  USE oscdft_functions,     ONLY : oscdft_electrons,&
                                   oscdft_scf_energy,&
                                   oscdft_print_energies,&
                                   oscdft_print_ns
#endif
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT (IN) :: printout
  !! * If printout>0, prints on output the total energy;
  !! * if printout>1, also prints decomposition into energy contributions.
  REAL(DP),INTENT (IN) :: exxen
  !! current estimate of the exchange energy
  !
  ! ... local variables
  !
  REAL(DP) :: dr2
  !! the norm of the diffence between potential
  REAL(DP) :: charge
  !! the total charge
  REAL(DP) :: deband_hwf
  !! deband for the Harris-Weinert-Foulkes functional
  REAL(DP) :: mag
  !! local magnetization
  INTEGER :: i
  !! counter on polarization
  INTEGER :: idum
  !! dummy counter on iterations
  INTEGER :: iter
  !! counter on iterations
  INTEGER :: nt
  !! counter on atomic types
  INTEGER :: ios, kilobytes
  !
  REAL(DP) :: tr2_min
  !! estimated error on energy coming from diagonalization
  REAL(DP) :: descf
  !! correction for variational energy
  REAL(DP) :: en_el=0.0_DP
  !! electric field contribution to the total energy
  REAL(DP) :: eext=0.0_DP
  !! external forces contribution to the total energy
  LOGICAL :: first, exst
  !! auxiliary variables for calculating and storing temporary copies of
  !! the charge density and of the HXC-potential
  !
  TYPE(scf_type) :: rhoin
  !! used to store rho_in of current/next iteration
  !
  ! ... external functions
  !
  REAL(DP), EXTERNAL :: ewald, get_clock
  REAL(DP) :: etot_cmp_paw(nat,2,2)
  ! 
  REAL(DP), ALLOCATABLE :: taupbc(:,:)
  !! atomic positions centered around r=0 - for grimme-d3
  INTEGER:: atnum(1:nat), na
  !! auxiliary variables for grimme-d3
  LOGICAL :: lhb
  !! if .TRUE. then background states are present (DFT+U)
  !
  lhb = .FALSE.
  IF ( lda_plus_u )  THEN
     DO nt = 1, ntyp
        IF (is_hubbard_back(nt)) lhb = .TRUE.
     ENDDO
  ENDIF
  !
  iter = 0
  dr2  = 0.0_dp
  IF ( restart ) THEN
     CALL restart_in_electrons( iter, dr2, ethr, et )
     !$acc update device (et)
  END IF
  !
  WRITE( stdout, 9000 ) get_clock( 'PWSCF' )
  !
  CALL memstat( kilobytes )
  IF ( kilobytes > 0 ) WRITE( stdout, 9001 ) kilobytes/1000.0
  !
  CALL start_clock( 'electrons' )
  !
  FLUSH( stdout )
  !
  ! ... calculates the ewald contribution to total energy
  !
  IF ( do_comp_esm ) THEN
     ewld = esm_ewald()
  ELSE
     ewld = ewald( alat, nat, nsp, ityp, zv, at, bg, tau, &
                omega, g, gg, ngm, gcutm, gstart, gamma_only, strf )
  ENDIF
  !
  IF ( llondon ) THEN
     elondon = energy_london( alat , nat , ityp , at ,bg , tau )
  ELSE
     elondon = 0.d0
  ENDIF
  !
  ! Grimme-D3 correction to the energy
  !
  IF (ldftd3) THEN
     CALL start_clock('energy_dftd3')
     ! taupbc are atomic positions in alat units, centered around r=0
     ALLOCATE ( taupbc(3,nat) )
     taupbc(:,:) = tau(:,:)
     CALL cryst_to_cart( nat, taupbc, bg, -1 ) 
     taupbc(:,:) = taupbc(:,:) - NINT(taupbc(:,:))
     CALL cryst_to_cart( nat, taupbc, at,  1 ) 
     DO na = 1, nat
        atnum(na) = get_atomic_number(TRIM(atm(ityp(na))))
     ENDDO
     call dftd3_pbc_dispersion(dftd3, alat*taupbc, atnum, alat*at, edftd3)
     edftd3=edftd3*2.d0
     DEALLOCATE( taupbc)
     CALL stop_clock('energy_dftd3')
  ELSE
     edftd3= 0.0
  ENDIF
  !
  !
  CALL create_scf_type( rhoin )
  !
  WRITE( stdout, 9002 )
  FLUSH( stdout )
  !
  CALL open_mix_file( iunmix, 'mix', exst )
  !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%          iterate !          %%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
  DO idum = 1, niter
     !
     IF ( check_stop_now() ) THEN
        conv_elec=.FALSE.
        CALL save_in_electrons (iter, dr2, ethr, et )
        GO TO 10
     ENDIF
     iter = iter + 1
     IF(scissor) sci_iter = iter
     !
     IF (use_gpu .and. (iverbosity >= 1) ) CALL dev_buf%print_report(stdout)
     IF (use_gpu .and. (iverbosity >= 1) ) CALL pin_buf%print_report(stdout)
     !
     WRITE( stdout, 9010 ) iter, ecutwfc, mixing_beta
     !
     FLUSH( stdout )
     !
     ! ... Convergence threshold for iterative diagonalization is
     ! ... automatically updated during self consistency
     !
     IF ( iter > 1 ) THEN
        !
#if defined (__OSCDFT)
        IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==2)) THEN
           ! do nothing
        ELSE
#endif
        IF ( iter == 2 ) THEN
           !
           IF ( lgcscf ) THEN
              ethr = 1.D-5
           ELSE
              ethr = 1.D-2
           ENDIF
           !
        ENDIF
#if defined (__OSCDFT)
        ENDIF
#endif
        !
        IF ( orbital_resolved ) THEN
           ethr = MIN( ethr, 0.05D0*dr2 / MAX( 1.D0, nelec ) )
           ! ... tighten diagonalization in orbital-resolved DFT+U
           ! ... for faster convergence
        ELSE
           !
           ethr = MIN( ethr, 0.1D0*dr2 / MAX( 1.D0, nelec ) )
        ENDIF
        ! ... do not allow convergence threshold to become too small:
        ! ... iterative diagonalization may become unstable
        ethr = MAX( ethr, 1.D-13 )
        !
     ENDIF
     !
     first = ( iter == 1 )
     !
     ! ... deband = - \sum_v <\psi_v | V_h + V_xc |\psi_v> is calculated a
     ! ... first time here using the input density and potential ( to be
     ! ... used to calculate the Harris-Weinert-Foulkes energy )
     !
     deband_hwf = delta_e()
     !
     ! ... save input density to rhoin
     !
     call scf_type_COPY( rho, rhoin )
     !
     scf_step: DO
        !
        ! ... tr2_min is set to an estimate of the error on the energy
        ! ... due to diagonalization - used only for the first scf iteration
        !
        tr2_min = 0.D0
        !
        IF ( first ) tr2_min = ethr*MAX( 1.D0, nelec ) 
        !
        ! ... diagonalization of the KS hamiltonian
        !
#if defined (__OSCDFT)
        IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==1)) THEN
           CALL oscdft_electrons(oscdft_ctx, iter, et, nbnd, nkstot, nks)
        ELSE
#endif
        IF ( lelfield ) THEN
           CALL c_bands_efield( iter )
        ELSE
           CALL c_bands( iter )
        ENDIF
#if defined (__OSCDFT)
        END IF
#endif
        !
        IF ( stopped_by_user ) THEN
           conv_elec=.FALSE.
           CALL save_in_electrons( iter-1, dr2, ethr, et )
           GO TO 10
        ENDIF
        !
        IF (one_atom_occupations) CALL new_evc()
        !
        ! ... xk, wk, isk, et, wg are distributed across pools;
        ! ... the first node has a complete copy of xk, wk, isk,
        ! ... while eigenvalues et and weights wg must be
        ! ... explicitly collected to the first node
        ! ... this is done here for et, in sum_band for wg
        !
        CALL poolrecover( et, nbnd, nkstot, nks )
        !
        ! ... the new density is computed here. For PAW:
        ! ... sum_band computes new becsum (stored in uspp modules)
        ! ... and a subtly different copy in rho%bec (scf module)
        !
        ! ... for charge self-consistent calculations call update
        ! ... of Fermi weights from DMFT here
        !
        IF ( dmft .AND. idum == 1) THEN
            CALL weights()
            CALL dmft_update()
        ELSEIF ( dmft .AND. idum > 1) THEN
            WRITE( stdout, '(5X,"WARNING: electron_maxstep > 1 not recommended for dmft = .true.")')
        END IF
        !
        CALL sum_band()
        !
        ! ... if DMFT update was made, make sure it is only done in the first iteration
        ! ... (generally in this mode it should only run a single iteration, but just to make sure!)
        !
        IF ( dmft .AND. idum == 1) THEN
            dmft_updated = .TRUE.
            DEALLOCATE(v_dmft)
        END IF
        !
        ! ... the Harris-Weinert-Foulkes energy is computed here using only
        ! ... quantities obtained from the input density
        !
        hwf_energy = eband + deband_hwf + (etxc - etxcc) + ewld + ehart + demet
        If ( okpaw ) hwf_energy = hwf_energy + epaw
        IF ( lda_plus_u ) hwf_energy = hwf_energy + eth
        IF ( lrism ) hwf_energy = hwf_energy + esol + vsol
        IF ( sic ) hwf_energy = hwf_energy + esic
        IF ( scissor ) hwf_energy = hwf_energy + esci
        !
        IF ( lda_plus_u )  THEN
           !
           ! Write the occupation matrices
           !
           IF ( iverbosity > 0 .OR. first ) THEN
              IF (lda_plus_u_kind.EQ.0) THEN
                IF (noncolin) THEN
                    CALL write_ns_nc()
                ELSE        
                    CALL write_ns()
                ENDIF 
              ELSEIF (lda_plus_u_kind.EQ.1) THEN
                 IF (noncolin) THEN
                    CALL write_ns_nc()
                 ELSE
                    CALL write_ns()
                 ENDIF
              ELSEIF (lda_plus_u_kind.EQ.2) THEN
                 IF (noncolin) THEN
                    CALL write_nsg_nc()   
                 ELSE 
                    CALL write_nsg()
                 ENDIF
              ENDIF
           ENDIF
           !
           ! Keep the Hubbard potential fixed, i.e. keep the 
           ! occupation matrix equal to the ground-state one.
           ! This is needed for the calculation of Hubbard parameters
           ! in a self-consistent way.
           !
           IF (hub_pot_fix) THEN
             IF (lda_plus_u_kind.EQ.0) THEN
                IF (noncolin) THEN
                   ! call occupation counting routine before resetting ns
                   IF (orbital_resolved) CALL alpha_m_nc_trace(rho%ns_nc) 
                   rho%ns_nc = rhoin%ns_nc
                ELSE
                   ! call occupation counting routine before resetting ns
                   IF (orbital_resolved) CALL alpha_m_trace(rho%ns) 
                   rho%ns = rhoin%ns ! back to input values
                ENDIF
                !
                IF (lhb) rho%nsb = rhoin%nsb
             ELSEIF (lda_plus_u_kind.EQ.1) THEN
                CALL errore('electrons_scf', &
                  & 'hub_pot_fix is not implemented for lda_plus_u_kind=1',1)
             ELSEIF (lda_plus_u_kind.EQ.2) THEN
                IF (noncolin) CALL errore('electrons_scf', &
                & 'hub_pot_fix is not implemented for (lda_plus_u_kind=2 .AND. noncolin)',1)
                nsgnew = nsg
             ENDIF
           ENDIF
           !
           IF ( first .AND. starting_pot == 'atomic' ) THEN
              IF (lda_plus_u_kind.EQ.0) THEN
                 CALL ns_adj()     
                 IF (noncolin) THEN
                    rhoin%ns_nc = rho%ns_nc                    
                 ELSE
                    rhoin%ns = rho%ns
                    IF (lhb) rhoin%nsb = rho%nsb
                 ENDIF   
              ELSEIF (lda_plus_u_kind.EQ.1) THEN
                 CALL ns_adj()
                 IF (noncolin) THEN
                    rhoin%ns_nc = rho%ns_nc
                 ELSE
                    rhoin%ns = rho%ns
                 ENDIF
              ELSEIF (lda_plus_u_kind.EQ.2) THEN
                 CALL nsg_adj()
              ENDIF
           ENDIF
           IF ( iter <= niter_with_fixed_ns ) THEN
              WRITE( stdout, '(/,5X,"RESET ns to initial values (iter <= mixing_fixed_ns)",/)')
              IF (lda_plus_u_kind.EQ.0) THEN
                 IF (noncolin) THEN   
                    rho%ns_nc = rhoin%ns_nc  
                 ELSE        
                    rho%ns = rhoin%ns
                    IF (lhb) rhoin%nsb = rho%nsb
                 ENDIF
              ELSEIF (lda_plus_u_kind.EQ.1) THEN
                 IF (noncolin) THEN
                    rho%ns_nc = rhoin%ns_nc
                 ELSE
                    rho%ns = rhoin%ns
                 ENDIF
              ELSEIF (lda_plus_u_kind.EQ.2) THEN
                 nsgnew = nsg
              ENDIF
           ENDIF
           !
        ENDIF
        !
        ! ... calculate total and absolute magnetization
        !
        IF ( lsda .OR. noncolin ) CALL compute_magnetization()
        !
        ! ... eband  = \sum_v \epsilon_v    is calculated by sum_band
        ! ... deband = - \sum_v <\psi_v | V_h + V_xc |\psi_v>
        ! ... eband + deband = \sum_v <\psi_v | T + Vion |\psi_v>
        !
        deband = delta_e()
        !
        ! ... mix_rho mixes several quantities: rho in g-space, tauk (for
        ! ... meta-gga), ns and ns_nc (for lda+u) and becsum (for paw)
        ! ... The mixing is done on pool 0 only (image parallelization
        ! ... inside mix_rho => rho_ddot => PAW_ddot is no longer there)
        !
        IF ( my_pool_id == root_pool ) CALL mix_rho( rho, rhoin, &
                mixing_beta, dr2, tr2_min, iter, nmix, iunmix, conv_elec )
        !
        ! ... Results are broadcast from pool 0 to others to prevent trouble
        ! ... on machines unable to yield the same results for the same 
        ! ... calculations on the same data, performed on different procs
        !
        IF ( lda_plus_u )  THEN
           ! ... For DFT+U, ns and ns_nc are also broadcast inside each pool
           ! ... to ensure consistency on all processors of all pools
           IF (noncolin) THEN
              IF (ALLOCATED(rhoin%ns_nc)) CALL mp_bcast( rhoin%ns_nc, root_pool, intra_pool_comm )
           ELSE
              IF (ALLOCATED(rhoin%ns)) CALL mp_bcast( rhoin%ns, root_pool, intra_pool_comm )
           ENDIF
           ! DFT+U+V: this variable is not in "mix-type" variable rhoin
           IF (lda_plus_u_kind.EQ.2) THEN
              IF (ALLOCATED(nsg) ) CALL mp_bcast ( nsg, root_pool, inter_pool_comm)
           ENDIF
        ENDIF
        !
        CALL bcast_scf_type( rhoin, root_pool, inter_pool_comm )
        CALL mp_bcast( dr2, root_pool, inter_pool_comm )
        CALL mp_bcast( conv_elec, root_pool, inter_pool_comm )
        !
        IF (.NOT. scf_must_converge .AND. idum == niter) conv_elec = .TRUE.
        !
        ! ... if convergence is achieved or if the self-consistency error
        ! ... (dr2) is smaller than the estimated error due to diagonalization
        ! ... (tr2_min), rhoin and rho are unchanged: rhoin contains the input
        ! ... density and rho contains the output density.
        ! ... In all other cases, rhoin contains the mixed charge density 
        ! ... (the new input density) while rho is unchanged
        !
        IF ( first .and. nat > 0) THEN
           !
           ! ... first scf iteration: check if the threshold on diagonalization
           ! ... (ethr) was small enough wrt the error in self-consistency (dr2)
           ! ... if not, perform a new diagonalization with reduced threshold
           !
           first = .FALSE.
           !
           IF ( dr2 < tr2_min ) THEN
              !
              WRITE( stdout, '(/,5X,"Threshold (ethr) on eigenvalues was ", &
                               &    "too large:",/,5X,                      &
                               & "Diagonalizing with lowered threshold",/)' )
              !
              ethr = MAX( 0.1D0*dr2 / MAX( 1.D0, nelec ), 1.D-13 )
              !
              CYCLE scf_step
              !
           ENDIF
           !
        ENDIF
        !
        IF ( .NOT. conv_elec ) THEN
           !
           ! ... no convergence yet: calculate new potential from mixed
           ! ... charge density (i.e. the new estimate)
           !
           IF ( orbital_resolved .AND. (.NOT. apply_U) ) THEN
              !
              IF ( ethr .LE. 1.D-4 .OR. iter .GT. 10 ) THEN
                 !
                 ! ... activate the orbital-resolved Hubbard corrections
                 ! ... once the eigenstates (and thus, the eigenvectors) 
                 ! ... have stablilized. Also turn them on after 10 iterations
                 ! ... but print a warning message in this case.
                 WRITE( stdout, '(/,5X,47("="))')
                 WRITE( stdout, '(/,5X,"Switching ON", &
                    & " orbital-resolved Hubbard corrections")' )
                 IF ( iter .GT. 10 ) &
                    WRITE( stdout, '(/,5X,"WARNING: check convergence of eigenstates")')
                 WRITE( stdout, '(/,5X,47("="))')
                 !
                 apply_U = .TRUE.
              ENDIF
           ENDIF
           !
           !
           CALL v_of_rho( rhoin, rho_core, rhog_core, &
                          ehart, etxc, vtxc, eth, etotefield, charge, v )
           !
           IF (lrism) THEN
              CALL rism_calc3d(rhoin%of_g(:, 1), esol, vsol, v%of_r, dr2)
           ENDIF
           !
           IF (okpaw) THEN
              CALL PAW_potential( rhoin%bec, ddd_paw, epaw,etot_cmp_paw )
              CALL PAW_symmetrize_ddd( ddd_paw )
           ENDIF
           !
           ! ... estimate correction needed to have variational energy:
           ! ... T + E_ion (eband + deband) are calculated in sum_band
           ! ... and delta_e using the output charge density rho;
           ! ... E_H (ehart) and E_xc (etxc) are calculated in v_of_rho
           ! ... above, using the mixed charge density rhoin%of_r.
           ! ... delta_escf corrects for this difference at first order
           !
           descf = delta_escf()
           !
           ! ... now copy the mixed charge density in R- and G-space in rho
           !
           CALL scf_type_COPY( rhoin, rho )
           !
#if defined (__OSCDFT)
           IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==2)) THEN
              IF (lda_plus_u .AND. .NOT.oscdft_ctx%conv) THEN
                 IF (lda_plus_u_kind.EQ.0) THEN
                    CALL write_ns()
                 ELSEIF (lda_plus_u_kind.EQ.2) THEN
                    CALL write_nsg()
                 ENDIF
              ENDIF
           ENDIF
#endif
           !
           IF (lda_plus_u .AND. lda_plus_u_kind.EQ.2) nsgnew = nsg
           !
           IF ( lgcscf ) THEN
              CALL gcscf_set_nelec( charge )
           END IF
           !
        ELSE 
           !
           ! ... convergence reached:
           ! ... 1) the output HXC-potential is saved in v
           ! ... 2) vnew contains V(out)-V(in) ( used to correct the forces ).
           !
           vnew%of_r(:,:) = v%of_r(:,:)
           !
           IF (lda_plus_u .AND. lda_plus_u_kind.EQ.2) nsg = nsgnew
           !
           CALL v_of_rho( rho,rho_core,rhog_core, &
                          ehart, etxc, vtxc, eth, etotefield, charge, v )
           !
           IF (lrism) THEN
              CALL rism_calc3d(rho%of_g(:, 1), esol, vsol, v%of_r, tr2)
           END IF
           !
           vnew%of_r(:,:) = v%of_r(:,:) - vnew%of_r(:,:)
           !
           IF (okpaw) THEN
              CALL PAW_potential( rho%bec, ddd_paw, epaw, etot_cmp_paw )
              CALL PAW_symmetrize_ddd( ddd_paw )
           ENDIF
           !
           ! ... note that rho is here the output, not mixed, charge density
           ! ... so correction for variational energy is no longer needed
           !
           descf = 0._dp
           !
        ENDIF 
        !
        ! ... if we didn't cycle before we can exit the do-loop
        !
        EXIT scf_step
        !
     ENDDO scf_step
     !
     plugin_etot = 0.0_dp
     !
#if defined (__LEGACY_PLUGINS) 
     CALL plugin_scf_energy (plugin_etot, rhoin) 
     !
     CALL plugin_scf_potential(rhoin, conv_elec, dr2, vltot)
#endif 
#if defined (__ENVIRON)
     IF (use_environ) THEN
        CALL calc_environ_energy(plugin_etot, .TRUE.)
        CALL calc_environ_potential(rhoin, conv_elec, dr2, vltot)
     END IF
#endif
#if defined (__OSCDFT)
     IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==1)) THEN
        CALL oscdft_scf_energy(oscdft_ctx, plugin_etot)
     END IF
#endif
     !
     ! ... define the total local potential (external + scf)
     !
     CALL set_vrs( vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, &
                   nspin, doublegrid )
     !
     ! ... in the US case we have to recompute the self-consistent
     ! ... term in the nonlocal potential
     ! ... PAW: newd contains PAW updates of NL coefficients
     !
     CALL newd()
     !
     IF ( lelfield ) en_el =  calc_pol ( )
     !
     IF ( report > 0 ) THEN
        IF ( conv_elec .OR.  MOD(iter,report) == 0 ) CALL report_mag()
     ELSE IF ( report < 0 ) THEN
        IF ( conv_elec ) CALL report_mag(SAVE_LOCALS=.TRUE.)
     END IF
     !
     WRITE( stdout, 9000 ) get_clock( 'PWSCF' )
     !
     IF ( conv_elec ) WRITE( stdout, 9101 )
 
     IF ( conv_elec ) THEN 
           scf_error = dr2
           n_scf_steps = iter
     ENDIF  

     !
     IF ( conv_elec .OR. dmft_updated ) THEN
        !
        ! iverbosity == 0 for the PW code
        ! iverbosity >  2 for the HP code
        !
        IF ( lda_plus_u .AND. (iverbosity == 0 .OR. iverbosity > 2) ) THEN
           !
           ! Recompute the occupation matrix:
           ! needed when computing U (and V) in a SCF way
           IF (hub_pot_fix) THEN
              IF (lda_plus_u_kind.EQ.0) THEN
                 IF (noncolin) THEN
                    CALL new_ns_nc(rho%ns_nc)
                 ELSE
                    CALL new_ns(rho%ns)
                 ENDIF
                 IF (lhb) CALL new_nsb(rho%nsb)
              ELSEIF (lda_plus_u_kind.EQ.2) THEN
                 CALL new_nsg()
              ENDIF
           ENDIF
           !
           ! Write the occupation matrices
           IF (lda_plus_u_kind == 0) THEN
              IF (noncolin) THEN
                 IF ( orbital_resolved .AND. hub_pot_fix ) &
                    CALL alpha_m_nc_trace(rho%ns_nc)      
                 CALL write_ns_nc()
              ELSE
                 IF ( orbital_resolved .AND. hub_pot_fix ) &
                    CALL alpha_m_trace(rho%ns)
                 CALL write_ns()
              ENDIF
           ELSEIF (lda_plus_u_kind == 1) THEN
              IF (noncolin) THEN
                 CALL write_ns_nc()
              ELSE
                 CALL write_ns()
              ENDIF
           ELSEIF (lda_plus_u_kind == 2) THEN
              IF (noncolin) THEN 
                 CALL write_nsg_nc()
              ELSE
                 CALL write_nsg()
              ENDIF
           ENDIF
           !
        ENDIF
#if defined (__OSCDFT)
        IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==1)) THEN
           CALL oscdft_print_ns(oscdft_ctx)
        END IF
#endif
        CALL print_ks_energies()
        !
     ENDIF
     !
     IF ( ABS( charge - nelec ) / charge > 1.D-7 ) THEN
        WRITE( stdout, 9050 ) charge, nelec
        IF ( ABS( charge - nelec ) / charge > 1.D-3 ) THEN
           IF (.not.lgauss) THEN
              CALL errore( 'electrons', 'charge is wrong: smearing is needed', 1 )
           ELSE
              CALL errore( 'electrons', 'charge is wrong', 1 )
           ENDIF
        ENDIF
     ENDIF
     !
     etot = eband + ( etxc - etxcc ) + ewld + ehart + deband + demet + descf
     ! for hybrid calculations, add the current estimate of exchange energy
     ! (it will subtracted later if exx_is_active to be replaced with a better estimate)
     etot = etot - exxen
     hwf_energy = hwf_energy - exxen ! [LP]
     !
     IF (okpaw) etot = etot + epaw
     IF ( lda_plus_u ) etot = etot + eth
     !
     IF ( lelfield ) etot = etot + en_el
     ! not sure about the HWF functional in the above case
     IF( textfor ) THEN
        eext = alat*compute_eextfor()
        etot = etot + eext
        hwf_energy = hwf_energy + eext
     ENDIF
     IF (llondon) THEN
        etot = etot + elondon
        hwf_energy = hwf_energy + elondon
     ENDIF
     !
     ! grimme-d3 dispersion energy
     IF (ldftd3) THEN
        etot = etot + edftd3
        hwf_energy = hwf_energy + edftd3
     ENDIF
     !
     ! calculate the xdm energy contribution with converged density
     IF (lxdm .and. conv_elec) THEN
        exdm = energy_xdm()  
        etot = etot + exdm
        hwf_energy = hwf_energy + exdm
     ENDIF
     IF (mbd_vdw) THEN
        ! factor 2 converts from Hartree to Ry units
        etot = etot + 2.0d0*EmbdvdW
        hwf_energy = hwf_energy + 2.0d0*EmbdvdW
     ELSE IF (ts_vdw) THEN
        ! factor 2 converts from Hartree to Ry units
        etot = etot + 2.0d0*EtsvdW
        hwf_energy = hwf_energy + 2.0d0*EtsvdW
     ENDIF
     !
     IF ( tefield ) THEN
        etot = etot + etotefield
        hwf_energy = hwf_energy + etotefield
     ENDIF
     ! TB gate energy
     IF ( gate ) THEN
        etot = etot + etotgatefield
        hwf_energy = hwf_energy + etotgatefield
     ENDIF
     !
     IF ( lgcscf .AND. (.NOT. gcscf_ignore_mun) ) THEN
        etot = etot + egrand
        hwf_energy = hwf_energy + egrand
     END IF
     !
     IF ( lfcp ) THEN
        etot = etot + fcp_mu * tot_charge
        hwf_energy = hwf_energy + fcp_mu * tot_charge
     END IF
     !
     IF ( lrism ) THEN
        etot = etot + esol + vsol
        hwf_energy = hwf_energy + esol + vsol
     END IF
     !
     IF (sic) etot = etot + esic
     IF (scissor) etot = etot + esci
     !
     ! ... adds possible external contribution from plugins to the energy
     !
     etot = etot + plugin_etot 
     !
     CALL print_energies ( printout )
     !
     IF ( conv_elec ) THEN
        !
        ! ... if system is charged add a Makov-Payne correction to the energy
        ! ... (not in case of hybrid functionals: it is added at the end)
        !
        IF ( do_makov_payne .AND. printout/= 0 ) CALL makov_payne( etot, tau, &
                rho%of_r(:,1), rho%of_g(:,1), strf, vloc, gamma_only,&
                etot_in_hartree=.false., output_in_hartree = .false., &
                vacuum_level = .true. )
        !
        ! ... print out ESM potentials if desired
        !
        IF ( do_comp_esm ) CALL esm_printpot( rho%of_g )
        !
        ! ... print out 3D-RISM potentials if desired
        !
        IF ( lrism ) CALL rism_printpot()
        !
        WRITE( stdout, 9110 ) iter
        !
        ! ... jump to the end
        !
        GO TO 10
        !
     ENDIF
     !
     ! ... uncomment the following line if you wish to monitor the evolution
     ! ... of the force calculation during self-consistency
     !
     !CALL forces()
     !
     ! ... it can be very useful to track internal clocks during
     ! ... self-consistency for benchmarking purposes
#if defined(__PW_TRACK_ELECTRON_STEPS)
     CALL print_clock_pw()
#endif
     !
  ENDDO
  n_scf_steps = iter
  scf_error = dr2
  !
  WRITE( stdout, 9101 )
  WRITE( stdout, 9120 ) iter
  !
10  FLUSH( stdout )
  !
  ! ... exiting: unless disabled by user (disk_io='minimal' or 'none'),
  ! ... write charge density (also ldaU ns coefficients and PAW becsum) 
  !
  IF ( io_level > -2 ) CALL write_scf( rho, nspin )
  !
  ! ... keep mixing info if no converged achieved, delete it otherwise
  !
  IF ( conv_elec .or. io_level < 0 ) THEN
     CALL close_mix_file( iunmix, 'delete' )
  ELSE
     CALL close_mix_file( iunmix, 'keep' )
  ENDIF
  !
  IF(scissor) sci_iter = 0
  !
  IF ( output_drho /= ' ' ) CALL remove_atomic_rho()
  call destroy_scf_type ( rhoin )
  CALL stop_clock( 'electrons' )
  !
  RETURN
  !
  ! ... formats
  !
9000 FORMAT(/'     total cpu time spent up to now is ',F10.1,' secs' )
9001 FORMAT(/'     per-process dynamical memory: ',f7.1,' Mb' )
9002 FORMAT(/'     Self-consistent Calculation' )
9010 FORMAT(/'     iteration #',I3,'     ecut=', F9.2,' Ry',5X,'beta=',F5.2 )
9050 FORMAT(/'     WARNING: integrated charge=',F15.8,', expected=',F15.8 )
9101 FORMAT(/'     End of self-consistent calculation' )
9110 FORMAT(/'     convergence has been achieved in ',i3,' iterations' )
9120 FORMAT(/'     convergence NOT achieved after ',i3,' iterations: stopping' )
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     SUBROUTINE compute_magnetization()
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       INTEGER :: ir
       !
       !
       IF ( lsda ) THEN
          !
          magtot = 0.D0
          absmag = 0.D0
          !
          DO ir = 1, dfftp%nnr
             !
             mag = rho%of_r(ir,2)
             !
             magtot = magtot + mag
             absmag = absmag + ABS( mag )
             !
          ENDDO
          !
          magtot = magtot * omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
          absmag = absmag * omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
          !
          CALL mp_sum( magtot, intra_bgrp_comm )
          CALL mp_sum( absmag, intra_bgrp_comm )
          !
          IF (two_fermi_energies.and.lgauss) bfield(3)=0.5D0*(ef_up-ef_dw)
          !
       ELSEIF ( noncolin ) THEN
          !
          magtot_nc = 0.D0
          absmag    = 0.D0
          !
          DO ir = 1,dfftp%nnr
             !
             mag = SQRT( rho%of_r(ir,2)**2 + &
                         rho%of_r(ir,3)**2 + &
                         rho%of_r(ir,4)**2 )
             !
             DO i = 1, 3
                !
                magtot_nc(i) = magtot_nc(i) + rho%of_r(ir,i+1)
                !
             ENDDO
             !
             absmag = absmag + ABS( mag )
             !
          ENDDO
          !
          CALL mp_sum( magtot_nc, intra_bgrp_comm )
          CALL mp_sum( absmag, intra_bgrp_comm )
          !
          DO i = 1, 3
             !
             magtot_nc(i) = magtot_nc(i) * omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
             !
          ENDDO
          !
          absmag = absmag * omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
          !
       ENDIF
       !
       RETURN
       !
     END SUBROUTINE compute_magnetization
     !
     !-----------------------------------------------------------------------
     FUNCTION delta_e()
       !-----------------------------------------------------------------------
       ! This function computes delta_e, where:
       !
       ! ... delta_e =  - \int rho%of_r(r)  v%of_r(r)
       !                - \int rho%kin_r(r) v%kin_r(r) [for Meta-GGA]
       !                - \sum rho%ns       v%ns       [for DFT+Hubbard]
       !                - \sum becsum       D1_Hxc     [for PAW]
       !
       USE xc_lib,  ONLY : xclib_dft_is
       !
       IMPLICIT NONE
       !
       REAL(DP) :: delta_e
       REAL(DP) :: delta_e_hub
       INTEGER  :: ir, na1, nt1, na2, nt2, m1, m2, equiv_na2, viz, is, is1, i, j
       !
       delta_e = 0._DP
       IF ( nspin==2 ) THEN
          !
          DO ir = 1,dfftp%nnr
            delta_e = delta_e - ( rho%of_r(ir,1) + rho%of_r(ir,2) ) * v%of_r(ir,1) &  ! up
                              - ( rho%of_r(ir,1) - rho%of_r(ir,2) ) * v%of_r(ir,2)    ! dw
          ENDDO 
          delta_e = 0.5_DP*delta_e
          !
       ELSE
          delta_e = - SUM( rho%of_r(:,1:nspin_mag)*v%of_r(:,1:nspin_mag) )
       ENDIF
       !
       IF ( xclib_dft_is('meta') ) &
          delta_e = delta_e - SUM( rho%kin_r(:,1:nspin_mag)*v%kin_r(:,1:nspin_mag) )
       !
       delta_e = omega * delta_e / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
       !
       CALL mp_sum( delta_e, intra_bgrp_comm )
       !
       IF (lda_plus_u .AND. (.NOT.hub_pot_fix)) THEN
         IF (lda_plus_u_kind.EQ.0) THEN
            IF (noncolin) THEN
              delta_e_hub = - SUM( rho%ns_nc(:,:,:,:)*v%ns_nc(:,:,:,:) )
              delta_e = delta_e + delta_e_hub              
            ELSE        
               delta_e_hub = - SUM( rho%ns(:,:,:,:)*v%ns(:,:,:,:) )
               IF (lhb) delta_e_hub = delta_e_hub - SUM(rho%nsb(:,:,:,:)*v%nsb(:,:,:,:))
               IF (nspin==1) delta_e_hub = 2.d0 * delta_e_hub
               delta_e = delta_e + delta_e_hub
            ENDIF  
         ELSEIF (lda_plus_u_kind.EQ.1) THEN
            IF (noncolin) THEN
              delta_e_hub = - SUM( rho%ns_nc(:,:,:,:)*v%ns_nc(:,:,:,:) )
              delta_e = delta_e + delta_e_hub
            ELSE
              delta_e_hub = - SUM( rho%ns(:,:,:,:)*v%ns(:,:,:,:) )
              IF (nspin==1) delta_e_hub = 2.d0 * delta_e_hub
              delta_e = delta_e + delta_e_hub
            ENDIF
         ELSEIF (lda_plus_u_kind.EQ.2) THEN
            delta_e_hub = 0.0_dp
            IF (noncolin) THEN
               DO is = 1, npol
                  DO is1 = 1, npol
                     i = npol * (is1-1) + is
                     DO na1 = 1, nat
                        nt1 = ityp(na1)
                        DO viz = 1, neighood(na1)%num_neigh
                           na2 = neighood(na1)%neigh(viz)
                           equiv_na2 = at_sc(na2)%at
                           nt2 = ityp(equiv_na2)
                           IF ( ANY(v_nsg(:,:,viz,na1,i).NE.0.0d0) ) THEN
                              DO m1 = 1, ldim_u(nt1)
                                 DO m2 = 1, ldim_u(nt2)
                                    delta_e_hub = delta_e_hub - &
                                       nsgnew(m2,m1,viz,na1,i)*v_nsg(m2,m1,viz,na1,i)
                                 ENDDO
                              ENDDO
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ELSE
               DO is = 1, nspin
                  DO na1 = 1, nat
                     nt1 = ityp(na1)
                     DO viz = 1, neighood(na1)%num_neigh
                        na2 = neighood(na1)%neigh(viz)
                        equiv_na2 = at_sc(na2)%at
                        nt2 = ityp(equiv_na2)
                        IF ( ANY(v_nsg(:,:,viz,na1,is).NE.0.0d0) ) THEN
                           DO m1 = 1, ldim_u(nt1)
                              DO m2 = 1, ldim_u(nt2)
                                 delta_e_hub = delta_e_hub - &
                                     nsgnew(m2,m1,viz,na1,is)*v_nsg(m2,m1,viz,na1,is)
                              ENDDO
                           ENDDO
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
            IF (nspin==1) delta_e_hub = 2.d0 * delta_e_hub
            delta_e = delta_e + delta_e_hub
         ENDIF
       ENDIF
       !
       IF (okpaw) delta_e = delta_e - SUM( ddd_paw(:,:,:)*rho%bec(:,:,:) )
       !
       RETURN
       !
     END FUNCTION delta_e
     !
     !-----------------------------------------------------------------------
     FUNCTION delta_escf()
       !-----------------------------------------------------------------------
       ! This function calculates the difference between the Hartree and XC energy
       ! at first order in the charge density difference delta_rho(r).
       !
       ! ... delta_escf = - \int \delta rho%of_r(r)  v%of_r(r)
       !                  - \int \delta rho%kin_r(r) v%kin_r(r) [for Meta-GGA]
       !                  - \sum \delta rho%ns       v%ns       [for DFT+Hubbard]
       !                  - \sum \delta becsum       D1         [for PAW] 
       !
       USE xc_lib, ONLY : xclib_dft_is
       !
       IMPLICIT NONE
       REAL(DP) :: delta_escf, delta_escf_hub, rho_dif(2)
       INTEGER  :: ir, na1, nt1, na2, nt2, m1, m2, equiv_na2, viz, is, is1, i, j
       !
       delta_escf=0._dp
       IF ( nspin==2 ) THEN
          !
          DO ir=1, dfftp%nnr
             !
             rho_dif = rhoin%of_r(ir,:) - rho%of_r(ir,:)
             !
             delta_escf = delta_escf - ( rho_dif(1) + rho_dif(2) ) * v%of_r(ir,1) &  !up
                                     - ( rho_dif(1) - rho_dif(2) ) * v%of_r(ir,2)    !dw
          ENDDO
          delta_escf = 0.5_dp*delta_escf
          !
       ELSE
         delta_escf = -SUM( ( rhoin%of_r(:,1:nspin_mag)-rho%of_r(:,1:nspin_mag) )*v%of_r(:,1:nspin_mag) )
       ENDIF
       !
       IF ( xclib_dft_is('meta') ) &
          delta_escf = delta_escf - &
                       SUM( (rhoin%kin_r(:,:)-rho%kin_r(:,:) )*v%kin_r(:,:))
       !
       delta_escf = omega * delta_escf / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
       !
       CALL mp_sum( delta_escf, intra_bgrp_comm )
       !
       IF (lda_plus_u .AND. (.NOT.hub_pot_fix)) THEN
          IF (lda_plus_u_kind.EQ.0) THEN
             IF (noncolin) THEN 
                delta_escf_hub = -SUM((rhoin%ns_nc(:,:,:,:)-rho%ns_nc(:,:,:,:))*v%ns_nc(:,:,:,:))
                delta_escf = delta_escf + delta_escf_hub
             ELSE        
                delta_escf_hub = -SUM((rhoin%ns(:,:,:,:)-rho%ns(:,:,:,:))*v%ns(:,:,:,:))
                IF (lhb) delta_escf_hub = delta_escf_hub - &
                            SUM((rhoin%nsb(:,:,:,:)-rho%nsb(:,:,:,:))*v%nsb(:,:,:,:))
                IF ( nspin==1 ) delta_escf_hub = 2.d0 * delta_escf_hub
                delta_escf = delta_escf + delta_escf_hub
             ENDIF   
          ELSEIF (lda_plus_u_kind.EQ.1) THEN
             IF (noncolin) THEN
                delta_escf_hub = -SUM((rhoin%ns_nc(:,:,:,:)-rho%ns_nc(:,:,:,:))*v%ns_nc(:,:,:,:))
                delta_escf = delta_escf + delta_escf_hub
             ELSE
                delta_escf_hub = -SUM((rhoin%ns(:,:,:,:)-rho%ns(:,:,:,:))*v%ns(:,:,:,:))
                IF ( nspin==1 ) delta_escf_hub = 2.d0 * delta_escf_hub
                delta_escf = delta_escf + delta_escf_hub
             ENDIF
          ELSEIF (lda_plus_u_kind.EQ.2) THEN
             delta_escf_hub = 0.0_dp
             IF (noncolin) THEN
               DO is = 1, npol
                  DO is1 = 1, npol
                     i =  npol*(is1-1) + is
                     DO na1 = 1, nat
                        nt1 = ityp(na1)
                        DO viz = 1, neighood(na1)%num_neigh
                           na2 = neighood(na1)%neigh(viz)
                           equiv_na2 = at_sc(na2)%at
                           nt2 = ityp(equiv_na2)
                           IF ( ANY(v_nsg(:,:,viz,na1,i).NE.0.0d0) ) THEN
                              DO m1 = 1, ldim_u(nt1)
                                 DO m2 = 1, ldim_u(nt2)
                                    delta_escf_hub = delta_escf_hub - &
                                          (nsg(m2,m1,viz,na1,i)-nsgnew(m2,m1,viz,na1,i)) * &
                                          v_nsg(m2,m1,viz,na1,i)
                                 ENDDO
                              ENDDO
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ELSE
               DO is = 1, nspin
                  DO na1 = 1, nat
                     nt1 = ityp(na1)
                     DO viz = 1, neighood(na1)%num_neigh
                        na2 = neighood(na1)%neigh(viz)
                        equiv_na2 = at_sc(na2)%at
                        nt2 = ityp(equiv_na2)
                        IF ( ANY(v_nsg(:,:,viz,na1,is).NE.0.0d0) ) THEN
                           DO m1 = 1, ldim_u(nt1)
                              DO m2 = 1, ldim_u(nt2)
                                 delta_escf_hub = delta_escf_hub - &
                                      (nsg(m2,m1,viz,na1,is)-nsgnew(m2,m1,viz,na1,is)) * &
                                       v_nsg(m2,m1,viz,na1,is)
                              ENDDO
                           ENDDO
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
             IF ( nspin==1 ) delta_escf_hub = 2.d0 * delta_escf_hub
             delta_escf = delta_escf + delta_escf_hub
          ENDIF
       ENDIF

       IF ( okpaw ) delta_escf = delta_escf - &
                                 SUM(ddd_paw(:,:,:)*(rhoin%bec(:,:,:)-rho%bec(:,:,:)))

       RETURN
       !
     END FUNCTION delta_escf
     !
     !-----------------------------------------------------------------------
     FUNCTION calc_pol( ) RESULT ( en_el )
       !-----------------------------------------------------------------------
       !
       USE kinds,     ONLY : DP
       USE constants, ONLY : pi
       USE bp,        ONLY : lelfield, ion_pol, el_pol, fc_pol, l_el_pol_old, &
                             el_pol_acc, el_pol_old, efield, l3dstring, gdir, &
                             transform_el, efield_cart
       !
       IMPLICIT NONE
       REAL (DP) :: en_el
       !
       INTEGER :: i, j 
       REAL(DP):: sca, el_pol_cart(3),  el_pol_acc_cart(3)
       !
       IF (.NOT.l3dstring) THEN
          !
          CALL c_phase_field( el_pol(gdir), ion_pol(gdir), fc_pol(gdir), gdir )
          !
          IF (.NOT.l_el_pol_old) THEN
             l_el_pol_old = .TRUE.
             el_pol_old(gdir) = el_pol(gdir)
             en_el = -efield*(el_pol(gdir)+ion_pol(gdir))
             el_pol_acc(gdir) = 0.d0
          ELSE
             sca = (el_pol(gdir)-el_pol_old(gdir))/fc_pol(gdir)
             IF (sca < - pi) THEN
                el_pol_acc(gdir) = el_pol_acc(gdir)+2.d0*pi*fc_pol(gdir)
             ELSEIF (sca > pi) THEN
                el_pol_acc(gdir) = el_pol_acc(gdir)-2.d0*pi*fc_pol(gdir)
             ENDIF
             en_el = -efield*(el_pol(gdir)+ion_pol(gdir)+el_pol_acc(gdir))
             el_pol_old = el_pol
          ENDIF
          !
       ELSE
          !
          DO i = 1, 3
            CALL c_phase_field( el_pol(i), ion_pol(i), fc_pol(i), i )
          ENDDO
          el_pol_cart(:) = 0.d0
          DO i = 1, 3
             DO j = 1, 3
                !el_pol_cart(i)=el_pol_cart(i)+transform_el(j,i)*el_pol(j)
                el_pol_cart(i) = el_pol_cart(i)+at(i,j)*el_pol(j) / &
                                 (SQRT(at(1,j)**2.d0+at(2,j)**2.d0+at(3,j)**2.d0))
             ENDDO
          ENDDO
          !
          WRITE( stdout,'( "Electronic Dipole on Cartesian axes" )' )
          DO i = 1, 3
             WRITE(stdout,*) i, el_pol_cart(i)
          ENDDO
          !
          WRITE( stdout,'( "Ionic Dipole on Cartesian axes" )' )
          DO i = 1, 3
             WRITE(stdout,*) i, ion_pol(i)
          ENDDO
          !
          IF (.NOT.l_el_pol_old) THEN
             l_el_pol_old = .TRUE.
             el_pol_old(:) = el_pol(:)
             en_el = 0.d0
             DO i = 1, 3
                en_el = en_el-efield_cart(i)*(el_pol_cart(i)+ion_pol(i))
             ENDDO
             el_pol_acc(:) = 0.d0
          ELSE
             DO i = 1, 3
                sca = (el_pol(i)-el_pol_old(i))/fc_pol(i)
                IF (sca < - pi) THEN
                   el_pol_acc(i) = el_pol_acc(i)+2.d0*pi*fc_pol(i)
                ELSEIF (sca > pi) THEN
                   el_pol_acc(i) = el_pol_acc(i)-2.d0*pi*fc_pol(i)
                ENDIF
             ENDDO
             el_pol_acc_cart(:) = 0.d0
             DO i = 1, 3
                DO j = 1, 3
                   el_pol_acc_cart(i) = el_pol_acc_cart(i)+transform_el(j,i)*el_pol_acc(j)
                ENDDO
             ENDDO
             en_el = 0.d0
             DO i = 1, 3
                en_el = en_el-efield_cart(i)*(el_pol_cart(i)+ion_pol(i)+el_pol_acc_cart(i))
             ENDDO
             el_pol_old(:) = el_pol(:)
          ENDIF
          !
       ENDIF
       !
     END FUNCTION calc_pol
     !
     !-----------------------------------------------------------------------
     SUBROUTINE print_energies ( printout )
       !-----------------------------------------------------------------------
       !
       USE constants, ONLY : eps8, RYTOEV
       INTEGER, INTENT (IN) :: printout
       !
   
       IF ( printout == 0 ) RETURN
       IF ( conv_elec .AND. printout > 1 ) THEN
          !
          WRITE( stdout, 9081 ) etot
          IF ( only_paw ) WRITE( stdout, 9085 ) etot+total_core_energy
          IF ( iverbosity > 1 ) WRITE( stdout, 9082 ) hwf_energy
          IF ( dr2 > eps8 ) THEN
             WRITE( stdout, 9083 ) dr2
          ELSE
             WRITE( stdout, 9084 ) dr2
          END IF
          !
          IF ( lgcscf )  WRITE( stdout, 9181 ) tot_charge
          !
          IF ( lgauss ) then
             WRITE( stdout, 9070 ) demet
             WRITE( stdout, 9170 ) etot-demet
          END IF
          !
          IF ( lgauss ) then
             WRITE( stdout, 9061 )
          ELSE
             WRITE( stdout, 9060 )
          END IF
          !
          WRITE( stdout, 9062 ) (eband + deband), ehart, ( etxc - etxcc ), ewld
          !
          IF ( llondon ) WRITE ( stdout , 9074 ) elondon
          IF ( ldftd3 )  WRITE ( stdout , 9078 ) edftd3
          IF ( lxdm )    WRITE ( stdout , 9075 ) exdm
          IF ( mbd_vdw ) THEN
             WRITE ( stdout , 9076 ) 2.0d0*Embdvdw
          ELSEIF ( ts_vdw ) THEN
             WRITE ( stdout , 9076 ) 2.0d0*EtsvdW
          ENDIF
          IF ( textfor)  WRITE ( stdout , 9077 ) eext
          IF ( tefield )            WRITE( stdout, 9064 ) etotefield
          IF ( gate )               WRITE( stdout, 9065 ) etotgatefield
          IF ( lda_plus_u )         WRITE( stdout, 9066 ) eth
          IF ( ABS (descf) > eps8 ) WRITE( stdout, 9069 ) descf
          IF ( okpaw ) THEN
            WRITE( stdout, 9067 ) epaw
            ! Detailed printout of PAW energy components, if verbosity is high
            IF(iverbosity>0)THEN
            WRITE( stdout, 9068) SUM(etot_cmp_paw(:,1,1)), &
                                 SUM(etot_cmp_paw(:,1,2)), &
                                 SUM(etot_cmp_paw(:,2,1)), &
                                 SUM(etot_cmp_paw(:,2,2)), &
            SUM(etot_cmp_paw(:,1,1))+SUM(etot_cmp_paw(:,1,2))+ehart, &
            SUM(etot_cmp_paw(:,2,1))+SUM(etot_cmp_paw(:,2,2))+etxc-etxcc
            ENDIF
          ENDIF
          !
          IF ( lrism ) THEN
             WRITE( stdout, 9902 ) esol
             IF ( ABS( vsol ) > eps8 ) WRITE( stdout, 9903 ) vsol
          END IF
          !
          ! ... With Grandcanonical-SCF (GC-SCF) or Fictitious charge particle (FCP),
          ! ... etot is the grand potential energy Omega = E - muN, where -muN is
          ! ... the potentiostat contribution.
          !
          IF ( lgcscf .AND. (.NOT. gcscf_ignore_mun) ) THEN
             WRITE( stdout, 9072 ) egrand
          END IF
          !
          IF ( lfcp ) WRITE( stdout, 9072 ) fcp_mu*tot_charge
          !
       ELSE IF ( conv_elec ) THEN
          !
          WRITE( stdout, 9081 ) etot
          IF ( iverbosity > 1 ) WRITE( stdout, 9082 ) hwf_energy
          IF ( dr2 > eps8 ) THEN
             WRITE( stdout, 9083 ) dr2
          ELSE
             WRITE( stdout, 9084 ) dr2
          END IF
          IF ( lgauss ) then
             WRITE( stdout, 9070 ) demet
             WRITE( stdout, 9170 ) etot-demet
          ENDIF
          !
          IF( lgcscf ) WRITE( stdout, 9181 ) tot_charge
          !
       ELSE
          !
          WRITE( stdout, 9080 ) etot
          IF (dmft) THEN
            WRITE( stdout, *) "    DMFT detected, writing all energy contributions"
            WRITE( stdout, 9062 ) (eband + deband), ehart, ( etxc - etxcc ), ewld
            WRITE( stdout, 9301 ) eband
          ENDIF
          IF ( iverbosity > 1 ) WRITE( stdout, 9082 ) hwf_energy
          IF ( dr2 > eps8 ) THEN
             WRITE( stdout, 9083 ) dr2
          ELSE
             WRITE( stdout, 9084 ) dr2
          END IF
          !
          IF( lgcscf ) WRITE( stdout, 9182 ) tot_charge, ef*RYTOEV, &
                                             ABS( ef - gcscf_mu ) * RYTOEV
          !
       ENDIF
       !
#if defined(__LEGACY_PLUGINS)
       CALL plugin_print_energies()
#endif 
#if defined (__ENVIRON)
       IF (use_environ) CALL print_environ_energies('PW')
#endif
#if defined (__OSCDFT)
       IF (use_oscdft .AND. conv_elec .AND. (oscdft_ctx%inp%oscdft_type==1)) &
          CALL oscdft_print_energies(oscdft_ctx)
#endif
       !
       IF ( lsda ) WRITE( stdout, 9017 ) magtot, absmag
       !
       IF ( noncolin .AND. domag ) &
            WRITE( stdout, 9018 ) magtot_nc(1:3), absmag
       !
       IF ( i_cons == 3 .OR. i_cons == 4 )  &
            WRITE( stdout, 9071 ) bfield(1), bfield(2), bfield(3)
       IF ( i_cons /= 0 .AND. i_cons < 4 ) &
            WRITE( stdout, 9073 ) lambda
       !
       FLUSH( stdout )
       !
       RETURN
       !
9017 FORMAT(/'     total magnetization       =', F9.2,' Bohr mag/cell', &
            /'     absolute magnetization    =', F9.2,' Bohr mag/cell' )
9018 FORMAT(/'     total magnetization       =',3F9.2,' Bohr mag/cell' &
       &   ,/'     absolute magnetization    =', F9.2,' Bohr mag/cell' )
9060 FORMAT(/'     The total energy is the sum of the following terms:' )
9061 FORMAT(/'     The total energy is F=E-TS. E is the sum of the following terms:' )
9062 FORMAT( '     one-electron contribution =',F17.8,' Ry' &
            /'     hartree contribution      =',F17.8,' Ry' &
            /'     xc contribution           =',F17.8,' Ry' &
            /'     ewald contribution        =',F17.8,' Ry' )
9064 FORMAT( '     electric field correction =',F17.8,' Ry' )
9065 FORMAT( '     gate field correction     =',F17.8,' Ry' ) ! TB
9066 FORMAT( '     Hubbard energy            =',F17.8,' Ry' )
9067 FORMAT( '     one-center paw contrib.   =',F17.8,' Ry' )
9068 FORMAT( '      -> PAW hartree energy AE =',F17.8,' Ry' &
            /'      -> PAW hartree energy PS =',F17.8,' Ry' &
            /'      -> PAW xc energy AE      =',F17.8,' Ry' &
            /'      -> PAW xc energy PS      =',F17.8,' Ry' &
            /'      -> total E_H with PAW    =',F17.8,' Ry' &
            /'      -> total E_XC with PAW   =',F17.8,' Ry' )
9069 FORMAT( '     scf correction            =',F17.8,' Ry' )
9070 FORMAT( '     smearing contrib. (-TS)   =',F17.8,' Ry' )
9071 FORMAT( '     Magnetic field            =',3F12.7,' Ry' )
9072 FORMAT( '     pot.stat. contrib. (-muN) =',F17.8,' Ry' )
9073 FORMAT( '     lambda                    =',F11.2,' Ry' )
9074 FORMAT( '     Dispersion Correction     =',F17.8,' Ry' )
9075 FORMAT( '     Dispersion XDM Correction =',F17.8,' Ry' )
9076 FORMAT( '     Dispersion Correction     =',F17.8,' Ry' )
9077 FORMAT( '     External forces energy    =',F17.8,' Ry' )
9078 FORMAT( '     DFT-D3 Dispersion         =',F17.8,' Ry' )
9080 FORMAT(/'     total energy              =',0PF17.8,' Ry' )
9081 FORMAT(/'!    total energy              =',0PF17.8,' Ry' )
9082 FORMAT( '     Harris-Foulkes estimate   =',0PF17.8,' Ry' )
9083 FORMAT( '     estimated scf accuracy    <',0PF17.8,' Ry' )
9084 FORMAT( '     estimated scf accuracy    <',1PE17.1,' Ry' )
9085 FORMAT( '     total all-electron energy =',0PF17.6,' Ry' )
9170 FORMAT( '     internal energy E=F+TS    =',0PF17.8,' Ry' )
9181 FORMAT(                                                  &
            /'!    total charge of GC-SCF    =',0PF17.8,' e' )
9182 FORMAT(                                                  &
            /'     total charge of GC-SCF    =',0PF17.8,' e'  &
            /'     the Fermi energy          =',0PF17.8,' eV' &
            /'                        (error :',0PF17.8,' eV)')
9301 FORMAT( '     band energy (sum(wg*et))  =',F17.8,' Ry' )

9902 FORMAT( '     solvation energy (RISM)   =',F17.8,' Ry' )
9903 FORMAT( '     level-shifting contrib.   =',F17.8,' Ry' )
  END SUBROUTINE print_energies
  !
END SUBROUTINE electrons_scf
!
!----------------------------------------------------------------------------
FUNCTION exxenergyace( )
  !--------------------------------------------------------------------------
  !! Compute exchange energy using ACE
  !
  USE kinds,              ONLY : DP
  USE buffers,            ONLY : get_buffer
  USE exx,                ONLY : vexxace_gamma, vexxace_k, domat, &
                                 vexxace_gamma_gpu, vexxace_k_gpu
  USE klist,              ONLY : nks, ngk
  USE wvfct,              ONLY : nbnd, npwx, current_k
  USE lsda_mod,           ONLY : lsda, isk, current_spin
  USE io_files,           ONLY : iunwfc, nwordwfc
  USE mp_pools,           ONLY : inter_pool_comm
  USE mp_bands,           ONLY : intra_bgrp_comm
  USE mp,                 ONLY : mp_sum
  USE control_flags,      ONLY : gamma_only, use_gpu
  USE wavefunctions,      ONLY : evc
  !
  IMPLICIT NONE
  !
  REAL(DP) :: exxenergyace
  !! computed energy
  !
  ! ... local variables
  !
  REAL(DP) :: ex
  INTEGER :: ik, npw
  !
  domat = .TRUE.
  exxenergyace=0.0_dp
  !
  DO ik = 1, nks
     npw = ngk (ik)
     !
     current_k = ik
     IF ( lsda ) current_spin = isk(ik)
     !
     IF (nks > 1) THEN
        CALL get_buffer( evc, nwordwfc, iunwfc, ik )
        !$acc update device(evc)
     ENDIF
     !
     IF (gamma_only) THEN
        IF (use_gpu) THEN
          !$acc host_data use_device(evc)
          CALL vexxace_gamma_gpu( npw, nbnd, evc, ex )
          !$acc end host_data
        ELSE
          CALL vexxace_gamma( npw, nbnd, evc, ex )
        END IF
     ELSE
        IF (use_gpu) THEN
          !$acc host_data use_device(evc)
          CALL vexxace_k_gpu( npw, nbnd, evc, ex )
          !$acc end host_data
        ELSE
          CALL vexxace_k( npw, nbnd, evc, ex )
        ENDIF
     ENDIF
     exxenergyace = exxenergyace + ex
  ENDDO
  !
  CALL mp_sum( exxenergyace, inter_pool_comm )
  !
  domat = .FALSE.
  !
END FUNCTION exxenergyace
