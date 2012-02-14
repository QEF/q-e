!
! Copyright (C) 2001-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE electrons()
  !----------------------------------------------------------------------------
  !
  ! ... This routine is a driver of the self-consistent cycle.
  ! ... It uses the routine c_bands for computing the bands at fixed
  ! ... Hamiltonian, the routine sum_band to compute the charge density,
  ! ... the routine v_of_rho to compute the new potential and the routine
  ! ... mix_rho to mix input and output charge densities.
  ! ... It prints on output the total energy and its decomposition in
  ! ... the separate contributions.
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : eps8, pi
  USE io_global,            ONLY : stdout, ionode
  USE cell_base,            ONLY : at, bg, alat, omega, tpiba2
  USE ions_base,            ONLY : zv, nat, nsp, ityp, tau, compute_eextfor
  USE basis,                ONLY : starting_pot, starting_wfc
  USE bp,                   ONLY : lelfield
  USE fft_base,             ONLY : dfftp
  USE gvect,                ONLY : ngm, gstart, nl, nlm, g, gg, gcutm
  USE gvecs,                ONLY : doublegrid, ngms
  USE klist,                ONLY : xk, wk, nelec, ngk, nks, nkstot, lgauss
  USE lsda_mod,             ONLY : lsda, nspin, magtot, absmag, isk
  USE vlocal,               ONLY : strf
  USE wvfct,                ONLY : nbnd, et, npwx, ecutwfc
  USE ener,                 ONLY : etot, hwf_energy, eband, deband, ehart, &
                                   vtxc, etxc, etxcc, ewld, demet, epaw, &
                                   elondon
  USE scf,                  ONLY : scf_type, scf_type_COPY, &
                                   create_scf_type, destroy_scf_type, &
                                   rho, rho_core, rhog_core, &
                                   v, vltot, vrs, kedtau, vnew
  USE control_flags,        ONLY : mixing_beta, tr2, ethr, niter, nmix, &
                                   iprint, istep, lscf, lmd, conv_elec, &
                                   restart, io_level, do_makov_payne,  &
                                   gamma_only, iverbosity, textfor,     &
                                   llondon
  USE io_files,             ONLY : iunwfc, iunocc, nwordwfc, output_drho, &
                                   iunefield, iunpaw
  USE buffers,              ONLY : save_buffer
  USE ldaU,                 ONLY : eth, Hubbard_U, Hubbard_lmax, &
                                   niter_with_fixed_ns, lda_plus_u
  USE extfield,             ONLY : tefield, etotefield
  USE wavefunctions_module, ONLY : evc, psic
  USE noncollin_module,     ONLY : noncolin, magtot_nc, i_cons,  bfield, &
                                   lambda, report
  USE spin_orb,             ONLY : domag
  USE io_rho_xml,           ONLY : write_rho
  USE uspp,                 ONLY : okvan
  USE exx,                  ONLY : exxinit, exxenergy, exxenergy2, &
                                   fock0, fock1, fock2, dexx, exx_restart
  USE funct,                ONLY : dft_is_hybrid, exx_is_active
  USE control_flags,        ONLY : adapt_thr, tr2_init, tr2_multi
  USE funct,                ONLY : dft_is_meta
  USE mp_global,            ONLY : intra_pool_comm, npool, intra_bgrp_comm, nbgrp, mpime, &
                                   inter_bgrp_comm, my_bgrp_id
  USE mp,                   ONLY : mp_sum
  !
  USE london_module,        ONLY : energy_london
  !
  USE paw_variables,        ONLY : okpaw, ddd_paw, total_core_energy, only_paw
  USE paw_onecenter,        ONLY : PAW_potential
  USE paw_symmetry,         ONLY : PAW_symmetrize_ddd
  USE uspp_param,           ONLY : nh, nhm ! used for PAW
#ifdef __ENVIRON
  USE environ_base,         ONLY : do_environ, update_venviron,             &
                                   vltot_zero, environ_thr,                 &
                                   env_static_permittivity,                 & 
                                   env_surface_tension, env_pressure,       &
                                   deenviron, esolvent, ecavity, epressure
#endif
  USE dfunct,                 only : newd
  USE esm,                  ONLY : do_comp_esm, esm_printpot
  !
  !
  IMPLICIT NONE
  !
  ! ... a few local variables
  !
  REAL(DP) :: &
      dr2,          &! the norm of the diffence between potential
      charge,       &! the total charge
      deband_hwf,   &! deband for the Harris-Weinert-Foulkes functional
      mag           ! local magnetization
  INTEGER :: &
      i,            &! counter on polarization
      idum,         &! dummy counter on iterations
      iter,         &! counter on iterations
      ik_,          &! used to read ik from restart file
      kilobytes
  REAL(DP) :: &
      tr2_min,     &! estimated error on energy coming from diagonalization
      descf,       &! correction for variational energy
      en_el=0.0_DP,&! electric field contribution to the total energy
      eext=0.0_DP   ! external forces contribution to the total energy
  LOGICAL :: &
      first
  !
  ! ... auxiliary variables for calculating and storing temporary copies of
  ! ... the charge density and of the HXC-potential
  !
  type (scf_type), save :: rhoin ! used to store rho_in of current/next iteration
  !
  ! ... external functions
  !
  REAL(DP), EXTERNAL :: ewald, get_clock
  REAL(DP) :: tr2_final ! final threshold for exx minimization 
                        ! when using adaptive thresholds.
  iter = 0
  ik_  = 0
  dr2  = 0.0_dp
  tr2_final = tr2
  IF (dft_is_hybrid() .AND. adapt_thr ) THEN
     tr2= tr2_init
  ENDIF
  !
  IF ( restart ) THEN
     !
     CALL restart_in_electrons( iter, ik_, dr2 )
     !
     IF ( ik_ == -1000 ) THEN
        !
        conv_elec = .TRUE.
        !
        IF ( output_drho /= ' ' ) CALL remove_atomic_rho ()
        !
        RETURN
        !
     END IF
     !
     IF( exx_is_active() ) THEN
       iter = 0
       call save_in_electrons( iter, dr2 )
       WRITE( stdout, '(5x,"EXX: now go back to refine exchange calculation")')
     ELSE IF ( dft_is_hybrid() .AND. TRIM(starting_wfc) == 'file' ) THEN
       !
       ! suggested by Hannu Komsa: useful when several calculations with 
       ! different values of alpha have to be performed
       !
       call exx_restart(.true.)
       WRITE( stdout, '(5x,"EXX: now go back to refine exchange calculation")')
     ENDIF
  END IF
  !
  WRITE( stdout, 9000 ) get_clock( 'PWSCF' )
  !
  CALL memstat( kilobytes )
  !
  IF ( kilobytes > 0 ) WRITE( stdout, 9001 ) kilobytes/1000.0
  !
  CALL flush_unit( stdout )
  !
  IF ( .NOT. lscf ) THEN
     !
     CALL non_scf (ik_)
     !
     conv_elec = .TRUE.
     !
     RETURN
     !
  END IF
  !
  CALL start_clock( 'electrons' )
  !
  if ( exx_is_active())  then
     CALL v_of_rho( rho, rho_core, rhog_core, &
                    ehart, etxc, vtxc, eth, etotefield, charge, v)
     CALL set_vrs( vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid )
  end if
  !  
#ifdef __ENVIRON
  IF ( do_environ ) THEN
    vltot_zero = vltot
    CALL environ_initions( nat, nsp, ityp, zv, tau ) 
    CALL environ_initcell( dfftp%nr1*dfftp%nr2*dfftp%nr3, omega ) 
  END IF
#endif
  !  
  CALL flush_unit( stdout )
  !
  ! ... calculates the ewald contribution to total energy
  !
  ewld = ewald( alat, nat, nsp, ityp, zv, at, bg, tau, &
                omega, g, gg, ngm, gcutm, gstart, gamma_only, strf )
  !
  !
  !
  elondon = 0.d0
  !
  IF ( llondon ) THEN
  !
  elondon = energy_london ( alat , nat , ityp , at ,&
                                         bg , tau )
  END IF
  !
  call create_scf_type ( rhoin )
  !
10 CONTINUE
  !
  ! ... Convergence threshold for iterative diagonalization
  !
  ! ... for the first scf iteration of each ionic step (after the first),
  ! ... the threshold is fixed to a default value of 1.D-6
  !
  IF ( istep > 0 ) ethr = 1.D-6
  !
  WRITE( stdout, 9002 )
  !
  CALL flush_unit( stdout )

  !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%          iterate !          %%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
  DO idum = 1, niter
     !
     IF ( check_stop_now() ) RETURN
     !
     iter = iter + 1
     !
     WRITE( stdout, 9010 ) iter, ecutwfc, mixing_beta
     !
     CALL flush_unit( stdout )
     !
     ! ... Convergence threshold for iterative diagonalization is
     ! ... automatically updated during self consistency
     !
     IF ( iter > 1 .AND. ik_ == 0 ) THEN
        !
        IF ( iter == 2 ) ethr = 1.D-2
        ethr = MIN( ethr, 0.1D0*dr2 / MAX( 1.D0, nelec ) )
        ! ... do not allow convergence threshold to become too small:
        ! ... iterative diagonalization may become unstable
        ethr = MAX( ethr, 1.D-13 )
        !
     END IF
     !
     first = ( iter == 1 )
     !
     ! ... deband = - \sum_v <\psi_v | V_h + V_xc |\psi_v> is calculated a
     ! ... first time here using the input density and potential ( to be
     ! ... used to calculate the Harris-Weinert-Foulkes energy )
     !
     deband_hwf = delta_e()
     !
     ! save input current density in rhoin
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
        IF ( lelfield ) THEN
           !
           CALL c_bands_efield ( iter, ik_, dr2 )
           !
        ELSE
           !
           CALL c_bands( iter, ik_, dr2 )
           !
        END IF
        !
        IF ( check_stop_now() ) RETURN
        !
        ! ... xk, wk, isk, et, wg are distributed across pools;
        ! ... the first node has a complete copy of xk, wk, isk,
        ! ... while eigenvalues et and weights wg must be
        ! ... explicitely collected to the first node
        ! ... this is done here for et, in sum_band for wg
        !
        CALL poolrecover( et, nbnd, nkstot, nks )
        !
        ! ... the new density is computed here
        CALL sum_band()
        ! PAW : sum_band computes new becsum (stored in uspp modules) and a
        ! subtly different copy in rho%bec (scf module)
        !
        ! ... the Harris-Weinert-Foulkes energy is computed here using only
        ! ... quantities obtained from the input density
        !
        hwf_energy = eband + deband_hwf + (etxc - etxcc) + ewld + ehart + demet
        !
        IF ( lda_plus_u )  THEN
           !
           hwf_energy = hwf_energy + eth
           !
           IF ( iverbosity > 0 .OR. first ) CALL write_ns()
           !
           IF ( first .AND. istep == 0 .AND. starting_pot == 'atomic' ) THEN
              CALL ns_adj()
               rhoin%ns = rho%ns
           END IF
           IF ( iter <= niter_with_fixed_ns ) THEN
              WRITE( stdout, '(/,5X,"RESET ns to initial values (iter <= mixing_fixed_ns)",/)')
              rho%ns = rhoin%ns
           END IF
           !
        END IF
        IF (okpaw) hwf_energy = hwf_energy + epaw
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
        ! ... mix_rho mixes several quantities: rho in g-space, tauk (for meta-gga)
        ! ... ns (for lda+u) and becsum (for paw)
        !
        CALL mix_rho( rho, rhoin, mixing_beta, dr2, tr2_min, iter, nmix, conv_elec )
        !
        ! ... if convergence is achieved or if the self-consistency error
        ! ... (dr2) is smaller than the estimated error due to diagonalization
        ! ... (tr2_min), rhoin and rho are unchanged: rhoin contains the input
        ! ...  density and rho contains the output density
        ! ... In the other cases rhoin contains the mixed charge density 
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
              ethr = 0.1D0*dr2 / MAX( 1.D0, nelec )
              !
              CYCLE scf_step
              !
           END IF
           !
        END IF
        !
        not_converged_electrons : &
        IF ( .NOT. conv_elec ) THEN
           ! ... no convergence yet: calculate new potential from mixed
           ! ... charge density (i.e. the new estimate)
           !
           CALL v_of_rho( rhoin, rho_core, rhog_core, &
                          ehart, etxc, vtxc, eth, etotefield, charge, v)
           IF (okpaw) THEN
              CALL PAW_potential(rhoin%bec, ddd_paw, epaw)
              CALL PAW_symmetrize_ddd(ddd_paw)
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
           ! ... write the charge density to file
           ! ... also write ldaU ns coeffs and PAW becsum
           CALL write_rho( rho, nspin )
           !
        ELSE not_converged_electrons
           !
           ! ... convergence reached:
           ! ... 1) the output HXC-potential is saved in vr
           ! ... 2) vnew contains V(out)-V(in) ( used to correct the forces ).
           !
           vnew%of_r(:,:) = v%of_r(:,:)
           !
           CALL v_of_rho( rho,rho_core,rhog_core, &
                          ehart, etxc, vtxc, eth, etotefield, charge, v)
           !
           IF (okpaw) THEN
              CALL PAW_potential(rho%bec, ddd_paw, epaw)
              CALL PAW_symmetrize_ddd(ddd_paw)
           ENDIF

           !
           vnew%of_r(:,:) = v%of_r(:,:) - vnew%of_r(:,:)
           !
           ! ... note that rho is here the output, not mixed, charge density
           ! ... so correction for variational energy is no longer needed
           !
           descf = 0._dp
           !
        END IF not_converged_electrons
        !
        IF ( exx_is_active() ) THEN
           !
           fock1 = exxenergy2()
           fock2 = fock0
           !
        ELSE
           !
           fock0 = 0.D0
           fock1 = 0.D0
           fock2 = 0.D0
           !
        END IF
        !
        ! ... if we didn't cycle before we can exit the do-loop
        !
        EXIT scf_step
        !
     END DO scf_step
     !
#ifdef __ENVIRON
     ! ... computes the external environment contribution to energy and potential
     !
     IF ( do_environ  )  THEN
       !
       vltot = vltot_zero
       !
       CALL calc_eenviron( dfftp%nnr, nspin, rhoin%of_r, vltot_zero, &
                           deenviron, esolvent, ecavity, epressure )
       !
       update_venviron = .NOT. conv_elec .AND. dr2 .LT. environ_thr
       !
       IF ( update_venviron ) WRITE( stdout, 9200 )
       !
       CALL calc_venviron( update_venviron, dfftp%nnr, nspin, rhoin%of_r, vltot )
       !
     END IF
#endif
     !
     ! ... define the total local potential (external + scf)
     !
     CALL set_vrs( vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid )
     !
     ! ... in the US case we have to recompute the self-consistent
     ! ... term in the nonlocal potential
     !
     ! ... PAW: newd contains PAW updates of NL coefficients
     CALL newd()
     !
     ! ... save converged wfc if they have not been written previously
     !
     IF ( nks == 1 .AND. (io_level < 2) ) &
        CALL save_buffer ( evc, nwordwfc, iunwfc, nks )
     !
     ! ... calculate the polarization
     !
     IF ( lelfield ) THEN
        CALL calc_pol (en_el)
     ELSE
        en_el=0.d0 
     ENDIF
     !
     ! ... write recover file
     !
     CALL save_in_electrons( iter, dr2 )
     !
     IF ( ( MOD( iter, report ) == 0 ) .OR. &
          ( report /= 0 .AND. conv_elec ) ) THEN
        !
        IF ( noncolin .AND. domag .or. i_cons==1) CALL report_mag()
        !
     END IF
     !
     WRITE( stdout, 9000 ) get_clock( 'PWSCF' )
     !
     IF ( conv_elec ) WRITE( stdout, 9101 )
     !
     IF ( conv_elec .OR. MOD( iter, iprint ) == 0 ) THEN
        !
        IF ( lda_plus_U .AND. iverbosity < 1 ) CALL write_ns ( )
        CALL print_ks_energies()
        !
     END IF
     !
     IF ( ABS( charge - nelec ) / charge > 1.D-7 ) THEN
        WRITE( stdout, 9050 ) charge, nelec
        IF ( ABS( charge - nelec ) / charge > 1.D-3 ) THEN
           IF (.not.lgauss) THEN
              CALL errore( 'electrons', 'charge is wrong: smearing is needed', 1 )
           ELSE
              CALL errore( 'electrons', 'charge is wrong', 1 )
           END IF
        END IF
     END IF
     !
     etot = eband + ( etxc - etxcc ) + ewld + ehart + deband + demet + descf +en_el
     IF (okpaw) etot = etot + epaw
     IF( textfor ) THEN
        eext =  compute_eextfor()
        etot = etot + eext
     END IF
     IF (llondon) THEN
        etot = etot + elondon
        hwf_energy = hwf_energy + elondon
     END IF
     !
     etot = etot - 0.5D0*fock0
     hwf_energy = hwf_energy -0.5D0*fock0
     !
     IF ( dft_is_hybrid() .AND. conv_elec ) THEN
        !
        first = .NOT. exx_is_active()
        !
        CALL exxinit()
        !
        IF ( first ) THEN
           !
           fock0 = exxenergy2()
           !
           CALL v_of_rho( rho, rho_core,rhog_core, &
                          ehart, etxc, vtxc, eth, etotefield, charge, v)
           IF (okpaw) CALL PAW_potential(rho%bec, ddd_PAW, epaw)
           !
           CALL set_vrs( vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid )
           !
           conv_elec = .false.
           iter = 0
           CALL save_in_electrons( iter, dr2 )
           WRITE( stdout,'(5x,"EXX: now go back to refine exchange calculation")')
           WRITE( stdout, * ) fock0
           !
           GO TO 10
           !
        END IF
        !
        fock2 = exxenergy2()
        !
        dexx = fock1 - 0.5D0*( fock0 + fock2 )
        !
        etot = etot  - dexx
        hwf_energy = hwf_energy - dexx
        !
        WRITE( stdout, * ) fock0, fock1, fock2
        WRITE( stdout, 9066 ) dexx
        !
        fock0 = fock2
        !
     END IF
     !
     IF ( lda_plus_u ) etot = etot + eth
     IF ( tefield ) THEN
        etot = etot + etotefield
        hwf_energy = hwf_energy + etotefield
     END IF
     !
#ifdef __ENVIRON
     !
     ! ... adds the external environment contribution to the energy
     !
     IF ( do_environ ) etot = etot + deenviron + esolvent + ecavity + epressure
#endif
     !
     IF ( ( conv_elec .OR. MOD( iter, iprint ) == 0 ) .AND. .NOT. lmd ) THEN
        !
        IF ( dr2 > eps8 ) THEN
           WRITE( stdout, 9081 ) etot, hwf_energy, dr2
        ELSE
           WRITE( stdout, 9083 ) etot, hwf_energy, dr2
        END IF
        IF ( only_paw ) WRITE( stdout, 9085 ) etot+total_core_energy
        !
        WRITE( stdout, 9060 ) &
            ( eband + deband ), ehart, ( etxc - etxcc ), ewld
        !
        IF ( llondon ) WRITE ( stdout , 9074 ) elondon
        !
        IF ( dft_is_hybrid()) THEN
           WRITE( stdout, 9062 ) - fock1
           WRITE( stdout, 9064 ) 0.5D0*fock2
        ENDIF
        !
        IF ( textfor)             WRITE( stdout, &
            '(/5x,"Energy of the external Forces = ", F18.8)' ) eext
        IF ( tefield )            WRITE( stdout, 9061 ) etotefield
        IF ( lda_plus_u )         WRITE( stdout, 9065 ) eth
        IF ( ABS (descf) > eps8 ) WRITE( stdout, 9069 ) descf
        IF ( okpaw )              WRITE( stdout, 9067 ) epaw
        !
        ! ... With Fermi-Dirac population factor, etot is the electronic
        ! ... free energy F = E - TS , demet is the -TS contribution
        !
        IF ( lgauss ) WRITE( stdout, 9070 ) demet
        !
     ELSE IF ( conv_elec .AND. lmd ) THEN
        !
        IF ( dr2 > eps8 ) THEN
           WRITE( stdout, 9081 ) etot, hwf_energy, dr2
        ELSE
           WRITE( stdout, 9083 ) etot, hwf_energy, dr2
        END IF
        !
     ELSE
        !
        IF ( dr2 > eps8 ) THEN
           WRITE( stdout, 9080 ) etot, hwf_energy, dr2
        ELSE
           WRITE( stdout, 9082 ) etot, hwf_energy, dr2
        END IF
     END IF
     !
#ifdef __ENVIRON
     IF ( do_environ )  THEN
        IF ( env_static_permittivity .GT. 1.D0 ) WRITE( stdout, 9201 ) esolvent
        IF ( env_surface_tension .GT. 0.D0 ) WRITE( stdout, 9202 ) ecavity
        IF ( env_pressure .NE. 0.D0 ) WRITE( stdout, 9203 ) epressure
     ENDIF
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
     CALL flush_unit( stdout )
     !
     IF ( conv_elec ) THEN
        !
        IF ( dft_is_hybrid() .AND. dexx > tr2_final ) THEN  
           !
           conv_elec = .false.
           iter = 0

           CALL save_in_electrons( iter, dr2 )
           !
           WRITE (stdout,*) " NOW GO BACK TO REFINE HYBRID CALCULATION"

           IF ( adapt_thr ) THEN
              tr2 = MAX(tr2_multi * dexx, tr2_final)
              WRITE( stdout, 9121 ) tr2
           ENDIF
           !
           GO TO 10
           !
        END IF
        !
        ! ... if system is charged add a Makov-Payne correction to the energy
        !
        IF ( do_makov_payne ) CALL makov_payne( etot )
        !
        ! ... print out ESM potentials if desired
        !
        IF ( do_comp_esm ) CALL esm_printpot()
        !
        WRITE( stdout, 9110 ) iter
        !
        ! ... jump to the end
        !
        IF ( output_drho /= ' ' ) CALL remove_atomic_rho()
        !
        CALL stop_clock( 'electrons' )
        !
        call destroy_scf_type ( rhoin )
        !
        RETURN
        !
     END IF
     !
     ! ... uncomment the following line if you wish to monitor the evolution
     ! ... of the force calculation during self-consistency
     !
     !CALL forces()
     !
  END DO
  !
  WRITE( stdout, 9101 )
  WRITE( stdout, 9120 ) iter
  !
  CALL flush_unit( stdout )
  !
  IF ( output_drho /= ' ' ) CALL remove_atomic_rho()
  !
  CALL stop_clock( 'electrons' )
  !
  RETURN
  !
  ! ... formats
  !
9000 FORMAT(/'     total cpu time spent up to now is ',F10.1,' secs' )
9001 FORMAT(/'     per-process dynamical memory: ',f7.1,' Mb' )
9002 FORMAT(/'     Self-consistent Calculation' )
9010 FORMAT(/'     iteration #',I3,'     ecut=', F9.2,' Ry',5X,'beta=',F4.2 )
9017 FORMAT(/'     total magnetization       =', F9.2,' Bohr mag/cell', &
            /'     absolute magnetization    =', F9.2,' Bohr mag/cell' )
9018 FORMAT(/'     total magnetization       =',3F9.2,' Bohr mag/cell' &
       &   ,/'     absolute magnetization    =', F9.2,' Bohr mag/cell' )
9050 FORMAT(/'     WARNING: integrated charge=',F15.8,', expected=',F15.8 )
9060 FORMAT(/'     The total energy is the sum of the following terms:',/,&
            /'     one-electron contribution =',F17.8,' Ry' &
            /'     hartree contribution      =',F17.8,' Ry' &
            /'     xc contribution           =',F17.8,' Ry' &
            /'     ewald contribution        =',F17.8,' Ry' )
9061 FORMAT( '     electric field correction =',F17.8,' Ry' )
9062 FORMAT( '     - averaged Fock potential =',F17.8,' Ry' )
9064 FORMAT( '     + Fock energy             =',F17.8,' Ry' )
9065 FORMAT( '     Hubbard energy            =',F17.8,' Ry' )
9066 FORMAT( '     est. exchange err (dexx)  =',F17.8,' Ry' )
9067 FORMAT( '     one-center paw contrib.   =',F17.8,' Ry' )
9069 FORMAT( '     scf correction            =',F17.8,' Ry' )
9070 FORMAT( '     smearing contrib. (-TS)   =',F17.8,' Ry' )
9071 FORMAT( '     Magnetic field            =',3F12.7,' Ry' )
9072 FORMAT( '     Magnetic field            =',F12.7, ' Ry' )
9073 FORMAT( '     lambda                    =',F11.2,' Ry' )
9074 FORMAT( '     Dispersion Correction     =',F17.8,' Ry' )
9080 FORMAT(/'     total energy              =',0PF17.8,' Ry' &
            /'     Harris-Foulkes estimate   =',0PF17.8,' Ry' &
            /'     estimated scf accuracy    <',0PF17.8,' Ry' )
9081 FORMAT(/'!    total energy              =',0PF17.8,' Ry' &
            /'     Harris-Foulkes estimate   =',0PF17.8,' Ry' &
            /'     estimated scf accuracy    <',0PF17.8,' Ry' )
9082 FORMAT(/'     total energy              =',0PF17.8,' Ry' &
            /'     Harris-Foulkes estimate   =',0PF17.8,' Ry' &
            /'     estimated scf accuracy    <',1PE17.1,' Ry' )
9083 FORMAT(/'!    total energy              =',0PF17.8,' Ry' &
            /'     Harris-Foulkes estimate   =',0PF17.8,' Ry' &
            /'     estimated scf accuracy    <',1PE17.1,' Ry' )
9085 FORMAT(/'     total all-electron energy =',0PF17.6,' Ry' )
9101 FORMAT(/'     End of self-consistent calculation' )
9110 FORMAT(/'     convergence has been achieved in ',i3,' iterations' )
9120 FORMAT(/'     convergence NOT achieved after ',i3,' iterations: stopping' )
9121 FORMAT(/'     scf convergence threshold =',1PE17.1,' Ry' )
#ifdef __ENVIRON
9200 FORMAT(/'     add environment contribution to local potential')
9201 FORMAT( '     solvation energy          =',F17.8,' Ry' ) 
9202 FORMAT( '     cavitation energy         =',F17.8,' Ry' ) 
9203 FORMAT( '     PV energy                 =',F17.8,' Ry' ) 
#endif
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
             mag = rho%of_r(ir,1) - rho%of_r(ir,2)
             !
             magtot = magtot + mag
             absmag = absmag + ABS( mag )
             !
          END DO
          !
          magtot = magtot * omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
          absmag = absmag * omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
          !
#ifdef __BANDS
          CALL mp_sum( magtot, intra_bgrp_comm )
          CALL mp_sum( absmag, intra_bgrp_comm )
#else
          CALL mp_sum( magtot, intra_pool_comm )
          CALL mp_sum( absmag, intra_pool_comm )
#endif
          !
       ELSE IF ( noncolin ) THEN
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
             END DO
             !
             absmag = absmag + ABS( mag )
             !
          END DO
          !
#ifdef __BANDS
          CALL mp_sum( magtot, intra_bgrp_comm )
          CALL mp_sum( absmag, intra_bgrp_comm )
#else
          CALL mp_sum( magtot, intra_pool_comm )
          CALL mp_sum( absmag, intra_pool_comm )
#endif
          !
          DO i = 1, 3
             !
             magtot_nc(i) = magtot_nc(i) * omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
             !
          END DO
          !
          absmag = absmag * omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
          !
       END IF
       !
       RETURN
       !
     END SUBROUTINE compute_magnetization
     !
     !-----------------------------------------------------------------------
     FUNCTION check_stop_now()
       !-----------------------------------------------------------------------
       !
       USE check_stop,    ONLY : global_check_stop_now => check_stop_now
       !
       IMPLICIT NONE
       !
       LOGICAL :: check_stop_now
       INTEGER :: unit
       !
       unit = stdout
       !
       check_stop_now = global_check_stop_now( unit )
       !
       IF ( check_stop_now ) conv_elec = .FALSE.
       !
       RETURN
       !
     END FUNCTION check_stop_now
     !
     !-----------------------------------------------------------------------
     FUNCTION delta_e()
       !-----------------------------------------------------------------------
       ! ... delta_e = - \int rho%of_r(r)  v%of_r(r)
       !               - \int rho%kin_r(r) v%kin_r(r) [for Meta-GGA]
       !               - \sum rho%ns       v%ns       [for LDA+U]
       !               - \sum becsum       D1_Hxc     [for PAW]
       IMPLICIT NONE
       REAL(DP) :: delta_e, delta_e_hub
       !
       delta_e = - SUM( rho%of_r(:,:)*v%of_r(:,:) )
       !
       IF ( dft_is_meta() ) &
          delta_e = delta_e - SUM( rho%kin_r(:,:)*v%kin_r(:,:) )
       !
       delta_e = omega * delta_e / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
       !
#ifdef __BANDS
       CALL mp_sum( delta_e, intra_bgrp_comm )
#else
       CALL mp_sum( delta_e, intra_pool_comm )
#endif
       !
       if (lda_plus_u) then
          delta_e_hub = - SUM (rho%ns(:,:,:,:)*v%ns(:,:,:,:))
          if (nspin==1) delta_e_hub = 2.d0 * delta_e_hub
          delta_e = delta_e + delta_e_hub
       end if
       !
       IF (okpaw) delta_e = delta_e - SUM(ddd_paw(:,:,:)*rho%bec(:,:,:))
       !
       RETURN
       !
     END FUNCTION delta_e
     !
     !-----------------------------------------------------------------------
     FUNCTION delta_escf()
       !-----------------------------------------------------------------------
       !
       ! ... delta_escf = - \int \delta rho%of_r(r)  v%of_r(r)
       !                  - \int \delta rho%kin_r(r) v%kin_r(r) [for Meta-GGA]
       !                  - \sum \delta rho%ns       v%ns       [for LDA+U]
       !                  - \sum \delta becsum       D1         [for PAW]
       ! ... calculates the difference between the Hartree and XC energy
       ! ... at first order in the charge density difference \delta rho(r)
       IMPLICIT NONE
       !
       REAL(DP) :: delta_escf, delta_escf_hub
       !
       delta_escf = - SUM( ( rhoin%of_r(:,:)-rho%of_r(:,:) )*v%of_r(:,:) )
       !
       IF ( dft_is_meta() ) &
          delta_escf = delta_escf - &
                       SUM( (rhoin%kin_r(:,:)-rho%kin_r(:,:) )*v%kin_r(:,:))
       !
       delta_escf = omega * delta_escf / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
       !
#ifdef __BANDS
       CALL mp_sum( delta_escf, intra_bgrp_comm )
#else
       CALL mp_sum( delta_escf, intra_pool_comm )
#endif
       !
       if (lda_plus_u) then
          delta_escf_hub = - SUM((rhoin%ns(:,:,:,:)-rho%ns(:,:,:,:))*v%ns(:,:,:,:))
          if (nspin==1) delta_escf_hub = 2.d0 * delta_escf_hub
          delta_escf = delta_escf + delta_escf_hub
       end if

       IF (okpaw) delta_escf = delta_escf - &
                               SUM(ddd_paw(:,:,:)*(rhoin%bec(:,:,:)-rho%bec(:,:,:)))

       RETURN
       !
     END FUNCTION delta_escf
     !
     !-----------------------------------------------------------------------
     SUBROUTINE calc_pol ( en_el )
       !-----------------------------------------------------------------------
       !
       USE kinds,     ONLY : DP
       USE constants, ONLY : pi
       USE bp,        ONLY : lelfield, ion_pol, el_pol, fc_pol, l_el_pol_old, &
                             el_pol_acc, el_pol_old, efield, l3dstring, gdir, &
                             transform_el, efield_cart
       !
       IMPLICIT NONE
       REAL (DP), INTENT(out) :: en_el
       !
       INTEGER :: i, j 
       REAL(DP):: sca, el_pol_cart(3),  el_pol_acc_cart(3)
       !
       IF (.not.l3dstring) THEN
          CALL c_phase_field(el_pol(gdir),ion_pol(gdir),fc_pol(gdir),gdir)
          if (.not.l_el_pol_old) then
             l_el_pol_old=.true.
             el_pol_old(gdir)=el_pol(gdir)
             en_el=-efield*(el_pol(gdir)+ion_pol(gdir))
             el_pol_acc(gdir)=0.d0
          else
             sca=(el_pol(gdir)-el_pol_old(gdir))/fc_pol(gdir)
             if(sca < - pi) then
                el_pol_acc(gdir)=el_pol_acc(gdir)+2.d0*pi*fc_pol(gdir)
             else if(sca > pi) then
                el_pol_acc(gdir)=el_pol_acc(gdir)-2.d0*pi*fc_pol(gdir)
             endif
             en_el=-efield*(el_pol(gdir)+ion_pol(gdir)+el_pol_acc(gdir))
             el_pol_old=el_pol
          endif
       ELSE
          do i=1,3
            CALL c_phase_field(el_pol(i),ion_pol(i),fc_pol(i),i)
          enddo
          el_pol_cart(:)=0.d0
          do i=1,3
             do j=1,3
                !el_pol_cart(i)=el_pol_cart(i)+transform_el(j,i)*el_pol(j)
                el_pol_cart(i)=el_pol_cart(i)+at(i,j)*el_pol(j)/(dsqrt(at(1,j)**2.d0+at(2,j)**2.d0+at(3,j)**2.d0))
             enddo
          enddo

          write(stdout,'( "Electronic Dipole on Cartesian axes" )')
          do i=1,3
             write(stdout,*) i, el_pol_cart(i)
          enddo

          write(stdout,'( "Ionic Dipole on Cartesian axes" )')
          do i=1,3
             write(stdout,*) i, ion_pol(i)
          enddo

          if(.not.l_el_pol_old) then
             l_el_pol_old=.true.
             el_pol_old(:)=el_pol(:)
             en_el=0.d0
             do i=1,3
                en_el=en_el-efield_cart(i)*(el_pol_cart(i)+ion_pol(i))
             enddo
             el_pol_acc(:)=0.d0
          else
             do i=1,3
                sca=(el_pol(i)-el_pol_old(i))/fc_pol(i)
                if(sca < - pi) then
                   el_pol_acc(i)=el_pol_acc(i)+2.d0*pi*fc_pol(i)
                else if(sca > pi) then
                   el_pol_acc(i)=el_pol_acc(i)-2.d0*pi*fc_pol(i)
                endif
             enddo
             el_pol_acc_cart(:)=0.d0
             do i=1,3
                do j=1,3
                   el_pol_acc_cart(i)=el_pol_acc_cart(i)+transform_el(j,i)*el_pol_acc(j)
                enddo
             enddo
             en_el=0.d0
             do i=1,3
                en_el=en_el-efield_cart(i)*(el_pol_cart(i)+ion_pol(i)+el_pol_acc_cart(i))
             enddo
             el_pol_old(:)=el_pol(:)
          endif
       ENDIF
       !
     END SUBROUTINE calc_pol
     !
END SUBROUTINE electrons
