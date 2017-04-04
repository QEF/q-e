!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
! TB
! included monopole related energy
!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------
SUBROUTINE electrons()
  !----------------------------------------------------------------------------
  !
  ! ... General self-consistency loop, also for hybrid functionals
  ! ... For non-hybrid functionals it just calls "electron_scf"
  !
  USE kinds,                ONLY : DP
  USE check_stop,           ONLY : check_stop_now, stopped_by_user
  USE io_global,            ONLY : stdout, ionode
  USE fft_base,             ONLY : dfftp
  USE gvecs,                ONLY : doublegrid
  USE gvect,                ONLY : ecutrho
  USE lsda_mod,             ONLY : nspin, magtot, absmag
  USE ener,                 ONLY : etot, hwf_energy, eband, deband, ehart, &
                                   vtxc, etxc, etxcc, ewld, demet, epaw, &
                                   elondon, ef_up, ef_dw
  USE scf,                  ONLY : rho, rho_core, rhog_core, v, vltot, vrs, &
                                   kedtau, vnew
  USE control_flags,        ONLY : tr2, niter, conv_elec, restart, lmd, &
                                   do_makov_payne
  USE io_files,             ONLY : iunmix, output_drho, &
                                   iunres, iunefield, seqopn
  USE ldaU,                 ONLY : eth
  USE extfield,             ONLY : tefield, etotefield
  USE wvfct,                ONLY : nbnd, wg, et
  USE klist,                ONLY : nks
  USE noncollin_module,     ONLY : noncolin, magtot_nc, i_cons,  bfield, &
                                   lambda, report
  USE uspp,                 ONLY : okvan
  USE exx,                  ONLY : exxinit, exxenergy2, exxenergy, exxbuff, &
                                   fock0, fock1, fock2, dexx, use_ace
  USE funct,                ONLY : dft_is_hybrid, exx_is_active
  USE control_flags,        ONLY : adapt_thr, tr2_init, tr2_multi, gamma_only
  !
  USE paw_variables,        ONLY : okpaw, ddd_paw, total_core_energy, only_paw
  USE paw_onecenter,        ONLY : PAW_potential
  USE paw_symmetry,         ONLY : PAW_symmetrize_ddd
  USE ions_base,            ONLY : nat
  !
  !
  IMPLICIT NONE
  !
  ! ... a few local variables
  !
  REAL(DP) :: &
      charge,       &! the total charge
      ee, exxen      ! used to compute exchange energy
  REAL(dp), EXTERNAL :: exxenergyace
  INTEGER :: &
      idum,         &! dummy counter on iterations
      iter,         &! counter on iterations
      printout,     &
      ik, ios
  REAL(DP) :: &
      tr2_min,     &! estimated error on energy coming from diagonalization
      tr2_final     ! final threshold for exx minimization 
                    ! when using adaptive thresholds.
  LOGICAL :: first, exst
  REAL(DP) :: etot_cmp_paw(nat,2,2)
  !
  !
  exxen = 0.0d0
  iter = 0
  first = .true.
  tr2_final = tr2
  IF ( dft_is_hybrid() ) THEN
     !printout = 0  ! do not print etot and energy components at each scf step
     printout = 1  ! print etot, not energy components at each scf step
  ELSE IF ( lmd ) THEN
     printout = 1  ! print etot, not energy components at each scf step
  ELSE
     printout = 2  ! print etot and energy components at each scf step
  END IF
  IF (dft_is_hybrid() .AND. adapt_thr ) tr2= tr2_init
  fock0 = 0.D0
  fock1 = 0.D0
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
        ELSE IF ( iter < 0 .OR. iter > niter ) THEN
           iter = 0
        ELSE 
           READ (iunres, *) fock0, fock1, fock2
           ! FIXME: et and wg should be read from xml file
           READ (iunres, *) (wg(1:nbnd,ik),ik=1,nks)
           READ (iunres, *) (et(1:nbnd,ik),ik=1,nks)
           CLOSE ( unit=iunres, status='delete')
           ! ... if restarting here, exx was already active
           ! ... initialize stuff for exx
           first = .false.
           CALL exxinit()
           ! FIXME: ugly hack, overwrites exxbuffer from exxinit
           CALL seqopn (iunres, 'restart_exx', 'unformatted', exst)
           IF (exst) READ (iunres, iostat=ios) exxbuff
           IF (ios /= 0) WRITE(stdout,'(5x,"Error in EXX restart!")')
           !
           CALL v_of_rho( rho, rho_core, rhog_core, &
               ehart, etxc, vtxc, eth, etotefield, charge, v)
           IF (okpaw) CALL PAW_potential(rho%bec, ddd_PAW, epaw,etot_cmp_paw)
           CALL set_vrs( vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, &
                         nspin, doublegrid )
           !
           WRITE(stdout,'(5x,"Calculation (EXX) restarted from iteration #", &
                        & i6)') iter 
        END IF
     END IF
     CLOSE ( unit=iunres, status='delete')
  END IF
  !
  DO idum=1,niter
     !
     iter = iter + 1
     !
     ! ... Self-consistency loop. For hybrid functionals the exchange potential
     ! ... is calculated with the orbitals at previous step (none at first step)
     !
     CALL electrons_scf ( printout, exxen )
     !
     IF ( .NOT. dft_is_hybrid() ) RETURN
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
           WRITE (iunres, *) fock0, fock1, fock2
           WRITE (iunres, *) (wg(1:nbnd,ik),ik=1,nks)
           WRITE (iunres, *) (et(1:nbnd,ik),ik=1,nks)
           CLOSE (unit=iunres, status='keep')
           CALL seqopn (iunres, 'restart_exx', 'unformatted', exst)
           WRITE (iunres) exxbuff
           CLOSE (unit=iunres, status='keep')
        END IF
        RETURN
     END IF
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
        CALL exxinit()
        IF ( use_ace) THEN
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
        IF (okpaw) CALL PAW_potential(rho%bec, ddd_PAW, epaw,etot_cmp_paw)
        CALL set_vrs( vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, &
             nspin, doublegrid )
        !
     ELSE
        !
        ! fock1 is the exchange energy calculated for orbitals at step n,
        !       using orbitals at step n-1 in the expression of exchange
        !
        IF ( use_ace) THEN
           fock1 = exxenergyace()
        ELSE
           fock1 = exxenergy2()
        ENDIF
        !
        ! Set new orbitals for the calculation of the exchange term
        !
        CALL exxinit()
        !
        ! fock2 is the exchange energy calculated for orbitals at step n,
        !       using orbitals at step n in the expression of exchange 
        ! fock0 is fock2 at previous step
        !
        fock0 = fock2
        IF ( use_ace) THEN
           fock2 = exxenergyace()
        ELSE
           fock2 = exxenergy2()
        ENDIF
        !
        ! check for convergence. dexx is positive definite: if it isn't,
        ! there is some numerical problem. One such cause could be that
        ! the treatment of the divergence in exact exchange has failed. 
        !
        dexx = fock1 - 0.5D0*(fock0+fock2)
        IF ( dexx < 0d0 ) THEN
!           WRITE(stdout,'(5x,a,1e12.3)') "BEWARE: negative dexx:", dexx
!           dexx = ABS(dexx)
          WRITE( stdout, * ) "dexx:", dexx
          CALL errore( 'electrons', 'dexx is negative! &
           & Check that exxdiv_treatment is appropriate for the system,&
           & or ecutfock may be too low', 1 )
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
        END IF
        IF(dexx>1.d-8)THEN
          WRITE( stdout, 9067 ) dexx
        ELSE
          WRITE( stdout, 9068 ) dexx
        ENDIF
        
        WRITE( stdout, 9062 ) - fock1
        WRITE( stdout, 9064 ) 0.5D0*fock2
        !
        IF ( dexx < tr2_final ) THEN
           IF ( do_makov_payne ) CALL makov_payne( etot )
           WRITE( stdout, 9101 )
           RETURN
        END IF
        !
        IF ( adapt_thr ) THEN
           tr2 = MAX(tr2_multi * dexx, tr2_final)
           WRITE( stdout, 9121 ) tr2
        ENDIF
     ENDIF
     !
     WRITE( stdout,'(/5x,"EXX: now go back to refine exchange calculation")')
     !
     IF ( check_stop_now() ) THEN
        WRITE(stdout,'(5x,"Calculation (EXX) stopped after iteration #", &
                        & i6)') iter
        conv_elec=.FALSE.
        CALL seqopn (iunres, 'restart_e', 'formatted', exst)
        WRITE (iunres, *) iter, tr2, dexx
        WRITE (iunres, *) fock0, fock1, fock2
        ! FIXME: et and wg are written to xml file
        WRITE (iunres, *) (wg(1:nbnd,ik),ik=1,nks)
        WRITE (iunres, *) (et(1:nbnd,ik),ik=1,nks)
        CLOSE (unit=iunres, status='keep')
        RETURN
     END IF
     !
  END DO
  !
  WRITE( stdout, 9120 ) iter
  FLUSH( stdout )
  !
  RETURN
  !
  ! ... formats
  !
9062 FORMAT( '     - averaged Fock potential =',0PF17.8,' Ry' )
9064 FORMAT( '     + Fock energy             =',0PF17.8,' Ry' )
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
  !
  ! ... This routine is a driver of the self-consistent cycle.
  ! ... It uses the routine c_bands for computing the bands at fixed
  ! ... Hamiltonian, the routine sum_band to compute the charge density,
  ! ... the routine v_of_rho to compute the new potential and the routine
  ! ... mix_rho to mix input and output charge densities.
  ! ... If printout > 0, prints on output the total energy;
  ! ... if printout > 1, also prints decomposition into energy contributions
  !
  USE kinds,                ONLY : DP
  USE check_stop,           ONLY : check_stop_now, stopped_by_user
  USE io_global,            ONLY : stdout, ionode
  USE cell_base,            ONLY : at, bg, alat, omega, tpiba2
  USE ions_base,            ONLY : zv, nat, nsp, ityp, tau, compute_eextfor
  USE basis,                ONLY : starting_pot
  USE bp,                   ONLY : lelfield
  USE fft_base,             ONLY : dfftp
  USE gvect,                ONLY : ngm, gstart, nl, nlm, g, gg, gcutm
  USE gvecs,                ONLY : doublegrid, ngms
  USE klist,                ONLY : xk, wk, nelec, ngk, nks, nkstot, lgauss, &
                                   two_fermi_energies, tot_charge
  USE lsda_mod,             ONLY : lsda, nspin, magtot, absmag, isk
  USE vlocal,               ONLY : strf
  USE wvfct,                ONLY : nbnd, et
  USE gvecw,                ONLY : ecutwfc
  USE ener,                 ONLY : etot, hwf_energy, eband, deband, ehart, &
                                   vtxc, etxc, etxcc, ewld, demet, epaw, &
                                   elondon, ef_up, ef_dw, exdm, ef
  USE scf,                  ONLY : scf_type, scf_type_COPY, bcast_scf_type,&
                                   create_scf_type, destroy_scf_type, &
                                   open_mix_file, close_mix_file, &
                                   rho, rho_core, rhog_core, v, vltot, vrs, &
                                   kedtau, vnew
  USE control_flags,        ONLY : mixing_beta, tr2, ethr, niter, nmix, &
                                   iprint, conv_elec, &
                                   restart, io_level, do_makov_payne,  &
                                   gamma_only, iverbosity, textfor,     &
                                   llondon, scf_must_converge, lxdm, ts_vdw
  USE control_flags,        ONLY : n_scf_steps, scf_error

  USE io_files,             ONLY : iunmix, output_drho, &
                                   iunres, iunefield, seqopn
  USE ldaU,                 ONLY : eth, Hubbard_U, Hubbard_lmax, &
                                   niter_with_fixed_ns, lda_plus_u
  USE extfield,             ONLY : tefield, etotefield, monopole, etotmonofield !TB
  USE noncollin_module,     ONLY : noncolin, magtot_nc, i_cons,  bfield, &
                                   lambda, report
  USE spin_orb,             ONLY : domag
  USE io_rho_xml,           ONLY : write_rho
  USE uspp,                 ONLY : okvan
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp_pools,             ONLY : root_pool, my_pool_id, inter_pool_comm
  USE mp,                   ONLY : mp_sum, mp_bcast
  !
  USE london_module,        ONLY : energy_london
  USE xdm_module,           ONLY : energy_xdm
  USE tsvdw_module,         ONLY : EtsvdW
  !
  USE paw_variables,        ONLY : okpaw, ddd_paw, total_core_energy, only_paw
  USE paw_onecenter,        ONLY : PAW_potential
  USE paw_symmetry,         ONLY : PAW_symmetrize_ddd
  USE dfunct,               ONLY : newd
  USE esm,                  ONLY : do_comp_esm, esm_printpot, esm_ewald
  USE fcp_variables,        ONLY : lfcpopt, lfcpdyn
  USE iso_c_binding,        ONLY : c_int
  !
  USE plugin_variables,     ONLY : plugin_etot
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT (IN) :: printout
  REAL(DP),INTENT (IN) :: exxen    ! current estimate of the exchange energy
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
      ios
  INTEGER(kind=c_int) :: kilobytes
  REAL(DP) :: &
      tr2_min,     &! estimated error on energy coming from diagonalization
      descf,       &! correction for variational energy
      en_el=0.0_DP,&! electric field contribution to the total energy
      eext=0.0_DP   ! external forces contribution to the total energy
  LOGICAL :: &
      first, exst
  !
  ! ... auxiliary variables for calculating and storing temporary copies of
  ! ... the charge density and of the HXC-potential
  !
  type (scf_type) :: rhoin ! used to store rho_in of current/next iteration
  !
  ! ... external functions
  !
  REAL(DP), EXTERNAL :: ewald, get_clock
  REAL(DP) :: etot_cmp_paw(nat,2,2)
  !
  iter = 0
  dr2  = 0.0_dp
  IF ( restart ) CALL restart_in_electrons (iter, dr2, ethr, et )
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
  END IF
  !
  IF ( llondon ) THEN
     elondon = energy_london ( alat , nat , ityp , at ,bg , tau )
  ELSE
     elondon = 0.d0
  END IF
  !
  call create_scf_type ( rhoin )
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
     END IF
     iter = iter + 1
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
        IF ( lelfield ) THEN
           CALL c_bands_efield ( iter )
        ELSE
           CALL c_bands( iter )
        END IF
        !
        IF ( stopped_by_user ) THEN
           conv_elec=.FALSE.
           CALL save_in_electrons (iter-1, dr2, ethr, et )
           GO TO 10
        END IF
        !
        ! ... xk, wk, isk, et, wg are distributed across pools;
        ! ... the first node has a complete copy of xk, wk, isk,
        ! ... while eigenvalues et and weights wg must be
        ! ... explicitely collected to the first node
        ! ... this is done here for et, in sum_band for wg
        !
        CALL poolrecover( et, nbnd, nkstot, nks )
        !
        ! ... the new density is computed here. For PAW:
        ! ... sum_band computes new becsum (stored in uspp modules)
        ! ... and a subtly different copy in rho%bec (scf module)
        !
        CALL sum_band()
        !
        ! ... the Harris-Weinert-Foulkes energy is computed here using only
        ! ... quantities obtained from the input density
        !
        hwf_energy = eband + deband_hwf + (etxc - etxcc) + ewld + ehart + demet
        If ( okpaw ) hwf_energy = hwf_energy + epaw
        IF ( lda_plus_u ) hwf_energy = hwf_energy + eth
        !
        IF ( lda_plus_u )  THEN
           !
           IF ( iverbosity > 0 .OR. first ) THEN
              IF (noncolin) THEN
                 CALL write_ns_nc()
              ELSE
                 CALL write_ns()
              ENDIF
           ENDIF
           !
           IF ( first .AND. starting_pot == 'atomic' ) THEN
              CALL ns_adj()
              IF (noncolin) THEN
                 rhoin%ns_nc = rho%ns_nc
              ELSE
                 rhoin%ns = rho%ns
              ENDIF
           END IF
           IF ( iter <= niter_with_fixed_ns ) THEN
              WRITE( stdout, '(/,5X,"RESET ns to initial values (iter <= mixing_fixed_ns)",/)')
              IF (noncolin) THEN
                 rho%ns_nc = rhoin%ns_nc
              ELSE
                 rho%ns = rhoin%ns
              ENDIF
           END IF
           !
        END IF
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
        ! ... Results are broadcast from pool 0 to others to prevent trouble
        ! ... on machines unable to yield the same results from the same 
        ! ... calculation on same data, performed on different procs
        ! ... The mixing should be done on pool 0 only as well, but inside
        ! ... mix_rho there is a call to rho_ddot that in the PAW case 
        ! ... contains a hidden parallelization level on the entire image
        !
        ! IF ( my_pool_id == root_pool ) 
        CALL mix_rho ( rho, rhoin, mixing_beta, dr2, tr2_min, iter, nmix, &
                       iunmix, conv_elec )
        CALL bcast_scf_type ( rhoin, root_pool, inter_pool_comm )
        CALL mp_bcast ( dr2, root_pool, inter_pool_comm )
        CALL mp_bcast ( conv_elec, root_pool, inter_pool_comm )
        !
        if (.not. scf_must_converge .and. idum == niter) conv_elec = .true.
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
        IF ( .NOT. conv_elec ) THEN
           !
           ! ... no convergence yet: calculate new potential from mixed
           ! ... charge density (i.e. the new estimate)
           !
           CALL v_of_rho( rhoin, rho_core, rhog_core, &
                          ehart, etxc, vtxc, eth, etotefield, charge, v)
           IF (okpaw) THEN
              CALL PAW_potential(rhoin%bec, ddd_paw, epaw,etot_cmp_paw)
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
        ELSE 
           !
           ! ... convergence reached:
           ! ... 1) the output HXC-potential is saved in v
           ! ... 2) vnew contains V(out)-V(in) ( used to correct the forces ).
           !
           vnew%of_r(:,:) = v%of_r(:,:)
           CALL v_of_rho( rho,rho_core,rhog_core, &
                          ehart, etxc, vtxc, eth, etotefield, charge, v)
           vnew%of_r(:,:) = v%of_r(:,:) - vnew%of_r(:,:)
           !
           IF (okpaw) THEN
              CALL PAW_potential(rho%bec, ddd_paw, epaw,etot_cmp_paw)
              CALL PAW_symmetrize_ddd(ddd_paw)
           ENDIF
           !
           ! ... note that rho is here the output, not mixed, charge density
           ! ... so correction for variational energy is no longer needed
           !
           descf = 0._dp
           !
        END IF 
        !
        ! ... if we didn't cycle before we can exit the do-loop
        !
        EXIT scf_step
        !
     END DO scf_step
     !
     plugin_etot = 0.0_dp
     !
     CALL plugin_scf_energy(plugin_etot,rhoin)
     !
     CALL plugin_scf_potential(rhoin,conv_elec,dr2)
     !
     ! ... define the total local potential (external + scf)
     !
     CALL sum_vrs( dfftp%nnr, nspin, vltot, v%of_r, vrs )
     !
     ! ... interpolate the total local potential
     !
     CALL interpolate_vrs( dfftp%nnr, nspin, doublegrid, kedtau, v%kin_r, vrs )
     !
     ! ... in the US case we have to recompute the self-consistent
     ! ... term in the nonlocal potential
     ! ... PAW: newd contains PAW updates of NL coefficients
     !
     CALL newd()
     !
     IF ( lelfield ) en_el =  calc_pol ( )
     !
     IF ( ( MOD(iter,report) == 0 ) .OR. ( report /= 0 .AND. conv_elec ) ) THEN
        !
        IF ( (noncolin .AND. domag) .OR. i_cons==1 .OR. nspin==2) CALL report_mag()
        !
     END IF
     !
     WRITE( stdout, 9000 ) get_clock( 'PWSCF' )
     !
     IF ( conv_elec ) WRITE( stdout, 9101 )
!  these values are assigned to global variables  because these information are needed for XML  printout 
!  P.D. 
     IF ( conv_elec ) THEN 
           scf_error = dr2
           n_scf_steps = iter
     END IF  

     !
     IF ( conv_elec .OR. MOD( iter, iprint ) == 0 ) THEN
        !
        IF ( lda_plus_U .AND. iverbosity == 0 ) THEN
           IF (noncolin) THEN
              CALL write_ns_nc()
           ELSE
              CALL write_ns()
           ENDIF
        ENDIF
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
     END IF
     IF (llondon) THEN
        etot = etot + elondon
        hwf_energy = hwf_energy + elondon
     END IF
     ! calculate the xdm energy contribution with converged density
     if (lxdm .and. conv_elec) then
        exdm = energy_xdm()
        etot = etot + exdm
        hwf_energy = hwf_energy + exdm
     end if
     IF (ts_vdw) THEN
        ! factor 2 converts from Ha to Ry units
        etot = etot + 2.0d0*EtsvdW
        hwf_energy = hwf_energy + 2.0d0*EtsvdW
     END IF
     !
     IF ( tefield ) THEN
        etot = etot + etotefield
        hwf_energy = hwf_energy + etotefield
     END IF
     ! TB monopole energy
     IF ( monopole) THEN
        etot = etot + etotmonofield
        hwf_energy = hwf_energy + etotmonofield
     END IF
     !
     IF ( lfcpopt .or. lfcpdyn ) THEN
        etot = etot + ef * tot_charge
        hwf_energy = hwf_energy + ef * tot_charge
     ENDIF
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
        IF ( do_makov_payne .AND. printout/= 0 ) CALL makov_payne( etot )
        !
        ! ... print out ESM potentials if desired
        !
        IF ( do_comp_esm ) CALL esm_printpot()
        !
        WRITE( stdout, 9110 ) iter
        !
        ! ... jump to the end
        !
        GO TO 10
        !
     END IF
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
  END DO
  !
  WRITE( stdout, 9101 )
  WRITE( stdout, 9120 ) iter
  !
10  FLUSH( stdout )
  !
  ! ... exiting: write (unless disabled) the charge density to file
  ! ... (also write ldaU ns coefficients and PAW becsum)
  !
  IF ( io_level > -1 ) CALL write_rho( rho, nspin )
  !
  ! ... delete mixing info if converged, keep it if not
  !
  IF ( conv_elec ) THEN
     CALL close_mix_file( iunmix, 'delete' )
  ELSE
     CALL close_mix_file( iunmix, 'keep' )
  END IF
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
          CALL mp_sum( magtot, intra_bgrp_comm )
          CALL mp_sum( absmag, intra_bgrp_comm )
          !
          IF (two_fermi_energies.and.lgauss) bfield(3)=0.5D0*(ef_up-ef_dw)
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
          CALL mp_sum( magtot_nc, intra_bgrp_comm )
          CALL mp_sum( absmag, intra_bgrp_comm )
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
     FUNCTION delta_e()
       !-----------------------------------------------------------------------
       ! ... delta_e = - \int rho%of_r(r)  v%of_r(r)
       !               - \int rho%kin_r(r) v%kin_r(r) [for Meta-GGA]
       !               - \sum rho%ns       v%ns       [for LDA+U]
       !               - \sum becsum       D1_Hxc     [for PAW]
       USE funct,  ONLY : dft_is_meta
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
       CALL mp_sum( delta_e, intra_bgrp_comm )
       !
       if (lda_plus_u) then
         if (noncolin) then
           delta_e_hub = - SUM (rho%ns_nc(:,:,:,:)*v%ns_nc(:,:,:,:))
           delta_e = delta_e + delta_e_hub
         else
           delta_e_hub = - SUM (rho%ns(:,:,:,:)*v%ns(:,:,:,:))
           if (nspin==1) delta_e_hub = 2.d0 * delta_e_hub
           delta_e = delta_e + delta_e_hub
         endif
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
       !
       USE funct,  ONLY : dft_is_meta
       IMPLICIT NONE
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
       CALL mp_sum( delta_escf, intra_bgrp_comm )
       !
       if (lda_plus_u) then
         if (noncolin) then
           delta_escf_hub = - SUM((rhoin%ns_nc(:,:,:,:)-rho%ns_nc(:,:,:,:))*v%ns_nc(:,:,:,:))
           delta_escf = delta_escf + delta_escf_hub
         else
           delta_escf_hub = - SUM((rhoin%ns(:,:,:,:)-rho%ns(:,:,:,:))*v%ns(:,:,:,:))
           if (nspin==1) delta_escf_hub = 2.d0 * delta_escf_hub
           delta_escf = delta_escf + delta_escf_hub
         endif
       end if

       IF (okpaw) delta_escf = delta_escf - &
                               SUM(ddd_paw(:,:,:)*(rhoin%bec(:,:,:)-rho%bec(:,:,:)))

       RETURN
       !
     END FUNCTION delta_escf
     !
     !-----------------------------------------------------------------------
     FUNCTION calc_pol ( ) RESULT ( en_el )
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
                el_pol_cart(i)=el_pol_cart(i)+at(i,j)*el_pol(j) / &
                               (sqrt(at(1,j)**2.d0+at(2,j)**2.d0+at(3,j)**2.d0))
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
     END FUNCTION calc_pol
     !
     !-----------------------------------------------------------------------
     SUBROUTINE print_energies ( printout )
       !-----------------------------------------------------------------------
       !
       USE constants, ONLY : eps8
       INTEGER, INTENT (IN) :: printout
       !
   
       IF ( printout == 0 ) RETURN
       IF ( ( conv_elec .OR. MOD(iter,iprint) == 0 ) .AND. printout > 1 ) THEN
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
          IF ( lxdm )    WRITE ( stdout , 9075 ) exdm
          IF ( ts_vdw )  WRITE ( stdout , 9076 ) 2.0d0*EtsvdW
          IF ( textfor)  WRITE ( stdout , 9077 ) eext
          IF ( tefield )            WRITE( stdout, 9061 ) etotefield
          IF ( monopole )           WRITE( stdout, 9062 ) etotmonofield ! TB
          IF ( lda_plus_u )         WRITE( stdout, 9065 ) eth
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
          ! ... With Fermi-Dirac population factor, etot is the electronic
          ! ... free energy F = E - TS , demet is the -TS contribution
          !
          IF ( lgauss ) WRITE( stdout, 9070 ) demet
          !
          ! ... With Fictitious charge particle (FCP), etot is the grand
          ! ... potential energy Omega = E - muN, -muN is the potentiostat
          ! ... contribution.
          !
          IF ( lfcpopt .or. lfcpdyn ) WRITE( stdout, 9072 ) ef*tot_charge
          !
       ELSE IF ( conv_elec ) THEN
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
       CALL plugin_print_energies()
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
9060 FORMAT(/'     The total energy is the sum of the following terms:',/,&
            /'     one-electron contribution =',F17.8,' Ry' &
            /'     hartree contribution      =',F17.8,' Ry' &
            /'     xc contribution           =',F17.8,' Ry' &
            /'     ewald contribution        =',F17.8,' Ry' )
9061 FORMAT( '     electric field correction =',F17.8,' Ry' )
9062 FORMAT( '     monopole field correction =',F17.8,' Ry' ) ! TB
9065 FORMAT( '     Hubbard energy            =',F17.8,' Ry' )
9067 FORMAT( '     one-center paw contrib.   =',F17.8,' Ry' )
9068 FORMAT( '      -> PAW hartree energy AE =',F17.8,' Ry' &
            /'      -> PAW hartree energy PS =',F17.8,' Ry' &
            /'      -> PAW xc energy AE      =',F17.8,' Ry' &
            /'      -> PAW xc energy PS      =',F17.8,' Ry' &
            /'      -> total E_H with PAW    =',F17.8,' Ry'& 
            /'      -> total E_XC with PAW   =',F17.8,' Ry' )
9069 FORMAT( '     scf correction            =',F17.8,' Ry' )
9070 FORMAT( '     smearing contrib. (-TS)   =',F17.8,' Ry' )
9071 FORMAT( '     Magnetic field            =',3F12.7,' Ry' )
9072 FORMAT( '     pot.stat. contrib. (-muN) =',F17.8,' Ry' )
9073 FORMAT( '     lambda                    =',F11.2,' Ry' )
9074 FORMAT( '     Dispersion Correction     =',F17.8,' Ry' )
9075 FORMAT( '     Dispersion XDM Correction =',F17.8,' Ry' )
9076 FORMAT( '     Dispersion T-S Correction =',F17.8,' Ry' )
9077 FORMAT( '     External forces energy    =',F17.8,' Ry' )
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

  END SUBROUTINE print_energies
  !
END SUBROUTINE electrons_scf
!
!----------------------------------------------------------------------------
FUNCTION exxenergyace ( )
  !--------------------------------------------------------------------------
  !
  ! ... Compute exchange energy using ACE
  !
  USE kinds,    ONLY : DP
  USE buffers,  ONLY : get_buffer
  USE exx,      ONLY : vexxace_gamma, vexxace_k, domat
  USE klist,    ONLY : nks, ngk, igk_k
  USE wvfct,    ONLY : nbnd, npwx, current_k
  USE lsda_mod, ONLY : lsda, isk, current_spin
  USE io_files, ONLY : iunwfc, nwordwfc
  USE mp_pools, ONLY : inter_pool_comm
  USE mp_bands, ONLY : intra_bgrp_comm
  USE mp,       ONLY : mp_sum
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions_module, ONLY : evc
  !
  IMPLICIT NONE
  !
  REAL (dp) :: exxenergyace  ! computed energy
  !
  REAL (dp) :: ex
  INTEGER :: ik, npw
  !
  domat = .true.
  exxenergyace=0.0_dp
  DO ik = 1, nks
     npw = ngk (ik)
     current_k = ik
     IF ( lsda ) current_spin = isk(ik)
     IF (nks > 1) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
     IF (gamma_only) THEN
        call vexxace_gamma ( npw, nbnd, evc, ex )
     ELSE
        call vexxace_k ( npw, nbnd, evc, ex )
     END IF
     exxenergyace = exxenergyace + ex
  END DO
  CALL mp_sum( exxenergyace, inter_pool_comm )
  domat = .false.
  !
END FUNCTION exxenergyace
