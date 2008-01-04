!
! Copyright (C) 2001-2006 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
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
  USE constants,            ONLY : eps8
  USE io_global,            ONLY : stdout, ionode
  USE cell_base,            ONLY : at, bg, alat, omega, tpiba2
  USE ions_base,            ONLY : zv, nat, nsp, ityp, tau
  USE basis,                ONLY : startingpot
  USE gvect,                ONLY : ngm, gstart, nr1, nr2, nr3, nrx1, nrx2, &
                                   nrx3, nrxx, nl, nlm, g, gg, ecutwfc, gcutm
  USE gsmooth,              ONLY : doublegrid, ngms
  USE klist,                ONLY : xk, wk, nelec, ngk, nks, nkstot, lgauss
  USE lsda_mod,             ONLY : lsda, nspin, magtot, absmag, isk
  USE vlocal,               ONLY : strf
  USE wvfct,                ONLY : nbnd, et, npwx
  USE ener,                 ONLY : etot, hwf_energy, eband, deband, ehart, &
                                   vtxc, etxc, etxcc, ewld, demet
  USE scf,                  ONLY : scf_type, scf_type_COPY, &
                                   create_scf_type, destroy_scf_type, &
                                   rho, rho_core, rhog_core, &
                                   v, vltot, vrs, kedtau, vnew
  USE control_flags,        ONLY : mixing_beta, tr2, ethr, niter, nmix, &
                                   iprint, istep, lscf, lmd, conv_elec, &
                                   restart, io_level, assume_isolated,  &
                                   gamma_only, iverbosity
  USE io_files,             ONLY : iunwfc, iunocc, nwordwfc, output_drho, &
                                   iunefield
  USE buffers,              ONLY : save_buffer
  USE ldaU,                 ONLY : eth, Hubbard_U, Hubbard_lmax, &
                                   niter_with_fixed_ns, lda_plus_u
  USE extfield,             ONLY : tefield, etotefield
  USE wavefunctions_module, ONLY : evc, psic
  USE noncollin_module,     ONLY : noncolin, magtot_nc, i_cons,  bfield, &
                                   lambda, report
  USE spin_orb,             ONLY : domag
  USE bp,                   ONLY : lelfield
  USE io_rho_xml,           ONLY : write_rho
  USE uspp,                 ONLY : okvan
#if defined (EXX)
  USE exx,                  ONLY : exxinit, init_h_wfc, exxenergy, exxenergy2
  USE funct,                ONLY : dft_is_hybrid, exx_is_active
#endif
  USE funct,                ONLY : dft_is_meta
  USE mp_global,            ONLY : intra_pool_comm, npool
  USE mp,                   ONLY : mp_sum
  !
  USE paw_variables,        ONLY : okpaw, ddd_paw
  USE paw_onecenter,        ONLY : PAW_potential
  USE uspp_param,           ONLY : nh, nhm ! used for PAW
  !
  IMPLICIT NONE
  !
  ! ... a few local variables
  !
#if defined (EXX)
  REAL(DP) :: dexx
  REAL(DP) :: fock0, fock1, fock2
#endif
  REAL(DP) :: &
      dr2,          &! the norm of the diffence between potential
      charge,       &! the total charge
      deband_hwf,   &! deband for the Harris-Weinert-Foulkes functional
      mag           ! local magnetization
  INTEGER :: &
      i,            &! counter on polarization
      is,           &! counter on spins
      ik,           &! counter on k points
      idum,         &! dummy counter on iterations
      iter,         &! counter on iterations
      ik_,          &! used to read ik from restart file
      kilobytes
  REAL(DP) :: &
      tr2_min,     &! estimated error on energy coming from diagonalization
      descf         ! correction for variational energy
  LOGICAL :: &
      exst, first
  !
  ! ... auxiliary variables for calculating and storing temporary copies of
  ! ... the charge density and of the HXC-potential
  !
  type (scf_type), save :: rhoin ! used to store rho_in of current/next iteration
  !
  ! ... external functions
  !
  REAL(DP), EXTERNAL :: ewald, get_clock
  !
  ! ... additional variables for PAW
  REAL (DP) :: e_PAW ! deband, descf and E corrections from PAW
  REAL (DP) :: correction1c                 ! total PAW correction
  !
  iter = 0
  ik_  = 0
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
     !
     RETURN
     !
  END IF
  !
  CALL start_clock( 'electrons' )
  !
  ! ... calculates the ewald contribution to total energy
  !
  ewld = ewald( alat, nat, nsp, ityp, zv, at, bg, tau, &
                omega, g, gg, ngm, gcutm, gstart, gamma_only, strf )
  !
  call create_scf_type ( rhoin )

#if defined (EXX)
10 CONTINUE
#endif
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
        !
        ethr = MIN( ethr, 0.1D0*dr2 / MAX( 1.D0, nelec ) )
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
        ! PAW : sum_band computes new becsum (in uspp modules)
        CALL sum_band()
        !
        ! ... the Harris-Weinert-Foulkes energy is computed here using only
        ! ... quantities obtained from the input density
        !
        hwf_energy = eband + deband_hwf + (etxc - etxcc) + ewld + ehart + demet
        IF ( lda_plus_u ) hwf_energy = hwf_energy + eth
        !
        IF ( lda_plus_u )  THEN
           !
           IF ( iverbosity > 0 .OR. first ) CALL write_ns()
           !
           IF ( first .AND. istep == 0 .AND. &
                startingpot == 'atomic' ) CALL ns_adj()
           IF ( iter <= niter_with_fixed_ns ) rho%ns = rhoin%ns
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
           IF (okpaw) CALL PAW_potential(rhoin%bec, ddd_paw, e_PAW)
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
           !
           CALL write_rho( rho%of_r, nspin )
           IF ( lda_plus_u ) THEN
              IF ( ionode ) THEN
                 CALL seqopn( iunocc, 'occup', 'FORMATTED', exst )
                 WRITE( iunocc, * ) rho%ns
                 CLOSE( UNIT = iunocc, STATUS = 'KEEP' )
              END IF
           END IF
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
           if (okpaw) CALL PAW_potential(rho%bec, ddd_PAW, e_PAW)
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
#if defined (EXX)
        IF ( exx_is_active() ) THEN
           !
           fock1 = exxenergy2()
           fock2 = fock0
           !
        ELSE
           !
           fock0 = 0.D0
           !
        END IF
#endif
        !
        ! ... if we didn't cycle before we can exit the do-loop
        !
        EXIT scf_step
        !
     END DO scf_step
     !
     ! ... define the total local potential (external + scf)
     !
     CALL set_vrs( vrs, vltot, v%of_r, kedtau, v%kin_r, nrxx, nspin, doublegrid )
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
     IF ( lelfield ) CALL c_phase_field()
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
        IF ( lda_plus_U .AND. iverbosity == 0 ) CALL write_ns ( )
        CALL print_ks_energies()
        !
     END IF
     !
     IF ( ABS( charge - nelec ) / charge > 1.D-7 ) THEN
        WRITE( stdout, 9050 ) charge, nelec
        IF ( ABS( charge - nelec ) / charge > 1.D-3 ) &
           CALL errore( 'electrons', 'charge is wrong', 1 )
     END IF
     !
     etot = eband + ( etxc - etxcc ) + ewld + ehart + deband + demet + descf
     !
     IF (okpaw) THEN
        correction1c = e_PAW
        etot       = etot       + correction1c
        hwf_energy = hwf_energy + correction1c
     END IF
     !
#if defined (EXX)
     !
     etot = etot - 0.5D0*fock0
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
           CALL v_of_rho( rho, rho_core,rhog_core, &
                          ehart, etxc, vtxc, eth, etotefield, charge, v)
           IF (okpaw) CALL PAW_potential(rho%bec, ddd_PAW, e_PAW)
           !
           CALL set_vrs( vrs, vltot, v%of_r, kedtau, v%kin_r, nrxx, nspin, doublegrid )
           !
           WRITE( stdout, * ) " NOW GO BACK TO REFINE HYBRID CALCULATION"
           WRITE( stdout, * ) fock0
           !
           iter = 0
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
        !
        WRITE( stdout, * ) fock0, fock1, fock2
        WRITE( stdout, 9066 ) dexx
        !
        fock0 = fock2
        !
     END IF
     !
#endif
     !
     IF ( lda_plus_u ) etot = etot + eth
     IF ( tefield ) THEN
        etot = etot + etotefield
        hwf_energy = hwf_energy + etotefield
     END IF
     !
     IF ( ( conv_elec .OR. MOD( iter, iprint ) == 0 ) .AND. .NOT. lmd ) THEN
        !
        IF ( dr2 > eps8 ) THEN
           WRITE( stdout, 9081 ) etot, hwf_energy, dr2
        ELSE
           WRITE( stdout, 9083 ) etot, hwf_energy, dr2
        END IF
        !
        WRITE( stdout, 9060 ) &
            ( eband + deband ), ehart, ( etxc - etxcc ), ewld
        !
#if defined (EXX)
        !
        WRITE( stdout, 9062 ) fock1
        WRITE( stdout, 9063 ) fock2
        WRITE( stdout, 9064 ) 0.5D0*fock2
        !
#endif
        !
        IF ( tefield ) WRITE( stdout, 9061 ) etotefield
        IF ( lda_plus_u ) WRITE( stdout, 9065 ) eth
        IF ( ABS (descf) > eps8 ) WRITE( stdout, 9069 ) descf
        IF( okpaw) WRITE( stdout, 9067 ) correction1c
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
        !
     END IF
     !
     IF ( lsda ) WRITE( stdout, 9017 ) magtot, absmag
     !
     IF ( noncolin .AND. domag ) &
        WRITE( stdout, 9018 ) magtot_nc(1:3), absmag
     !
     IF ( i_cons == 3 .OR. i_cons == 4 )  &
        WRITE( stdout, 9071 ) bfield(1), bfield(2), bfield(3)
     IF ( i_cons == 5 ) &
        WRITE( stdout, 9072 ) bfield(3)
     IF ( i_cons /= 0 .AND. i_cons < 4 ) &
        WRITE( stdout, 9073 ) lambda
     !
     CALL flush_unit( stdout )
     !
     IF ( conv_elec ) THEN
        !
#if defined (EXX)
        !
        IF ( dft_is_hybrid() .AND. dexx > tr2 ) THEN
           !
           WRITE (stdout,*) " NOW GO BACK TO REFINE HYBRID CALCULATION"
           !
           iter = 0
           !
           GO TO 10
           !
        END IF
#endif
        !
        ! ... if system is charged add a Makov-Payne correction to the energy
        !
        IF ( assume_isolated ) CALL makov_payne( etot )
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
9000 FORMAT(/'     total cpu time spent up to now is ',F9.2,' secs' )
9001 FORMAT(/'     per-process dynamical memory: ',f7.1,' Mb' )
9002 FORMAT(/'     Self-consistent Calculation' )
9010 FORMAT(/'     iteration #',I3,'     ecut=', F9.2,' Ry',5X,'beta=',F4.2 )
9017 FORMAT(/'     total magnetization       =', F9.2,' Bohr mag/cell', &
            /'     absolute magnetization    =', F9.2,' Bohr mag/cell' )
9018 FORMAT(/'     total magnetization       =',3F9.2,' Bohr mag/cell' &
       &   ,/'     absolute magnetization    =', F9.2,' Bohr mag/cell' )
9050 FORMAT(/'     WARNING: integrated charge=',F15.8,', expected=',F15.8 )
9060 FORMAT(/'     The total energy is the sum of the following terms:',/,&
            /'     one-electron contribution =',F15.8,' Ry' &
            /'     hartree contribution      =',F15.8,' Ry' &
            /'     xc contribution           =',F15.8,' Ry' &
            /'     ewald contribution        =',F15.8,' Ry' )
9061 FORMAT( '     electric field correction =',F15.8,' Ry' )
9062 FORMAT( '     Fock energy 1             =',F15.8,' Ry' )
9063 FORMAT( '     Fock energy 2             =',F15.8,' Ry' )
9064 FORMAT( '     Half Fock energy 2        =',F15.8,' Ry' )
9065 FORMAT( '     Hubbard energy            =',F15.8,' Ry' )
9066 FORMAT( '     dexx                      =',F15.8,' Ry' )
9067 FORMAT( '     one-center paw contrib.   =',F15.8,' Ry' )
9069 FORMAT( '     scf correction            =',F15.8,' Ry' )
9070 FORMAT( '     smearing contrib. (-TS)   =',F15.8,' Ry' )
9071 FORMAT( '     Magnetic field            =',3F12.7,' Ry' )
9072 FORMAT( '     Magnetic field            =',F12.7, ' Ry' )
9073 FORMAT( '     lambda                    =',F11.2,' Ry' )
9080 FORMAT(/'     total energy              =',0PF15.8,' Ry' &
            /'     Harris-Foulkes estimate   =',0PF15.8,' Ry' &
            /'     estimated scf accuracy    <',0PF15.8,' Ry' )
9081 FORMAT(/'!    total energy              =',0PF15.8,' Ry' &
            /'     Harris-Foulkes estimate   =',0PF15.8,' Ry' &
            /'     estimated scf accuracy    <',0PF15.8,' Ry' )
9082 FORMAT(/'     total energy              =',0PF15.8,' Ry' &
            /'     Harris-Foulkes estimate   =',0PF15.8,' Ry' &
            /'     estimated scf accuracy    <',1PE15.1,' Ry' )
9083 FORMAT(/'!    total energy              =',0PF15.8,' Ry' &
            /'     Harris-Foulkes estimate   =',0PF15.8,' Ry' &
            /'     estimated scf accuracy    <',1PE15.1,' Ry' )
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
          DO ir = 1, nrxx
             !
             mag = rho%of_r(ir,1) - rho%of_r(ir,2)
             !
             magtot = magtot + mag
             absmag = absmag + ABS( mag )
             !
          END DO
          !
          magtot = magtot * omega / ( nr1*nr2*nr3 )
          absmag = absmag * omega / ( nr1*nr2*nr3 )
          !
          CALL mp_sum( magtot, intra_pool_comm )
          CALL mp_sum( absmag, intra_pool_comm )
          !
       ELSE IF ( noncolin ) THEN
          !
          magtot_nc = 0.D0
          absmag    = 0.D0
          !
          DO ir = 1,nrxx
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
          CALL mp_sum( magtot_nc, intra_pool_comm )
          CALL mp_sum( absmag,    intra_pool_comm )
          !
          DO i = 1, 3
             !
             magtot_nc(i) = magtot_nc(i) * omega / ( nr1*nr2*nr3 )
             !
          END DO
          !
          absmag = absmag * omega / ( nr1*nr2*nr3 )
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
       USE control_flags, ONLY : lpath
       USE check_stop,    ONLY : global_check_stop_now => check_stop_now
       USE io_files,      ONLY : iunpath
       !
       IMPLICIT NONE
       !
       LOGICAL :: check_stop_now
       INTEGER :: unit
       !
       unit = stdout
       !
       IF ( lpath ) unit = iunpath
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
       REAl(DP):: delta_e_paw
       !
       delta_e = - SUM( rho%of_r(:,:)*v%of_r(:,:) )
       !
       IF ( dft_is_meta() ) &
          delta_e = delta_e - SUM( rho%kin_r(:,:)*v%kin_r(:,:) )
       !
       delta_e = omega * delta_e / ( nr1*nr2*nr3 )
       !
       CALL mp_sum( delta_e, intra_pool_comm )
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
       delta_escf = omega * delta_escf / ( nr1*nr2*nr3 )
       !
       CALL mp_sum( delta_escf, intra_pool_comm )
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
END SUBROUTINE electrons
