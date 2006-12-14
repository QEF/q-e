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
  USE constants,            ONLY : eps8, rytoev
  USE io_global,            ONLY : stdout, ionode
  USE cell_base,            ONLY : at, bg, alat, omega, tpiba2
  USE ions_base,            ONLY : zv, nat, nsp, ityp, tau
  USE basis,                ONLY : startingpot
  USE gvect,                ONLY : ngm, gstart, nr1, nr2, nr3, nrx1, nrx2, &
                                   nrx3, nrxx, nl, nlm, g, gg, ecutwfc, gcutm
  USE gsmooth,              ONLY : doublegrid, ngms
  USE klist,                ONLY : xk, wk, degauss, nelec, ngk, nks, nkstot, &
                                   lgauss, ngauss, two_fermi_energies, &
                                   nelup, neldw
  USE lsda_mod,             ONLY : lsda, nspin, magtot, absmag, isk
  USE ktetra,               ONLY : ltetra, ntetra, tetra  
  USE vlocal,               ONLY : strf, vnew  
  USE wvfct,                ONLY : nbnd, et, gamma_only, wg, npwx
  USE fixed_occ,            ONLY : f_inp, tfixed_occ
  USE ener,                 ONLY : etot, hwf_energy, eband, deband, ehart, &
                                   vtxc, etxc, etxcc, ewld, demet, ef,     &
                                   ef_up, ef_dw 
  USE scf,                  ONLY : rho, rhog, rho_core, rhog_core, &
                                   vr, vltot, vrs 
  USE control_flags,        ONLY : mixing_beta, tr2, ethr, niter, nmix, &
                                   iprint, istep, lscf, lmd, conv_elec, &
                                   lbands, restart, reduce_io, iverbosity
  USE io_files,             ONLY : prefix, iunwfc, iunocc, nwordwfc, iunpath, &
                                   output_drho, iunefield
  USE ldaU,                 ONLY : ns, nsnew, eth, Hubbard_U, Hubbard_lmax, &
                                   niter_with_fixed_ns, lda_plus_u  
  USE extfield,             ONLY : tefield, etotefield  
  USE wavefunctions_module, ONLY : evc, psic
  USE noncollin_module,     ONLY : noncolin, npol, magtot_nc, factlist, &
                                   pointlist, pointnum, mcons, i_cons,  &
                                   bfield, lambda, vtcon, report
  USE spin_orb,             ONLY : domag
  USE bp,                   ONLY : lelfield, lberry, nberrycyc
  USE io_rho_xml,           ONLY : write_rho
  USE uspp,                 ONLY : okvan
  USE realus,               ONLY : tqr
#if defined (EXX)
  USE exx,                  ONLY : exxinit, init_h_wfc, exxenergy, exxenergy2 
  USE funct,                ONLY : dft_is_hybrid, exx_is_active
#endif
  USE mp_global,            ONLY : root_pool, intra_pool_comm, inter_pool_comm
  USE mp,                   ONLY : mp_sum, mp_bcast
  !
  IMPLICIT NONE
  !
  ! ... a few local variables
  !  
#if defined (EXX)
  REAL(DP) :: dexx
  REAL(DP) :: fock0, fock1, fock2
#endif
  INTEGER, ALLOCATABLE :: &
      ngk_g(:)       ! number of plane waves summed on all nodes
  CHARACTER(LEN=256) :: &
      flmix          !
  REAL(DP) :: &
      dr2,          &! the norm of the diffence between potential
      charge,       &! the total charge
      mag,          &! local magnetization
      ehomo, elumo   ! highest occupied and lowest onuccupied levels
  INTEGER :: &
      i,            &! counter on polarization
      is,           &! counter on spins
      ik,           &! counter on k points
      kbnd,         &! counter on bands
      ibnd_up,      &! counter on bands
      ibnd_dw,    &! counter on bands
      ibnd,         &! counter on bands
      idum,         &! dummy counter on iterations
      iter,         &! counter on iterations
      ik_            ! used to read ik from restart file
  REAL(DP) :: &
       tr2_min,     &! estimated error on energy coming from diagonalization
       descf         ! correction for variational energy
  REAL(DP), ALLOCATABLE :: &
      wg_g(:,:)      ! temporary array used to collect array wg from pools
                     ! and then print occupations on stdout
  LOGICAL :: &
      exst, first
  !
  ! ... auxiliary variables for calculating and storing temporary copies of
  ! ... the charge density and of the HXC-potential
  !
  COMPLEX(DP), ALLOCATABLE :: rhognew(:,:)
  REAL(DP),    ALLOCATABLE :: rhonew(:,:)
  REAL(DP),    ALLOCATABLE :: vaux(:,:)
  !
  ! ... variables needed for electric field calculation
  !
  INTEGER                  :: inberry
  COMPLEX(DP), ALLOCATABLE :: psi(:,:)
  !
  ! ... external functions
  !
  REAL(DP), EXTERNAL :: ewald, get_clock, DNRM2
  !
  !
  CALL start_clock( 'electrons' )
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
        IF ( output_drho /= ' ' ) CALL remove_atomic_rho
        !
        CALL stop_clock( 'electrons' )
        !
        RETURN
        !
     END IF
     !
  END IF
  !
  WRITE( stdout, 9000 ) get_clock( 'PWSCF' )
  !
  CALL flush_unit( stdout )
  !
  IF ( .NOT. lscf ) THEN
     !
     CALL non_scf()
     !
     RETURN
     !
  END IF
  !
  ! ... calculates the ewald contribution to total energy
  !
  ewld = ewald( alat, nat, nsp, ityp, zv, at, bg, tau, &
                omega, g, gg, ngm, gcutm, gstart, gamma_only, strf )
  !               
  IF ( reduce_io ) THEN
     !
     flmix = ' '
     !
  ELSE
     !
     flmix = 'mix'
     !
  END IF
  !
  ! ... Convergence threshold for iterative diagonalization
  !
  ! ... for the first scf iteration of each ionic step (after the first),
  ! ... the threshold is fixed to a default value of 1.D-6
  !
#if defined (EXX)
10 CONTINUE
#endif
  !
  IF ( istep > 0 ) ethr = 1.D-6
  !
  WRITE( stdout, 9001 )
  !
  CALL flush_unit( stdout )
  !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%          iterate !          %%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
  DO idum = 1, niter
     !
     IF ( idum > 1 .AND. check_stop_now() ) RETURN
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
        ethr = MAX( MIN( ethr , ( dr2 / nelec * 0.1D0 ) ) , &
                       ( tr2 / nelec * 0.01D0 ) )
        !
     END IF
     !
     first = ( iter == 1 )
     !
     scf_step: DO 
        !
        ! ... tr2_min is set to an estimate of the error on the energy
        ! ... due to diagonalization - used only in the first scf iteration
        !
        tr2_min = 0.D0
        !
        IF ( first ) tr2_min = nelec * ethr
        !
        ! ... diagonalization of the KS hamiltonian
        !
        IF ( lelfield ) THEN
           !
           DO inberry = 1, nberrycyc
              !
              ALLOCATE( psi( npwx, nbnd ) )
              !
              DO ik = 1, nks
                 !
                 CALL davcio( psi, nwordwfc, iunwfc,    ik, -1 )
                 CALL davcio( psi, nwordwfc, iunefield, ik,  1 )
                 !
              END DO
              !
              DEALLOCATE( psi )
              !
              CALL c_bands( iter, ik_, dr2 )
              !
           END DO
           !
        ELSE
           !
           CALL c_bands( iter, ik_, dr2 )
           !
        END IF
        !
        IF ( iter > 1 .AND. check_stop_now() ) RETURN
        !
        ! ... deband = - \sum_v <\psi_v | V_h + V_xc |\psi_v> is calculated a
        ! ... first time here using the input density and potential ( to be
        ! ... used to calculate the Harris-Weinert-Foulkes energy )
        !
        deband = delta_e()
        !
        ! ... the new density is computed here
        !
        CALL sum_band()
        !
        ! ... bring output charge density (now stored in rho) to G-space
        ! ... (rhognew) for mixing
        !
        ALLOCATE( rhognew( ngm, nspin ) )
        !
        DO is = 1, nspin
           !
           psic(:) = rho(:,is)
           !
           CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1 )
           !
           rhognew(:,is) = psic(nl(:))
           !
           IF ( okvan .AND. tqr ) THEN
              !
              ! ... in case the augmentation charges are computed in real space
              ! ... we apply an FFT filter to the density in real space to
              ! ... remove features that are not compatible with the FFT grid.
              !
              psic(:) = ( 0.D0, 0.D0 )
              !
              psic(nl(:)) = rhognew(:,is)
              !
              IF ( gamma_only ) psic(nlm(:)) = CONJG( rhognew(:,is) )
              !
              CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1 )
              !
              rho(:,is) = psic(:)
              !
           END IF
           !
        END DO
        !
        ! ... the Harris-Weinert-Foulkes energy is computed here using only
        ! ... quantities obtained from the input density
        !
        hwf_energy = eband + deband + ( etxc - etxcc ) + ewld + ehart + demet
        !
        IF ( lda_plus_u )  THEN
           !
           CALL write_ns()
           !
           IF ( first .AND. istep == 0 .AND. &
                startingpot == 'atomic' ) CALL ns_adj()
           !
           IF ( iter <= niter_with_fixed_ns ) nsnew = ns 
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
        CALL mix_rho( rhognew, rhog, nsnew, ns, mixing_beta, &
                      dr2, tr2_min, iter, nmix, flmix, conv_elec )
        !
        ! ... if convergence is achieved or if the self-consistency error
        ! ... (dr2) is smaller than the estimated error due to diagonalization
        ! ... (tr2_min), rhog and rhognew are unchanged: rhog contains the
        ! ... input density and rhognew contains the output density, both in
        ! ... G-space.
        ! ... in the other cases rhog now contains mixed charge density in
        ! ... G-space.
        !
        IF ( conv_elec ) THEN
           !
           ! ... if convergence is achieved, rhognew is copied into rhog so
           ! ... that rho and rhog contain the same charge density, one in
           ! ... R-space, the other in G-space
           !
           rhog(:,:) = rhognew(:,:)
           !
        END IF
        !
        DEALLOCATE( rhognew )
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
              ethr = dr2 / nelec
              !
              CYCLE scf_step
              !
           END IF
           !
        END IF
        !
        IF ( .NOT. conv_elec ) THEN
           !
           ! ... bring mixed charge density (rhog) from G- to R-space (rhonew)
           !
           ALLOCATE( rhonew( nrxx, nspin ) )
           !
           DO is = 1, nspin
              !
              psic(:) = ( 0.D0, 0.D0 )
              !
              psic(nl(:)) = rhog(:,is)
              !
              IF ( gamma_only ) psic(nlm(:)) = CONJG( rhog(:,is) )
              !
              CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1 )
              !
              rhonew(:,is) = psic(:)
              !
           END DO
           !
           ! ... no convergence yet: calculate new potential from mixed
           ! ... charge density (i.e. the new estimate) 
           !
           CALL v_of_rho( rhonew, rhog, rho_core, rhog_core, &
                          ehart, etxc, vtxc, etotefield, charge, vr )
           !
           ! ... estimate correction needed to have variational energy:
           ! ... T + E_ion (eband + deband) are calculated in sum_band
           ! ... and delta_e using the output charge density rho;
           ! ... E_H (ehart) and E_xc (etxc) are calculated in v_of_rho
           ! ... above, using the mixed charge density rhonew.
           ! ... delta_escf corrects for this difference at first order
           !
           descf = delta_escf()
           !
           ! ... now rho contains the mixed charge density in R-space
           !
           rho(:,:) = rhonew(:,:)
           !
           DEALLOCATE( rhonew )
           !
           ! ... write the charge density to file
           !
           CALL write_rho( rho, nspin )
           !
        ELSE
           !
           ALLOCATE( vaux( nrxx, nspin ) )
           !
           ! ... convergence reached:
           ! ... 1) the new HXC-potential is temporarily stored in vaux and
           ! ...    then copied into vr.
           ! ... 2) vnew contains V(out)-V(in) ( used to correct the forces ).
           !
           CALL v_of_rho( rho, rhog, rho_core, rhog_core, &
                          ehart, etxc, vtxc, etotefield, charge, vaux )
           !
           vnew(:,:) = vaux(:,:) - vr(:,:)
           !
           vr(:,:) = vaux(:,:)
           !
           DEALLOCATE( vaux )
           !
           ! ... note that rho is the output, not mixed, charge density
           ! ... so correction for variational energy is no longer needed
           !
           descf = 0.D0
           !
        END IF
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
        EXIT scf_step
        !
     END DO scf_step
     !
     ! ... define the total local potential (external + scf)
     !
     CALL set_vrs( vrs, vltot, vr, nrxx, nspin, doublegrid )
     !
     IF ( lda_plus_u ) THEN  
        !
        IF ( ionode ) THEN
           !
           CALL seqopn( iunocc, 'occup', 'FORMATTED', exst )
           !
           WRITE( iunocc, * ) ns
           !
           CLOSE( UNIT = iunocc, STATUS = 'KEEP' )
           !
        END IF
        !
     END IF
     !
     ! ... in the US case we have to recompute the self-consistent
     ! ... term in the nonlocal potential
     !
     CALL newd()
     !
     ! ... save converged wfc if they have not been written previously
     !     
     IF ( nks == 1 .AND. reduce_io ) &
        CALL davcio( evc, nwordwfc, iunwfc, nks, 1 )
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
        IF ( noncolin .AND. domag ) CALL report_mag()
        !
     END IF
     !
     WRITE( stdout, 9000 ) get_clock( 'PWSCF' )
     !
     IF ( conv_elec ) WRITE( stdout, 9101 )
     !
     CALL flush_unit( stdout )
     !
     IF ( conv_elec .OR. MOD( iter, iprint ) == 0 ) THEN
        !
        ALLOCATE ( ngk_g (nkstot) ) 
        !
        ngk_g(1:nks) = ngk(:)
        !
        CALL mp_sum( ngk_g(1:nks), intra_pool_comm )
        !
        CALL ipoolrecover( ngk_g, 1, nkstot, nks )
        !
        CALL mp_bcast( ngk_g, root_pool, intra_pool_comm )
        CALL mp_bcast( ngk_g, root_pool, inter_pool_comm )
        !
        CALL poolrecover( et, nbnd, nkstot, nks )
        !
        DO ik = 1, nkstot
           !
           IF ( lsda ) THEN
              !
              IF ( ik == 1 ) WRITE( stdout, 9015)
              IF ( ik == ( 1 + nkstot / 2 ) ) WRITE( stdout, 9016)
              !
           END IF
           !
           IF ( conv_elec ) THEN
              WRITE( stdout, 9021 ) ( xk(i,ik), i = 1, 3 ), ngk_g(ik)
           ELSE
              WRITE( stdout, 9020 ) ( xk(i,ik), i = 1, 3 )
           END IF
           !
           WRITE( stdout, 9030 ) ( et(ibnd,ik) * rytoev, ibnd = 1, nbnd )
           !
           IF( iverbosity > 0 ) THEN
               !
               ALLOCATE( wg_g( nbnd, nkstot ) )
               !
               wg_g = wg
               CALL poolrecover( wg_g, nbnd, nkstot, nks )
               !
               WRITE( stdout, 9032 )
               WRITE( stdout, 9030 ) ( wg_g(ibnd,ik), ibnd = 1, nbnd )
               !
               DEALLOCATE( wg_g )
               !
           END IF
           !
        END DO
        !
        DEALLOCATE ( ngk_g )
        !
        IF ( lgauss .OR. ltetra ) THEN
           !
           IF ( two_fermi_energies ) THEN
              WRITE( stdout, 9041 ) ef_up*rytoev, ef_dw*rytoev
           ELSE
              WRITE( stdout, 9040 ) ef*rytoev
           END IF
           !
        ELSE
           !
           IF ( tfixed_occ ) THEN
              ibnd = 0
              DO kbnd = 1, nbnd
                 IF ( nspin == 1 .OR. nspin == 4 ) THEN
                    IF ( f_inp(kbnd,1) > 0.D0 ) ibnd = kbnd
                 ELSE
                    IF ( f_inp(kbnd,1) > 0.D0 ) ibnd_up   = kbnd
                    IF ( f_inp(kbnd,2) > 0.D0 ) ibnd_dw = kbnd
                 END IF
              END DO
           ELSE
              IF ( nspin == 1 ) THEN
                 ibnd = NINT( nelec ) / 2
              ELSE
                 ibnd    = NINT( nelec )
                 ibnd_up = NINT( nelup )
                 ibnd_dw = NINT( neldw )
              END IF
           END IF
           !
           IF ( ionode .AND. nbnd > ibnd ) THEN
              !
              IF ( nspin == 1 .OR. nspin == 4 ) THEN
                 ehomo = MAXVAL( et(ibnd,  1:nkstot) )
                 elumo = MINVAL( et(ibnd+1,1:nkstot) )
              ELSE
                 IF ( ibnd_up == 0 ) THEN
                    !
                    ehomo = MAXVAL( et(ibnd_dw,1:nkstot/2) )
                    elumo = MINVAL( et(ibnd_up+1,nkstot/2+1:nkstot) )
                    !
                 ELSE IF ( ibnd_dw == 0 ) THEN
                    !
                    ehomo = MAXVAL( et(ibnd_up,1:nkstot/2) )
                    elumo = MINVAL( et(ibnd_dw+1,nkstot/2+1:nkstot) )
                    !
                 ELSE
                    !
                    ehomo = MAX( MAXVAL( et(ibnd_up,1:nkstot/2) ), &
                                 MAXVAL( et(ibnd_dw,nkstot/2+1:nkstot) ) )
                    elumo = MIN( MINVAL( et(ibnd_up+1,1:nkstot/2) ), &
                                 MINVAL( et(ibnd_dw+1,nkstot/2+1:nkstot) ) )
                    !
                 END IF
              END IF
              !
              WRITE( stdout, 9042 ) ehomo*rytoev, elumo*rytoev
              !
           END IF
        END IF
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
           CALL v_of_rho( rho, rhog, rho_core, rhog_core, &
                          ehart, etxc, vtxc, etotefield, charge, vr )
           !
           CALL set_vrs( vrs, vltot, vr, nrxx, nspin, doublegrid )
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
        !
        !   With Fermi-Dirac population factor, etot is the electronic
        !   free energy F = E - TS , demet is the -TS contribution
        !
        IF ( degauss /= 0.0d0 ) WRITE( stdout, 9070 ) demet
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
        WRITE( stdout, 9110 )
        !
        IF ( output_drho /= ' ' ) CALL remove_atomic_rho()
        !
        CALL stop_clock( 'electrons' )
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
  WRITE( stdout, 9120 )
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
9001 FORMAT(/'     Self-consistent Calculation' )
9010 FORMAT(/'     iteration #',I3,'     ecut=',F9.2,' Ry',5X,'beta=',F4.2 )
9015 FORMAT(/' ------ SPIN UP ------------'/ )
9016 FORMAT(/' ------ SPIN DOWN ----------'/ )
9017 FORMAT(/'     total magnetization       =', F9.2,' Bohr mag/cell', &
            /'     absolute magnetization    =', F9.2,' Bohr mag/cell' )
9018 FORMAT(/'     total magnetization       =',3f9.2,' Bohr mag/cell' &
       &   ,/'     absolute magnetization    =', f9.2,' Bohr mag/cell' )
9020 FORMAT(/'          k =',3F7.4,'     band energies (ev):'/ )
9021 FORMAT(/'          k =',3F7.4,' (',I6,' PWs)   bands (ev):'/ )
9030 FORMAT( '  ',8F9.4 )
9032 FORMAT(/'     occupation numbers ' )
9042 FORMAT(/'     highest occupied, lowest unoccupied level (ev): ',2F10.4 )
9041 FORMAT(/'     the spin up/dw Fermi energies are ',2F10.4,' ev' )
9040 FORMAT(/'     the Fermi energy is ',F10.4,' ev' )
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
9066 FORMAT( '     dexx                      =',F15.8,' Ry' )
9065 FORMAT( '     Hubbard energy            =',F15.8,' Ry' )
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
9110 FORMAT(/'     convergence has been achieved' )
9120 FORMAT(/'     convergence NOT achieved, stopping' )
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     SUBROUTINE non_scf()
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       REAL(DP), EXTERNAL :: efermit, efermig
       !
       !
       WRITE( stdout, 9002 )
       !
       CALL flush_unit( stdout )
       !
       iter = 1
       !
       ! ... diagonalization of the KS hamiltonian
       !
       IF ( lelfield ) THEN
          !
          DO inberry = 1, nberrycyc
             !
             ALLOCATE( psi( npwx, nbnd ) )
             !
             DO ik=1,nks
                !
                CALL davcio( psi, nwordwfc, iunwfc,    ik, -1 )
                CALL davcio( psi, nwordwfc, iunefield, ik,  1 )
                !
             END DO
             !
             DEALLOCATE( psi )
             !
             CALL c_bands( iter, ik_, dr2 )
             !
          END DO
          !
       ELSE
          !
          CALL c_bands( iter, ik_, dr2 )
          !
       END IF
       !
       conv_elec = .TRUE.
       !
       CALL poolrecover( et, nbnd, nkstot, nks )
       !
       WRITE( stdout, 9000 ) get_clock( 'PWSCF' )
       !
       WRITE( stdout, 9102 )
       !
       ! ... write band eigenvalues
       !
       DO ik = 1, nkstot
          !
          IF ( lsda ) THEN
             !   
             IF ( ik == 1 ) WRITE( stdout, 9015 )
             IF ( ik == ( 1 + nkstot / 2 ) ) WRITE( stdout, 9016 )
             !
          END IF
          !
          WRITE( stdout, 9020 ) ( xk(i,ik), i = 1, 3 )
          WRITE( stdout, 9030 ) ( et(ibnd,ik)*rytoev, ibnd = 1, nbnd )
          !
       END DO
       !
       IF ( lgauss .AND. .NOT. lbands ) THEN
          !
          ef = efermig( et, nbnd, nks, nelec, wk, degauss, ngauss, 0, isk )
          !
          WRITE( stdout, 9040 ) ef*rytoev
          !
       ELSE IF ( ltetra .AND. .NOT. lbands ) THEN
          !
          ef = efermit( et, nbnd, nks, nelec, nspin, ntetra, tetra, 0, isk )
          !
          WRITE( stdout, 9040 ) ef*rytoev
          !
       ELSE IF (.NOT. lbands) THEN
          !
          IF ( tfixed_occ ) THEN
              ibnd = 0
              DO kbnd = 1, nbnd
                 IF ( nspin == 1 .OR. nspin == 4 ) THEN
                    IF ( f_inp(kbnd,1) > 0.D0 ) ibnd = kbnd
                 ELSE
                    IF ( f_inp(kbnd,1) > 0.D0 ) ibnd_up   = kbnd
                    IF ( f_inp(kbnd,2) > 0.D0 ) ibnd_dw = kbnd
                 END IF
              END DO
          ELSE
             IF ( nspin == 1 ) THEN
                ibnd = NINT( nelec ) / 2
             ELSE
                ibnd = NINT( nelec )
             END IF
          END IF
          !
          IF ( ionode .AND. nbnd > ibnd ) THEN
             !
             IF ( nspin == 1 .OR. nspin == 4 ) THEN
                ehomo = MAXVAL( et(ibnd  ,1:nkstot) )
                elumo = MINVAL( et(ibnd+1,1:nkstot) )
             !
             ELSE
                ehomo = MAX( MAXVAL( et(ibnd_up,1:nkstot/2) ), &
                             MAXVAL( et(ibnd_dw,nkstot/2+1:nkstot) ) )
                elumo = MIN( MINVAL( et(ibnd_up+1,1:nkstot/2) ), &
                             MINVAL( et(ibnd_dw+1,nkstot/2+1:nkstot) ) )
             END IF
             IF ( ehomo < elumo ) &
                WRITE( stdout, 9042 ) ehomo*rytoev, elumo*rytoev
             !
          END IF
          !
       END IF
       !
       CALL flush_unit( stdout )
       !
       ! ... do a Berry phase polarization calculation if required
       !
       IF ( lberry ) CALL c_phase()
       !
       CALL stop_clock( 'electrons' )
       !
9000 FORMAT(/'     total cpu time spent up to now is ',F9.2,' secs' )
9002 FORMAT(/'     Band Structure Calculation' )
9015 FORMAT(/' ------ SPIN UP ------------'/ )
9016 FORMAT(/' ------ SPIN DOWN ----------'/ )
9020 FORMAT(/'          k =',3F7.4,'     band energies (ev):'/ )
9030 FORMAT( '  ',8F9.4 )
9040 FORMAT(/'     the Fermi energy is ',F10.4,' ev' )
9042 FORMAT(/'     Highest occupied, lowest unoccupied level (ev): ',2F10.4 )
9102 FORMAT(/'     End of band structure calculation' )
       !
     END SUBROUTINE non_scf
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
             mag = rho(ir,1) - rho(ir,2)
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
             mag = SQRT( rho(ir,2)**2 + rho(ir,3)**2 + rho(ir,4)**2 )
             !
             DO i = 1, 3
                !
                magtot_nc(i) = magtot_nc(i) + rho(ir,i+1)
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
       ENDIF
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
       !
       IMPLICIT NONE
       !
       LOGICAL :: check_stop_now
       INTEGER :: unit
       !
       !
       IF ( lpath ) THEN  
          !
          unit = iunpath
          !  
       ELSE
          !
          unit = stdout
          !   
       END IF
       !
       check_stop_now = global_check_stop_now( unit )
       !
       IF ( check_stop_now ) THEN
          !  
          conv_elec = .FALSE.
          !
          RETURN          
          !
       END IF              
       !
     END FUNCTION check_stop_now
     !
     !-----------------------------------------------------------------------
     FUNCTION delta_e()
       !-----------------------------------------------------------------------
       !
       ! ... delta_e = - \int rho(r) V_scf(r)
       !
       USE kinds, ONLY : DP
       !
       IMPLICIT NONE
       !   
       REAL(DP) :: delta_e
       !
       INTEGER :: ipol
       !
       !
       delta_e = 0.D0
       !
       DO ipol = 1, nspin
          !
          delta_e = delta_e - SUM( rho(:,ipol)*vr(:,ipol) )
          !
       END DO
       !
       delta_e = omega * delta_e / ( nr1*nr2*nr3 )
       !
       CALL mp_sum( delta_e, intra_pool_comm )
       !
       RETURN
       !
     END FUNCTION delta_e
     !
     !-----------------------------------------------------------------------
     FUNCTION delta_escf()
       !-----------------------------------------------------------------------
       !
       ! ... delta_escf = - \int \delta rho(r) V_scf(r)
       ! ... calculates the difference between the Hartree and XC energy
       ! ... at first order in the charge density difference \delta rho(r) 
       !
       USE kinds, ONLY : DP
       !
       IMPLICIT NONE
       !   
       REAL(DP) :: delta_escf
       !
       INTEGER :: ipol
       !
       !
       delta_escf = 0.D0
       !
       DO ipol = 1, nspin
          !
          delta_escf = delta_escf - &
                       SUM( ( rhonew(:,ipol) - rho(:,ipol) )*vr(:,ipol) )
          !
       END DO
       !
       delta_escf = omega * delta_escf / ( nr1*nr2*nr3 )
       !
       CALL mp_sum( delta_escf, intra_pool_comm )
       !
       RETURN
       !
     END FUNCTION delta_escf
     !
END SUBROUTINE electrons
