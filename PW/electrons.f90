!
! Copyright (C) 2001-2004 PWSCF group
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
  ! ... Hamiltonian, the routine sum_band to compute the charge
  ! ... density, the routine v_of_rho to compute the new potential
  ! ... and the routine mix_rho to mix input and output charge densities
  ! ... It prints on output the total energy and its decomposition in
  ! ... the separate contributions.
  !
  USE kinds,                ONLY : DP
  USE parameters,           ONLY : npk 
  USE constants,            ONLY : eps8, rytoev
  USE io_global,            ONLY : stdout, ionode
  USE cell_base,            ONLY : at, bg, alat, omega, tpiba2
  USE ions_base,            ONLY : zv, nat, ntyp => nsp, ityp, tau
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
  USE wvfct,                ONLY : nbnd, et, gamma_only, wg,npwx
  USE fixed_occ,            ONLY : f_inp, tfixed_occ
  USE ener,                 ONLY : etot, eband, deband, ehart, vtxc, etxc, &
                                   etxcc, ewld, demet, ef, ef_up, ef_dw 
  USE scf,                  ONLY : rho, vr, vltot, vrs, rho_core
  USE control_flags,        ONLY : mixing_beta, tr2, ethr, niter, nmix,       &
                                   iprint, istep, lscf, lmd, conv_elec,       &
                                   lbands, restart, reduce_io, iverbosity
  USE io_files,             ONLY : prefix, iunwfc, iunocc, nwordwfc, iunpath, &
                                   output_drho, iunefield
  USE ldaU,                 ONLY : ns, nsnew, eth, Hubbard_U, &
                                   niter_with_fixed_ns, Hubbard_lmax, &
                                   lda_plus_u  
  USE extfield,             ONLY : tefield, etotefield  
  USE wavefunctions_module, ONLY : evc, evc_nc, psic
  USE noncollin_module,     ONLY : noncolin, npol, magtot_nc
  USE noncollin_module,     ONLY : factlist, pointlist, pointnum, mcons,&
                                   i_cons, bfield, lambda, vtcon, report
  USE spin_orb,             ONLY : domag
  USE bp,                   ONLY : lelfield, lberry, nberrycyc
#if defined (EXX)
  USE exx,                  ONLY : exxinit, init_h_wfc, exxenergy, exxenergy2 
  USE funct,                ONLY : dft_is_hybrid, exx_is_active
#endif
  !
  IMPLICIT NONE
  !
  ! ... a few local variables
  !  
#if defined (EXX)
  REAL(DP) :: dexx
  REAL(DP) :: fock0,  fock1,  fock2
#endif
  INTEGER :: &
      ngkp(npk)        !  number of plane waves summed on all nodes
  CHARACTER (LEN=256) :: &
      flmix            !
  REAL(DP) :: &
      dr2,            &!  the norm of the diffence between potential
      charge,         &!  the total charge
      mag,            &!  local magnetization
      ehomo, elumo,   &!  highest occupied and lowest onuccupied levels
      tcpu             !  cpu time
   INTEGER :: &
      i,              &!  counter on polarization
      is,             &!  counter on spins
      ig,             &!  counter on G-vectors
      ik,             &!  counter on k points
      kbnd,           &!  counter on bands
      ibnd_up,ibnd_down, &!  counter on bands
      ibnd,           &!  counter on bands
      idum,           &!  dummy counter on iterations
      iter,           &!  counter on iterations
      ik_              !  used to read ik from restart file
  INTEGER :: &
      ldim2           !
  REAL(DP) :: &
       tr2_min,      &! estimated error on energy coming from diagonalization
       descf          ! correction for variational energy

  REAL(DP), ALLOCATABLE :: &
      wg_g(:,:)        ! temporary array used to collect array wg from pools
                       ! and then print occupations on stdout
  LOGICAL :: &
      exst, first
  !
  ! ... external functions
  !
  REAL(DP), EXTERNAL :: ewald, get_clock
  !
  ! ... auxiliary variables for calculating and storing rho in G-space
  !
  COMPLEX(DP), ALLOCATABLE :: rhog(:,:)
  COMPLEX(DP), ALLOCATABLE :: rhognew(:,:)
  REAL(DP),    ALLOCATABLE :: rhonew(:,:)
  !
  ! ... variables needed for electric field calculation
  !
  COMPLEX(DP), ALLOCATABLE  :: psi(:,:)
  INTEGER                   :: inberry
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
  tcpu = get_clock( 'PWSCF' )
  WRITE( stdout, 9000 ) tcpu
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
  ewld = ewald( alat, nat, ntyp, ityp, zv, at, bg, tau, omega, &
                g, gg, ngm, gcutm, gstart, gamma_only, strf )
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
  ! ... for the first scf iteration of each ionic step after the first,
  ! ... the threshold is fixed to a default value of 1.D-5
  !
#if defined (EXX)
10 continue
#endif

  IF ( istep > 1 ) ethr = 1.D-5
  !
  WRITE( stdout, 9001 )
  !
  CALL flush_unit( stdout )
  !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%          iterate !          %%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
  ! ... bring starting rho to G-space
  !
  IF ( .not. ALLOCATED(rhog) ) ALLOCATE (rhog(ngm, nspin))
  do is = 1, nspin
     psic(:) = rho (:, is)
     call cft3 (psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1)
     rhog(:, is) = psic ( nl(:) )
  end do
  !
  DO idum = 1, niter
     !
     IF ( check_stop_now() ) RETURN
     !
     !  
     iter = iter + 1
     !
     WRITE( stdout, 9010 ) iter, ecutwfc, mixing_beta
     !
     CALL flush_unit( stdout )
     !
     ! ... Convergence threshold for iterative diagonalization
     ! ... is automatically updated during self consistency
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
        IF ( first ) THEN
           !
           tr2_min = nelec * ethr
           !
        ELSE
           !
           tr2_min = 0.D0
           !
        END IF
        !
        ! ... diagonalization of the KS hamiltonian
        !
        if(lelfield) then
          do inberry=1,nberrycyc
            ALLOCATE (psi(npwx,nbnd))
            do ik=1,nks
               call davcio(psi,nwordwfc,iunwfc,ik,-1)
               call davcio(psi,nwordwfc,iunefield,ik,1)
             enddo
             DEALLOCATE( psi)
            CALL c_bands( iter, ik_, dr2 )
          enddo
        else
          CALL c_bands( iter, ik_, dr2 )
         endif
        !
        IF ( check_stop_now() ) RETURN
        !
        CALL sum_band()
        !
        IF ( lda_plus_u )  THEN
           !
           ldim2 = ( 2 * Hubbard_lmax + 1 )**2
           !
           CALL write_ns()
           !
           IF ( first .AND. istep == 1 .AND. &
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
        ! ... delta_e = - int rho(r) (V_H + V_xc)(r) dr
        !
        deband = delta_e()
        !
        ! ... bring newly calculated (in sum_band) rho to G-space for mixing
        !
        ALLOCATE (rhognew(ngm, nspin))
        do is = 1, nspin
           psic(:) = rho (:, is)
           call cft3 (psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)
           rhognew (:, is) = psic ( nl(:) )
        end do
        !
        CALL mix_rho( rhognew, rhog, nsnew, ns, mixing_beta, &
             dr2, tr2_min, iter, nmix, flmix, conv_elec )
        !
        DEALLOCATE (rhognew)
        !
        ! ... for the first scf iteration it is controlled that the 
        ! ... threshold is small enough for the diagonalization to 
        ! ... be adequate
        !
        IF ( first ) THEN
           !
           first = .FALSE.
           !
           IF ( dr2 < tr2_min ) THEN
              !
              ! ... a new diagonalization is needed       
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
           ! ... bring mixed rho from G- to R-space
           !
           ALLOCATE (rhonew (nrxx, nspin) )
           do is = 1, nspin
              psic( :) = (0.d0, 0.d0)
              psic( nl(:) ) = rhog (:, is)
              if (gamma_only) psic( nlm(:) ) = CONJG (rhog (:, is))
              call cft3 (psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, +1)
              rhonew (:, is) = psic (:)
           end do
           !
           ! ... no convergence yet: calculate new potential from 
           ! ... new estimate of the charge density 
           !
           CALL v_of_rho( rhonew, rho_core, nr1, nr2, nr3, nrx1, nrx2,   &
                          nrx3, nrxx, nl, ngm, gstart, nspin, g, gg, alat, &
                          omega, ehart, etxc, vtxc, etotefield, charge, vr )
           !
           ! ... estimate correction needed to have variational energy 
           !
           descf = delta_escf()
           !
           ! ... write the charge density to file
           !
           CALL io_pot( 1, 'rho', rhonew, nspin )
           !
           DEALLOCATE (rhonew )
           !
        ELSE
           !
           ! ... convergence reached: store V(out)-V(in) in vnew
           ! ... Used to correct the forces
           !
           CALL v_of_rho( rho, rho_core, nr1, nr2, nr3, nrx1, nrx2, nrx3,   &
                          nrxx, nl, ngm, gstart, nspin, g, gg, alat, omega, &
                          ehart, etxc, vtxc, etotefield, charge, vnew )
           !
           vnew = vnew - vr
           !
           ! ... correction for variational energy no longer needed
           !
           descf = 0.D0
           !
           DEALLOCATE (rhog)
           !
        END IF
#if defined (EXX)
        if (exx_is_active()) then
           fock1 = exxenergy2()
           fock2 = fock0
        else
           fock0 = 0.d0
        end if
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
     ! ... In the US case we need to recompute the self consistent term
     ! ... in the nonlocal potential.
     !
     CALL newd()
     !
     ! ... write the potential to file
     !
     CALL io_pot( 1, 'pot', vr, nspin )     
     !
     ! ... save converged wfc if they have not been written previously
     !     
     IF ( noncolin ) THEN
        !
        IF ( nks == 1 .AND. reduce_io ) &
           CALL davcio( evc_nc, nwordwfc, iunwfc, nks, 1 )
        !
     ELSE
        !
        IF ( nks == 1 .AND. reduce_io ) &
           CALL davcio( evc, nwordwfc, iunwfc, nks, 1 )
        !
     END IF
     !
     ! ... calculate the polarization
     !
     IF ( lelfield ) CALL c_phase_field( )
     !
     ! ... write recover file
     !
     CALL save_in_electrons( iter, dr2 )
     !
     IF ( ( MOD(iter,report) == 0 ).OR. &
          ( report /= 0 .AND. conv_elec ) ) THEN
        !
        IF ( noncolin .and. domag ) CALL report_mag()
        !
     END IF
     !
     tcpu = get_clock( 'PWSCF' )
     WRITE( stdout, 9000 ) tcpu
     !
     IF ( conv_elec ) WRITE( stdout, 9101 )
     !
     CALL flush_unit( stdout )
     !
     IF ( conv_elec .OR. MOD( iter, iprint ) == 0 ) THEN
        !
#if defined (__PARA)
        !
        ngkp(1:nks) = ngk(1:nks)
        !
        CALL ireduce( nks, ngkp )
        CALL ipoolrecover( ngkp, 1, nkstot, nks )
        CALL poolrecover( et, nbnd, nkstot, nks )
        !
#endif
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
#if defined (__PARA)
              WRITE( stdout, 9021 ) ( xk(i,ik), i = 1, 3 ), ngkp(ik)
#else
              WRITE( stdout, 9021 ) ( xk(i,ik), i = 1, 3 ), ngk(ik)
#endif
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
        IF ( lgauss .OR. ltetra ) THEN
           IF ( two_fermi_energies ) THEN
              WRITE( stdout, 9041 ) ef_up * rytoev, ef_dw * rytoev
           ELSE
              WRITE( stdout, 9040 ) ef * rytoev
           END IF
        ELSE
           !
           IF ( tfixed_occ ) THEN
              ibnd = 0
              DO kbnd=1,nbnd
                 IF (nspin==1.OR.nspin==4) THEN
                    IF (f_inp(kbnd,1)>0.d0) ibnd=kbnd
                 ELSE
                    IF (f_inp(kbnd,1)>0.d0) ibnd_up=kbnd
                    IF (f_inp(kbnd,2)>0.d0) ibnd_down=kbnd
                 END IF
              END DO
           ELSE
              IF ( nspin == 1 ) THEN
                 ibnd =  NINT( nelec ) / 2.D0
              ELSE
                 ibnd =  NINT( nelec )
              END IF
           END IF
           !
           IF ( ionode .AND. nbnd > ibnd ) THEN
              !
              IF (nspin==1.OR.nspin==4) THEN
                 ehomo = MAXVAL ( et( ibnd  , 1:nkstot) )
                 elumo = MINVAL ( et( ibnd+1, 1:nkstot) )
              ELSE 
                 ehomo = MAX(MAXVAL(et(ibnd_up, 1:nkstot/2)), &
                             MAXVAL(et(ibnd_down, nkstot/2+1:nkstot)))
                 elumo = MIN(MINVAL(et(ibnd_up+1, 1:nkstot/2) ), &
                             MINVAL(et(ibnd_down+1, nkstot/2+1:nkstot)))
              ENDIF
              !
              WRITE( stdout, 9042 ) ehomo * rytoev, elumo * rytoev
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
     etot = etot - 0.5d0 * fock0
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
           CALL v_of_rho( rho, rho_core, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
                          nrxx, nl, ngm, gstart, nspin, g, gg, alat, omega, &
                          ehart, etxc, vtxc, etotefield, charge, vr )
           !
           CALL set_vrs( vrs, vltot, vr, nrxx, nspin, doublegrid )
           !
           WRITE (stdout,*) " NOW GO BACK TO REFINE HYBRID CALCULATION"
           WRITE (stdout,*) fock0
           !
           iter = 0
           !
           GO TO 10
           !
        END IF
        !
        fock2 = exxenergy2()
        !
        dexx = fock1 - 0.5D0 * ( fock0 + fock2 )
        !
        etot = etot  - dexx
        !
        WRITE( stdout,*) fock0, fock1, fock2
        WRITE( stdout, 9066 ) dexx
        !
        fock0 = fock2
        !
     END IF
     !
#endif
     !
     IF ( lda_plus_u ) etot = etot + eth
     IF ( tefield )    etot = etot + etotefield
     !
     IF ( ( conv_elec .OR. MOD( iter, iprint ) == 0 ) .AND. &
          ( .NOT. lmd ) ) THEN
        !  
        IF ( dr2 > eps8 ) THEN
           WRITE( stdout, 9081 ) etot, dr2
        ELSE
           WRITE( stdout, 9083 ) etot, dr2
        END IF
        !
        WRITE( stdout, 9060 ) &
            eband, ( eband + deband ), ehart, ( etxc - etxcc ), ewld
        !
#if defined (EXX)
        !
        WRITE( stdout, 9062 ) fock1
        WRITE( stdout, 9063 ) fock2
        WRITE( stdout, 9064 ) 0.5D0 * fock2
        !
#endif
        !
        IF ( tefield ) WRITE( stdout, 9061 ) etotefield
        IF ( lda_plus_u ) WRITE( stdout, 9065 ) eth
        IF ( degauss /= 0.0 ) WRITE( stdout, 9070 ) demet
        !
     ELSE IF ( conv_elec .AND. lmd ) THEN
        !
        IF ( dr2 > eps8 ) THEN
           WRITE( stdout, 9081 ) etot, dr2
        ELSE
           WRITE( stdout, 9083 ) etot, dr2
        END IF
        !
     ELSE
        !
        IF ( dr2 > eps8 ) THEN
           WRITE( stdout, 9080 ) etot, dr2
        ELSE
           WRITE( stdout, 9082 ) etot, dr2
        END IF
        !
     END IF
     !
     IF ( lsda ) WRITE( stdout, 9017 ) magtot, absmag
     !
     IF ( noncolin .AND. domag ) &
        WRITE( stdout, 9018 ) ( magtot_nc(i), i = 1, 3 ), absmag
     !
     IF ( i_cons == 3 .OR. i_cons == 4 )  &
        WRITE( stdout, 9071 ) bfield(1), bfield(2),bfield(3)
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
9010 FORMAT(/'     iteration #',I3,'     ecut=',F9.2,' ryd',5X,'beta=',F4.2 )
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
9060 FORMAT(/'     band energy sum           =',  F15.8,' ryd' &
            /'     one-electron contribution =',  F15.8,' ryd' &
            /'     hartree contribution      =',  F15.8,' ryd' &
            /'     xc contribution           =',  F15.8,' ryd' &
            /'     ewald contribution        =',  F15.8,' ryd' )
9061 FORMAT( '     electric field correction =',  F15.8,' ryd' )
9062 FORMAT( '     Fock energy 1             =',  F15.8,' ryd' )
9063 FORMAT( '     Fock energy 2             =',  F15.8,' ryd' )
9064 FORMAT( '     Half Fock energy 2        =',  F15.8,' ryd' )
9066 FORMAT( '     dexx                      =',  F15.8,' ryd' )
9065 FORMAT( '     Hubbard energy            =',F15.8,' ryd' )
9070 FORMAT( '     correction for metals     =',F15.8,' ryd' )
9071 FORMAT( '     Magnetic field            =',3F12.7,' ryd' )
9072 FORMAT( '     Magnetic field            =', F12.7,' ryd' )
9073 FORMAT( '     lambda                    =', F11.2,' ryd' )
9080 FORMAT(/'     total energy              =',0PF15.8,' ryd' &
            /'     estimated scf accuracy    <',0PF15.8,' ryd' )
9081 FORMAT(/'!    total energy              =',0PF15.8,' ryd' &
            /'     estimated scf accuracy    <',0PF15.8,' ryd' )
9082 FORMAT(/'     total energy              =',0PF15.8,' ryd' &
            /'     estimated scf accuracy    <',1PE15.1,' ryd' )
9083 FORMAT(/'!    total energy              =',0PF15.8,' ryd' &
            /'     estimated scf accuracy    <',1PE15.1,' ryd' )
9085 FORMAT(/'     total energy              =',0PF15.8,' ryd' &
            /'     potential mean squ. error =',1PE15.1,' ryd^2' )
9086 FORMAT(/'!    total energy              =',0PF15.8,' ryd' &
            /'     potential mean squ. error =',1PE15.1,' ryd^2' )
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
       tcpu = get_clock( 'PWSCF' )
       WRITE( stdout, 9000 ) tcpu
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
          WRITE( stdout, 9030 ) ( et(ibnd,ik) * rytoev, ibnd = 1, nbnd )
          !
       END DO
       !
       IF ( lgauss .AND. .NOT. lbands ) THEN
          !
          ef = efermig( et, nbnd, nks, nelec, wk, degauss, ngauss, 0, isk )
          !
          WRITE( stdout, 9040 ) ef * rytoev
          !
       ELSE IF ( ltetra .AND. .NOT. lbands ) THEN
          !
          ef = efermit( et, nbnd, nks, nelec, nspin, ntetra, tetra, 0, isk )
          !
          WRITE( stdout, 9040 ) ef * rytoev
          !
       ELSE 
          !
          IF ( tfixed_occ ) THEN
              ibnd = 0
              DO kbnd=1,nbnd
                 IF (nspin==1.OR.nspin==4) THEN
                    IF (f_inp(kbnd,1)>0.d0) ibnd=kbnd
                 ELSE
                    IF (f_inp(kbnd,1)>0.d0) ibnd_up=kbnd
                    IF (f_inp(kbnd,2)>0.d0) ibnd_down=kbnd
                 END IF
              END DO
          ELSE
             IF ( nspin == 1 ) THEN
                ibnd =  nint (nelec) / 2.d0
             ELSE
                ibnd =  nint (nelec)
             END IF
          END IF
          !
          IF ( ionode .AND. nbnd > ibnd ) THEN
             !
             IF (nspin==1.OR.nspin==4) THEN
                ehomo = MAXVAL ( et( ibnd  , 1:nkstot) )
                elumo = MINVAL ( et( ibnd+1, 1:nkstot) )
             !
             ELSE
                ehomo = MAX(MAXVAL(et(ibnd_up, 1:nkstot/2)), &
                            MAXVAL(et(ibnd_down, nkstot/2+1:nkstot)))
                elumo = MIN(MINVAL(et(ibnd_up+1, 1:nkstot/2) ), &
                            MINVAL(et(ibnd_down+1, nkstot/2+1:nkstot)))
             END IF
             IF ( ehomo < elumo ) &
                  WRITE( stdout, 9042 ) ehomo * rytoev, elumo * rytoev
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
       IF ( output_drho /= ' ' ) CALL remove_atomic_rho()
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
          magtot = magtot * omega / ( nr1 * nr2 * nr3 )
          absmag = absmag * omega / ( nr1 * nr2 * nr3 )
          !
          CALL reduce( 1, magtot )
          CALL reduce( 1, absmag )
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
          CALL reduce( 3, magtot_nc )
          CALL reduce( 1, absmag )
          !
          DO i = 1, 3
             !
             magtot_nc(i) = magtot_nc(i) * omega / ( nr1 * nr2 * nr3 )
             !
          END DO
          !
          absmag = absmag * omega / ( nr1 * nr2 * nr3 )
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
     FUNCTION delta_e ( )
       !-----------------------------------------------------------------------
       !
       ! ... delta_e = - \int rho(r) V_scf(r)
       !
       USE kinds
       !
       IMPLICIT NONE
       !   
       REAL (DP) :: delta_e
       !
       INTEGER :: ipol
       !
       !
       delta_e = 0.D0
       !
       DO ipol = 1, nspin
          !
          delta_e = delta_e - SUM( rho(:,ipol) * vr(:,ipol) )
          !
       END DO
       !
       delta_e = omega * delta_e / ( nr1 * nr2 * nr3 )
       !
       CALL reduce( 1, delta_e )
       !
       RETURN
       !
     END FUNCTION delta_e
     !
     !-----------------------------------------------------------------------
     FUNCTION delta_escf ( )
       !-----------------------------------------------------------------------
       !
       ! ... delta_escf = - \int \delta rho(r) V_scf(r)
       ! ... this is the correction needed to have variational energy
       !
       USE kinds
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
                       SUM( ( rhonew(:,ipol) - rho(:,ipol) ) * vr(:,ipol) )
          !
       END DO
       !
       delta_escf = omega * delta_escf / ( nr1 * nr2 * nr3 )
       !
       CALL reduce( 1, delta_escf )
       !
       RETURN
       !
     END FUNCTION delta_escf
     !
END SUBROUTINE electrons
