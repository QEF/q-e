!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!  
#include "f_defs.h"
!
#if defined (__PARA)
#define NRXX ncplane*npp(me_pool+1)
  ! This is needed in mix_pot whenever nproc is not a divisor of nr3.
#else
#define NRXX nrxx
#endif
!
!----------------------------------------------------------------------------
SUBROUTINE electrons()
  !----------------------------------------------------------------------------
  !
  ! ... This routine is a driver of the self-consistent cycle.
  ! ... It uses the routine c_bands for computing the bands at fixed
  ! ... Hamiltonian, the routine sum_bands to compute the charge
  ! ... density, the routine v_of_rho to compute the new potential
  ! ... and the routine mix_potential to mix input and output
  ! ... potentials.
  !
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
                                   nrx3, nrxx, nl, g, gg, ecutwfc, gcutm
  USE gsmooth,              ONLY : doublegrid, ngms
  USE klist,                ONLY : xk, wk, degauss, nelec, ngk, nks, nkstot, &
                                   lgauss, ngauss, two_fermi_energies
  USE lsda_mod,             ONLY : lsda, nspin, magtot, absmag, isk
  USE ktetra,               ONLY : ltetra, ntetra, tetra  
  USE vlocal,               ONLY : strf, vnew  
  USE wvfct,                ONLY : nbnd, et, gamma_only, wg  
  USE ener,                 ONLY : etot, eband, deband, ehart, vtxc, etxc, &
                                   etxcc, ewld, demet, ef, ef_up, ef_dw 
  USE scf,                  ONLY : rho, rho_save, vr, vltot, vrs, rho_core
  USE control_flags,        ONLY : mixing_beta, tr2, ethr, ngm0, &
                                   niter, nmix, imix, iprint, istep, &
                                   lscf, lpath, lmd, conv_elec, restart, &
                                   reduce_io, iverbosity
  USE io_files,             ONLY : prefix, iunwfc, iunocc, nwordwfc, iunpath, &
                                   output_drho
  USE ldaU,                 ONLY : ns, nsnew, eth, Hubbard_U, &
                                   niter_with_fixed_ns, Hubbard_lmax, &
                                   lda_plus_u  
  USE extfield,             ONLY : tefield, etotefield  
  USE bp,                   ONLY : lberry  
  USE wavefunctions_module, ONLY : evc, evc_nc
  USE noncollin_module,     ONLY : noncolin, npol, magtot_nc
  USE noncollin_module,     ONLY : factlist, pointlist, pointnum, mcons,&
                                   i_cons, bfield, lambda, vtcon, report
  USE spin_orb,             ONLY : domag
  USE mp_global,            ONLY : me_pool
  USE pfft,                 ONLY : npp, ncplane
  !
  IMPLICIT NONE
  !
  ! ... a few local variables
  !  
  INTEGER :: &
      ngkp(npk)        !  number of plane waves summed on all nodes
  CHARACTER (LEN=42) :: &
      flmix            !
  REAL(KIND=DP) :: &
      dr2,            &!  the norm of the diffence between potential
      charge,         &!  the total charge
      mag,            &!  local magnetization
      tcpu             !  cpu time
   INTEGER :: &
      i,              &!  counter on polarization
      ir,             &!  counter on the mesh points
      ig,             &!
      ik,             &!  counter on k points
      ibnd,           &!  counter on bands
      idum,           &!  dummy counter on iterations
      iter,           &!  counter on iterations
      ik_              !  used to read ik from restart file
  INTEGER :: &
      ldim2           !
  REAL (KIND=DP) :: &
      ehart_new,      &!
      etxc_new,       &!
      vtxc_new,       &!
      etotefield_new, &!
      charge_new      !
  REAL (KIND=DP) :: &
      ethr_min         ! minimal threshold for diagonalization at the first scf
                       ! iteration 
  REAL (KIND=DP), ALLOCATABLE :: &
      wg_g(:,:)        ! temporary array used to recover from pools array wg,
                       ! and then print occupations on stdout
  LOGICAL :: &
      exst
  REAL (KIND=DP) :: dr2old, lambda0
  !
  ! ... external functions
  !
  REAL (KIND=DP), EXTERNAL :: ewald, get_clock
  !
  !
  CALL start_clock( 'electrons' )
  !
  iter = 0
  ik_  = 0
  if (i_cons==3) then
     dr2old=0.d0
     lambda0=lambda
     lambda=0.003
  endif
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
#if defined (FLUSH)
  CALL flush( stdout )
#endif
  !
  IF ( .NOT. lscf ) THEN
     !
     WRITE( stdout, 9002 )
     !
#if defined (FLUSH)
     CALL flush( stdout )
#endif
     !
     IF ( imix >= 0 ) rho_save = rho
     !  
     iter = iter + 1
     !
     ! ... diagonalization of the KS hamiltonian
     !
     CALL c_bands( iter, ik_, dr2 )
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
     IF ( lgauss ) THEN
        !
        call efermig (et, nbnd, nks, nelec, wk, degauss, ngauss, ef, 0, isk)
        WRITE( stdout, 9040 ) ef * rytoev
        !
     ELSE IF ( ltetra ) THEN
        !
        CALL efermit (et, nbnd, nks, nelec, nspin, ntetra, tetra, ef, 0, isk)
        WRITE( stdout, 9040 ) ef * rytoev
        !
     END IF
     !
#if defined (FLUSH)
     CALL flush( stdout )
#endif
     !
     ! ... do a Berry phase polarization calculation if required
     !
     IF ( lberry ) CALL c_phase()
     !
     IF ( output_drho /= ' ' ) CALL remove_atomic_rho()
     !
     CALL stop_clock( 'electrons' )
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
     flmix = ' '
  ELSE
     flmix = 'flmix'
  END IF
  !
  ! ... Convergence threshold for iterative diagonalization
  !
  ! ... for the first scf iteration of each ionic step (except than for the
  ! ... first) the threshold is fixed to a default value of 1.D-5
  !
  IF ( istep > 1 ) ethr = 1.D-5
  !
  IF ( imix >= 0 ) ngm0 = ngm
  !
  WRITE( stdout, 9001 )
  !
#if defined (FLUSH)
  CALL flush( stdout )
#endif
  !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%          iterate !          %%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
  DO idum = 1, niter
     !
     IF ( imix >= 0 ) rho_save = rho
     !  
     iter = iter + 1
     !
     WRITE( stdout, 9010 ) iter, ecutwfc, mixing_beta
     !
     IF ( lpath ) THEN
        !
        CALL flush( stdout )
        !
     ELSE
        !
#if defined (FLUSH)
        CALL flush( stdout )
#endif
        !
     END IF
     !
     ! ... Convergence threshold for iterative diagonalization
     ! ... is automatically updated during self consistency
     !
     IF ( iter > 1 .AND. ik_ == 0 ) THEN
        !
        IF ( imix >= 0 ) THEN
           !
           IF ( iter == 2 ) ethr = 1.D-2
           !
           ethr = MAX( MIN( ethr , ( dr2 / nelec * 0.1D0 ) ) , &
                       ( tr2 / nelec * 0.01D0 ) )
           !
        ELSE
           !
           ethr = MAX( MIN( ( ethr * 0.5D0 ) , &
                            ( SQRT( dr2 ) * 0.001D0 ) ) , 1.D-12 )
           !
        END IF
        !
     END IF
     !
     ! ... diagonalization of the KS hamiltonian
     !
     CALL c_bands( iter, ik_, dr2 )
     !
     ! ... the program checks if the maximum CPU time has been exceeded
     ! ... or if the user has required a soft exit
     !
     IF ( check_stop_now() ) RETURN
     !
     CALL sum_band()
     !
     IF ( lda_plus_u ) CALL write_ns()
     !
     IF ( iter == 1 .AND. lda_plus_u .AND. &
          startingpot == 'atomic' .AND. istep == 1 ) CALL ns_adj()
     !
     ! ... calculate total and absolute magnetization
     !
     IF ( lsda .OR. noncolin ) CALL compute_magnetization()
     !
     CALL v_of_rho( rho, rho_core, nr1, nr2, nr3, nrx1, nrx2, nrx3,   &
                    nrxx, nl, ngm, gstart, nspin, g, gg, alat, omega, &
                    ehart, etxc, vtxc, etotefield, charge, vnew )
     !
     deband = delta_e ( )
     !
     IF ( imix >= 0 ) THEN
        !
        IF ( lda_plus_u .AND. iter <= niter_with_fixed_ns ) THEN
           !
           ldim2 = ( 2 * Hubbard_lmax + 1 )**2
           nsnew = ns 
           !
        END IF
        !
        IF ( iter == 1 ) THEN
           !
           ! ... for the first scf iteration ethr_min is set for a check 
           ! ... in mix_rho ( in mix_rho ethr_min = dr2 * ethr_min )
           !
           ethr_min = 1.D0 / nelec
           !
        ELSE
           !
           ! ... otherwise ethr_min is set to a negative number: 
           ! ... no check is needed
           !
           ethr_min = - 1.D0
           !
        END IF      
        !
        IF ( noncolin ) THEN
           !
           CALL mix_rho_nc( rho, rho_save, nsnew, ns, mixing_beta, dr2, iter, &
                            nmix, flmix, conv_elec )
           !
        ELSE
           !
           CALL mix_rho( rho, rho_save, nsnew, ns, mixing_beta, dr2, ethr, &
                         ethr_min, iter, nmix, flmix, conv_elec )
           !
        END IF
        !
        ! ... for the first scf iteration it is controlled that the threshold 
        ! ... is small enough for the diagonalization to be adequate
        !
        IF ( iter == 1 .AND. ethr >= ethr_min ) THEN
           !
           ! ... a new diagonalization is needed       
           !
           WRITE( stdout, '(/,5X,"Threshold (ethr) on eigenvalues was ", &
                            &    "too large:",/,                         &
                            & 5X,"Diagonalizing with lowered threshold",/)' )
           !
           CALL c_bands( iter, ik_, dr2 )
           !
           ! ... check if the maximum CPU time has been exceeded
           ! ... or if the user has required a soft exit
           !
           IF ( check_stop_now() ) RETURN
           !
           CALL sum_band()
           !
           IF ( lda_plus_u ) CALL write_ns()
           !
           ! ... calculate total and absolute magnetization
           !
           IF ( lsda .OR. noncolin ) CALL compute_magnetization()
           !
           CALL v_of_rho( rho, rho_core, nr1, nr2, nr3, nrx1, nrx2, nrx3,   &
                          nrxx, nl, ngm, gstart, nspin, g, gg, alat, omega, &
                          ehart, etxc, vtxc, etotefield, charge, vnew )
           !
           deband = delta_e ( )
           !
           IF ( lda_plus_u .AND. iter <= niter_with_fixed_ns ) THEN
              !
              ldim2 = ( 2 * Hubbard_lmax + 1 )**2
              nsnew = ns             
              !
           END IF
           !
           ! ... ethr_min is set to a negative number: no check is needed
           !
           IF (  noncolin ) THEN
              !
              CALL mix_rho_nc( rho, rho_save, nsnew, ns, mixing_beta, dr2, &
                               iter, nmix, flmix, conv_elec )
              !
           ELSE
              !
              CALL mix_rho( rho, rho_save, nsnew, ns, mixing_beta, dr2, ethr, &
                           -1.D0, iter, nmix, flmix, conv_elec )
              !
           END IF
           !
        END IF             
        !
        vnew =  vnew - vr
        !
        CALL v_of_rho( rho_save, rho_core, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
                       nrxx, nl, ngm, gstart, nspin, g, gg, alat, omega,    &
                       ehart_new, etxc_new, vtxc_new, etotefield_new,       &
                       charge_new, vr )
        !
     ELSE 
        !
        ! ... old style potential mixing
        !
        CALL vpack( NRXX, nrxx, nspin, vnew, vr, + 1 )
        !
        CALL mix_potential( ( nspin * NRXX ), vnew, vr, mixing_beta, dr2, &
                            tr2, iter, nmix, flmix, conv_elec )
        !
        CALL vpack( NRXX, nrxx, nspin, vnew, vr, -1 )
        !
        if (i_cons.eq.3) then
           if (dr2*5.d0.lt.dr2old.or.dr2<tr2*1.d2) then
              lambda=min(lambda*2.d0,lambda0)
           endif
           dr2old=dr2
        endif
     END IF
     !
     ! ... On output vnew contains V(out)-V(in). Used to correct the forces
     ! ... define the total local potential (external + scf)
     !
     CALL set_vrs( vrs, vltot, vr, nrxx, nspin, doublegrid )
     !
     IF ( lda_plus_u ) THEN  
        !
        ldim2 = ( 2 * Hubbard_lmax + 1 )**2
        !
        IF ( iter > niter_with_fixed_ns .AND. imix < 0 ) ns = nsnew
        !
        IF ( ionode ) THEN
           !
           CALL seqopn( iunocc, TRIM( prefix )//'.occup', 'FORMATTED', exst )
           !
           WRITE( iunocc, * ) ns
           !
           CLOSE( UNIT = iunocc, STATUS = 'KEEP' )
           !
        END IF
        !
     END IF
     !
     ! ... In the US case we need to recompute the self consistent term in
     ! ... the nonlocal potential.
     !
     CALL newd()
     !
     ! ... write the potential (and rho) on file
     !     
     IF ( imix >= 0 ) CALL io_pot( 1, TRIM( prefix )//'.rho', rho_save, nspin )
     !
     CALL io_pot( 1, TRIM( prefix )//'.pot', vr, nspin )     
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
     IF ( lpath ) THEN
        !
        CALL flush( stdout )
        !
     ELSE
        !
#if defined (FLUSH)
        CALL flush( stdout )
#endif
        !
     END IF
     !
     IF ( ( conv_elec .OR. MOD( iter, iprint ) == 0 ) .AND. &
          (.NOT. lmd) ) THEN
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
        IF ( lgauss .OR. ltetra ) then
           IF (two_fermi_energies) then
              WRITE( stdout, 9041 ) ef_up * rytoev, ef_dw * rytoev
           ELSE
              WRITE( stdout, 9040 ) ef * rytoev
           END IF
        END IF
        !
     END IF
     !
     IF ( ( ABS( charge - nelec ) / charge ) > 1.D-7 ) &
        WRITE( stdout, 9050 ) charge
     !
     IF ( ( imix >= 0 ) .AND. &
          ( ABS( charge_new - nelec ) / charge_new ) > 1.D-7 ) &
        WRITE( stdout, 9051 ) charge_new
     !
     !
     etot = eband + ( etxc - etxcc ) + ewld + ehart + deband + demet
     !
     IF ( lda_plus_u ) etot = etot + eth
     IF ( tefield )    etot = etot + etotefield
     !
     IF ( ( conv_elec .OR. MOD( iter, iprint ) == 0 ) .AND. &
          ( .NOT. lmd ) ) THEN
        !  
        IF ( imix >= 0 ) THEN
           !
           IF ( dr2 > eps8 ) THEN
              WRITE( stdout, 9081 ) etot, dr2
           ELSE
              WRITE( stdout, 9083 ) etot, dr2
           END IF
           !
        ELSE
           !
           WRITE( stdout, 9086 ) etot, dr2
           !
        END IF
        !
        WRITE( stdout, 9060 ) &
            eband, ( eband + deband ), ehart, ( etxc - etxcc ), ewld
        !
        IF ( tefield ) WRITE( stdout, 9061 ) etotefield
        IF ( lda_plus_u ) WRITE( stdout, 9065 ) eth
        IF ( degauss /= 0.0 ) WRITE( stdout, 9070 ) demet
        !
     ELSE IF ( conv_elec .AND. lmd ) THEN
        !
        IF ( imix >= 0 ) THEN
           !   
           IF ( dr2 > eps8 ) THEN
              WRITE( stdout, 9081 ) etot, dr2
           ELSE
              WRITE( stdout, 9083 ) etot, dr2
           END IF
           !   
        ELSE
           !   
           WRITE( stdout, 9086 ) etot, dr2
           !   
        END IF
        !
     ELSE
        !
        IF ( imix >=  0 ) THEN
           !   
           IF ( dr2 > eps8 ) THEN
              WRITE( stdout, 9080 ) etot, dr2
           ELSE
              WRITE( stdout, 9082 ) etot, dr2
           END IF
           !   
        ELSE
           !   
           WRITE( stdout, 9085 ) etot, dr2
           !   
        END IF
        !
     END IF
     !
     IF ( lsda ) WRITE( stdout, 9017 ) magtot, absmag
     !
     IF ( noncolin .AND. domag ) &
        WRITE( stdout, 9018 ) ( magtot_nc(i), i = 1, 3 ), absmag
     !
     IF ( i_cons == 3 .OR. i_cons == 4 )  WRITE(stdout, 9071) bfield(1), &
                                                    bfield(2),bfield(3)
     IF ( i_cons == 5 )  WRITE(stdout, 9072) bfield(3)
     IF (i_cons.ne.0.and.i_cons.lt.4) WRITE( stdout,9073) lambda


     !
     IF ( lpath ) THEN
        !
        CALL flush( stdout )
        !
     ELSE
        !
#if defined (FLUSH)
        CALL flush( stdout )
#endif
        !
     END IF
     !
     IF ( conv_elec ) THEN
        !
        WRITE( stdout, 9110 )
        !
        ! ... jump to the end
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
     IF ( imix >= 0 ) rho = rho_save
     !
  END DO
  !
  WRITE( stdout, 9101 )
  WRITE( stdout, 9120 )
  !
  IF ( lpath ) THEN
     !
     CALL flush( stdout )
     !
  ELSE
     !
#if defined (FLUSH)
     CALL flush( stdout )
#endif
     !
  END IF
  !
  IF ( output_drho /= ' ' ) CALL remove_atomic_rho()
  !
  CALL stop_clock( 'electrons' )
  !
  RETURN
  !
  ! ... formats
  !
9000 FORMAT(/'     total cpu time spent up to now is ',F9.2,' secs')
9001 FORMAT(/'     Self-consistent Calculation')
9002 FORMAT(/'     Band Structure Calculation')
9010 FORMAT(/'     iteration #',I3,'     ecut=',F9.2,' ryd',5X,'beta=',F4.2)
9015 FORMAT(/' ------ SPIN UP ------------'/)
9016 FORMAT(/' ------ SPIN DOWN ----------'/)
9017 FORMAT(/'     total magnetization       =', F9.2,' Bohr mag/cell', &
            /'     absolute magnetization    =', F9.2,' Bohr mag/cell')
9018 FORMAT(/'     total magnetization       =',3f9.2,' Bohr mag/cell' &
       &   ,/'     absolute magnetization    =', f9.2,' Bohr mag/cell')
9020 FORMAT(/'          k =',3F7.4,'     band energies (ev):'/)
9021 FORMAT(/'          k =',3F7.4,' (',I6,' PWs)   bands (ev):'/)
9030 FORMAT( '  ',8F9.4)
9032 FORMAT(/'     occupation numbers ')
9041 FORMAT(/'     the spin up/dw Fermi energies are ',2F10.4,' ev')
9040 FORMAT(/'     the Fermi energy is ',F10.4,' ev')
9050 FORMAT(/'     integrated charge         =',F15.8)
9051 FORMAT(/'     integrated charge_new     =',F15.8)
9060 FORMAT(/'     band energy sum           =',  F15.8,' ryd' &
            /'     one-electron contribution =',  F15.8,' ryd' &
            /'     hartree contribution      =',  F15.8,' ryd' &
            /'     xc contribution           =',  F15.8,' ryd' &
            /'     ewald contribution        =',  F15.8,' ryd' )
9061 FORMAT( '     electric field correction =',  F15.8,' ryd' )
9065 FORMAT( '     Hubbard energy            =',F15.8,' ryd')
9070 FORMAT( '     correction for metals     =',F15.8,' ryd')
9071 FORMAT( '     Magnetic field            =',3F12.7,' ryd')
9072 FORMAT( '     Magnetic field            =', F12.7,' ryd')
9073 FORMAT( '     lambda                    =', F11.2,' ryd')

9080 FORMAT(/'     total energy              =',0PF15.8,' ryd' &
            /'     estimated scf accuracy    <',0PF15.8,' ryd')
9081 FORMAT(/'!    total energy              =',0PF15.8,' ryd' &
            /'     estimated scf accuracy    <',0PF15.8,' ryd')
9082 FORMAT(/'     total energy              =',0PF15.8,' ryd' &
            /'     estimated scf accuracy    <',1PE15.1,' ryd')
9083 FORMAT(/'!    total energy              =',0PF15.8,' ryd' &
            /'     estimated scf accuracy    <',1PE15.1,' ryd')
9085 FORMAT(/'     total energy              =',0PF15.8,' ryd' &
            /'     potential mean squ. error =',1PE15.1,' ryd^2')
9086 FORMAT(/'!    total energy              =',0PF15.8,' ryd' &
            /'     potential mean squ. error =',1PE15.1,' ryd^2')
! 9090 FORMAT(/'     the final potential is written on file ',A14)
! 9100 FORMAT(/'     this iteration took ',F9.2,' cpu secs')
9101 FORMAT(/'     End of self-consistent calculation')
9102 FORMAT(/'     End of band structure calculation')
9110 FORMAT(/'     convergence has been achieved')
9120 FORMAT(/'     convergence NOT achieved, stopping')
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     SUBROUTINE compute_magnetization()
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
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
       USE check_stop, ONLY : global_check_stop_now => check_stop_now
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
       USE kinds
       !
       IMPLICIT NONE
       !   
       REAL(kind=DP) :: delta_e
       !
       INTEGER :: i, ipol
       !
       delta_e = 0.d0
       !
       DO ipol=1, nspin
          DO i = 1, nrxx
             delta_e = delta_e - rho(i,ipol) * vr(i,ipol)
          END DO
       END DO
       !
       delta_e = omega * delta_e / (nr1 * nr2 * nr3)
       !
#ifdef __PARA
       CALL reduce (1, delta_e)
#endif
       !
       RETURN
       !
     END FUNCTION delta_e
     !
END SUBROUTINE electrons
