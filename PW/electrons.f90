!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!  
#include "machine.h"
!
! ... uncomment the following line to use the "old ethr" 
! ... for the first iteration  
!
!#define OLDSTYLE
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
  USE io_global,            ONLY : stdout
  USE brilz,                ONLY : at, bg, alat, omega, tpiba2
  USE basis,                ONLY : nat, ntyp, ityp, tau   
  USE gvect,                ONLY : ngm, gstart, nr1, nr2, nr3, nrx1, nrx2, &
                                   nrx3, nrxx, nl, g, gg, ecutwfc, gcutm
  USE gsmooth,              ONLY : doublegrid  
  USE klist,                ONLY : xk, degauss, nelec, ngk, nks, nkstot, &
                                   lgauss    
  USE lsda_mod,             ONLY : lsda, nspin  
  USE ktetra,               ONLY : ltetra  
  USE pseud,                ONLY : zv    
  USE vlocal,               ONLY : strf, vnew  
  USE wvfct,                ONLY : nbnd, et, gamma_only  
  USE ener,                 ONLY : etot, eband, deband, ehart, vtxc, etxc, &
                                   etxcc, ewld, demet, ef  
  USE scf,                  ONLY : rho, rho_save, vr, vltot, vrs, rho_core
  USE control_flags,                ONLY : mixing_beta, tr2, time_max, ethr, ngm0, &
                                   niter, nmix, imix, iprint, istep, iswitch, &
                                   lscf, lneb, lmd, conv_elec, restart, &
                                   reduce_io  
  USE io_files,             ONLY : prefix, iunwfc, iunocc, nwordwfc, iunneb, &
                                   output_drho, iunexit, exit_file
  USE ldaU,                 ONLY : ns, nsnew, eth, Hubbard_U, &
                                   niter_with_fixed_ns, Hubbard_lmax, &
                                   lda_plus_u  
  USE extfield,             ONLY : tefield, etotefield  
  USE bp,                   ONLY : lberry  
  USE wavefunctions_module, ONLY : evc
#ifdef __PARA
  USE para,                 ONLY : me, mypool, npp, ncplane
  USE mp,                   ONLY : mp_barrier
#endif
  !
  IMPLICIT NONE
  !
  ! ... a few local variables
  !  
#ifdef __PARA
  INTEGER :: &
      ngkp(npk)       !  number of plane waves summed on all nodes
#define NRXX ncplane*npp(me)
  ! This is needed in mix_pot whenever nproc is not a divisor of nr3.
#else
#define NRXX nrxx
#endif
  CHARACTER :: &
      flmix * 42      !
  REAL(KIND=DP) :: &
      de,            &!  the correction energy
      dr2,           &!  the norm of the diffence between potential
      charge,        &!  the total charge
      mag,           &!  local magnetization
      magtot,        &!  total magnetization
      absmag,        &!  total absolute magnetization
      tcpu          !
  INTEGER :: &
      i,             &!  counter on polarization
      ir,            &!  counter on the mesh points
      ig,            &!
      ik,            &!  counter on k points
      ibnd,          &!  counter on bands
      idum,          &!  dummy counter on iterations
      iter,          &!  counter on iterations
      ik_             !  used to read ik from restart file
  INTEGER :: &
      ldim2           !
  REAL (KIND=DP) :: &
      ehart_new,     &!
      etxc_new,      &!
      vtxc_new,      &!
      charge_new      !
  REAL (KIND=DP) :: &
      ethr_min        ! minimal threshold for diagonalization at the first scf
                      ! iteration of a MD calculation 
  REAL (KIND=DP), EXTERNAL :: ewald, get_clock
  LOGICAL :: &
      exst,          &!
      file_exists     !  .TRUE. if a soft exit has been required
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
     IF ( ik_ == - 1000 ) THEN
        !
        conv_elec = .TRUE.
        !
        ! ...jump to the end 
        !
        IF ( output_drho /= ' ' ) CALL remove_atomic_rho
        CALL stop_clock( 'electrons' )
        !
        RETURN
        !
     END IF
     !
  END IF
  !
  IF ( lscf ) THEN
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
  END IF 
  !
  ! ... Convergence threshold for iterative diagonalization
  !
  ! ... for the first scf iteration of each ionic step (except than for the
  ! ... first) the threshold is fixed to a default value of 1.D-5
  !
#if defined (OLDSTYLE)
  IF ( .FALSE. ) THEN
#else
  IF ( istep > 1 ) THEN
#endif     
     !
     ethr = 1.D-5
     !
  END IF  
  !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%          iterate !          %%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
  IF ( imix >= 0 ) ngm0 = ngm
  !
  DO idum = 1, niter
     !
     tcpu = get_clock( 'PWSCF' )
     WRITE( stdout, 9000 ) tcpu
     !
     IF ( imix >= 0 ) CALL DCOPY( ( nspin * nrxx), rho, 1, rho_save, 1 )
     !
     iter = iter + 1
     !
     IF ( lscf ) THEN
        WRITE( stdout, 9010 ) iter, ecutwfc, mixing_beta
     ELSE
        WRITE( stdout, 9009 )
     END IF
#ifdef FLUSH
     CALL flush( 6 )
#endif
     !
     ! ... Convergence threshold for iterative diagonalization
     ! ... is automatically updated during self consistency
     !
     IF ( lscf .AND. iter > 1 .AND. ik_ == 0 ) THEN
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
     CALL c_bands( iter, ik_, dr2 )
     !
     ! ... skip all the rest if not lscf
     !
     IF ( .NOT. lscf ) THEN
        !
        conv_elec = .TRUE.
        !
#ifdef __PARA
        CALL poolrecover( et, nbnd, nkstot, nks )
#endif
        !
        DO ik = 1, nkstot
           IF ( lsda ) THEN
              IF ( ik == 1 ) WRITE( stdout, 9015 )
              IF ( ik == ( 1 + nkstot / 2 ) ) WRITE( stdout, 9016 )
           END IF
           WRITE( stdout, 9020 ) ( xk(i,ik), i = 1, 3 )
           WRITE( stdout, 9030 ) ( et(ibnd,ik) * 13.6058, ibnd = 1, nbnd )
        END DO
        !
        ! ... do a Berry phase polarization calculation if required
        !
        IF ( lberry ) CALL c_phase()
        !
        ! ... jump to the end
        !
        IF ( output_drho /= ' ' ) CALL remove_atomic_rho()
        CALL stop_clock( 'electrons' )
        !
        RETURN
        !
     END IF
     !
     ! ... the program checks if the maximum CPU time has been exceeded
     !
     CALL check_cpu_time()
     !
     ! ... the program checks if the user has required a soft exit
     !
     CALL check_soft_exit()
     !
     CALL sum_band()
     !
     IF ( lda_plus_u ) CALL write_ns()
     !
     !IF ( iter == 1 .AND. lda_plus_u .AND. &
     !     input_pot /=  ' ' .AND. istep == 1 ) CALL ns_adj()
     !
     ! ... calculate total and absolute magnetization
     !
     IF ( lsda ) CALL compute_magnetization()
     !
     CALL v_of_rho( rho, rho_core, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
                    nrxx, nl, ngm, gstart, nspin, g, gg, alat, omega, &
                    ehart, etxc, vtxc, charge, vnew )
     !
     CALL delta_e( nr1, nr2, nr3, nrxx, rho, vr, vnew, omega, de, &
                   deband, nspin )
     !
     IF ( imix >= 0 ) THEN
        !
        IF ( lda_plus_u .AND. iter <= niter_with_fixed_ns ) THEN
           !
           ldim2 = ( 2 * Hubbard_lmax + 1 )**2
           CALL DCOPY( ( ldim2 * nspin * nat ), ns, 1, nsnew, 1 )
           !
        END IF
        !
#if defined (OLDSTYLE)
        IF ( .FALSE. ) THEN
#else
        IF ( iter == 1 ) THEN
#endif        
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
        CALL mix_rho( rho, rho_save, nsnew, ns, mixing_beta, dr2, ethr, &
                      ethr_min, iter, nmix, flmix, conv_elec )
        !
        ! ... for the first scf iteration it is controlled that the threshold 
        ! ... is small enought for the diagonalization to be adequate
        !
#if defined (OLDSTYLE)
        IF ( .FALSE. ) THEN
#else
        IF ( iter == 1 .AND. ethr >= ethr_min ) THEN
#endif
           !
           ! ... a new diagonalization is needed       
           !
           WRITE( stdout, '(/,5X,"WARNING: ethr not adequate")' )
           WRITE( stdout, '(5X,"         a new diagonalization is needed",/)' )
           !
           CALL c_bands( iter, ik_, dr2 )
           !
           ! ... the program checks if the maximum CPU time has been exceeded
           !
           CALL check_cpu_time()
           !
           ! ... the programs checks if the user has required a soft exit
           !
           CALL check_soft_exit()
           !
           CALL sum_band()
           !
           IF ( lda_plus_u ) CALL write_ns()
           !
           ! ... calculate total and absolute magnetization
           !
           IF ( lsda ) CALL compute_magnetization()
           !
           CALL v_of_rho( rho, rho_core, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
                          nrxx, nl, ngm, gstart, nspin, g, gg, alat, omega, &
                          ehart, etxc, vtxc, charge, vnew )
           !
           CALL delta_e( nr1, nr2, nr3, nrxx, rho, vr, vnew, omega, de, &
                         deband, nspin )              
           !
           IF ( lda_plus_u .AND. iter <= niter_with_fixed_ns ) THEN
              !
              ldim2 = ( 2 * Hubbard_lmax + 1 )**2
              CALL DCOPY( ( ldim2 * nspin * nat ), ns, 1, nsnew, 1 )                 
              !
           END IF
           !
           ! ... ethr_min is set to a negative number: no check is needed
           !
           CALL mix_rho( rho, rho_save, nsnew, ns, mixing_beta, dr2, ethr, &
                         -1.D0, iter, nmix, flmix, conv_elec )
           !
        END IF             
        !
        CALL DAXPY( ( nspin * nrxx ), -1.D0, vr, 1, vnew, 1 )
        !
        CALL v_of_rho( rho_save, rho_core, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
                       nrxx, nl, ngm, gstart, nspin, g, gg, alat, omega, &
                       ehart_new, etxc_new, vtxc_new, charge_new, vr )
        !
     ELSE 
        !
        ! ... old style potential mixing
        !
        CALL vpack( NRXX, nrxx, nspin, vnew, vr, + 1 )
        !
        CALL mix_potential( ( nspin * NRXX ), vnew, vr, mixing_beta, dr2, tr2, &
                            iter, nmix, flmix, conv_elec )
        !
        CALL vpack( NRXX, nrxx, nspin, vnew, vr, - 1 )
        !
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
        IF ( iter > niter_with_fixed_ns .AND. imix < 0 ) &
           CALL DCOPY( ( ldim2 * nspin * nat ), nsnew, 1, ns, 1 )    
#ifdef __PARA
        IF ( me == 1 .AND. mypool == 1 ) THEN
#endif
           CALL seqopn( iunocc, TRIM( prefix )//'.occup', 'formatted', exst )
           WRITE( iunocc, * ) ns
           CLOSE( UNIT = iunocc, STATUS = 'KEEP' )
#ifdef __PARA
        END IF
#endif
     END IF
     !
     ! ... In the US case we need to recompute the self consistent term in
     ! ... the nonlocal potential.
     !
     CALL newd()
     !
     ! ... write the potential (and rho) on file
     !     
     IF ( imix >= 0 ) CALL io_pot( +1, TRIM( prefix )//'.rho', rho_save, nspin )
     CALL io_pot( +1, TRIM( prefix )//'.pot', vr, nspin )     
     !
     ! ... save converged wfc if they have not been written previously
     !     
     IF ( nks == 1 .AND. reduce_io ) &
        CALL davcio( evc, nwordwfc, iunwfc, nks, 1 )
     !
     ! ... write recover file
     !
     CALL save_in_electrons( iter, dr2 )
     !
     IF ( ( conv_elec .OR. MOD( iter, iprint )  == 0 ) .AND. &
          iswitch <= 2 ) THEN
        !  
        !IF ( lda_plus_u ) CALL write_ns()
        !
#ifdef __PARA
        DO ik = 1, nks
           ngkp(ik) = ngk(ik)
        END DO
        !
        CALL ireduce( nks, ngkp )
        CALL ipoolrecover( ngkp, 1, nkstot, nks )
        CALL poolrecover( et, nbnd, nkstot, nks )
#endif
        !
        DO ik = 1, nkstot
           !
           IF ( lsda ) THEN
              IF ( ik == 1 ) WRITE( stdout, 9015)
              IF ( ik == ( 1 + nkstot / 2 ) ) WRITE( stdout, 9016)
           END IF
           !
           IF ( conv_elec ) THEN
#ifdef __PARA
              WRITE( stdout, 9021 ) ( xk(i,ik), i = 1, 3 ), ngkp(ik)
#else
              WRITE( stdout, 9021 ) ( xk(i,ik), i = 1, 3 ), ngk(ik)
#endif
           ELSE
              WRITE( stdout, 9020 ) ( xk(i,ik), i = 1, 3 )
           END IF
           !
           WRITE( stdout, 9030 ) ( et(ibnd,ik) * 13.6058, ibnd = 1, nbnd )
           !
        END DO
        !
        IF ( lgauss .OR. ltetra ) WRITE( stdout, 9040 ) ef * 13.6058
        !
     END IF
     !
     IF ( ( ABS( charge - nelec ) / charge ) > 1.D-7 ) &
        WRITE( stdout, 9050 ) charge
     !
     etot = eband + ( etxc - etxcc ) + ewld + ehart + deband + demet
     !
     IF ( lda_plus_u ) etot = etot + eth
     !
     IF ( tefield ) etot = etot + etotefield
     !
     IF ( ( conv_elec .OR. MOD( iter, iprint ) == 0 ) .AND. &
          ( iswitch <= 2 ) ) THEN
        !  
        IF ( imix >= 0 ) THEN
           WRITE( stdout, 9081 ) etot, dr2
        ELSE
           WRITE( stdout, 9086 ) etot, dr2
        END IF
        !
        WRITE( stdout, 9060 ) &
            eband, ( eband + deband ), ehart, ( etxc - etxcc ), ewld
        !
        IF ( tefield ) WRITE( stdout, 9061 ) etotefield
        IF ( lda_plus_u ) WRITE( stdout, 9065 ) eth
        IF ( degauss /= 0.0 ) WRITE( stdout, 9070 ) demet
        !
     ELSE IF ( conv_elec .AND. iswitch > 2 ) THEN
        !
        IF ( imix >= 0 ) THEN
           WRITE( stdout, 9081 ) etot, dr2
        ELSE
           WRITE( stdout, 9086 ) etot, dr2
        END IF
        !
     ELSE
        !
        IF ( imix >=  0 ) THEN
           WRITE( stdout, 9080 ) etot, dr2
        ELSE
           WRITE( stdout, 9085 ) etot, dr2
        END IF
        !
     END IF
     !
     IF ( lsda ) WRITE( stdout, 9017 ) magtot, absmag
     !
#ifdef FLUSH
     CALL flush( 6 )
#endif
     IF ( conv_elec ) THEN
        !
        WRITE( stdout, 9110 )
        !
        ! ... jump to the end
        !
        IF ( output_drho /= ' ' ) CALL remove_atomic_rho
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
     IF ( imix >= 0 ) CALL DCOPY( ( nspin * nrxx), rho_save, 1, rho, 1 )
     !
  END DO
  !
  WRITE( stdout, 9120 )
  !
  ! <------- jump here if not scf
  !
  IF ( output_drho /= ' ' ) CALL remove_atomic_rho
  !
  CALL stop_clock( 'electrons' )
  !
  RETURN
  !
  ! ... formats
  !
9000 FORMAT(/'     total cpu time spent up to now is ',F9.2,' secs')
9009 FORMAT(/'     Band Structure Calculation')
9010 FORMAT(/'     iteration #',I3,'     ecut=',F9.2,' ryd',5X, &
             'beta=',F4.2)
9015 FORMAT(/' ------ SPIN UP ------------'/)
9016 FORMAT(/' ------ SPIN DOWN ----------'/)
9017 FORMAT(/'     total magnetization       =',F9.2,' Bohr mag/cell', &
            /'     absolute magnetization    =',F9.2,' Bohr mag/cell')
9020 FORMAT(/'          k =',3F7.4,'     band energies (ev):'/)
9021 FORMAT(/'          k =',3F7.4,' (',I5,' PWs)   bands (ev):'/)
9030 FORMAT( '  ',8F9.4)
9040 FORMAT(/'     the Fermi energy is ',F10.4,' ev')
9050 FORMAT(/'     integrated charge         =',F15.8)
9060 FORMAT(/'     band energy sum           =',  F15.8,' ryd' &
            /'     one-electron contribution =',  F15.8,' ryd' &
            /'     hartree contribution      =',  F15.8,' ryd' &
            /'     xc contribution           =',  F15.8,' ryd' &
            /'     ewald contribution        =',  F15.8,' ryd' )
9061 FORMAT( '     electric field correction =',  F15.8,' ryd' )
9065 FORMAT( '     Hubbard energy            =',F15.8,' ryd')
9070 FORMAT( '     correction for metals     =',F15.8,' ryd')
9080 FORMAT(/'     total energy              =',0PF15.8,' ryd' &
            /'     estimated scf accuracy    <',0PF15.8,' ryd')
9081 FORMAT(/'!    total energy              =',0PF15.8,' ryd' &
            /'     estimated scf accuracy    <',0PF15.8,' ryd')
9085 FORMAT(/'     total energy              =',0PF15.8,' ryd' &
            /'     potential mean squ. error =',1PE15.1,' ryd^2')
9086 FORMAT(/'!    total energy              =',0PF15.8,' ryd' &
            /'     potential mean squ. error =',1PE15.1,' ryd^2')
9090 FORMAT(/'     the final potential is written on file ',A14)
9100 FORMAT(/'     this iteration took ',F9.2,' cpu secs')
9110 FORMAT(/'     convergence has been achieved')
9120 FORMAT(/'     convergence NOT achieved. stopping ...')
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
       magtot = 0.D0
       absmag = 0.D0
       !
       DO ir = 1, nrxx
          !   
          mag    = rho(ir,1) - rho(ir,2)
          magtot = magtot + mag
          absmag = absmag + ABS( mag )
          !
       END DO
       !
       magtot = magtot * omega / ( nr1 * nr2 * nr3 )
       absmag = absmag * omega / ( nr1 * nr2 * nr3 )
       !
#ifdef __PARA
       CALL reduce( 1, magtot )
       CALL reduce( 1, absmag )
#endif
       !
       RETURN
       !
     END SUBROUTINE compute_magnetization
     !
     !-----------------------------------------------------------------------
     SUBROUTINE check_cpu_time()
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       !     
       tcpu = get_clock( 'PWSCF' )
       !
       IF ( tcpu > time_max ) THEN
          !
          IF ( lneb ) THEN  
             WRITE( iunneb, '(5X,"Maximum CPU time exceeded",2F15.2)' ) &
                 tcpu, time_max
          ELSE
             WRITE( stdout, '(5X,"Maximum CPU time exceeded",2F15.2)' ) &
                 tcpu, time_max
          END IF       
          !
          CALL stop_pw( .FALSE. )
          !
       END IF
       !
     END SUBROUTINE check_cpu_time
     !       
     !
     !-----------------------------------------------------------------------
     SUBROUTINE check_soft_exit()
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       !
       INQUIRE( FILE = TRIM( exit_file ), EXIST = file_exists )
       !
       IF ( file_exists ) THEN
          !
#ifdef __PARA
          !
          ! ... all jobs are syncronized
          !
          CALL mp_barrier()
          !
          IF ( me == 1 .AND. mypool == 1 ) THEN
#endif
             OPEN( UNIT = iunexit, FILE = TRIM( exit_file ), STATUS = "OLD" )
             CLOSE( UNIT = iunexit, STATUS = "DELETE" )
#ifdef __PARA
          END IF
#endif
          !
          IF ( lneb ) THEN  
             WRITE( iunneb, '(/,5X,"WARNING :  soft exit required",/, &
                              & 5X,"stopping in electrons ...",/)' )  
          ELSE
             WRITE( stdout, '(/,5X,"WARNING :  soft exit required",/, &
                              & 5X,"stopping in electrons ...",/)' )  
          END IF
          !  
          CALL stop_pw( .FALSE. )
          !
       END IF              
       !
     END SUBROUTINE check_soft_exit 
     !
END SUBROUTINE electrons
