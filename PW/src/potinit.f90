!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE potinit()
  !----------------------------------------------------------------------------
  !
  ! ... This routine initializes the self consistent potential in the array
  ! ... vr. There are three possible cases:
  !
  ! ... a) the code is restarting from a broken run:
  ! ...    read rho from data stored during the previous run
  ! ... b) the code is performing a non-scf calculation following a scf one:
  ! ...    read rho from the file produced by the scf calculation
  ! ... c) the code starts a new calculation:
  ! ...    calculate rho as a sum of atomic charges
  ! 
  ! ... In all cases the scf potential is recalculated and saved in vr
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : pi
  USE io_global,            ONLY : stdout
  USE cell_base,            ONLY : alat, omega
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE starting_scf,         ONLY : starting_pot
  USE klist,                ONLY : nelec
  USE lsda_mod,             ONLY : lsda, nspin
  USE fft_base,             ONLY : dfftp
  USE gvect,                ONLY : ngm, gstart, g, gg, ig_l2g
  USE gvecs,                ONLY : doublegrid
  USE control_flags,        ONLY : lscf, gamma_only, restart, sic
  USE scf,                  ONLY : rho, rho_core, rhog_core, &
                                   vltot, v, vrs, kedtau
  USE xc_lib,               ONLY : xclib_dft_is
  USE ener,                 ONLY : ehart, etxc, vtxc, epaw, esol, vsol
  USE ldaU,                 ONLY : lda_plus_u, Hubbard_lmax, eth, &
                                   niter_with_fixed_ns, lda_plus_u_kind, &
                                   nsg, nsgnew, apply_U, hub_pot_fix, &
                                   orbital_resolved
  USE noncollin_module,     ONLY : noncolin, domag, report, lforcet
  USE io_files,             ONLY : restart_dir, input_drho, check_file_exist
  USE mp,                   ONLY : mp_sum
  USE mp_bands ,            ONLY : intra_bgrp_comm, root_bgrp
  USE io_global,            ONLY : ionode, ionode_id
  USE io_rho_xml,           ONLY : read_scf
  USE io_base,              ONLY : read_rhog
  USE fft_rho,              ONLY : rho_g2r, rho_r2g
  !
  USE uspp,                 ONLY : becsum
  USE paw_variables,        ONLY : okpaw, ddd_paw
  USE paw_init,             ONLY : PAW_atomic_becsum
  USE paw_onecenter,        ONLY : PAW_potential
  !
  USE pwcom,                ONLY : report_mag 
  USE rism_module,          ONLY : lrism, rism_init3d, rism_calc3d
  !
#if defined (__ENVIRON)
  USE plugin_flags,         ONLY : use_environ
  USE environ_pw_module,    ONLY : calc_environ_potential
#endif
  !
  IMPLICIT NONE
  !
  REAL(DP)                  :: charge           ! the starting charge
  REAL(DP)                  :: etotefield       !
  REAL(DP)                  :: fact
  INTEGER                   :: is
  LOGICAL                   :: exst 
  CHARACTER(LEN=320)        :: filename
  COMPLEX (DP), ALLOCATABLE :: work(:,:)
  !
  CALL start_clock('potinit')
  !
  filename = TRIM (restart_dir( )) // 'charge-density'
#if defined __HDF5
  exst     =  check_file_exist( TRIM(filename) // '.hdf5' )
#else 
  exst     =  check_file_exist( TRIM(filename) // '.dat' )
#endif
  !
  IF ( starting_pot == 'file' .AND. exst ) THEN
     !
     ! ... Cases a) and b): the charge density is read from file
     ! ... this also reads rho%ns if DFT+U, rho%bec if PAW, rho%kin if metaGGA
     !
     ! ... if we restart from a preexisting charge density, the eigenstates
     ! ... are considered stable and we can apply orbital-resolved Hubbard 
     ! ... corrections starting from the first iteration
     IF ( orbital_resolved ) apply_U = .TRUE.
     !
     IF ( .NOT.lforcet ) THEN
        CALL read_scf ( rho, nspin, gamma_only )
     ELSE
        IF ( okpaw )  CALL errore( 'potinit', &
                                   'force theorem with PAW not implemented', 1 )
        !
        ! ... 'force theorem' calculation of MAE: read rho only from previous
        ! ... lsda calculation, set noncolinear magnetization from angles
        ! ... (not if restarting! the charge density saved to file in that
        ! ...  case has already the required magnetization direction)
        !
        CALL read_rhog ( filename, root_bgrp, intra_bgrp_comm, &
             ig_l2g, nspin, rho%of_g, gamma_only )
        IF ( .NOT. restart ) &
           CALL nc_magnetization_from_lsda ( dfftp%ngm, nspin, rho%of_g )
     END IF
     !
     IF ( lscf ) THEN
        !
        WRITE( stdout, '(/5X, &
             & "The initial density is read from file :"/5X,A,/)' ) &
            TRIM( filename )
        !
     ELSE
        !
        WRITE( stdout, '(/5X, &
             & "The potential is recalculated from file :"/5X,A,/)' ) &
            TRIM( filename )
        !
     END IF
     !
     IF ( input_drho /= ' ' ) THEN
        !
        filename = TRIM( restart_dir( )) // input_drho
        CALL read_rhog ( filename, root_bgrp, intra_bgrp_comm, &
             ig_l2g, nspin, v%of_g, gamma_only )
        ! 
        WRITE( UNIT = stdout, &
               FMT = '(/5X,"a scf correction to at. rho is read from",A)' ) &
            TRIM( filename )
        !
        ALLOCATE( work( dfftp%ngm, nspin ) )
        CALL atomic_rho_g( work, nspin )
        rho%of_g(:,1) = work(:,1) + v%of_g(:,1)
        DEALLOCATE(work)
        !
     END IF
     !
  ELSE
     !
     ! ... Case c): the potential is built from a superposition 
     ! ... of atomic charges contained in the array rho_at
     !
     IF ( starting_pot == 'file' .AND. .NOT. exst ) &
        WRITE( stdout, '(5X,"Cannot read rho : file not found")' )
     !
     WRITE( UNIT = stdout, &
            FMT = '(/5X,"Initial potential from superposition of free atoms")' )
     !
     CALL atomic_rho_g( rho%of_g, nspin )

     ! ... in the DFT+U(+V) case set the initial value of ns (or nsg)
     !
     IF (lda_plus_u) THEN
        !
        IF (lda_plus_u_kind == 0) THEN
           IF ( hub_pot_fix ) &
              CALL errore( 'potinit', &
                     'cannot apply Hubbard alpha without &
                     &restarting from a converged potential', 1 )
           IF ( orbital_resolved .AND. (.NOT. apply_U) ) THEN
              WRITE( stdout, '(/,5X,47("="))')
              WRITE( stdout, '(/,5X,"Not restarting from a converged ", &
                                &    "potential:",/,5X,             &
              & "Orbital-resolved Hubbard corrections not yet active")')
              WRITE( stdout, '(/,5X,47("="))')
              !
           ENDIF
           IF (noncolin) THEN
              CALL init_ns_nc() 
           ELSE 
              CALL init_ns()
           ENDIF   
        ELSEIF (lda_plus_u_kind == 1) THEN
           IF (noncolin) THEN
              CALL init_ns_nc()
           ELSE
              CALL init_ns()
           ENDIF
        ELSEIF (lda_plus_u_kind == 2) THEN
           CALL init_nsg()
        ENDIF
        !
     ENDIF

     ! ... in the paw case uses atomic becsum
     IF ( okpaw )      CALL PAW_atomic_becsum()
     !
     IF ( input_drho /= ' ' ) THEN
        !
        filename = TRIM( restart_dir( )) // input_drho
        CALL read_rhog ( filename, root_bgrp, intra_bgrp_comm, &
             ig_l2g, nspin, v%of_g, gamma_only )
        !
        WRITE( UNIT = stdout, &
               FMT = '(/5X,"a scf correction to at. rho is read from",A)' ) &
            TRIM( filename )
        !
        rho%of_g = rho%of_g + v%of_g
        !
     END IF
     !
  END IF
  !
  ! ... check the integral of the starting charge, renormalize if needed
  !
  charge = 0.D0
  IF ( gstart == 2 ) THEN
     charge = omega*REAL( rho%of_g(1,1) )
  END IF
  CALL mp_sum(  charge , intra_bgrp_comm )
  !
  IF ( lscf .AND. ABS( charge - nelec ) > ( 1.D-7 * charge ) ) THEN
     !
     IF ( charge > 1.D-8 .AND. nat > 0 ) THEN
        WRITE( stdout, '(/,5X,"starting charge ",F12.4, &
                         & ", renormalised to ",F12.4)') charge, nelec
        rho%of_g = rho%of_g / charge * nelec
     ELSE 
        WRITE( stdout, '(/,5X,"Starting from uniform charge")')
        rho%of_g(:,1:nspin) = (0.0_dp,0.0_dp)
        IF ( gstart == 2 ) rho%of_g(1,1) = nelec / omega
     ENDIF
     !
  ELSE IF ( .NOT. lscf .AND. ABS( charge - nelec ) > (1.D-3 * charge ) ) THEN
     !
     CALL errore( 'potinit', 'starting and expected charges differ', 1 )
     !
  END IF
  !
  ! ... bring starting rho from G- to R-space
  !
  CALL rho_g2r (dfftp, rho%of_g, rho%of_r)
  !
  ! .... initialize polaron density as spin-density
  !
  IF(sic) THEN
     rho%pol_g(:,1) = rho%of_g(:,2)
     rho%pol_g(:,2) = rho%of_g(:,2)
     CALL rho_g2r (dfftp, rho%pol_g, rho%pol_r)
  END IF   
  !
  IF  ( xclib_dft_is('meta') ) THEN
     !
     IF (starting_pot /= 'file') THEN
        ! ... define a starting (TF) guess for rho%kin_r from rho%of_r
        fact = (3.d0/5.d0)*(3.d0*pi*pi)**(2.0/3.0)
        IF ( nspin == 1) THEN
           rho%kin_r(:,1) = fact * abs(rho%of_r(:,1))**(5.0/3.0)
        ELSE ! IF ( nspin == 2) THEN 
           ! ... NB: for LSDA rho is (tot,magn), rho_kin is (up,down) 
           rho%kin_r(:,1) = ( rho%of_r(:,1) + rho%of_r(:,2) ) / 2.0_dp
           rho%kin_r(:,2) = ( rho%of_r(:,1) - rho%of_r(:,2) ) / 2.0_dp
           ! multiplication by nspin: see Eq.2.9 of 10.1103/PhysRevA.20.397
           DO is = 1, nspin
              rho%kin_r(:,is) = fact * abs(rho%kin_r(:,is)*nspin)**(5.0/3.0)/nspin
           END DO
        END IF
        ! ... bring it to g-space
        CALL rho_r2g (dfftp, rho%kin_r, rho%kin_g)
     ELSE
        ! ... rho%kin was read from file in G-space, bring it to R-space
        CALL rho_g2r (dfftp, rho%kin_g, rho%kin_r)
     ENDIF
     !
  END IF
  !
  ! ... initialize 3D-RISM
  !
  IF (lrism) CALL rism_init3d()
  !
  ! ... plugin contribution to local potential
  !
#if defined(__LEGACY_PLUGINS)
  CALL plugin_scf_potential(rho, .FALSE., -1.d0, vltot)
#endif 
#if defined (__ENVIRON)
  IF (use_environ) CALL calc_environ_potential(rho, .FALSE., -1.D0, vltot)
#endif
  !
  ! ... compute the potential and store it in v
  !
  CALL v_of_rho( rho, rho_core, rhog_core, &
                 ehart, etxc, vtxc, eth, etotefield, charge, v )
  IF (okpaw) CALL PAW_potential(rho%bec, ddd_paw, epaw)
  !
  ! ... calculate 3D-RISM to get the solvation potential
  !
  IF (lrism) CALL rism_calc3d(rho%of_g(:, 1), esol, vsol, v%of_r, -1.0_DP)
  !
  ! ... define the total local potential (external+scf)
  !
  CALL set_vrs( vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid )
  ! ... write on output the parameters used in the DFT+U(+V) calculation
  !
  IF ( lda_plus_u ) THEN
     !
     IF (niter_with_fixed_ns>0) &
     WRITE( stdout, '(5X,"Number of Hubbard iterations with fixed ns =",I3)') &
         niter_with_fixed_ns
     !
     ! ... info about starting occupations
     WRITE( stdout, '(/5X,"STARTING HUBBARD OCCUPATIONS:")')
     !
     IF (lda_plus_u_kind == 0) THEN
        IF (noncolin) THEN
           CALL write_ns_nc() 
        ELSE   
           CALL write_ns()
        ENDIF
     ELSEIF (lda_plus_u_kind == 1) THEN
        IF (noncolin) THEN
           CALL write_ns_nc()
        ELSE
           CALL write_ns()
        ENDIF
     ELSEIF (lda_plus_u_kind == 2) THEN
        nsgnew = nsg
        IF(noncolin) THEN
           CALL write_nsg_nc()
        ELSE
           CALL write_nsg()
        ENDIF
     ENDIF
     !
  END IF
  !
  IF ( report /= 0 .AND. &
       noncolin .AND. domag .AND. lscf ) CALL report_mag()
  !
  CALL stop_clock('potinit')
  !
  RETURN
  !
END SUBROUTINE potinit
!
!-------------
SUBROUTINE nc_magnetization_from_lsda ( ngm, nspin, rho )
  !-------------
  !
  USE kinds,     ONLY: dp
  USE constants, ONLY: pi
  USE io_global, ONLY: stdout
  USE noncollin_module, ONLY: angle1, angle2
  !
  IMPLICIT NONE
  INTEGER, INTENT (in):: ngm, nspin
  COMPLEX(dp), INTENT (inout):: rho(ngm,nspin)
  !
  IF ( nspin < 4 ) RETURN
  !---  
  !  set up noncollinear m_x,y,z from collinear m_z (AlexS) 
  !
  WRITE(stdout,*)
  WRITE(stdout,*) '-----------'
  WRITE(stdout,'("Spin angles Theta, Phi (degree) = ",2f8.4)') &
       angle1(1)/PI*180.d0, angle2(1)/PI*180.d0 
  WRITE(stdout,*) '-----------'
  !
  ! now set rho(2)=magn*sin(theta)*cos(phi)   x
  !         rho(3)=magn*sin(theta)*sin(phi)   y
  !         rho(4)=magn*cos(theta)            z
  !
  rho(:,4) = rho(:,2)*cos(angle1(1))
  rho(:,2) = rho(:,2)*sin(angle1(1))
  rho(:,3) = rho(:,2)*sin(angle2(1))
  rho(:,2) = rho(:,2)*cos(angle2(1))
  !
  RETURN
  !
END SUBROUTINE nc_magnetization_from_lsda
