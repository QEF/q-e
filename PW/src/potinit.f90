!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
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
  USE basis,                ONLY : starting_pot
  USE klist,                ONLY : nelec
  USE lsda_mod,             ONLY : lsda, nspin
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE gvect,                ONLY : ngm, gstart, g, gg, ig_l2g
  USE gvecs,                ONLY : doublegrid
  USE control_flags,        ONLY : lscf, gamma_only
  USE scf,                  ONLY : rho, rho_core, rhog_core, &
                                   vltot, v, vrs, kedtau
  USE funct,                ONLY : dft_is_meta
  USE ener,                 ONLY : ehart, etxc, vtxc, epaw
  USE ldaU,                 ONLY : lda_plus_u, Hubbard_lmax, eth, &
                                   niter_with_fixed_ns
  USE noncollin_module,     ONLY : noncolin, report
  USE io_files,             ONLY : tmp_dir, prefix, postfix, input_drho, check_file_exist
  USE spin_orb,             ONLY : domag, lforcet
  USE mp,                   ONLY : mp_sum
  USE mp_bands ,            ONLY : intra_bgrp_comm, root_bgrp
  USE io_global,            ONLY : ionode, ionode_id
  USE io_rho_xml,           ONLY : read_scf
  USE io_base,              ONLY : read_rhog
  USE fft_rho,              ONLY : rho_g2r, rho_r2g
  !
  USE uspp,                 ONLY : becsum
  USE paw_variables,        ONLY : okpaw, ddd_PAW
  USE paw_init,             ONLY : PAW_atomic_becsum
  USE paw_onecenter,        ONLY : PAW_potential
  !
  IMPLICIT NONE
  !
  REAL(DP)              :: charge           ! the starting charge
  REAL(DP)              :: etotefield       !
  REAL(DP)              :: fact
  INTEGER               :: is
  LOGICAL               :: exst 
  CHARACTER(LEN=320)    :: filename
  !
  CALL start_clock('potinit')
  !
  filename = TRIM(tmp_dir) // TRIM (prefix) // postfix // 'charge-density'
#if defined __HDF5
  exst     =  check_file_exist( TRIM(filename) // '.hdf5' )
#else 
  exst     =  check_file_exist( TRIM(filename) // '.dat' )
#endif
  !
  IF ( starting_pot == 'file' .AND. exst ) THEN
     !
     ! ... Cases a) and b): the charge density is read from file
     ! ... this also reads rho%ns if lda+U, rho%bec if PAW, rho%kin if metaGGA
     !
     IF ( .NOT.lforcet ) THEN
        CALL read_scf ( rho, nspin, gamma_only )
     ELSE
        !
        ! ... 'force theorem' calculation of MAE: read rho only from previous
        ! ... lsda calculation, set noncolinear magnetization from angles
        !
        CALL read_rhog ( filename, root_bgrp, intra_bgrp_comm, &
             ig_l2g, nspin, rho%of_g, gamma_only )
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

     ! ... in the lda+U case set the initial value of ns
     IF (lda_plus_u) THEN
        !
        IF (noncolin) THEN
           CALL init_ns_nc()
        ELSE
           CALL init_ns()
        ENDIF
        !
     ENDIF

     ! ... in the paw case uses atomic becsum
     IF ( okpaw )      CALL PAW_atomic_becsum()
     !
     IF ( input_drho /= ' ' ) THEN
        !
        IF ( nspin > 1 ) CALL errore &
             ( 'potinit', 'spin polarization not allowed in drho', 1 )
        !
        filename = TRIM(tmp_dir) // TRIM (prefix) // postfix // input_drho
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
        WRITE( stdout, '(/,5X,"starting charge ",F10.5, &
                         & ", renormalised to ",F10.5)') charge, nelec
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
  IF  ( dft_is_meta() ) THEN
     IF (starting_pot /= 'file') THEN
        ! ... define a starting (TF) guess for rho%kin_r from rho%of_r
        ! ... to be verified for LSDA: rho is (tot,magn), rho_kin is (up,down)
        fact = (3.d0/5.d0)*(3.d0*pi*pi)**(2.0/3.0)
        DO is = 1, nspin
           rho%kin_r(:,is) = fact * abs(rho%of_r(:,is)*nspin)**(5.0/3.0)/nspin
        END DO
        !if (nspin==2) then
        !     rho%kin_r(:,1) = fact * abs(rho%of_r(:,1)+rho%of_r(:,2))**(5.0/3.0)/2.0
        !     rho%kin_r(:,2) = fact * abs(rho%of_r(:,1)-rho%of_r(:,2))**(5.0/3.0)/2.0
        !endif
        ! ... bring it to g-space
        CALL rho_r2g (dfftp, rho%kin_r, rho%kin_g)
     ELSE
        ! ... rho%kin was read from file in G-space, bring it to R-space
        CALL rho_g2r (dfftp, rho%kin_g, rho%kin_r)
     ENDIF
     !
  END IF
  !
  ! ... plugin contribution to local potential
  !
  CALL plugin_scf_potential(rho,.FALSE.,-1.d0,vltot)
  !
  ! ... compute the potential and store it in v
  !
  CALL v_of_rho( rho, rho_core, rhog_core, &
                 ehart, etxc, vtxc, eth, etotefield, charge, v )
  IF (okpaw) CALL PAW_potential(rho%bec, ddd_PAW, epaw)
  !
  ! ... define the total local potential (external+scf)
  !
  CALL set_vrs( vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid )
  !
  ! ... write on output the parameters used in the lda+U calculation
  !
  IF ( lda_plus_u ) THEN
     !
     WRITE( stdout, '(5X,"Number of +U iterations with fixed ns =",I3)') &
         niter_with_fixed_ns
     WRITE( stdout, '(5X,"Starting occupations:")')
     !
     IF (noncolin) THEN
       CALL write_ns_nc()
     ELSE
       CALL write_ns()
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
  rho(:,2) = rho(:,4)*sin(angle1(1))
  rho(:,3) = rho(:,2)*sin(angle2(1))
  rho(:,4) = rho(:,4)*cos(angle1(1))
  rho(:,2) = rho(:,2)*cos(angle2(1))
  !
  RETURN
  !
END SUBROUTINE nc_magnetization_from_lsda

