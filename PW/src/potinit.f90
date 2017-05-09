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
  USE gvect,                ONLY : ngm, gstart, nl, g, gg
  USE gvecs,                ONLY : doublegrid
  USE control_flags,        ONLY : lscf
  USE scf,                  ONLY : rho, rho_core, rhog_core, &
                                   vltot, v, vrs, kedtau
  USE funct,                ONLY : dft_is_meta
  USE wavefunctions_module, ONLY : psic
  USE ener,                 ONLY : ehart, etxc, vtxc, epaw
  USE ldaU,                 ONLY : lda_plus_u, Hubbard_lmax, eth, &
                                   niter_with_fixed_ns
  USE noncollin_module,     ONLY : noncolin, report
  USE io_files,             ONLY : tmp_dir, prefix, input_drho
  USE spin_orb,             ONLY : domag, lforcet
  USE mp,                   ONLY : mp_sum
  USE mp_bands ,            ONLY : intra_bgrp_comm
  USE io_global,            ONLY : ionode, ionode_id
  USE io_rho_xml,           ONLY : read_scf
  USE xml_io_base,          ONLY : read_rho, check_file_exst
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
  CHARACTER(LEN=256)    :: filename
  !
  CALL start_clock('potinit')
  !
#if defined __HDF5
  filename = TRIM(tmp_dir) // TRIM (prefix) // '.save/charge-density.hdf5'
  exst = check_file_exst( TRIM(filename))
#else 
  ! check for both .dat/ and .xml extensions (compatibility reasons) 
  !
  filename =  TRIM( tmp_dir ) // TRIM( prefix ) // '.save/charge-density.dat'
  exst     =  check_file_exst( TRIM(filename) )
  !
  IF ( .NOT. exst ) THEN
      !
      filename =  TRIM( tmp_dir ) // TRIM( prefix ) // '.save/charge-density.xml'
      exst     =  check_file_exst( TRIM(filename) )
      !
  ENDIF
#endif
  !
  !
  IF ( starting_pot == 'file' .AND. exst ) THEN
     !
     ! ... Cases a) and b): the charge density is read from file
     ! ... this also reads rho%ns if lda+U and rho%bec if PAW
     !
     IF ( .NOT.lforcet ) THEN
        CALL read_scf ( rho, nspin )
     ELSE
        !
        ! ... 'force theorem' calculation of MAE: read rho only from previous
        ! ... lsda calculation, set noncolinear magnetization from angles
        !
        CALL read_rho ( rho%of_r, 2 )
        CALL nc_magnetization_from_lsda ( dfftp%nnr, nspin, rho%of_r )
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
     CALL atomic_rho( rho%of_r, nspin )

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
        CALL read_rho ( v%of_r, 1, input_drho )
        !
        WRITE( UNIT = stdout, &
               FMT = '(/5X,"a scf correction to at. rho is read from",A)' ) &
            TRIM( input_drho )
        !
        rho%of_r = rho%of_r + v%of_r
        !
     END IF
     !
  END IF
  !
  ! ... check the integral of the starting charge
  !
  IF ( nspin == 2 ) THEN
     !
     charge = SUM ( rho%of_r(:,1:nspin) )*omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
     !
  ELSE
     !
     charge = SUM ( rho%of_r(:,1) )*omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
     !
  END IF
  !
  CALL mp_sum(  charge , intra_bgrp_comm )
  !
  IF ( lscf .AND. ABS( charge - nelec ) > ( 1.D-7 * charge ) ) THEN
     !
     IF ( charge > 1.D-8 .AND. nat > 0 ) THEN
        WRITE( stdout, '(/,5X,"starting charge ",F10.5, &
                         & ", renormalised to ",F10.5)') charge, nelec
        rho%of_r = rho%of_r / charge * nelec
     ELSE 
        WRITE( stdout, '(/,5X,"Starting from uniform charge")')
        IF ( nspin == 2 ) THEN
           rho%of_r(:,1:nspin) = nelec / omega / nspin
        ELSE
           rho%of_r(:,1) = nelec / omega
        END IF
     ENDIF
     !
  ELSE IF ( .NOT. lscf .AND. ABS( charge - nelec ) > (1.D-3 * charge ) ) THEN
     !
     CALL errore( 'potinit', 'starting and expected charges differ', 1 )
     !
  END IF
  !
  ! ... bring starting rho to G-space
  !
  DO is = 1, nspin
     !
     psic(:) = rho%of_r(:,is)
     !
     CALL fwfft ('Dense', psic, dfftp)
     !
     rho%of_g(:,is) = psic(nl(:))
     !
  END DO
  !
  if ( dft_is_meta()) then
     ! ... define a starting (TF) guess for rho%kin_r and rho%kin_g
     fact = (3.d0/5.d0)*(3.d0*pi*pi)**(2.0/3.0)
     !
     ! ... for obscure reasons this starting guess doesn't seem much better
     ! ... (and sometimes it is much worse) than starting from zero
     !
     !!! fact = 0.0_dp
     DO is = 1, nspin
        if (starting_pot /= 'file') rho%kin_r(:,is) = fact * abs(rho%of_r(:,is)*nspin)**(5.0/3.0)/nspin
        psic(:) = rho%kin_r(:,is)
        CALL fwfft ('Dense', psic, dfftp)
        rho%kin_g(:,is) = psic(nl(:))
     END DO
     !
  end if
  !
  ! ... plugin contribution to local potential
  !
  CALL plugin_scf_potential(rho,.FALSE.,-1.d0)
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
SUBROUTINE nc_magnetization_from_lsda ( nnr, nspin, rho )
  !-------------
  !
  USE kinds,     ONLY: dp
  USE constants, ONLY: pi
  USE io_global, ONLY: stdout
  USE noncollin_module, ONLY: angle1, angle2
  !
  IMPLICIT NONE
  INTEGER, INTENT (in):: nnr, nspin
  REAL(dp), INTENT (inout):: rho(nnr,nspin)
  !---  
  !  set up noncollinear m_x,y,z from collinear m_z (AlexS) 
  !
  WRITE(stdout,*)
  WRITE(stdout,*) '-----------'
  WRITE(stdout,'("Spin angles Theta, Phi (degree) = ",2f8.4)') &
       angle1(1)/PI*180.d0, angle2(1)/PI*180.d0 
  WRITE(stdout,*) '-----------'
  !
  ! On input, rho(1)=rho_up, rho(2)=rho_down
  ! Set rho(1)=rho_tot, rho(3)=rho_up-rho_down=magnetization
  ! 
  rho(:,3) = rho(:,1)-rho(:,2)
  rho(:,1) = rho(:,1)+rho(:,2)
  !
  ! now set rho(2)=magn*sin(theta)*cos(phi)   x
  !         rho(3)=magn*sin(theta)*sin(phi)   y
  !         rho(4)=magn*cos(theta)            z
  !
  rho(:,4) = rho(:,3)*cos(angle1(1))
  rho(:,2) = rho(:,3)*sin(angle1(1))
  rho(:,3) = rho(:,2)*sin(angle2(1))
  rho(:,2) = rho(:,2)*cos(angle2(1))
  !
  RETURN
  !
END SUBROUTINE nc_magnetization_from_lsda
