!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
#define ONE  (1.D0,0.D0)
#define ZERO (0.D0,0.D0)
!
!----------------------------------------------------------------------------
MODULE extrapolation
  !---------------------------------------------------------------------------
  !! Wfc and rho extrapolation.
  !
  USE kinds, ONLY: DP
  !
  REAL(DP) :: alpha0
  !! the mixing parameters for the extrapolation
  REAL(DP) :: beta0
  !! of the starting potential
  INTEGER :: history
  !! number of old steps available for potential updating
  INTEGER :: pot_order = 0
  !! type of potential updating ( see \(\texttt{update_pot}\) )
  INTEGER :: wfc_order = 0
  !! type of wavefunctions updating ( see \(\texttt{update_pot}\) )
  !
  PRIVATE
  PUBLIC :: pot_order, wfc_order
  PUBLIC :: update_file, update_neb, update_pot, extrapolate_charge
  !
  CONTAINS
!
!----------------------------------------------------------------------------
SUBROUTINE update_file()
  !----------------------------------------------------------------------------
  !! Reads, updates and rewrites the file containing atomic positions at
  !! two previous steps, used by potential and wavefunction extrapolation.  
  !! Requires the number of atoms (nat), current atomic positions (tau).  
  !! Produces the length of history and tau at current and two previous steps
  !! written to file $prefix.update.
  !
  USE io_global, ONLY : ionode
  USE io_files,  ONLY : iunupdate, seqopn
  USE ions_base, ONLY : nat, tau
  !
  IMPLICIT NONE
  !
  REAL(DP), ALLOCATABLE :: tauold(:,:,:)
  LOGICAL :: exst
  !
  IF ( ionode ) THEN
     !
     ALLOCATE( tauold( 3, nat, 3 ) )
     CALL seqopn( iunupdate, 'update', 'FORMATTED', exst ) 
     IF ( .NOT. exst ) THEN
        !
        ! ... file not present, start the procedure
        !
        history = 1
        tauold  = 0.D0
     ELSE
        READ( UNIT = iunupdate, FMT = * ) history
        READ( UNIT = iunupdate, FMT = * ) tauold
        REWIND( UNIT = iunupdate ) 
        !
        ! ... read and save the previous two steps ( three steps are saved )
        !
        tauold(:,:,3) = tauold(:,:,2)
        tauold(:,:,2) = tauold(:,:,1)
        tauold(:,:,1) = tau(:,:)
        !
        ! ... history is updated (a new ionic step has been done)
        !
        history = MIN( 3, ( history + 1 ) )
        !
     END IF
     !
     ! ... history and positions are written on file, file is closed
     !
     WRITE( UNIT = iunupdate, FMT = * ) history
     WRITE( UNIT = iunupdate, FMT = * ) tauold
     CLOSE( UNIT = iunupdate, STATUS = 'KEEP' )
     DEALLOCATE( tauold )
     !
  END IF
  !
END SUBROUTINE update_file
!
!----------------------------------------------------------------------------
SUBROUTINE update_neb( )
  !----------------------------------------------------------------------------
  !! Potential and wavefunction extrapolation for NEB.  
  !! Prepares file with previous steps for usage by \(\texttt{update_pot}\).
  !
  ! ... Must be merged soon with update_file for MD in PWscf.
  !
  USE io_global, ONLY : ionode, ionode_id
  USE io_files,  ONLY : iunupdate, seqopn
  USE mp,        ONLY : mp_bcast
  USE mp_images, ONLY : intra_image_comm
  USE ions_base, ONLY : nat, tau, nsp, ityp
  USE gvect,     ONLY : ngm, g, eigts1, eigts2, eigts3
#if defined(__CUDA)
  USE gvect,     ONLY : eigts1_d, eigts2_d, eigts3_d
#endif
  USE vlocal,    ONLY : strf
  USE cell_base, ONLY : bg
  USE fft_base,  ONLY : dfftp
  !
  IMPLICIT NONE
  !
  REAL(DP), ALLOCATABLE :: tauold(:,:,:)
  LOGICAL :: exst
  !
  ALLOCATE( tauold( 3, nat, 3 ) )
  !
  IF ( ionode ) THEN
     !
     CALL seqopn( iunupdate, 'update', 'FORMATTED', exst )
     IF ( exst ) THEN
        !
        READ( UNIT = iunupdate, FMT = * ) history
        READ( UNIT = iunupdate, FMT = * ) tauold
        !
     ELSE
        !
        ! ... file not present: create one (update_pot needs it)
        !
        history = 0
        tauold  = 0.D0
        WRITE( UNIT = iunupdate, FMT = * ) history
        WRITE( UNIT = iunupdate, FMT = * ) tauold
        !
     END IF
     !
     CLOSE( UNIT = iunupdate, STATUS = 'KEEP' )
     !
   END IF
   !
   CALL mp_bcast( history, ionode_id, intra_image_comm )
   CALL mp_bcast( tauold,  ionode_id, intra_image_comm )
   !
   IF ( history > 0 ) THEN
      !
      ! ... potential and wavefunctions are extrapolated only if
      ! ... we are starting a new self-consistency ( scf on the
      ! ... previous image was achieved )
      !
      IF ( pot_order > 0 ) THEN
         !
         ! ... structure factors of the old positions are computed
         ! ... (needed for the old atomic charge; update_pot will then
         ! ...  overwrite them with structure factors at curret positions)
         !
         CALL struc_fact( nat, tauold(:,:,1), nsp, ityp, ngm, g, bg, &
                          dfftp%nr1, dfftp%nr2, dfftp%nr3, strf,     &
                          eigts1, eigts2, eigts3 )
         ! sync duplicated version
#if defined(__CUDA)
         eigts1_d = eigts1
         eigts2_d = eigts2
         eigts3_d = eigts3
#endif
         !$acc update device(eigts1, eigts2, eigts3) 
         !
      END IF
      !
      CALL update_pot()
      !
   END IF
   !
   IF ( ionode ) THEN
      !
      ! ... save the previous two steps, for usage in next scf
      ! ... ( a total of three ionic steps is saved )
      !
      tauold(:,:,3) = tauold(:,:,2)
      tauold(:,:,2) = tauold(:,:,1)
      tauold(:,:,1) = tau(:,:)
      !
      ! ... update history count (will be used at next step)
      !
      history = MIN( 3, ( history + 1 ) )
      !
      ! ... update history file (must be deleted if scf conv. not reached)
      !
      CALL seqopn( iunupdate, 'update', 'FORMATTED', exst )
      !
      WRITE( UNIT = iunupdate, FMT = * ) history
      WRITE( UNIT = iunupdate, FMT = * ) tauold
      !
      CLOSE( UNIT = iunupdate, STATUS = 'KEEP' )
      !
   END IF
   !
   DEALLOCATE ( tauold )
   !
END SUBROUTINE update_neb
 
!
!----------------------------------------------------------------------------
SUBROUTINE update_pot()
  !----------------------------------------------------------------------------
  !! Update the potential by extrapolating the charge density and extrapolates
  !! the wave-functions.
  !
  !! Charge density extrapolation:
  !
  !! * pot_order=0 \(\rightarrow\) copy the old potential (nothing is done);
  !! * pot_order=1 \(\rightarrow\) subtract old atomic charge density and sum
  !!                               the new if dynamics is done the routine 
  !!                               extrapolates also the difference between
  !!                               the scf charge and the atomic one;
  !! * pot_order=2 \(\rightarrow\) first order extrapolation:
  !!                               \[ \rho(t+dt) = 
  !!                                  2\ \rho(t)-\rho(t-dt); \]
  !! * pot_order=3 \(\rightarrow\) second order extrapolation:
  !!                               \[ \rho(t+dt) = \rho(t) + \alpha_0\ (\rho(t)
  !!                                  -\rho(t-dt)) + \beta_0\ (\rho(t-dt)-
  !!                                   \rho(t-2 dt)). \]
  !
  !! Wave-function extrapolation:
  !
  !! * wfc_order = 0 \(\rightarrow\) nothing is done;
  !! * wfc_order = 2 \(\rightarrow\) first order extrapolation:
  !!                                 \[ |\psi(t+dt)\rangle = 2\ |\psi(t)\rangle-
  !!                                 |\psi(t-dt)\rangle; \]
  !! * wfc_order = 3 \(\rightarrow\) second order extrapolation:
  !!                                 \[ |\psi(t+dt)\rangle = |\psi(t)\rangle
  !!                               + \alpha_0\ ( |\psi(t)\rangle - |\psi(t-dt)\rangle)
  !!                               + \beta_0\ ( |\psi(t-dt)\rangle - |\psi(t-2 dt)\rangle). \]
  !
  !! The \(\alpha_0\) and \(\beta_0\) parameters are calculated in \(\texttt{find_alpha_and_beta}\)
  !! so that \(|\tau'-\tau(t+dt)|\) is minimum. \(\tau'\) and \(\tau(t+dt)\) are respectively
  !! the atomic positions at time t+dt and the extrapolated one:
  !! \[ \tau(t+dt) = \tau(t) + \alpha_0\ ( \tau(t)    - \tau(t-dt)   )
  !!                         + \beta_0\ ( \tau(t-dt) - \tau(t-2 dt) ). \]
  !
  USE io_files,      ONLY : iunupdate, nd_nmbr, seqopn, restart_dir
  USE io_files,      ONLY : wfc_dir, prefix
  USE io_global,     ONLY : ionode, ionode_id
  USE cell_base,     ONLY : bg
  USE ions_base,     ONLY : nat, tau, nsp, ityp
  USE gvect,         ONLY : ngm, g
  USE vlocal,        ONLY : strf
  USE mp,            ONLY : mp_bcast
  USE mp_images,     ONLY : intra_image_comm
  !
  IMPLICIT NONE
  !
  REAL(DP), ALLOCATABLE :: tauold(:,:,:)
  CHARACTER(LEN=256)    :: dirname
  INTEGER               :: rho_extr, wfc_extr
  LOGICAL               :: exists
  !
  !
  CALL start_clock( 'update_pot' )
  !
  ALLOCATE( tauold( 3, nat, 3 ) )
  !
  IF ( ionode ) THEN
     !
     CALL seqopn( iunupdate, 'update', 'FORMATTED', exists )
     !
     IF ( exists ) THEN
        !
        READ( UNIT = iunupdate, FMT = * ) history
        READ( UNIT = iunupdate, FMT = * ) tauold
        !
        ! ... find the best coefficients for the extrapolation
        ! ... of the charge density and of the wavefunctions
        ! ... (see Arias et al. PRB 45, 1538 (1992) )
        !
        CALL find_alpha_and_beta( nat, tau, tauold, alpha0, beta0 )
        !
        CLOSE( UNIT = iunupdate, STATUS = 'KEEP' )
        !
     ELSE
        !
        ! ... default values of extrapolation coefficients
        !
        alpha0 = 1.D0
        beta0  = 0.D0
        history = 0
        tauold = 0.0_dp
        !
        CLOSE( UNIT = iunupdate, STATUS = 'DELETE' )
        !
     END IF
     !
  END IF
  !
  CALL mp_bcast( alpha0, ionode_id, intra_image_comm )
  CALL mp_bcast( beta0,  ionode_id, intra_image_comm )
  CALL mp_bcast( history,ionode_id, intra_image_comm )
  CALL mp_bcast( tauold, ionode_id, intra_image_comm )
  !
  IF ( wfc_order > 0 ) THEN
     !
     ! ... determines the maximum effective order of the extrapolation on the
     ! ... basis of the files that are really available (for wavefunctions)
     !
     IF ( ionode ) THEN
        !
        wfc_extr = MIN( 1, history, wfc_order )
        !
        INQUIRE( FILE = TRIM( wfc_dir ) // &
            & TRIM( prefix ) // '.oldwfc' // nd_nmbr, EXIST = exists )
        !
        IF ( exists ) THEN
           !
           wfc_extr = MIN( 2, history, wfc_order  )
           !
           INQUIRE( FILE = TRIM( wfc_dir ) // &
               & TRIM( prefix ) // '.old2wfc' // nd_nmbr , EXIST = exists )
           !
           IF ( exists ) wfc_extr = MIN( 3, history, wfc_order )
           !
        END IF
        !
     END IF
     !
     CALL mp_bcast( wfc_extr, ionode_id, intra_image_comm )
     !
     !
     ! ... save tau(t+dt), replace with tau(t)
     ! ... extrapolate_wfcs needs tau(t) to evaluate S(t)
     ! ... note that structure factors have not yet been updated
     !
     tauold (:,:,2) = tau (:,:)
     tau (:,:) = tauold (:,:,1)
     !
     CALL extrapolate_wfcs( wfc_extr )
     !
     ! ... restore tau(t+dt)
     !
     tau (:,:) = tauold (:,:,2)
     !
  END IF
  !
  DEALLOCATE( tauold )
  !
  ! ... determines the maximum effective order of the extrapolation on the
  ! ... basis of the files that are really available (for the charge density)
  !
  rho_extr = MIN( 1, history, pot_order )
  !
  IF ( pot_order > 1 .AND. ionode ) THEN
     !
     ! FIXME: won't work for hdf5
     dirname =  restart_dir () 
     INQUIRE( FILE = TRIM( dirname ) // 'charge-old.dat', &
          EXIST = exists )
     !
     IF ( exists ) THEN
        !
        rho_extr = MIN( 2, history, pot_order )
        !
        INQUIRE( FILE = TRIM( dirname ) //'charge-old2.dat', &
             EXIST = exists )
        !
        IF ( exists ) rho_extr = MIN( 3, history, pot_order )
        !
     END IF
     !
  END IF
  !
  CALL mp_bcast( rho_extr, ionode_id, intra_image_comm )
  !
  CALL extrapolate_charge( dirname, rho_extr )
  !
  CALL stop_clock( 'update_pot' )
  !
  RETURN
  !
END SUBROUTINE update_pot
!
!----------------------------------------------------------------------------
SUBROUTINE extrapolate_charge( dirname, rho_extr )
  !----------------------------------------------------------------------------
  !! Charge density extrapolation.
  !
  USE io_global,            ONLY : stdout
  USE cell_base,            ONLY : omega, bg
  USE ions_base,            ONLY : nat, tau, nsp, ityp
  USE fft_base,             ONLY : dfftp, dffts
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE control_flags,        ONLY : gamma_only
  USE gvect,                ONLY : ngm, g, gg, gstart, eigts1, eigts2, eigts3, &
                                   mill, ig_l2g
#if defined(__CUDA)
  USE gvect,                ONLY : eigts1_d, eigts2_d, eigts3_d
#endif
  USE lsda_mod,             ONLY : lsda, nspin
  USE scf,                  ONLY : rho, rho_core, rhog_core, v
  USE ldaU,                 ONLY : eth
  USE ener,                 ONLY : ehart, etxc, vtxc, epaw
  USE extfield,             ONLY : etotefield
  USE cellmd,               ONLY : lmovecell, omega_old
  USE vlocal,               ONLY : strf
  USE noncollin_module,     ONLY : noncolin
  USE klist,                ONLY : nelec
  USE paw_variables,        ONLY : okpaw, ddd_paw
  USE paw_onecenter,        ONLY : PAW_potential
  USE mp_pools,             ONLY : my_pool_id
  USE mp_bands,             ONLY : my_bgrp_id, root_bgrp_id, root_bgrp, &
                                   intra_bgrp_comm, inter_bgrp_comm, nbgrp
  USE mp,                   ONLY : mp_bcast
  USE io_base,              ONLY : write_rhog, read_rhog
  USE fft_rho,              ONLY : rho_g2r, rho_r2g
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: rho_extr
  CHARACTER(LEN=*), INTENT(IN) :: dirname
  !
  ! ... local variables
  !
  COMPLEX (DP), ALLOCATABLE :: work(:,:), work1(:,:)
    ! work  is the difference between rho and atomic rho at time t
    ! work1 is the same thing at time t-dt
  REAL(DP) :: dum(3) = [0_dp, 0_dp, 0_dp]
  REAL(DP) :: charge
  !
  INTEGER :: is
  !
  IF ( rho_extr < 1 ) THEN
     !
     ! ... calculate structure factors for the new positions
     !
     IF ( lmovecell ) CALL scale_h()
     !
     CALL struc_fact( nat, tau, nsp, ityp, ngm, g, bg, &
          dfftp%nr1, dfftp%nr2, dfftp%nr3, strf, eigts1, eigts2, eigts3 )
     ! sync duplicated version
#if defined(__CUDA)
     eigts1_d = eigts1
     eigts2_d = eigts2
     eigts3_d = eigts3
#endif
     !$acc update device(eigts1, eigts2, eigts3) 
     !
     ! ... new charge density from extrapolated wfcs
     !
     IF ( rho_extr < 0 ) THEN
        !
        CALL sum_band ()
        !
        WRITE( UNIT = stdout, FMT = '(5X, &
             & "charge density from extrapolated wavefunctions")' )
     ELSE
        !
        WRITE( UNIT = stdout, FMT = '(5X, &
             & "charge density from previous step")' )
        !
     END IF
     !
     CALL set_rhoc()
     !
  ELSE
     ! 
     ALLOCATE( work( dfftp%ngm, 1 ) )
     !
     ! ... in the lsda case the magnetization will follow rigidly the density
     ! ... keeping fixed the value of zeta = mag / rho_tot.
     ! ... zeta is set here and put in rho%of_r(:,2) while rho%of_r(:,1) 
     ! ... contains the total valence charge
     ! FIXME: half extrapolation is done in G-space, half in real space
     !
     IF ( lsda .or. noncolin ) CALL rho2zeta( rho%of_r, rho_core, dfftp%nnr, nspin, 1 )
     !
     ! ... subtract the old atomic charge density
     !
     CALL atomic_rho_g( work, 1 )
     !
     rho%of_g(:,1) = rho%of_g(:,1) - work(:,1)
     !
     IF ( lmovecell ) rho%of_g(:,1) = rho%of_g(:,1) * omega_old
     !
     ! ... extrapolate the difference between the atomic charge and
     ! ... the self-consistent one
     !
     IF ( rho_extr == 1 ) THEN
        !
        ! ... if rho_extr = 1  update the potential subtracting to the charge
        ! ...                  density the "old" atomic charge and summing the
        ! ...                  new one
        !
        WRITE( UNIT = stdout, FMT = '(5X, &
             & "NEW-OLD atomic charge density approx. for the potential")' )
        !
        ! ... no need to save old charge density if no extrapolation is done
        !
        IF ( pot_order > 1 ) THEN
           ! Write G-space density
           IF ( my_pool_id == 0 .AND. my_bgrp_id == root_bgrp_id ) &
                CALL write_rhog( TRIM(dirname) // "charge-old", &
                root_bgrp, intra_bgrp_comm, dum, dum, dum, & 
                gamma_only, mill, ig_l2g, rho%of_g(:,1:1) )
           !
        END IF
        !
     ELSE IF ( rho_extr == 2 ) THEN
        !
        WRITE( UNIT = stdout, &
               FMT = '(5X,"first order charge density extrapolation")' )
        !
        ! ...   oldrho  ->  work
        !
        ! Read G-space density
        IF ( my_pool_id == 0 .AND. my_bgrp_id == root_bgrp_id ) &
             CALL read_rhog( TRIM(dirname) // "charge-old", &
             root_bgrp, intra_bgrp_comm, ig_l2g, 1, work(:,1:1) )
        IF( nbgrp > 1 ) CALL mp_bcast( work, root_bgrp_id, inter_bgrp_comm )
        !
        ! ...   rho%of_r   ->  oldrho
        ! ...   work  ->  oldrho2
        ! FIXME: rename files instead of rewriting them
        !
        IF ( my_pool_id == 0 .AND. my_bgrp_id == root_bgrp_id ) &
             CALL write_rhog( TRIM(dirname) // "charge-old", &
             root_bgrp, intra_bgrp_comm, dum, dum, dum, & 
             gamma_only, mill, ig_l2g, rho%of_g(:,1:1) )
        IF ( my_pool_id == 0 .AND. my_bgrp_id == root_bgrp_id ) &
             CALL write_rhog( TRIM(dirname) // "charge-old2", &
             root_bgrp, intra_bgrp_comm, dum, dum, dum, & 
             gamma_only, mill, ig_l2g, work(:,1:1) )
        !
        ! ... extrapolation
        !
        rho%of_g(:,1) = 2.D0*rho%of_g(:,1) - work(:,1)
        !
     ELSE IF ( rho_extr == 3 ) THEN
        !
        WRITE( UNIT = stdout, &
               FMT = '(5X,"second order charge density extrapolation")' )
        !
        ALLOCATE( work1( dfftp%ngm, 1 ) )
        !
        work1 = 0.D0
        !
        ! ...   oldrho2  ->  work1
        ! ...   oldrho   ->  work
        !
        IF ( my_pool_id == 0 .AND. my_bgrp_id == root_bgrp_id ) &
             CALL read_rhog( TRIM(dirname) // "charge-old2", &
             root_bgrp, intra_bgrp_comm, ig_l2g, 1, work1(:,1:1) )
        IF( nbgrp > 1 ) CALL mp_bcast( work1,root_bgrp_id, inter_bgrp_comm )
        IF ( my_pool_id == 0 .AND. my_bgrp_id == root_bgrp_id ) &
             CALL read_rhog( TRIM(dirname) // "charge-old", &
             root_bgrp, intra_bgrp_comm, ig_l2g, 1, work(:,1:1) )
        IF( nbgrp > 1 ) CALL mp_bcast( work, root_bgrp_id, inter_bgrp_comm )
        !
        ! ...   rho%of_r   ->  oldrho
        ! ...   work  ->  oldrho2
        !
        IF ( my_pool_id == 0 .AND. my_bgrp_id == root_bgrp_id ) &
             CALL write_rhog( TRIM(dirname) // "charge-old", &
             root_bgrp, intra_bgrp_comm, dum, dum, dum, & 
             gamma_only, mill, ig_l2g, rho%of_g(:,1:1) )
        IF ( my_pool_id == 0 .AND. my_bgrp_id == root_bgrp_id ) &
             CALL write_rhog( TRIM(dirname) // "charge-old2", &
             root_bgrp, intra_bgrp_comm, dum, dum, dum, & 
             gamma_only, mill, ig_l2g, work(:,1:1) )
        !
        rho%of_g(:,1) = rho%of_g(:,1) + alpha0*( rho%of_g(:,1) - work(:,1) ) + &
                               beta0*( work(:,1) - work1(:,1) )
        !
        DEALLOCATE( work1 )
        !
     END IF
     !
     IF ( lmovecell ) rho%of_g(:,1) = rho%of_g(:,1) / omega
     !
     ! ... calculate structure factors for the new positions
     !
     IF ( lmovecell ) CALL scale_h()
     !
     CALL struc_fact( nat, tau, nsp, ityp, ngm, g, bg, &
          dfftp%nr1, dfftp%nr2, dfftp%nr3, strf, eigts1, eigts2, eigts3 )
     ! sync duplicated version
#if defined(__CUDA)
     eigts1_d = eigts1
     eigts2_d = eigts2
     eigts3_d = eigts3
#endif
     !$acc update device(eigts1, eigts2, eigts3) 
     !
     CALL set_rhoc()
     !
     ! ... add atomic charges in the new positions
     !
     CALL atomic_rho_g( work, 1 )
     !
     rho%of_g(:,1) = rho%of_g(:,1) + work(:,1)
     !
     ! ... bring extrapolated rho to real space
     !
     CALL rho_g2r (dfftp, rho%of_g(:,1:1), rho%of_r(:,1:1))
     !
     ! ... reset up and down charge densities in the LSDA case
     !
     IF ( lsda .or. noncolin ) CALL rho2zeta( rho%of_r, rho_core, dfftp%nnr, nspin, -1 )
     !
     DEALLOCATE( work )
     !
  END IF
  !
  CALL rho_r2g (dfftp, rho%of_r, rho%of_g )
  !
  CALL v_of_rho( rho, rho_core, rhog_core, &
                 ehart, etxc, vtxc, eth, etotefield, charge, v )
  IF (okpaw) CALL PAW_potential(rho%bec, ddd_PAW, epaw)
  !
  IF ( ABS( charge - nelec ) / charge > 1.D-7 ) THEN
     !
     WRITE( stdout, &
            '(5X,"extrapolated charge ",F10.5,", renormalised to ",F10.5)') &
         charge, nelec
     !
     rho%of_r  = rho%of_r  / charge*nelec
     rho%of_g = rho%of_g / charge*nelec
     !
  END IF
  !
  RETURN
  !
END SUBROUTINE extrapolate_charge
!
!-----------------------------------------------------------------------
SUBROUTINE extrapolate_wfcs( wfc_extr )
  !-----------------------------------------------------------------------
  !! This routine extrapolate the wfc's after a "parallel alignment"
  !! of the basis of the \(t-dt\) and \(t\) time steps, according to a recipe
  !! by Mead, Rev. Mod. Phys., vol 64, pag. 51 (1992), Eqs. (3.20) and (3.29).
  !
  USE io_global,            ONLY : stdout
  USE klist,                ONLY : nks, ngk, xk, igk_k
  USE lsda_mod,             ONLY : lsda, current_spin, isk
  USE wvfct,                ONLY : nbnd, npwx
  USE ions_base,            ONLY : nat, tau
  USE io_files,             ONLY : nwordwfc, iunwfc, iunoldwfc, &
                                   iunoldwfc2, diropn
  USE buffers,              ONLY : get_buffer, save_buffer
  USE uspp,                 ONLY : nkb, vkb, okvan
  USE wavefunctions,        ONLY : evc
  USE noncollin_module,     ONLY : noncolin, npol
  USE control_flags,        ONLY : gamma_only
  USE becmod,               ONLY : allocate_bec_type, deallocate_bec_type, &
                                   bec_type, becp, calbec
  USE mp_images,            ONLY : intra_image_comm
  USE mp,                   ONLY : mp_barrier
  USE mp_bands,             ONLY : use_bgrp_in_hpsi
  USE wavefunctions_gpum,   ONLY : using_evc
  USE becmod_subs_gpum,     ONLY : using_becp_auto
  USE uspp_init,            ONLY : init_us_2
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: wfc_extr
  !
  ! ... local variables
  !
  INTEGER :: npw, ik, zero_ew, lwork, info
    ! do-loop variables
    ! counter on k-points
    ! number of zero 'eigenvalues' of the s_m matrix
    ! used by singular value decomposition (ZGESVD)
    ! flag returned by ZGESVD
  COMPLEX(DP), ALLOCATABLE :: sp_m(:,:), u_m(:,:), w_m(:,:), work(:)
    ! the overlap matrix s^+ (eq. 3.24)
    ! left unitary matrix in the SVD of sp_m
    ! right unitary matrix in the SVD of sp_m
    ! workspace for ZGESVD
  COMPLEX(DP), ALLOCATABLE :: evcold(:,:), aux(:,:)
    ! wavefunctions at previous iteration + workspace
  REAL(DP), ALLOCATABLE :: ew(:), rwork(:), rp_m(:,:)
    ! the eigenvalues of s_m
    ! workspace for ZGESVD
    ! real version of sp_m
  LOGICAL :: exst
  LOGICAL :: save_flag
  !
  CALL mp_barrier( intra_image_comm ) ! debug
  !
  save_flag = use_bgrp_in_hpsi ; use_bgrp_in_hpsi=.false.
  !
  CALL using_evc(0)
  !
  IF ( wfc_extr == 1 ) THEN
     !
     CALL diropn( iunoldwfc, 'oldwfc', 2*nwordwfc, exst )
     !
     DO ik = 1, nks
        !
        ! ... "now"  -> "old"
        !
        IF ( nks > 1 ) CALL get_buffer( evc, nwordwfc, iunwfc, ik )
        IF ( nks > 1 ) CALL using_evc(2)
        CALL davcio( evc, 2*nwordwfc, iunoldwfc, ik, +1 )
        !
     END DO
     !
     CLOSE( UNIT = iunoldwfc, STATUS = 'KEEP' )
     !
  ELSE 
     !
     CALL diropn( iunoldwfc, 'oldwfc', 2*nwordwfc, exst )
     IF ( wfc_extr > 2 .OR. wfc_order > 2 ) &
        CALL diropn( iunoldwfc2, 'old2wfc', 2*nwordwfc, exst )
     !
     IF ( wfc_extr == 2 ) THEN
        !
        WRITE( stdout, '(/5X,"first order wave-functions extrapolation")' )
        !
     ELSE
        !
        WRITE( stdout, '(/5X,"second order wave-functions extrapolation")' )
        !
     END IF
     !
     ALLOCATE( evcold( npwx*npol, nbnd ), aux( npwx*npol, nbnd ) )
     ALLOCATE( sp_m( nbnd, nbnd ), u_m( nbnd, nbnd ), w_m( nbnd, nbnd ), ew( nbnd ) )
     CALL allocate_bec_type ( nkb, nbnd, becp ) 
     !
     IF( SIZE( aux ) /= SIZE( evc ) ) &
        CALL errore('extrapolate_wfcs ', ' aux wrong size ', ABS( SIZE( aux ) - SIZE( evc ) ) ) 
     !
     ! query workspace
     !
     lwork = 5*nbnd
     !
     ALLOCATE( rwork( lwork ) )
     ALLOCATE( work( lwork ) )
     lwork = -1
     CALL ZGESVD( 'A', 'A', nbnd, nbnd, sp_m, nbnd, ew, u_m, &
                     nbnd, w_m, nbnd, work, lwork, rwork, info )
     !
     lwork = INT(work( 1 )) + 1
     !
     IF( lwork > SIZE( work ) ) THEN
        DEALLOCATE( work )
        ALLOCATE( work( lwork ) )
     END IF
     !
     zero_ew = 0
     !
     DO ik = 1, nks
        !
        ! ... read wavefcts as (t-dt), replace with wavefcts at (t)
        !
        CALL davcio( evcold, 2*nwordwfc, iunoldwfc, ik, -1 )
        IF ( nks > 1 ) CALL get_buffer( evc, nwordwfc, iunwfc, ik )
        IF ( nks > 1 ) CALL using_evc(2)
        CALL davcio(    evc, 2*nwordwfc, iunoldwfc, ik, +1 )
        !
        npw = ngk (ik)
        IF ( okvan ) THEN
           !
           ! ... Ultrasoft PP: calculate overlap matrix
           ! ... Required by s_psi:
           ! ... nonlocal pseudopotential projectors |beta>, <psi|beta>
           !
           IF ( nkb > 0 ) CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb )
           CALL using_becp_auto(2)
           CALL calbec( npw, vkb, evc, becp )
           CALL s_psi ( npwx, npw, nbnd, evc, aux )
           !
        ELSE
           !
           ! ... Norm-Conserving  PP: no overlap matrix
           !
           aux = evc
           !
        END IF
        !
        ! ... construct s^+_m = <psi(t)|S|psi(t-dt)>
        !
        IF ( gamma_only ) THEN
           ALLOCATE( rp_m ( nbnd, nbnd ) )
           CALL calbec ( npw, aux, evcold, rp_m )
           sp_m(:,:) = CMPLX(rp_m(:,:),0.0_DP,kind=DP)
           DEALLOCATE( rp_m )
        ELSE IF ( noncolin) THEN
           CALL calbec ( npwx*npol, aux, evcold, sp_m )
        ELSE
           CALL calbec ( npw, aux, evcold, sp_m )
        END IF
        !
        ! ... the unitary matrix [sp_m*s_m]^(-1/2)*sp_m (eq. 3.29) by means the
        ! ... singular value decomposition (SVD) of  sp_m = u_m*diag(ew)*w_m
        ! ... becomes u_m * w_m
        !
        CALL ZGESVD( 'A', 'A', nbnd, nbnd, sp_m, nbnd, ew, u_m, &
                     nbnd, w_m, nbnd, work, lwork, rwork, info )
        !
        ! ... check on eigenvalues
        !
        zero_ew = COUNT( ew(:) < 0.1D0 )
        !
        ! ... use sp_m to store u_m * w_m
        !
        CALL ZGEMM( 'N', 'N', nbnd, nbnd, nbnd, ONE, &
                    u_m, nbnd, w_m, nbnd, ZERO, sp_m, nbnd )
        !
        ! ... now use aux as workspace to calculate "aligned" wavefcts:
        !
        ! ... aux_i = sum_j evcold_j*s_m_ji (eq.3.21)
        !
        CALL ZGEMM( 'N', 'C', npw, nbnd, nbnd, ONE, &
                    evcold, npwx, sp_m, nbnd, ZERO, aux, npwx )
        !
        ! ... alpha0 and beta0 are calculated in "update_pot"
        ! ... for first-order interpolation, alpha0=1, beta0=0
        !
        CALL using_evc(1)
        IF ( wfc_extr == 3 ) THEN
           evc = ( 1.0_dp + alpha0 ) * evc + ( beta0 - alpha0 ) * aux
        ELSE
           evc = 2.0_dp * evc - aux
        END IF
        !
        IF ( wfc_order > 2 ) THEN
           !
           ! ... second-order interpolation:
           ! ... read wavefcts at (t-2dt), save aligned wavefcts at (t-dt)
           !
           IF ( wfc_extr == 3 ) &
               CALL davcio( evcold, 2*nwordwfc, iunoldwfc2, ik, -1 )
           !
           CALL davcio(    aux, 2*nwordwfc, iunoldwfc2, ik, +1 )
           !
           IF ( wfc_extr ==3 ) THEN
              !
              ! ... align wfcs at (t-2dt), add to interpolation formula 
              !
              CALL ZGEMM( 'N', 'C', npw, nbnd, nbnd, ONE, &
                          evcold, npwx, sp_m, nbnd, ZERO, aux, npwx )
              !
              evc = evc - beta0 * aux
              !
           END IF
           !
        END IF
        !
        ! ... save interpolated wavefunctions to file iunwfc
        !
        ! CALL using_evc(0) aldready done above
        IF ( nks > 1 ) CALL save_buffer( evc, nwordwfc, iunwfc, ik )
        !
     END DO
     !
     IF ( zero_ew > 0 ) &
        WRITE( stdout, '( 5X,"Message from extrapolate_wfcs: ",/,  &
                        & 5X,"the matrix <psi(t-dt)|psi(t)> has ", &
                        & I2," small (< 0.1) eigenvalues")' ) zero_ew
     !
     DEALLOCATE( u_m, w_m, ew, aux, evcold, sp_m )
     DEALLOCATE( work, rwork )
     CALL deallocate_bec_type ( becp ) 
     CALL using_becp_auto(2)
     !
     CLOSE( UNIT = iunoldwfc, STATUS = 'KEEP' )
     IF ( wfc_extr > 2 .OR. wfc_order > 2 ) &
        CLOSE( UNIT = iunoldwfc2, STATUS = 'KEEP' )
     !
  END IF
  !
  use_bgrp_in_hpsi = save_flag

  CALL mp_barrier( intra_image_comm ) ! debug
  !
  RETURN
  !
END SUBROUTINE extrapolate_wfcs
!
! ... this routine is used also by compute_scf (NEB) and compute_fes_grads
!
!----------------------------------------------------------------------------
SUBROUTINE find_alpha_and_beta( nat, tau, tauold, alpha0, beta0 )
  !----------------------------------------------------------------------------
  !! This routine finds the best coefficients \(\alpha_0\) and \(\beta_0\) so
  !! that \(|\tau(t+dt)-\tau'|\) is minimum, with:
  !! \[ \tau' = \tau(t) + \alpha_0\cdot( \tau(t)    - \tau(t-dt)   )
  !!                    + \beta_0 \cdot( \tau(t-dt) - \tau(t-2dt) ). \]
  !
  USE constants,     ONLY : eps16
  USE io_global,     ONLY : stdout
  !
  IMPLICIT NONE
  !
  INTEGER  :: nat
  REAL(DP) :: tau(3,nat)
  REAL(DP) :: tauold(3,nat,3)
  REAL(DP) :: alpha0
  REAL(DP) :: beta0
  !
  ! ... local variables
  !
  INTEGER :: na, ipol
  REAL(DP) :: a11, a12, a21, a22, b1, b2, c, det
  !
  !
  IF ( history <= 2 ) RETURN
  !
  ! ... solution of the linear system
  !
  a11 = 0.D0
  a12 = 0.D0
  a21 = 0.D0
  a22 = 0.D0
  b1  = 0.D0
  b2  = 0.D0
  c   = 0.D0
  !
  DO na = 1, nat
     !
     DO ipol = 1, 3
        !
        a11 = a11 + ( tauold(ipol,na,1) - tauold(ipol,na,2) )**2
        !
        a12 = a12 + ( tauold(ipol,na,1) - tauold(ipol,na,2) ) * &
                    ( tauold(ipol,na,2) - tauold(ipol,na,3) )
        !
        a22 = a22 + ( tauold(ipol,na,2) - tauold(ipol,na,3) )**2
        !
        b1 = b1 - ( tauold(ipol,na,1) - tau(ipol,na) ) * &
                  ( tauold(ipol,na,1) - tauold(ipol,na,2) )
        !
        b2 = b2 - ( tauold(ipol,na,1) - tau(ipol,na) ) * &
                  ( tauold(ipol,na,2) - tauold(ipol,na,3) )
        !
        c = c + ( tauold(ipol,na,1) - tau(ipol,na) )**2
        !
     END DO
     !
  END DO
  !
  a21 = a12
  !
  det = a11 * a22 - a12 * a21
  !
  IF ( det < - eps16 ) THEN
     !
     alpha0 = 0.D0
     beta0  = 0.D0
     !
     WRITE( UNIT = stdout, &
            FMT = '(5X,"WARNING: in find_alpha_and_beta  det = ",F10.6)' ) det
     !
  END IF
  !
  ! ... case det > 0:  a well defined minimum exists
  !
  IF ( det > eps16 ) THEN
     !
     alpha0 = ( b1 * a22 - b2 * a12 ) / det
     beta0  = ( a11 * b2 - a21 * b1 ) / det
     !
  ELSE
     !
     ! ... case det = 0 : the two increments are linearly dependent,
     ! ...                chose solution with alpha = b1 / a11 and beta = 0
     ! ...                ( discard oldest configuration )
     !
     alpha0 = 0.D0
     beta0  = 0.D0
     !
     IF ( a11 /= 0.D0 ) alpha0 = b1 / a11
     !
  END IF
  !
  RETURN
  !
END SUBROUTINE find_alpha_and_beta
  !
END MODULE extrapolation
