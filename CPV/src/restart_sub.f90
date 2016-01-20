!
! Copyright (C) 2002-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
SUBROUTINE from_restart( )
   !
   USE kinds,                 ONLY : DP
   USE control_flags,         ONLY : tbeg, taurdr, tfor, tsdp, iverbosity, &
                                     tsde, tzeroe, tzerop, nbeg, tranp, amprp,&
                                     thdyn, tzeroc, force_pairing, trhor, &
                                     ampre, trane, tpre, dt_old
   USE wavefunctions_module,  ONLY : c0_bgrp, cm_bgrp
   USE electrons_module,      ONLY : occn_info
   USE electrons_base,        ONLY : nspin, iupdwn, nupdwn, f, nbsp, nbsp_bgrp
   USE io_global,             ONLY : ionode, ionode_id, stdout
   USE cell_base,             ONLY : ainv, h, hold, deth, r_to_s, s_to_r, &
                                     velh, at, alat
   USE ions_base,             ONLY : na, nsp, iforce, vel_srt, nat, randpos
   USE time_step,             ONLY : tps, delt
   USE ions_positions,        ONLY : taus, tau0, tausm, taum, vels, fion, fionm, set_velocities
   USE ions_nose,             ONLY : xnhp0, xnhpm
   USE gvect,    ONLY : mill, eigts1, eigts2, eigts3 
   USE printout_base,         ONLY : printout_pos
   USE gvecs,                 ONLY : ngms
   USE gvecw,                 ONLY : ngw
   USE cp_interfaces,         ONLY : phfacs, strucf, prefor, calbec_bgrp, caldbec_bgrp
   USE energies,              ONLY : eself, dft_energy_type
   USE wave_base,             ONLY : rande_base
   USE efield_module,         ONLY : efield_berry_setup,  tefield, &
                                     efield_berry_setup2, tefield2
   USE uspp,                  ONLY : okvan, vkb, nkb, nlcc_any
   USE cp_main_variables,     ONLY : ht0, htm, lambdap, lambda, lambdam, eigr, &
                                     sfac, taub, irb, eigrb, edft, bec_bgrp, dbec, descla
   USE time_step,             ONLY : delt
   USE fft_base,              ONLY : dfftp
   USE matrix_inversion
   !
   IMPLICIT NONE

   INTEGER :: iss

   CALL start_clock( 'from_restart' )
   
   !
   ! ... We are restarting from file recompute ainv
   !
   CALL invmat( 3, h, ainv, deth )
   !
   ! ... Reset total time counter if the run is not strictly 'restart'
   !
   IF ( nbeg < 1 ) tps = 0.D0
   !
   IF ( taurdr ) THEN
      !
      ! ... Input positions read from input file and stored in tau0
      ! ... in readfile, only scaled positions are read
      !
      CALL r_to_s( tau0, taus, na, nsp, ainv )
      !
   END IF
   !
   IF ( ANY( tranp(1:nsp) ) ) THEN
      !
      ! ... Input positions are randomized
      !
      CALL randpos( taus, na, nsp, tranp, amprp, ainv, iforce )
      !
   END IF
   !
   IF ( tzerop .AND. tfor ) THEN
      !
      vel_srt(:,:) = 0.0_dp
      vels(:,:) = 0.0_dp
      CALL set_velocities( tausm, taus, vels, iforce, nat, delt )
      WRITE( stdout, '(" Ionic velocities set to zero")' )
      !
   END IF
   !
   CALL s_to_r( taus,  tau0, na, nsp, h )
   !
   !CALL s_to_r( tausm, taum, na, nsp, h )
   !BS: tausm to taum conversion should use hold in variable cell calculations...
   CALL s_to_r( tausm, taum, na, nsp, hold )
   !
   IF ( tzeroc ) THEN
      !
      hold = h
      velh = 0.D0
      !
      htm      = ht0
      ht0%hvel = 0.D0
      !
   END IF
   !
   fion = 0.D0
   !
   IF( force_pairing ) THEN
      cm_bgrp(:,iupdwn(2):nbsp) = cm_bgrp(:,1:nupdwn(2))
      c0_bgrp(:,iupdwn(2):nbsp) = c0_bgrp(:,1:nupdwn(2))
      lambdap( :, :, 2) =  lambdap( :, :, 1)
      lambda( :, :, 2) =  lambda( :, :, 1)
      lambdam( :, :, 2) = lambdam( :, :, 1)
   END IF 
   !
   IF ( tzeroe ) THEN
      !
      lambdam = lambda
      !
      cm_bgrp = c0_bgrp
      !
      WRITE( stdout, '(" Electronic velocities set to zero")' )
      !
   END IF
   !
   ! ... computes form factors and initializes nl-pseudop. according
   ! ... to starting cell (from ndr or again standard input)
   !
   IF ( okvan .or. nlcc_any ) THEN
      CALL initbox( tau0, alat, at, ainv, taub, irb )
      CALL phbox( taub, iverbosity, eigrb )
   END IF
   !
   CALL phfacs( eigts1, eigts2, eigts3, eigr, mill, taus, dfftp%nr1, dfftp%nr2, dfftp%nr3, nat )
   !
   CALL strucf( sfac, eigts1, eigts2, eigts3, mill, ngms )
   !
   CALL prefor( eigr, vkb )
   !
   CALL formf( .TRUE. , eself )
   !
   IF ( trane ) THEN
      !
      WRITE( stdout, 515 ) ampre
      !
515   FORMAT(   3X,'Initial random displacement of el. coordinates',/ &
                3X,'Amplitude = ',F10.6 )
      !
      CALL rande_base( c0_bgrp, ampre )
      !
      CALL gram_bgrp( vkb, bec_bgrp, nkb, c0_bgrp, ngw )
      !
      IF( force_pairing ) c0_bgrp(:,iupdwn(2):nbsp) = c0_bgrp(:,1:nupdwn(2))
      !
      cm_bgrp = c0_bgrp
      !
   END IF
   !
   CALL calbec_bgrp( 1, nsp, eigr, c0_bgrp, bec_bgrp )
   !
   IF ( tpre     ) CALL caldbec_bgrp( eigr, c0_bgrp, dbec, descla )
   !
   IF ( tefield  ) CALL efield_berry_setup( eigr, tau0 )
   IF ( tefield2 ) CALL efield_berry_setup2( eigr, tau0 )
   !
   CALL plugin_init_ions( tau0 )
   !
   edft%eself = eself
   !
   IF( tzerop .or. tzeroe .or. tzeroc ) THEN
      IF( .not. ( tzerop .and. tzeroe .and. ( tzeroc .or. .not. thdyn ) ) ) THEN
         IF( ionode ) THEN
            WRITE( stdout, * ) 'WARNING setting to ZERO ions, electrons and cell velocities without '
            WRITE( stdout, * ) 'setting to ZERO all velocities could generate meaningles trajectories '
         END IF
      END IF
   END IF
   !
   ! dt_old should be -1.0 here if untouched ...
   !
   if ( dt_old > 0.0d0 ) then
      tausm = taus - (taus-tausm)*delt/dt_old
      xnhpm = xnhp0 - (xnhp0-xnhpm)*delt/dt_old
      WRITE( stdout, '(" tausm & xnhpm were rescaled ")' )
   endif

   CALL stop_clock( 'from_restart' )
   !
   RETURN
   !
END SUBROUTINE from_restart
