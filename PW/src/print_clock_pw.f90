!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE print_clock_pw()
   !---------------------------------------------------------------------------
   !
   ! ... this routine prints out the clocks at the end of the run
   ! ... it tries to construct the calling tree of the program.
   !
   USE io_global,          ONLY : stdout
   USE control_flags,      ONLY : isolve, iverbosity, gamma_only, lxdm, &
        ts_vdw, ldftd3, llondon
   USE paw_variables,      ONLY : okpaw
   USE uspp,               ONLY : okvan
   USE realus,             ONLY : real_space
   USE noncollin_module,   ONLY : noncolin
   USE ldaU,               ONLY : lda_plus_u, lda_plus_u_kind, is_hubbard_back
   USE xc_lib,             ONLY : xclib_dft_is
   USE bp,                 ONLY : lelfield
   USE rism_module,        ONLY : rism_print_clock
   !
#if defined (__ENVIRON)
   USE plugin_flags,        ONLY : use_environ
   USE environ_base_module, ONLY : print_environ_clocks
#endif
#if defined (__OSCDFT)
   USE plugin_flags,        ONLY : use_oscdft
   USE oscdft_base,         ONLY : print_oscdft_clocks, oscdft_ctx
#endif
   !
   IMPLICIT NONE
   !
   !
   WRITE( stdout, * )
   !
   CALL print_clock( 'init_run' )
   CALL print_clock( 'electrons' )
   CALL print_clock( 'update_pot' )
   CALL print_clock( 'forces' )
   IF (ldftd3)  CALL print_clock('force_dftd3')
   IF (llondon) CALL print_clock('force_london')
   CALL print_clock( 'stress' )
   IF (ldftd3)  CALL print_clock('stres_dftd3')
   IF (llondon) CALL print_clock('stres_london')
   !
   WRITE( stdout, '(/5x,"Called by init_run:")' )
   CALL print_clock( 'aceinit0' )
   CALL print_clock( 'wfcinit' )
   IF ( iverbosity > 0 ) THEN
      CALL print_clock( 'wfcinit:atomic' )
      CALL print_clock( 'wfcinit:wfcrot' )
   END IF
   CALL print_clock( 'potinit' )
   CALL print_clock( 'realus' )
   IF ( iverbosity > 0 ) THEN
      CALL print_clock( 'realus:boxes' )
      CALL print_clock( 'realus:spher' )
      CALL print_clock( 'realus:tabp' )
   END IF
   CALL print_clock( 'hinit0' )
   IF (lxdm)   CALL print_clock('init_xdm')
   IF (ts_vdw) CALL print_clock('tsvdw_pair')
   !
   WRITE( stdout, '(/5x,"Called by electrons:")' )
   CALL print_clock( 'c_bands' )
   CALL print_clock( 'sum_band' )
   CALL print_clock( 'v_of_rho' )
   IF ( iverbosity > 0 ) THEN
      CALL print_clock( 'v_h' )
      CALL print_clock( 'v_xc' )
      CALL print_clock( 'v_xc_meta' )
   END IF
   CALL print_clock( 'newd' )
   CALL print_clock( 'PAW_pot')
   CALL print_clock( 'mix_rho' )

   CALL print_clock( 'vdW_energy' )
   CALL print_clock( 'vdW_ffts' )
   CALL print_clock( 'vdW_v' )
   CALL print_clock( 'vdW_kernel' ) 
   
   IF (lxdm) THEN
      CALL print_clock('energy_xdm')
      CALL print_clock('exdm:environ')
      CALL print_clock('exdm:paw_charge')
      CALL print_clock('exdm:rho')
   ELSE IF (ts_vdw) THEN
      CALL print_clock('tsvdw_rhotot')
      CALL print_clock('tsvdw_screen')
      CALL print_clock('tsvdw_veff')
      CALL print_clock('tsvdw_dveff')
      CALL print_clock('tsvdw_energy')
   ELSE IF (ldftd3) THEN
      CALL print_clock('energy_dftd3')
   ELSE IF (llondon) THEN
      CALL print_clock('energy_london')
   END IF
   !
   WRITE( stdout, '(/5x,"Called by c_bands:")' )
   CALL print_clock( 'init_us_2' )
   CALL print_clock( 'init_us_2:cpu' )
   CALL print_clock( 'init_us_2:gpu' )
   IF ( isolve == 0 ) THEN
      CALL print_clock( 'regterg' )    ; CALL print_clock( 'cegterg' )
   ELSE  IF (isolve == 1) THEN
      CALL print_clock( 'rcgdiagg' )   ; CALL print_clock( 'ccgdiagg' )
      CALL print_clock( 'wfcrot' )
   ELSE  IF (isolve == 3) THEN
      CALL print_clock( 'paro_gamma' ) ; CALL print_clock( 'paro_k' )
   ELSE IF ( isolve == 4 ) THEN
      CALL print_clock( 'rrmmdiagg' )  ; CALL print_clock( 'crmmdiagg' )
      CALL print_clock( 'wfcrot' )
      CALL print_clock( 'gsorth' )
   ENDIF
   !
   IF ( iverbosity > 0)  THEN
      WRITE( stdout, '(/5x,"Called by sum_band:")' )
      CALL print_clock( 'sum_band:weights' )
      CALL print_clock( 'sum_band:loop' )
      CALL print_clock( 'sum_band:buffer' )
      CALL print_clock( 'sum_band:init_us_2' )
      CALL print_clock( 'sum_band:calbec' )
      CALL print_clock( 'sum_band:becsum' )
      CALL print_clock( 'addusdens' )
      CALL print_clock( 'addusd:skk' )
      CALL print_clock( 'addusd:dgemm' )
      CALL print_clock( 'addusd:qvan2' )
   ENDIF
   !
   IF ( isolve == 0 ) THEN
      WRITE( stdout, '(/5x,"Called by *egterg:")' )
      IF ( gamma_only ) THEN
         CALL print_clock( 'rdiaghg' )
         IF ( iverbosity > 0 )  THEN
            CALL print_clock( 'regterg:overlap' )
            CALL print_clock( 'regterg:update' )
            CALL print_clock( 'regterg:last' )
            CALL print_clock( 'rdiaghg:choldc' )
            CALL print_clock( 'rdiaghg:inversion' )
            CALL print_clock( 'rdiaghg:paragemm' )
         ENDIF
      ELSE
         CALL print_clock( 'cdiaghg' )
         IF ( iverbosity > 0 )  THEN
            CALL print_clock( 'cegterg:overlap' )
            CALL print_clock( 'cegterg:update' )
            CALL print_clock( 'cegterg:last' )
            CALL print_clock( 'cdiaghg:choldc' )
            CALL print_clock( 'cdiaghg:inversion' )
            CALL print_clock( 'cdiaghg:paragemm' )
         END IF
      END IF
   ELSE IF ( isolve == 1 ) THEN
      WRITE( stdout, '(/5x,"Called by *cgdiagg:")' )
   ELSE IF ( isolve == 3 ) THEN
      WRITE( stdout, '(/5x,"Called by paro_*:")' )
      IF ( iverbosity > 0 )  THEN
         CALL print_clock( 'paro:init' )
         CALL print_clock( 'paro:pack' )
         CALL print_clock( 'paro:zero' )
         CALL print_clock( 'paro:mp_bar' )
         CALL print_clock( 'paro:mp_sum' )
         CALL print_clock( 'pcg' )
         CALL print_clock( 'pcg:hs_1psi' )
         CALL print_clock( 'pcg:ortho' )
         CALL print_clock( 'pcg:move1' )
         CALL print_clock( 'pcg:move2' )

         CALL print_clock( 'rotHSw' )
         CALL print_clock( 'rotHSw:move' )
         CALL print_clock( 'rotHSw:hc' )
         CALL print_clock( 'rotHSw:diag' )
         CALL print_clock( 'rotHSw:evc' )
         CALL print_clock( 'rotHSw:hc:b0' ) ; 
         CALL print_clock( 'rotHSw:hc:s1' ) ; call print_clock('rotHSw:hc:comp')
         CALL print_clock( 'rotHSw:hc:b1' ) ; 
         CALL print_clock( 'rotHSw:hc:s2' ) ; 
         CALL print_clock( 'rotHSw:hc:s3' ) ; call print_clock('rotHSw:hc:rs')
         CALL print_clock( 'rotHSw:hc:b2' ) ; call print_clock('rotHSw:hc:sy')
         CALL print_clock( 'rotHSw:hc:s4' ) ; CALL print_clock('rotHSw:hc:b3' ) 
         CALL print_clock( 'rotHSw:ev:b0' ) ; 
         CALL print_clock( 'rotHSw:ev:b3' ) ; call print_clock('rotHSw:ev:bc')
         CALL print_clock( 'rotHSw:ev:s5' ) ; 
         CALL print_clock( 'rotHSw:ev:b4' ) ; call print_clock('rotHSw:ev:comp')
         CALL print_clock( 'rotHSw:ev:s6' ) ;
         CALL print_clock( 'rotHSw:ev:b5' ) ; call print_clock('rotHSw:ev:sum')
         CALL print_clock( 'rotHSw:ev:s7' ) ; CALL print_clock('rotHSw:ev:b6' ) 
      END IF
   ELSE IF ( isolve == 4 ) THEN
      WRITE( stdout, '(/5x,"Called by *rmmdiagg:")' )
   END IF
   !
   CALL print_clock( 'h_psi' )
   CALL print_clock( 's_psi' )
   CALL print_clock( 'g_psi' )

   IF (real_space ) THEN
    WRITE( stdout, '(/5x,"Called by real space routines:")' )
    CALL print_clock ( 'realus' )
    CALL print_clock ( 'betapointlist' )
    CALL print_clock ( 'addusdens' )
    CALL print_clock ( 'calbec_rs' )
    CALL print_clock ( 's_psir' )
    CALL print_clock ( 'add_vuspsir' )
    CALL print_clock ( 'invfft_orbital' )
    CALL print_clock ( 'fwfft_orbital' )
    CALL print_clock ( 'v_loc_psir' )
   ENDIF
   !
   WRITE( stdout, '(/5x,"Called by h_psi:")' )
   CALL print_clock( 'h_psi:calbec' )
   CALL print_clock( 'vloc_psi' )
   CALL print_clock( 'vloc_psi:tg_gather' )
   CALL print_clock( 'v_loc_psir' )
   CALL print_clock( 'add_vuspsi' )
   CALL print_clock( 'add_vuspsir' )
   CALL print_clock( 'vhpsi' )
   CALL print_clock( 'h_psi_meta' )
   CALL print_clock( 'hs_1psi' )
   CALL print_clock( 's_1psi' )
   !
   WRITE( stdout, '(/5X,"General routines")' )
   !
   CALL print_clock( 'calbec' )
   CALL print_clock( 'fft' )
   CALL print_clock( 'ffts' )
   CALL print_clock( 'fftw' )
   CALL print_clock( 'fftc' )
   CALL print_clock( 'fftcw' )
   CALL print_clock( 'fftr' )
   CALL print_clock( 'interpolate' )
   CALL print_clock( 'davcio' )
   !    
   WRITE( stdout, * )
   !
#if defined (__MPI)
   WRITE( stdout, '(5X,"Parallel routines")' )
   !
   CALL print_clock( 'reduce' )
   CALL print_clock( 'fft_scatt_xy' )
   CALL print_clock( 'fft_scatt_yz' )
   CALL print_clock( 'fft_scatt_tg' )
   CALL print_clock( 'ALLTOALL' )
#endif
   CALL print_clock( 'localization' )
   CALL print_clock( 'measure' )
   !
   IF ( lda_plus_u ) THEN
      WRITE( stdout, '(/,5X,"Hubbard U routines")' )
      IF (lda_plus_u_kind.EQ.0) THEN
         CALL print_clock( 'new_ns' )
         IF (ANY(is_hubbard_back(:))) &
            CALL print_clock( 'new_nsb' )
      ELSEIF (lda_plus_u_kind.EQ.1) THEN
         IF (noncolin) THEN
            CALL print_clock( 'new_ns_nc' )
         ELSE
            CALL print_clock( 'new_ns' )
         ENDIF
      ELSEIF (lda_plus_u_kind.EQ.2) THEN
         CALL print_clock( 'new_nsg' )
         CALL print_clock( 'alloc_neigh' )
      ENDIF
      CALL print_clock( 'vhpsi' )
      CALL print_clock( 'force_hub' )
      CALL print_clock( 'stres_hub' )
   ENDIF
   !
   IF ( xclib_dft_is('hybrid') ) THEN
      WRITE( stdout, '(/,5X,"EXX routines")' )
      CALL print_clock( 'exx_grid' )
      CALL print_clock( 'exxinit' )
      CALL print_clock( 'vexx' )
      CALL print_clock( 'matcalc' )
      CALL print_clock( 'aceupdate' )
      CALL print_clock( 'vexxace' )
      CALL print_clock( 'vexxloc' )
      CALL print_clock( 'aceinit' )
      CALL print_clock( 'exxenergy' )
      IF( okvan) THEN
        WRITE( stdout, '(/,5X,"EXX+US  routines")' )
        CALL print_clock( 'becxx' )
        CALL print_clock( 'addusxx' )
        CALL print_clock( 'newdxx' )
        CALL print_clock( 'qvan_init' )
        CALL print_clock( 'nlxx_pot' )
      ENDIF
      IF ( okpaw ) THEN
        WRITE( stdout, '(/,5X,"EXX+PAW routines")' )
        CALL print_clock('PAW_newdxx')
        CALL print_clock('PAW_xx_nrg')
        CALL print_clock('PAW_keeq')
      ENDIF
   ENDIF
   !
   IF ( okpaw .AND. iverbosity > 0 ) THEN
      WRITE( stdout, '(/,5X,"PAW routines")' )
      ! radial routines:
      CALL print_clock ('PAW_pot')
      CALL print_clock ('PAW_newd')
      CALL print_clock ('PAW_int')
      CALL print_clock ('PAW_ddot')
      CALL print_clock ('PAW_rad_init')
      CALL print_clock ('PAW_energy')
      CALL print_clock ('PAW_symme')
      ! second level routines:
      CALL print_clock ('PAW_rho_lm')
      CALL print_clock ('PAW_h_pot')
      CALL print_clock ('PAW_xc_pot')
      CALL print_clock ('PAW_lm2rad')
      CALL print_clock ('PAW_rad2lm')
      ! third level, or deeper:
      CALL print_clock ('PAW_rad2lm3')
      CALL print_clock ('PAW_gcxc_v')
      CALL print_clock ('PAW_div')
      CALL print_clock ('PAW_grad')
   END IF

   IF ( lelfield ) THEN
      WRITE( stdout, '(/,5X,"Electric-field routines")' )
      call print_clock('h_epsi_set')
      call print_clock('h_epsi_apply')
      call print_clock('c_phase_field')
   END IF
   !
   CALL rism_print_clock()
   !
#if defined(__LEGACY_PLUGINS)
   CALL plugin_clock()
#endif 
#if defined (__ENVIRON)
   IF (use_environ) CALL print_environ_clocks()
#endif
#if defined (__OSCDFT)
   IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==1)) CALL print_oscdft_clocks(oscdft_ctx)
#endif
   !
   RETURN
   !
END SUBROUTINE print_clock_pw
