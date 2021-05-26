!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE report_clocks_pw()
   !---------------------------------------------------------------------------
   !
   ! ... this routine prints out the clocks at the end of the run
   ! ... it tries to construct the calling tree of the program.
   !
   USE io_global,          ONLY : stdout
   USE control_flags,      ONLY : isolve, iverbosity, gamma_only, lxdm
   USE paw_variables,      ONLY : okpaw
   USE uspp,               ONLY : okvan
   USE realus,             ONLY : real_space
   USE noncollin_module,   ONLY : noncolin
   USE ldaU,               ONLY : lda_plus_u, lda_plus_u_kind, is_hubbard_back
   USE xc_lib,             ONLY : xclib_dft_is
   USE bp,                 ONLY : lelfield
   USE qexsd_module,       ONLY : qexsd_allocate_clock_list
   !
   IMPLICIT NONE
   !
   !
   WRITE( stdout, * )
   CALL qexsd_allocate_clock_list('PW')
   !
   CALL report_clock( 'init_run' )
   CALL report_clock( 'electrons' )
   CALL report_clock( 'update_pot' )
   CALL report_clock( 'forces' )
   CALL report_clock( 'stress' )
   !
   CALL report_clock( 'wfcinit:atomic')
   CALL report_clock( 'wfcinit:wfcrot')
   CALL report_clock( 'potinit' )
   CALL report_clock( 'realus' )
   CALL report_clock( 'realus:boxes',iverbosity)
   CALL report_clock( 'realus:spher',iverbosity)
   CALL report_clock( 'realus:tabp' ,iverbosity)
   !
   CALL report_clock( 'hinit0' )
   IF (lxdm) THEN
      CALL report_clock('init_xdm')
   ENDIF
   !
   
   CALL report_clock( 'c_bands' )
   CALL report_clock( 'sum_band' )
   CALL report_clock( 'v_of_rho' )
   !
   CALL report_clock( 'v_h' ,iverbosity)
   CALL report_clock( 'v_xc',iverbosity )
   CALL report_clock( 'v_xc_meta',iverbosity )
   !
   CALL report_clock( 'newd' )
   CALL report_clock( 'PAW_pot')
   CALL report_clock( 'mix_rho')

   CALL report_clock( 'vdW_energy' )
   CALL report_clock( 'vdW_ffts' )
   CALL report_clock( 'vdW_v'  )
   CALL report_clock( 'vdW_kernel' ) 
   
   IF (lxdm) THEN
      CALL report_clock('energy_xdm')
      CALL report_clock('exdm:environ')
      CALL report_clock('exdm:paw_charge')
      CALL report_clock('exdm:rho')
   END IF

   !
   CALL report_clock( 'init_us_2' )
   IF ( isolve == 0 ) THEN
      CALL report_clock( 'regterg' )    ; CALL report_clock( 'cegterg' )
   ELSE  IF (isolve == 1) THEN
      CALL report_clock( 'rcgdiagg' )   ; CALL report_clock( 'ccgdiagg' )
      CALL report_clock( 'wfcrot' )
   ELSE  IF (isolve == 2) THEN
      CALL report_clock( 'ppcg_gamma' ) ; CALL report_clock( 'ppcg_k' )
      CALL report_clock( 'wfcrot' )
   ELSE  IF (isolve == 3) THEN
      CALL report_clock( 'paro_gamma' ) ; CALL report_clock( 'paro_k' )
   ENDIF
   !
   CALL report_clock( 'sum_band:weights' , iverbosity)
   CALL report_clock( 'sum_band:loop' ,iverbosity)
   CALL report_clock( 'sum_band:buffer' ,iverbosity)
   CALL report_clock( 'sum_band:init_us_2' ,iverbosity)
   CALL report_clock( 'sum_band:calbec' ,iverbosity)
   CALL report_clock( 'sum_band:becsum',iverbosity )
   CALL report_clock( 'addusdens' ,iverbosity)
   CALL report_clock( 'addusd:skk' ,iverbosity)
   CALL report_clock( 'addusd:dgemm' ,iverbosity)
   CALL report_clock( 'addusd:qvan2',iverbosity )
   !
   IF ( isolve == 0 ) THEN
      IF ( gamma_only ) THEN
         CALL report_clock( 'rdiaghg' )
         CALL report_clock( 'regterg:overlap', iverbosity)
         CALL report_clock( 'regterg:update' , iverbosity)
         CALL report_clock( 'regterg:last'   , iverbosity)
         CALL report_clock( 'rdiaghg:choldc' , iverbosity)
         CALL report_clock( 'rdiaghg:inversion', iverbosity )
         CALL report_clock( 'rdiaghg:paragemm' , iverbosity)
      ELSE
         CALL report_clock( 'cdiaghg' )
         !
         CALL report_clock( 'cegterg:overlap', iverbosity )
         CALL report_clock( 'cegterg:update',iverbosity )
         CALL report_clock( 'cegterg:last' ,iverbosity)
         CALL report_clock( 'cdiaghg:choldc' , iverbosity)
         CALL report_clock( 'cdiaghg:inversion' , iverbosity)
         CALL report_clock( 'cdiaghg:paragemm' , iverbosity)
      END IF
   ELSE IF ( isolve == 2 ) THEN
      CALL report_clock( 'ppcg:zgemm', iverbosity ) ; CALL report_clock( 'ppcg:dgemm',iverbosity )
      CALL report_clock( 'ppcg:hpsi' , iverbosity)
      CALL report_clock( 'ppcg:cholQR' ,iverbosity)
      CALL report_clock( 'ppcg:RR' , iverbosity)
      CALL report_clock( 'ppcg:ZTRSM' , iverbosity) ; CALL report_clock( 'ppcg:DTRSM', iverbosity )
      CALL report_clock( 'ppcg:lock', iverbosity )

   ELSE IF ( isolve == 3 ) THEN
      !
      CALL report_clock( 'paro:init' , iverbosity )
      CALL report_clock( 'paro:pack' , iverbosity )
      CALL report_clock( 'paro:zero' , iverbosity )
      CALL report_clock( 'paro:mp_bar' , iverbosity )
      CALL report_clock( 'paro:mp_sum' , iverbosity )
      CALL report_clock( 'pcg' , iverbosity )
      CALL report_clock( 'pcg:hs_1psi' , iverbosity )
      CALL report_clock( 'pcg:ortho' , iverbosity )
      CALL report_clock( 'pcg:move' , iverbosity )

      CALL report_clock( 'rotHSw' , iverbosity )
      CALL report_clock( 'rotHSw:move' , iverbosity )
      CALL report_clock( 'rotHSw:hc' , iverbosity )
      CALL report_clock( 'rotHSw:diag' , iverbosity )
      CALL report_clock( 'rotHSw:evc' , iverbosity )
      CALL report_clock( 'rotHSw:hc:b0' , iverbosity ) ; 
      CALL report_clock( 'rotHSw:hc:s1' , iverbosity ) ; call report_clock('rotHSw:hc:comp', iverbosity )
      CALL report_clock( 'rotHSw:hc:b1' , iverbosity ) ; 
      CALL report_clock( 'rotHSw:hc:s2' , iverbosity ) ; 
      CALL report_clock( 'rotHSw:hc:s3' , iverbosity ) ; call report_clock('rotHSw:hc:rs', iverbosity )
      CALL report_clock( 'rotHSw:hc:b2' , iverbosity ) ; call report_clock('rotHSw:hc:sy', iverbosity )
      CALL report_clock( 'rotHSw:hc:s4' , iverbosity ) ; CALL report_clock('rotHSw:hc:b3' , iverbosity ) 
      CALL report_clock( 'rotHSw:ev:b0' , iverbosity ) ; 
      CALL report_clock( 'rotHSw:ev:b3' , iverbosity ) ; call report_clock('rotHSw:ev:bc', iverbosity )
      CALL report_clock( 'rotHSw:ev:s5' , iverbosity ) ; 
      CALL report_clock( 'rotHSw:ev:b4' , iverbosity ) ; call report_clock('rotHSw:ev:comp', iverbosity )
      CALL report_clock( 'rotHSw:ev:s6' , iverbosity ) ;
      CALL report_clock( 'rotHSw:ev:b5' , iverbosity ) ; call report_clock('rotHSw:ev:sum', iverbosity )
      CALL report_clock( 'rotHSw:ev:s7' , iverbosity ) ; CALL report_clock('rotHSw:ev:b6' , iverbosity ) 
   END IF
   !
   CALL report_clock( 'h_psi' )
   CALL report_clock( 's_psi' )
   CALL report_clock( 'g_psi' )

   IF (real_space ) THEN
     CALL report_clock ( 'realus' )
     CALL report_clock ( 'betapointlist' )
     CALL report_clock ( 'addusdens' )
     CALL report_clock ( 'calbec_rs' )
     CALL report_clock ( 's_psir' )
     CALL report_clock ( 'add_vuspsir' )
     CALL report_clock ( 'invfft_orbital' )
     CALL report_clock ( 'fwfft_orbital' )
     CALL report_clock ( 'v_loc_psir' )
   ENDIF
   !
   CALL report_clock( 'h_psi:calbec' )
   CALL report_clock( 'vloc_psi' )
   CALL report_clock( 'vloc_psi:tg_gather' )
   CALL report_clock( 'v_loc_psir' )
   CALL report_clock( 'add_vuspsi' )
   CALL report_clock( 'add_vuspsir' )
   CALL report_clock( 'vhpsi' )
   CALL report_clock( 'h_psi_meta' )
   CALL report_clock( 'hs_1psi' )
   CALL report_clock( 's_1psi' )
   !
   !
   CALL report_clock( 'calbec' )
   CALL report_clock( 'fft' )
   CALL report_clock( 'ffts' )
   CALL report_clock( 'fftw' )
   CALL report_clock( 'fftc' )
   CALL report_clock( 'fftcw' )
   CALL report_clock( 'interpolate' )
   CALL report_clock( 'davcio' )
   !    
   !
#if defined (__MPI)
   !
   CALL report_clock( 'reduce' )
   CALL report_clock( 'fft_scatt_xy' )
   CALL report_clock( 'fft_scatt_yz' )
   CALL report_clock( 'fft_scatt_tg' )
   CALL report_clock( 'ALLTOALL' )
#endif
   CALL report_clock( 'localization' )
   CALL report_clock( 'measure' )
   !
   IF ( lda_plus_u ) THEN
      IF (lda_plus_u_kind.EQ.0) THEN
         CALL report_clock( 'new_ns' )
         IF (ANY(is_hubbard_back(:))) &
            CALL report_clock( 'new_nsb' )
      ELSEIF (lda_plus_u_kind.EQ.1) THEN
         IF (noncolin) THEN
            CALL report_clock( 'new_ns_nc' )
         ELSE
            CALL report_clock( 'new_ns' )
         ENDIF
      ELSEIF (lda_plus_u_kind.EQ.2) THEN
         CALL report_clock( 'new_nsg' )
         CALL report_clock( 'alloc_neigh' )
      ENDIF
      CALL report_clock( 'new_ns' )
      CALL report_clock( 'vhpsi' )
      CALL report_clock( 'force_hub' )
      CALL report_clock( 'stres_hub' )
   ENDIF
   !
   IF ( xclib_dft_is('hybrid') ) THEN
      CALL report_clock( 'exx_grid' )
      CALL report_clock( 'exxinit' )
      CALL report_clock( 'vexx' )
      CALL report_clock( 'matcalc' )
      CALL report_clock( 'aceupdate' )
      CALL report_clock( 'vexxace' )
      CALL report_clock( 'vexxloc' )
      CALL report_clock( 'aceinit' )
      CALL report_clock( 'exxenergy' )
      IF( okvan) THEN
        CALL report_clock( 'becxx' )
        CALL report_clock( 'addusxx' )
        CALL report_clock( 'newdxx' )
        CALL report_clock( 'qvan_init' )
        CALL report_clock( 'nlxx_pot' )
      ENDIF
      IF ( okpaw ) THEN
        CALL report_clock('PAW_newdxx')
        CALL report_clock('PAW_xx_nrg')
        CALL report_clock('PAW_keeq')
      ENDIF
   ENDIF
   !
   IF ( okpaw ) THEN
      ! radial routines:
      CALL report_clock ('PAW_pot', iverbosity)
      CALL report_clock ('PAW_newd', iverbosity)
      CALL report_clock ('PAW_int', iverbosity)
      CALL report_clock ('PAW_ddot', iverbosity)
      CALL report_clock ('PAW_rad_init', iverbosity)
      CALL report_clock ('PAW_energy', iverbosity)
      CALL report_clock ('PAW_symme', iverbosity)
      ! second level routines:
      CALL report_clock ('PAW_rho_lm', iverbosity)
      CALL report_clock ('PAW_h_pot', iverbosity)
      CALL report_clock ('PAW_xc_pot', iverbosity)
      CALL report_clock ('PAW_lm2rad', iverbosity)
      CALL report_clock ('PAW_rad2lm', iverbosity)
      ! third level, or deeper:
      CALL report_clock ('PAW_rad2lm3', iverbosity)
      CALL report_clock ('PAW_gcxc_v', iverbosity)
      CALL report_clock ('PAW_div', iverbosity)
      CALL report_clock ('PAW_grad', iverbosity)
      !END IF 
   END IF

   IF ( lelfield ) THEN
      call report_clock('h_epsi_set')
      call report_clock('h_epsi_apply')
      call report_clock('c_phase_field')
   END IF
   !
   !
   RETURN
   !
  CONTAINS
    SUBROUTINE report_clock(label, verbflag_)
      USE qexsd_module, ONLY: qexsd_add_label
      IMPLICIT NONE 
      CHARACTER(*),INTENT(IN) :: label
      INTEGER,OPTIONAL,INTENT(IN) :: verbflag_
      !
      INTEGER  :: verbflag
      IF (PRESENT(verbflag_)) THEN
         verbflag = verbflag_
      ELSE 
         verbflag = 42
      END IF
      CALL qexsd_add_label(label)
    END SUBROUTINE report_clock 
END SUBROUTINE report_clocks_pw
