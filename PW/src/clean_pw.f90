
! Copyright (C) 2001-2025 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE clean_pw( lflag )
  !----------------------------------------------------------------------
  !! This routine deallocates dynamically allocated arrays.
  !
  !! * If lflag=.TRUE.  all arrays are deallocated (end of calculation);
  !! * if lflag=.FALSE. ion-related variables and arrays allocated.
  !
  !! At the very beginning of the calculation (routines iosys, read_file,
  !! setup, read_pseudo) are not deallocated; all others arrays are.
  !! This is used when a new calculation has to be performed (e.g. in neb,
  !! phonon, vc-relax). Beware: the new calculation should not CALL any
  !! of the routines mentioned above.
  !
  USE cellmd,               ONLY : lmovecell
  USE ions_base,            ONLY : deallocate_ions_base
  USE fixed_occ,            ONLY : f_inp
  USE ktetra,               ONLY : deallocate_tetra
  USE klist,                ONLY : deallocate_igk
  USE gvect,                ONLY : deallocate_gvect
  USE vlocal,               ONLY : strf, vloc
  USE wvfct,                ONLY : g2kin, et, wg, btype
  USE force_mod,            ONLY : force
  USE scf,                  ONLY : rho, v, vltot, rho_core, rhog_core, &
                                   vrs, kedtau, destroy_scf_type, vnew
  USE symm_base,            ONLY : irt
  USE symme,                ONLY : sym_rho_deallocate
  USE wavefunctions,        ONLY : deallocate_wfc, psic, psic_nc
  USE uspp,                 ONLY : deallocate_uspp
  USE uspp_param,           ONLY : upf
  USE atwfc_mod,            ONLY : deallocate_tab_atwfc
  USE m_gth,                ONLY : deallocate_gth
  USE ldaU,                 ONLY : deallocate_hubbard, order_um
  USE extfield,             ONLY : forcefield, forcegate
  USE fft_base,             ONLY : dfftp, dffts  
  USE fft_base,             ONLY : pstickdealloc
  USE fft_types,            ONLY : fft_type_deallocate
  USE noncollin_module,     ONLY : deallocate_noncol
  USE dynamics_module,      ONLY : deallocate_dyn_vars
  USE paw_init,             ONLY : deallocate_paw_internals
  USE atom,                 ONLY : msh, rgrid
  USE radial_grids,         ONLY : deallocate_radial_grid
  USE wannier_new,          ONLY : use_wannier
  !
  USE london_module,        ONLY : dealloca_london
  USE xdm_module,           ONLY : cleanup_xdm
  USE constraints_module,   ONLY : deallocate_constraint
  USE realus,               ONLY : deallocate_realsp
  USE pseudo_types,         ONLY : deallocate_pseudo_upf
  USE bp,                   ONLY : deallocate_bp_efield
  USE exx,                  ONLY : deallocate_exx
  USE Coul_cut_2D,          ONLY : cutoff_2D, lr_Vloc 
  !
  USE control_flags,        ONLY : ts_vdw, mbd_vdw, use_gpu
  USE tsvdw_module,         ONLY : tsvdw_finalize
  USE libmbd_interface,     ONLY : clean_mbd
  USE dftd3_qe,             ONLY : dftd3_clean
  !
  USE control_flags,        ONLY : sic, scissor
  USE sic_mod,              ONLY : deallocate_sic
  USE sci_mod,              ONLY : deallocate_scissor
  !
  USE rism_module,          ONLY : deallocate_rism
#if defined (__ENVIRON)
  USE plugin_flags,         ONLY : use_environ
  USE environ_base_module,  ONLY : clean_environ
#endif
#if defined (__OSCDFT)
   USE plugin_flags,     ONLY : use_oscdft
   USE oscdft_base,      ONLY : oscdft_ctx
#endif
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: lflag
  !! see routine main comments.
  !
  ! ... local variables
  !
  INTEGER :: nt, nr1, nr2, nr3, istat
  !
  IF ( lflag ) THEN
     !
     ! ... arrays allocated at the very beginning of the calculation
     !
     IF( ALLOCATED( upf ) ) THEN
        DO nt = 1, SIZE( upf )
           CALL deallocate_pseudo_upf( upf( nt ) )
        ENDDO
        DEALLOCATE( upf )
     ENDIF
     !
     IF (ALLOCATED(msh)) DEALLOCATE( msh )
     !
     CALL deallocate_radial_grid( rgrid )
     !
     CALL deallocate_ions_base()
     !
     IF ( ALLOCATED( force )      ) DEALLOCATE( force      )
     IF ( ALLOCATED( forcefield ) ) DEALLOCATE( forcefield )
     IF ( ALLOCATED( forcegate )  ) DEALLOCATE( forcegate  )
     IF ( ALLOCATED( irt )        ) DEALLOCATE( irt        )
     !
     CALL dealloca_london()
     CALL cleanup_xdm()
     CALL dftd3_clean()
     CALL deallocate_constraint()
     CALL deallocate_tetra()
     !
  ENDIF
  !
  CALL deallocate_bp_efield()
  !
  CALL deallocate_hubbard( lflag )
  !
  IF ( ALLOCATED( f_inp ) .AND. lflag )  DEALLOCATE( f_inp )
  !
  ! ... arrays in gvect module
  !
  CALL deallocate_gvect( lmovecell )
  !
  CALL sym_rho_deallocate()
  !
  ! ... arrays allocated in allocate_fft.f90 ( and never deallocated )
  !
  CALL destroy_scf_type( rho  )
  CALL destroy_scf_type( v    )
  CALL destroy_scf_type( vnew )
  !
  IF ( ALLOCATED( order_um))     DEALLOCATE (order_um) 
  IF ( ALLOCATED( kedtau ) )     DEALLOCATE( kedtau )
  IF ( ALLOCATED( vltot  ) )     DEALLOCATE( vltot  )
  IF ( ALLOCATED( rho_core  ) )  DEALLOCATE( rho_core  )
  IF ( ALLOCATED( rhog_core ) )  DEALLOCATE( rhog_core )
  IF ( ALLOCATED( psic    ) )    DEALLOCATE( psic    )
  IF ( ALLOCATED( psic_nc ) )    DEALLOCATE( psic_nc )
  !$acc exit data delete(vrs)
  IF ( ALLOCATED( vrs     ) )    DEALLOCATE( vrs     )
  !
  ! ... arrays allocated in allocate_locpot.f90 ( and never deallocated )
  !
  IF ( ALLOCATED( vloc )      )  DEALLOCATE( vloc      )
  IF ( ALLOCATED( cutoff_2D ) )  DEALLOCATE( cutoff_2D )
  IF ( ALLOCATED( lr_Vloc )   )  DEALLOCATE( lr_Vloc   )
  IF ( ALLOCATED( strf )      )  DEALLOCATE( strf      )
  !
  CALL deallocate_tab_atwfc()
  CALL deallocate_uspp() 
  !
  CALL deallocate_gth( lflag ) 
  CALL deallocate_noncol()
  CALL deallocate_igk()
  !
  ! ... arrays allocated in init_run.f90 ( and never deallocated )
  !
  !$acc exit data delete(g2kin)
  IF ( ALLOCATED( g2kin ) )      DEALLOCATE( g2kin )
  !$acc exit data delete(et)
  IF ( ALLOCATED( et ) )         DEALLOCATE( et )
  IF ( ALLOCATED( wg ) )         DEALLOCATE( wg )
  IF ( ALLOCATED( btype ) )      DEALLOCATE( btype )
  !
  ! ... arrays allocated in allocate_wfc.f90 ( and never deallocated )
  !
  CALL deallocate_wfc ( )
  !
  ! ... fft structures allocated in data_structure.f90  
  !
  ! UGLY HACK WARNING: unlike previous versions, fft_type_deallocate
  ! removes all information about FFT grids, including FFT dimensions.
  ! If however FFT dimensions were set from input data, one may end
  ! up with a different grid if FFT grids are re-initialized later.
  ! The following workaround restores the previous functionality.
  ! TODO: replace clean_pw with more fine-grained cleaning routines.
  !
  nr1 = dfftp%nr1; nr2 = dfftp%nr2; nr3 = dfftp%nr3
  CALL fft_type_deallocate( dfftp )
  dfftp%nr1 = nr1; dfftp%nr2 = nr2; dfftp%nr3 = nr3
  !
  nr1 = dffts%nr1; nr2 = dffts%nr2; nr3 = dffts%nr3
  CALL fft_type_deallocate( dffts )
  !
  dffts%nr1 = nr1; dffts%nr2 = nr2; dffts%nr3 = nr3
  !
  ! ... stick-owner matrix allocated in sticks_base
  !
  CALL pstickdealloc()
  !
  ! ... arrays allocated for dynamics
  !
  CALL deallocate_dyn_vars()
  !
  ! ... additional arrays for PAW
  !
  CALL deallocate_paw_internals()
  !
  ! ... arrays for real-space algorithm
  !
  CALL deallocate_realsp()
  !
  ! for Wannier_ac
  IF (use_wannier) CALL wannier_clean()
  !
  CALL deallocate_exx() 
  !
  IF(sic) CALL deallocate_sic()
  IF(scissor) CALL deallocate_scissor()
  !
  IF (ts_vdw .or. mbd_vdw) CALL tsvdw_finalize()
  IF (mbd_vdw) CALL clean_mbd()
  !
  ! ... arrays for RISM
  !
  CALL deallocate_rism( lflag )
  !
#if defined (__LEGACY_PLUGINS) 
  CALL plugin_clean( 'PW', lflag )
#endif 
#if defined (__ENVIRON)
  IF (use_environ) CALL clean_environ('PW', lflag)
#endif
#if defined (__OSCDFT)
     IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==2)) THEN
        DEALLOCATE (oscdft_ctx%inp%occupation)
     ENDIF
#endif
  CALL   plugin_clean('PW', lflag) 
  !
  RETURN
  !
END SUBROUTINE clean_pw
