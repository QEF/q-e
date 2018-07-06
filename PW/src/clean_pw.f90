
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! TB
! included deallocation of forcefield of gate 'forcegate'
!
!----------------------------------------------------------------------
SUBROUTINE clean_pw( lflag )
  !----------------------------------------------------------------------
  !    
  ! ... This routine deallocates dynamically allocated arrays
  ! ... if lflag=.TRUE.  all arrays are deallocated (end of calculation)
  ! ... if lflag=.FALSE. ion-related variables and arrays allocated
  ! ... at the very beginning of the calculation (routines iosys, read_file,
  ! ... setup, read_pseudo) are not deallocated; all others arrays are.
  ! ... This is used when a new calculation has to be performed (e.g. in neb,
  ! ... phonon, vc-relax). Beware: the new calculation should not call any
  ! ... of the routines mentioned above
  !
  USE basis,                ONLY : swfcatom
  USE cellmd,               ONLY : lmovecell
  USE ions_base,            ONLY : deallocate_ions_base
  USE gvect,                ONLY : g, gg, gl, igtongl, mill, &
                                   eigts1, eigts2, eigts3
  USE fixed_occ,            ONLY : f_inp
  USE ktetra,               ONLY : deallocate_tetra
  USE klist,                ONLY : deallocate_igk
  USE gvect,                ONLY : ig_l2g
  USE vlocal,               ONLY : strf, vloc
  USE wvfct,                ONLY : g2kin, et, wg, btype
  USE force_mod,            ONLY : force
  USE scf,                  ONLY : rho, v, vltot, rho_core, rhog_core, &
                                   vrs, kedtau, destroy_scf_type, vnew
  USE symm_base,            ONLY : irt
  USE symme,                ONLY : sym_rho_deallocate
  USE wavefunctions_module, ONLY : evc, psic, psic_nc
  USE us,                   ONLY : qrad, tab, tab_at, tab_d2y, spline_ps
  USE uspp,                 ONLY : deallocate_uspp
  USE uspp_param,           ONLY : upf
  USE m_gth,                ONLY : deallocate_gth
  USE ldaU,                 ONLY : deallocate_ldaU
  USE extfield,             ONLY : forcefield, forcegate
  USE fft_base,             ONLY : dfftp, dffts  
  USE fft_base,             ONLY : pstickdealloc
  USE fft_types,            ONLY : fft_type_deallocate
  USE spin_orb,             ONLY : lspinorb, fcoef
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
  USE control_flags,        ONLY : ts_vdw
  USE tsvdw_module,         ONLY : tsvdw_finalize
  USE dftd3_qe,             ONLY : dftd3_clean
  !
  USE wavefunctions_module_gpum, ONLY : deallocate_wavefunctions_module_gpu
  USE wvfct_gpum,                ONLY : using_et, using_et_d
  USE wvfct_gpum,                ONLY : using_g2kin, using_g2kin_d
  USE gvect_gpum,                ONLY : deallocate_gvect_gpu !using_g, using_gg, using_g_d, using_gg_d
  USE scf_gpum,                  ONLY : using_vrs, using_vrs_d
  USE uspp_gpum,                 ONLY : deallocate_uspp_gpu
  USE us_gpum,                   ONLY : deallocate_us_gpu
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: lflag
  !
  INTEGER :: nt, nr1,nr2,nr3
  !
  IF ( lflag ) THEN
     !
     ! ... arrays allocated at the very beginning of the calculation
     !
     IF( ALLOCATED( upf ) ) THEN
        DO nt = 1, SIZE( upf )
           CALL deallocate_pseudo_upf( upf( nt ) )
        END DO
        DEALLOCATE( upf )
     END IF
     IF (ALLOCATED(msh)) DEALLOCATE (msh)
     CALL deallocate_radial_grid(rgrid)
     !
     CALL deallocate_ions_base()
     !
     IF ( ALLOCATED( force ) )      DEALLOCATE( force )
     IF ( ALLOCATED( forcefield ) ) DEALLOCATE( forcefield )
     IF ( ALLOCATED( forcegate ) )  DEALLOCATE( forcegate )
     IF ( ALLOCATED( irt ) )        DEALLOCATE( irt )
     !
     CALL dealloca_london()
     CALL cleanup_xdm()
     CALL dftd3_clean()
     CALL deallocate_constraint()
     CALL deallocate_tetra ( )
     !
  END IF
  !
  CALL deallocate_bp_efield()
  !
  CALL deallocate_ldaU ( lflag )
  !
  IF ( ALLOCATED( f_inp ) .and. lflag )      DEALLOCATE( f_inp )
  !
  ! ... arrays allocated in ggen.f90
  !
  IF ( ALLOCATED( ig_l2g ) )     DEALLOCATE( ig_l2g )
  IF ( .NOT. lmovecell ) THEN
     IF ( ASSOCIATED( gl ) )     DEALLOCATE ( gl )
  END IF
  !
  CALL sym_rho_deallocate ( )
  !
  ! ... arrays allocated in allocate_fft.f90 ( and never deallocated )
  !
  IF ( ALLOCATED( g ) )          DEALLOCATE( g )
  IF ( ALLOCATED( gg ) )         DEALLOCATE( gg )
  !
  CALL deallocate_gvect_gpu()
  !CALL using_g(2); CALL using_gg(2)       ! Trick to deallocate
  !CALL using_g_d(2); CALL using_gg_d(2)
  !
  IF ( ALLOCATED( igtongl ) )    DEALLOCATE( igtongl )  
  IF ( ALLOCATED( mill ) )       DEALLOCATE( mill )
  call destroy_scf_type(rho)
  call destroy_scf_type(v)
  call destroy_scf_type(vnew)
  IF ( ALLOCATED( kedtau ) )     DEALLOCATE( kedtau )
  IF ( ALLOCATED( vltot ) )      DEALLOCATE( vltot )
  IF ( ALLOCATED( rho_core ) )   DEALLOCATE( rho_core )
  IF ( ALLOCATED( rhog_core ) )  DEALLOCATE( rhog_core )
  IF ( ALLOCATED( psic ) )       DEALLOCATE( psic )
  IF ( ALLOCATED( psic_nc ) )    DEALLOCATE( psic_nc )
  IF ( ALLOCATED( vrs ) )        DEALLOCATE( vrs )
  CALL using_vrs(2); CALL using_vrs_d(2); ! Trick to deallocate
  if (spline_ps) then
    IF ( ALLOCATED( tab_d2y) )     DEALLOCATE( tab_d2y )
  endif
  !
  ! ... arrays allocated in allocate_locpot.f90 ( and never deallocated )
  !
  IF ( ALLOCATED( vloc ) )       DEALLOCATE( vloc ) 
  IF ( ALLOCATED( cutoff_2D ) )  DEALLOCATE( cutoff_2D )
  IF ( ALLOCATED( lr_Vloc ) )    DEALLOCATE( lr_Vloc )
  IF ( ALLOCATED( strf ) )       DEALLOCATE( strf )
  IF ( ALLOCATED( eigts1 ) )     DEALLOCATE( eigts1 )
  IF ( ALLOCATED( eigts2 ) )     DEALLOCATE( eigts2 )
  IF ( ALLOCATED( eigts3 ) )     DEALLOCATE( eigts3 )
  !
  ! ... arrays allocated in allocate_nlpot.f90 ( and never deallocated )
  !
  IF ( ALLOCATED( g2kin ) )      DEALLOCATE( g2kin )
  CALL using_g2kin(2) ;     CALL using_g2kin_d(2)
  IF ( ALLOCATED( qrad ) )       DEALLOCATE( qrad )
  IF ( ALLOCATED( tab ) )        DEALLOCATE( tab )
  IF ( ALLOCATED( tab_at ) )     DEALLOCATE( tab_at )
  IF ( lspinorb ) THEN
     IF ( ALLOCATED( fcoef ) )   DEALLOCATE( fcoef )
  END IF
  !
  CALL deallocate_igk ( )
  CALL deallocate_uspp() 
  CALL deallocate_uspp_gpu()
  CALL deallocate_us_gpu()
  !
  CALL deallocate_gth( lflag ) 
  CALL deallocate_noncol() 
  !
  ! ... arrays allocated in init_run.f90 ( and never deallocated )
  !
  IF ( ALLOCATED( et ) )         DEALLOCATE( et )
  CALL using_et(2); CALL using_et_d(2);
  IF ( ALLOCATED( wg ) )         DEALLOCATE( wg )
  IF ( ALLOCATED( btype ) )      DEALLOCATE( btype )
  !
  ! ... arrays allocated in allocate_wfc.f90 ( and never deallocated )
  !
  IF ( ALLOCATED( evc ) )        DEALLOCATE( evc )
  IF ( ALLOCATED( swfcatom ) )   DEALLOCATE( swfcatom )
  !
  CALL deallocate_wavefunctions_module_gpu()
  !
  ! ... fft structures allocated in data_structure.f90  
  !
  ! UGLY HACK WARNING: unlike previous versions, fft_type_deallocate
  ! removes all information about FFT grids, including FFT dimensions.
  ! If however FFT dimensions were set from input data, one may end
  ! up with a different grid if FFT grids are re-initialized later.
  ! The following workaround restores the previous functionality.
  ! TODO: replace clean_pw with more fine-grained cleaning routines.
  nr1 = dfftp%nr1; nr2 = dfftp%nr2; nr3 = dfftp%nr3
  CALL fft_type_deallocate( dfftp )
  dfftp%nr1 = nr1; dfftp%nr2 = nr2; dfftp%nr3 = nr3
  !
  nr1 = dffts%nr1; nr2 = dffts%nr2; nr3 = dffts%nr3
  CALL fft_type_deallocate( dffts )
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
  CALL  deallocate_realsp()
  !
  ! for Wannier_ac
  if (use_wannier) CALL wannier_clean()
  !
  CALL deallocate_exx ( ) 
  !
  IF (ts_vdw) CALL tsvdw_finalize()
  !
  CALL plugin_clean( 'PW', lflag )
  !
  RETURN
  !
END SUBROUTINE clean_pw
