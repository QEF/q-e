!
! Copyright (C) 2001-2006 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE clean_pw( lflag )
  !----------------------------------------------------------------------
  !    
  ! ... This routine deallocates all dynamically allocated arrays
  !
  USE ions_base,            ONLY : deallocate_ions_base
  USE input_parameters,     ONLY : deallocate_input_parameters
  USE gvect,                ONLY : g, gg, gl, nl, nlm, igtongl, ig1, ig2, ig3, &
                                   eigts1, eigts2, eigts3
  USE gsmooth,              ONLY : nls, nlsm, doublegrid
  USE fixed_occ,            ONLY : f_inp
  USE ktetra,               ONLY : tetra
  USE klist,                ONLY : ngk
  USE reciprocal_vectors,   ONLY : ig_l2g
  USE symme,                ONLY : irt
  USE vlocal,               ONLY : strf, vloc, vnew
  USE wvfct,                ONLY : igk, g2kin, et, wg, btype, gamma_only
  USE force_mod,            ONLY : force
  USE scf,                  ONLY : rho, rhog, vr, vltot, rho_core, rhog_core, &
                                   vrs, tauk, kedtau, kedtaur, taukg
  USE wavefunctions_module, ONLY : evc, psic, psic_nc
  USE us,                   ONLY : qrad, tab, tab_at, tab_d2y, spline_ps
  USE uspp,                 ONLY : deallocate_uspp
  USE ldaU,                 ONLY : ns, nsnew, swfcatom
  USE extfield,             ONLY : forcefield
  USE sticks,               ONLY : dfftp, dffts  
  USE stick_base,           ONLY : sticks_deallocate
  USE fft_types,            ONLY : fft_dlay_deallocate
  USE spin_orb,             ONLY : lspinorb, fcoef
  USE noncollin_module,     ONLY : deallocate_noncol
  USE dynamics_module,      ONLY : deallocate_dyn_vars
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: lflag
    ! if .TRUE. deallocate_ions_base is called 
    ! ( for instance neb and smd do not call it )
  !
  ! ... arrays allocated in input.f90, read_file.f90 or setup.f90
  !
  IF ( lflag ) THEN
     !
     CALL deallocate_ions_base()
     !
     IF ( ALLOCATED( force ) )      DEALLOCATE( force )
     IF ( ALLOCATED( forcefield ) ) DEALLOCATE( forcefield )
     !
  END IF
  !
  IF ( ALLOCATED( f_inp ) )      DEALLOCATE( f_inp )
  IF ( ALLOCATED( tetra ) )      DEALLOCATE( tetra )
  IF ( ALLOCATED( irt ) )        DEALLOCATE( irt )
  !
  ! ... arrays allocated in ggen.f90
  !
  IF ( ALLOCATED( ig_l2g ) )     DEALLOCATE( ig_l2g )
  IF ( ASSOCIATED( gl ) )        DEALLOCATE( gl )
  !
  ! ... arrays allocated in allocate_fft.f90 ( and never deallocated )
  !
  IF ( ALLOCATED( g ) )          DEALLOCATE( g )
  IF ( ALLOCATED( gg ) )         DEALLOCATE( gg )
  IF ( ALLOCATED( nl ) )         DEALLOCATE( nl )  
  IF ( gamma_only ) THEN
     IF ( ALLOCATED( nlm ) )     DEALLOCATE( nlm )
  END IF
  IF ( ALLOCATED( igtongl ) )    DEALLOCATE( igtongl )  
  IF ( ALLOCATED( ig1 ) )        DEALLOCATE( ig1 )
  IF ( ALLOCATED( ig2 ) )        DEALLOCATE( ig2 )
  IF ( ALLOCATED( ig3 ) )        DEALLOCATE( ig3 )
  IF ( ALLOCATED( rho ) )        DEALLOCATE( rho )
  IF ( ALLOCATED( rhog ) )       DEALLOCATE( rhog )
  IF ( ALLOCATED( tauk ) )       DEALLOCATE( tauk )
  IF ( ALLOCATED( kedtau ) )     DEALLOCATE( kedtau )
  IF ( ALLOCATED( kedtaur ) )    DEALLOCATE( kedtaur )
  IF ( ALLOCATED( taukg ) )      DEALLOCATE( taukg )
  IF ( ALLOCATED( vr ) )         DEALLOCATE( vr )
  IF ( ALLOCATED( vltot ) )      DEALLOCATE( vltot )
  IF ( ALLOCATED( vnew ) )       DEALLOCATE( vnew )
  IF ( ALLOCATED( rho_core ) )   DEALLOCATE( rho_core )
  IF ( ALLOCATED( rhog_core ) )  DEALLOCATE( rhog_core )
  IF ( ALLOCATED( psic ) )       DEALLOCATE( psic )
  IF ( ALLOCATED( psic_nc ) )    DEALLOCATE( psic_nc )
  IF ( ALLOCATED( vrs ) )        DEALLOCATE( vrs )
  if (spline_ps) then
    IF ( ALLOCATED( tab_d2y) )     DEALLOCATE( tab_d2y )
  endif
  IF ( doublegrid ) THEN
    IF ( ASSOCIATED( nls ) )     DEALLOCATE( nls )
  END IF
  IF ( doublegrid .AND. gamma_only ) THEN
     IF ( ASSOCIATED( nlsm ) )   DEALLOCATE( nlsm )
  END IF
  !
  ! ... arrays allocated in allocate_locpot.f90 ( and never deallocated )
  !
  IF ( ALLOCATED( vloc ) )       DEALLOCATE( vloc )
  IF ( ALLOCATED( strf ) )       DEALLOCATE( strf )
  IF ( ALLOCATED( eigts1 ) )     DEALLOCATE( eigts1 )
  IF ( ALLOCATED( eigts2 ) )     DEALLOCATE( eigts2 )
  IF ( ALLOCATED( eigts3 ) )     DEALLOCATE( eigts3 )
  !
  ! ... arrays allocated in allocate_nlpot.f90 ( and never deallocated )
  !
  IF ( ALLOCATED( ngk ) )        DEALLOCATE( ngk )
  IF ( ALLOCATED( igk ) )        DEALLOCATE( igk )
  IF ( ALLOCATED( g2kin ) )      DEALLOCATE( g2kin )
  IF ( ALLOCATED( qrad ) )       DEALLOCATE( qrad )
  IF ( ALLOCATED( ns ) )         DEALLOCATE( ns )
  IF ( ALLOCATED( nsnew ) )      DEALLOCATE( nsnew )
  IF ( ALLOCATED( tab ) )        DEALLOCATE( tab )
  IF ( ALLOCATED( tab_at ) )     DEALLOCATE( tab_at )
  IF ( lspinorb ) THEN
     IF ( ALLOCATED( fcoef ) )   DEALLOCATE( fcoef )
  END IF
  !
  CALL deallocate_uspp() 
  !
  CALL deallocate_noncol() 
  !
  ! ... arrays allocated in init_run.f90 ( and never deallocated )
  !
  IF ( ALLOCATED( et ) )         DEALLOCATE( et )
  IF ( ALLOCATED( wg ) )         DEALLOCATE( wg )
  IF ( ALLOCATED( btype ) )      DEALLOCATE( btype )
  !
  ! ... arrays allocated in allocate_wfc.f90 ( and never deallocated )
  !
  IF ( ALLOCATED( evc ) )        DEALLOCATE( evc )
  IF ( ALLOCATED( swfcatom ) )   DEALLOCATE( swfcatom )
  !
  ! ... fft structures allocated in data_structure.f90  
  !
  CALL fft_dlay_deallocate( dfftp )
  CALL fft_dlay_deallocate( dffts )
  !
  ! ... stick-owner matrix allocated in sticks_base
  !
  CALL sticks_deallocate()
  !
  ! ... arrays allocated for dynamics
  !
  CALL deallocate_dyn_vars()
  !
  CALL deallocate_input_parameters()
  !
  RETURN
  !
END SUBROUTINE clean_pw
