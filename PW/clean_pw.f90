SUBROUTINE clean_pw

  USE pwcom,		ONLY : tau, force, ityp, tetra, irt, ig_l2g, & 
                               gg, nl, igtongl, ig1, ig2, ig3, rho, rho_save, &
                               vr, vltot, vnew, rho_core, psic, vrs, nls, &
			       vloc, strf, eigts1, eigts2, eigts3, &
			       igk, igk_l2g, g2kin, indv, nhtol, nhtom, qq, &
			       dvan, deeq, qrad, vkb, qgm, becsum, ns, nsnew, &
			       tab, et, wg, swfcatom, doublegrid
  USE becmod,		ONLY : becp			       
  USE sticks,		ONLY : dfftp, dffts  
  USE stick_base,       ONLY : sticks_deallocate
  USE berry_phase,      ONLY : berry_closeup

#ifdef __SX6

  USE afftnec,		ONLY : auxp

#endif  

  USE fft_types,	ONLY : fft_dlay_deallocate

  IMPLICIT NONE
  
  
  !  arrays allocated in input.f90, read_file.f90 or setup.f90
  
  IF ( ALLOCATED( tau )	)	DEALLOCATE( tau )
  IF ( ALLOCATED( force ) )	DEALLOCATE( force )
  IF ( ALLOCATED( ityp ) )	DEALLOCATE( ityp )
  IF ( ALLOCATED( tetra ) )	DEALLOCATE( tetra )
  IF ( ALLOCATED( irt )	)	DEALLOCATE( irt )
  
  !  arrays allocated in ggen.f90
  
  IF ( ALLOCATED( ig_l2g ) )	DEALLOCATE( ig_l2g )
  
  !  arrays allocated in allocate_fft.f90 ( and never deallocated )
  
  IF ( ALLOCATED( gg ) )	DEALLOCATE( gg )
  IF ( ALLOCATED( nl ) )	DEALLOCATE( nl)  
  IF ( ALLOCATED( igtongl ) )	DEALLOCATE( igtongl )  
  IF ( ALLOCATED( ig1 ) )	DEALLOCATE( ig1 )
  IF ( ALLOCATED( ig2 )	)	DEALLOCATE( ig2 )
  IF ( ALLOCATED( ig3 )	)	DEALLOCATE( ig3 )
  IF ( ALLOCATED( rho )	)	DEALLOCATE( rho )
  IF ( ALLOCATED( rho_save ) )	DEALLOCATE( rho_save )
  IF ( ALLOCATED( vr ) )	DEALLOCATE( vr )
  IF ( ALLOCATED( vltot ) )	DEALLOCATE( vltot )
  IF ( ALLOCATED( vnew ) )	DEALLOCATE( vnew )
  IF ( ALLOCATED( rho_core ) )	DEALLOCATE( rho_core )
  IF ( ALLOCATED( psic ) )	DEALLOCATE( psic )
  IF ( ALLOCATED( vrs )	)	DEALLOCATE( vrs )
  IF ( doublegrid ) THEN
    IF ( ASSOCIATED( nls ) )	DEALLOCATE( nls )
  END IF
  
  !  arrays allocated in allocate_locpot.f90 ( and never deallocated )

  IF ( ALLOCATED( vloc ) )	DEALLOCATE( vloc )
  IF ( ALLOCATED( strf ) )	DEALLOCATE( strf )
  IF ( ALLOCATED( eigts1 ) )	DEALLOCATE( eigts1 )
  IF ( ALLOCATED( eigts2 ) )	DEALLOCATE( eigts2 )
  IF ( ALLOCATED( eigts3 ) )	DEALLOCATE( eigts3 )
  
  !  arrays allocated in allocate_nlpot.f90 ( and never deallocated )
  
  IF ( ALLOCATED( igk )	)	DEALLOCATE( igk )
  IF ( ALLOCATED( igk_l2g ) )	DEALLOCATE( igk_l2g )
  IF ( ALLOCATED( g2kin ) )	DEALLOCATE( g2kin )
  IF ( ALLOCATED( indv ) )	DEALLOCATE( indv )
  IF ( ALLOCATED( nhtol ) )	DEALLOCATE( nhtol )
  IF ( ALLOCATED( nhtom ) )	DEALLOCATE( nhtom )
  IF ( ALLOCATED( qq ) )	DEALLOCATE( qq )
  IF ( ALLOCATED( dvan ) )	DEALLOCATE( dvan )
  IF ( ALLOCATED( deeq ) )	DEALLOCATE( deeq )
  IF ( ALLOCATED( qrad ) )	DEALLOCATE( qrad )
  IF ( ALLOCATED( vkb ) )	DEALLOCATE( vkb )
  IF ( ALLOCATED( qgm ) )	DEALLOCATE( qgm )
  IF ( ALLOCATED( becsum ) )	DEALLOCATE( becsum )
  IF ( ALLOCATED( ns ) )	DEALLOCATE( ns )
  IF ( ALLOCATED( nsnew ) )	DEALLOCATE( nsnew )
  IF ( ALLOCATED( tab ) )	DEALLOCATE( tab )

  !  arrays allocated in allocate_wfc.f90 ( and never deallocated )
 
  IF ( ALLOCATED( et ) ) 	DEALLOCATE( et )
  IF ( ALLOCATED( wg ) )	DEALLOCATE( wg )
  IF ( ALLOCATED( becp ) )	DEALLOCATE( becp )
  IF ( ALLOCATED( swfcatom ) )	DEALLOCATE( swfcatom )

#ifdef __SX6

  !  arrays allocated in cft_3.f90 ( and never deallocated )

  IF ( ALLOCATED( auxp ) )	DEALLOCATE( auxp )

#endif 

  !  fft structures allocated in data_structure.f90  

  CALL fft_dlay_deallocate( dfftp )
  CALL fft_dlay_deallocate( dffts )

  !  sticks owners matrix allocates in sticks_base

  CALL sticks_deallocate( )

  !  deallocate indexes use in calculation of polarizability at gamma

  CALL berry_closeup( )

END SUBROUTINE clean_pw


