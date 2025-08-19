!
! Copyright (C) 2001-2022 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE stress( sigma )
  !----------------------------------------------------------------------
  !! Computes the total stress.
  !
  USE io_global,        ONLY : stdout
  USE kinds,            ONLY : DP
  USE cell_base,        ONLY : omega, alat, at, bg
  USE ions_base,        ONLY : nat, ntyp => nsp, ityp, tau, zv, atm
  USE constants,        ONLY : ry_kbar
  USE ener,             ONLY : etxc, vtxc
  USE gvect,            ONLY : ngm, gstart, g, gg, gcutm, gl
  USE fft_base,         ONLY : dfftp
  USE ldaU,             ONLY : lda_plus_u, Hubbard_projectors
  USE lsda_mod,         ONLY : nspin
  USE noncollin_module, ONLY : domag
  USE scf,              ONLY : rho, rho_core, rhog_core
  USE control_flags,    ONLY : iverbosity, gamma_only, llondon, ldftd3, lxdm, &
                               ts_vdw, mbd_vdw
  USE xc_lib,           ONLY : xclib_dft_is
  USE symme,            ONLY : symmatrix
  USE bp,               ONLY : lelfield
  USE uspp,             ONLY : okvan
  USE london_module,    ONLY : stres_london
  USE dftd3_api,        ONLY : get_atomic_number
  USE dftd3_qe,         ONLY : dftd3_pbc_gdisp, dftd3
  USE xdm_module,       ONLY : stress_xdm
  USE exx,              ONLY : exx_stress
  USE tsvdw_module,     ONLY : HtsvdW
  USE libmbd_interface, ONLY : HmbdvdW
  USE rism_module,      ONLY : lrism, stres_rism
  USE esm,              ONLY : do_comp_esm, esm_bc ! for ESM stress
  USE esm,              ONLY : esm_stres_har, esm_stres_ewa, esm_stres_loclong 
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: sigma(3,3)
  !! The stress tensor
  !
  ! ... local variables
  !
  REAL(DP) :: sigmakin(3,3),  sigmaloc(3,3), sigmahar(3,3),                &
              sigmaxc(3,3),   sigmaxcc(3,3), sigmaewa(3,3), sigmanlc(3,3), &
              sigmabare(3,3), sigmah(3,3),   sigmael(3,3),  sigmaion(3,3), &
              sigmad23(3,3),  sigmaxdm(3,3), sigma_ts(3,3), sigma_mbd(3,3),&
              sigma_nonloc_dft(3,3), sigmaexx(3,3)
  REAL(DP) :: sigmaloclong(3,3)  ! for ESM stress
  REAL(DP) :: sigmasol(3,3)      ! for RISM stress
  INTEGER  :: l, m
  !
  ! ... Auxiliary variables for Grimme-D3
  !
  INTEGER  :: atnum(1:nat)
  REAL(DP), ALLOCATABLE :: taupbc(:,:)
  REAL(DP), ALLOCATABLE :: force_d3(:,:)
  !
  WRITE( stdout, '(//5x,"Computing stress (Cartesian axis) and pressure"/)' )
  !
  IF ( lelfield .AND. okvan ) THEN
     CALL infomsg( 'stres', 'stress with USPP and electric fields (Berry) not implemented' )
     RETURN
  END IF
  !
  CALL start_clock( 'stress' )
  !
  !$acc update device( g, gg )
  !FIXME: I don't think the above line is needed
  !
  ! ... contribution from local potential
  !
  CALL stres_loc( sigmaloc )
  !
  IF ( do_comp_esm .AND. ( esm_bc /= 'pbc' ) ) THEN
     ! In ESM, sigmaloc has only short-range term: add long-range term
     CALL esm_stres_loclong( sigmaloclong, rho%of_g(:,1) )
     sigmaloc(:,:) = sigmaloc(:,:) + sigmaloclong(:,:)
  END IF
  !
  ! ... Hartree contribution
  !
  IF ( do_comp_esm .AND. ( esm_bc /= 'pbc' ) )  THEN ! for ESM stress
     CALL esm_stres_har( sigmahar, rho%of_g(:,1) )
  ELSE
     CALL stres_har( sigmahar )
  ENDIF
  !
  ! ... XC contribution (diagonal)
  !
  sigmaxc(:,:) = 0.d0
  DO l = 1, 3
     sigmaxc(l,l) = - (etxc - vtxc) / omega
  ENDDO
  !
  ! ... XC contribution: add gradient corrections (non diagonal)
  !
  CALL stres_gradcorr( rho, rho_core, rhog_core, nspin, domag, &
                       dfftp, g, alat, omega, sigmaxc )
  !
  ! ... meta-GGA contribution
  !
  CALL stres_mgga( sigmaxc )
  !
  ! ... core correction contribution
  !
  CALL stres_cc( sigmaxcc )
  !
  ! ... Ewald contribution
  !
  IF ( do_comp_esm .AND. ( esm_bc /= 'pbc' ) ) THEN ! for ESM stress
     CALL esm_stres_ewa( sigmaewa )
  ELSE
     CALL stres_ewa( alat, nat, ntyp, ityp, zv, at, bg, &
                     tau, omega, g, gg, ngm, gstart,    &
                     gamma_only, gcutm, sigmaewa )
  ENDIF
  !
  ! ... semi-empirical dispersion contribution: Grimme-D2 and D3
  !
  sigmad23( : , : ) = 0.d0
  IF ( llondon ) THEN
    sigmad23 = stres_london( alat , nat , ityp , at , bg , tau , omega )
  ELSE IF ( ldftd3 ) THEN
    CALL start_clock('stres_dftd3')
    ALLOCATE( force_d3(3,nat) )
    force_d3( : , : ) = 0.0_DP
    ! taupbc are atomic positions in alat units, centered around r=0
    ALLOCATE ( taupbc(3,nat) )
    taupbc(:,:) = tau(:,:)
    CALL cryst_to_cart( nat, taupbc, bg, -1 ) 
    taupbc(:,:) = taupbc(:,:) - NINT(taupbc(:,:))
    CALL cryst_to_cart( nat, taupbc, at,  1 ) 
    atnum(:) = get_atomic_number(atm(ityp(:)))
    CALL dftd3_pbc_gdisp( dftd3, alat*taupbc, atnum, alat*at, &
                         force_d3, sigmad23 )
    sigmad23 = 2.d0*sigmad23
    DEALLOCATE( taupbc )
    DEALLOCATE( force_d3 )
    CALL stop_clock('stres_dftd3')
  END IF
  !
  ! ... kinetic + nonlocal contribuition
  !
  CALL stres_knl( sigmanlc, sigmakin )
  !
  DO l = 1, 3
     DO m = 1, 3
        sigmabare(l, m) = sigmaloc(l, m) + sigmanlc(l, m)
     ENDDO
  ENDDO
  !
  ! ... Hubbard contribution
  !     (included by stres_knl if using beta as local projectors)
  !
  sigmah(:,:) = 0.d0
  IF ( lda_plus_u .AND. Hubbard_projectors /= 'pseudo' ) CALL stres_hub( sigmah )
  !
  ! ... Electric field contribution
  !
  sigmael(:,:)=0.d0
  sigmaion(:,:)=0.d0
  !the following is for calculating the improper stress tensor
!  call stress_bp_efield (sigmael )
!  call stress_ion_efield (sigmaion )
  !
  ! ... vdW dispersion contribution: xdm
  !
  sigmaxdm = 0._dp
  IF (lxdm) sigmaxdm = stress_xdm()
  !
  ! ... vdW dispersion contribution: Tkatchenko-Scheffler
  !
  sigma_ts = 0.0_DP
  ! ... vdW dispersion contribution: Many-Body Dispersion
  !
  sigma_mbd = 0.0_DP
  IF ( mbd_vdw ) THEN
    sigma_mbd = -2.0_DP*alat*MATMUL( HmbdvdW, TRANSPOSE(at) )/omega
  ELSE IF ( ts_vdw ) THEN ! Written like this as mbd_vdw needs ts_vdw
    sigma_ts = -2.0_DP*alat*MATMUL( HtsvdW, TRANSPOSE(at) )/omega
  ENDIF
  !
  ! ... DFT-non_local contribution
  !
  sigma_nonloc_dft(:,:) = 0.d0
  CALL stres_nonloc_dft( rho%of_r, rho_core, nspin, sigma_nonloc_dft )
  !
  ! ... The solvation contribution (3D-RISM)
  !
  sigmasol(:,:) = 0.d0
  IF (lrism) CALL stres_rism(sigmasol)
  !
  ! ... Sum all terms
  !
  sigma(:,:) = sigmakin(:,:) + sigmaloc(:,:) + sigmahar(:,:) +  &
               sigmaxc(:,:)  + sigmaxcc(:,:) + sigmaewa(:,:) +  &
               sigmanlc(:,:) + sigmah(:,:)   + sigmael(:,:)  +  &
               sigmaion(:,:) + sigmad23(:,:) + sigmaxdm(:,:) + &
               sigma_nonloc_dft(:,:) + sigma_ts(:,:) + sigma_mbd(:,:) + &
               sigmasol(:,:)
  !
  IF (xclib_dft_is('hybrid')) THEN
     sigmaexx = exx_stress()
     CALL symmatrix( sigmaexx )
     sigma(:,:) = sigma(:,:) + sigmaexx(:,:)
  ELSE
     sigmaexx = 0.d0
  ENDIF
  ! ... Resymmetrize the total stress. This should not be strictly necessary,
  !     but prevents loss of symmetry in long vc-bfgs runs

  CALL symmatrix( sigma )
  !
  ! ... write results in Ry/(a.u.)^3 and in kbar
  !
  IF ( do_comp_esm .AND. ( esm_bc /= 'pbc' ) ) THEN ! for ESM stress
     WRITE( stdout, 9000) (sigma(1,1) + sigma(2,2)) * ry_kbar/3d0, &
     sigma(1,1), sigma(1,2), 0.d0, sigma(1,1)*ry_kbar, sigma(1,2)*ry_kbar, 0.d0,&
     sigma(2,1), sigma(2,2), 0.d0, sigma(2,1)*ry_kbar, sigma(2,2)*ry_kbar, 0.d0,&
     0.d0      , 0.d0      , 0.d0, 0.d0              , 0.d0              , 0.d0
  ELSE
     WRITE( stdout, 9000) (sigma(1,1) + sigma(2,2) + sigma(3,3)) * ry_kbar/3d0, &
                          (sigma(l,1), sigma(l,2), sigma(l,3),                  &
                           sigma(l,1)*ry_kbar, sigma(l,2)*ry_kbar,              &
                           sigma(l,3)*ry_kbar, l=1,3)
  END IF

  IF ( iverbosity > 0 ) WRITE( stdout, 9005) &
     (sigmakin(l,1)*ry_kbar,sigmakin(l,2)*ry_kbar,sigmakin(l,3)*ry_kbar, l=1,3), &
     (sigmaloc(l,1)*ry_kbar,sigmaloc(l,2)*ry_kbar,sigmaloc(l,3)*ry_kbar, l=1,3), &
     (sigmanlc(l,1)*ry_kbar,sigmanlc(l,2)*ry_kbar,sigmanlc(l,3)*ry_kbar, l=1,3), &
     (sigmahar(l,1)*ry_kbar,sigmahar(l,2)*ry_kbar,sigmahar(l,3)*ry_kbar, l=1,3), &
     (sigmaxc (l,1)*ry_kbar,sigmaxc (l,2)*ry_kbar,sigmaxc (l,3)*ry_kbar, l=1,3), &
     (sigmaxcc(l,1)*ry_kbar,sigmaxcc(l,2)*ry_kbar,sigmaxcc(l,3)*ry_kbar, l=1,3), &
     (sigmaewa(l,1)*ry_kbar,sigmaewa(l,2)*ry_kbar,sigmaewa(l,3)*ry_kbar, l=1,3), &
     (sigmah  (l,1)*ry_kbar,sigmah  (l,2)*ry_kbar,sigmah  (l,3)*ry_kbar, l=1,3), &
     (sigmad23(l,1)*ry_kbar,sigmad23(l,2)*ry_kbar,sigmad23(l,3)*ry_kbar, l=1,3), &
     (sigmaxdm(l,1)*ry_kbar,sigmaxdm(l,2)*ry_kbar,sigmaxdm(l,3)*ry_kbar, l=1,3), &
     (sigma_nonloc_dft(l,1)*ry_kbar,sigma_nonloc_dft(l,2)*ry_kbar,sigma_nonloc_dft(l,3)*ry_kbar, l=1,3),&
     (sigma_ts(l,1)*ry_kbar,sigma_ts(l,2)*ry_kbar,sigma_ts(l,3)*ry_kbar, l=1,3), &
     (sigma_mbd(l,1)*ry_kbar,sigma_mbd(l,2)*ry_kbar,sigma_mbd(l,3)*ry_kbar, l=1,3), &
     (sigmasol(l,1)*ry_kbar,sigmasol(l,2)*ry_kbar,sigmasol(l,3)*ry_kbar, l=1,3)

  IF ( xclib_dft_is('hybrid') .AND. (iverbosity > 0) ) WRITE( stdout, 9006) &
     (sigmaexx(l,1)*ry_kbar,sigmaexx(l,2)*ry_kbar,sigmaexx(l,3)*ry_kbar, l=1,3)
9006 format (5x,'EXX     stress (kbar)',3f10.2/2(26x,3f10.2/)/ )
  !
  IF ( lelfield .AND. iverbosity > 0 ) THEN
     WRITE(stdout,*) "Stress tensor electronic el field part:"
     WRITE(stdout,*) (sigmael(l,1),sigmael(l,2),sigmael(l,3), l=1,3)
     WRITE(stdout,*) "Stress tensor electronic el field part:"
     WRITE(stdout,*) (sigmaion(l,1),sigmaion(l,2),sigmaion(l,3), l=1,3)
  ENDIF
  !
  CALL stop_clock( 'stress' )
  !
  RETURN
9000 format (10x,'total   stress  (Ry/bohr**3) ',18x,'(kbar)', &
             &5x,'P=',f12.2/3 (3f13.8,4x,3f12.2/))
9005 format &
         &  (5x,'kinetic stress (kbar)',3f10.2/2(26x,3f10.2/)/ &
         &   5x,'local   stress (kbar)',3f10.2/2(26x,3f10.2/)/ &
         &   5x,'nonloc. stress (kbar)',3f10.2/2(26x,3f10.2/)/ &
         &   5x,'hartree stress (kbar)',3f10.2/2(26x,3f10.2/)/ &
         &   5x,'exc-cor stress (kbar)',3f10.2/2(26x,3f10.2/)/ &
         &   5x,'corecor stress (kbar)',3f10.2/2(26x,3f10.2/)/ &
         &   5x,'ewald   stress (kbar)',3f10.2/2(26x,3f10.2/)/ &
         &   5x,'hubbard stress (kbar)',3f10.2/2(26x,3f10.2/)/ &
         &   5x,'DFT-D   stress (kbar)',3f10.2/2(26x,3f10.2/)/ &
         &   5x,'XDM     stress (kbar)',3f10.2/2(26x,3f10.2/)/ &
         &   5x,'dft-nl  stress (kbar)',3f10.2/2(26x,3f10.2/)/ &
         &   5x,'TS-vdW  stress (kbar)',3f10.2/2(26x,3f10.2/)/ &
         &   5x,'MDB     stress (kbar)',3f10.2/2(26x,3f10.2/)/ &
         &   5x,'3D-RISM stress (kbar)',3f10.2/2(26x,3f10.2/)) 
  !
END SUBROUTINE stress
