!
! Copyright (C) 2001-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine stress ( sigma )
  !----------------------------------------------------------------------
  !
  USE io_global,     ONLY : stdout
  USE kinds,         ONLY : DP
  USE cell_base,     ONLY : omega, alat, at, bg
  USE ions_base,     ONLY : nat, ntyp => nsp, ityp, tau, zv
  USE constants,     ONLY : ry_kbar
  USE ener,          ONLY : etxc, vtxc
  USE gvect,         ONLY : ngm, gstart, nl, g, gg, gcutm
  USE fft_base,      ONLY : dfftp
  USE ldaU,          ONLY : lda_plus_u, U_projection
  USE lsda_mod,      ONLY : nspin
  USE scf,           ONLY : rho, rho_core, rhog_core
  USE control_flags, ONLY : iverbosity, gamma_only, llondon, lxdm, ts_vdw
  USE noncollin_module, ONLY : noncolin
  USE funct,         ONLY : dft_is_meta, dft_is_gradient
  USE symme,         ONLY : symmatrix
  USE bp,            ONLY : lelfield
  USE uspp,          ONLY : okvan
  USE london_module, ONLY : stres_london
  USE xdm_module,    ONLY : stress_xdm
  USE exx,           ONLY : exx_stress
  USE funct,         ONLY : dft_is_hybrid
  use tsvdw_module,  only : HtsvdW
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: sigma(3,3)
  !
  real(DP) :: sigmakin (3, 3), sigmaloc (3, 3), sigmahar (3, 3), &
       sigmaxc (3, 3), sigmaxcc (3, 3), sigmaewa (3, 3), sigmanlc (3, 3), &
       sigmabare (3, 3), sigmah (3, 3), sigmael( 3, 3), sigmaion(3, 3), &
       sigmalon ( 3 , 3 ), sigmaxdm(3, 3), sigma_nonloc_dft (3 ,3), sigmaexx(3,3), sigma_ts(3,3)
  integer :: l, m
  !
  WRITE( stdout, '(//5x,"Computing stress (Cartesian axis) and pressure"/)')

  IF ( dft_is_meta() ) THEN
     CALL infomsg ('stress','Meta-GGA and stress not implemented')
     RETURN
  ELSE IF ( noncolin .AND. dft_is_gradient() ) then
     CALL infomsg('stres', 'noncollinear stress + GGA not implemented')
     RETURN
  ELSE IF ( lelfield .AND. okvan ) THEN
     CALL infomsg('stres', 'stress with USPP and electric fields (Berry) not implemented')
     RETURN
  END IF
  !
  call start_clock ('stress')
  !
  !   contribution from local  potential
  !
  call stres_loc (sigmaloc)
  !
  !  hartree contribution
  !
  call stres_har (sigmahar)
  !
  !  xc contribution (diagonal)
  !
  sigmaxc(:,:) = 0.d0
  do l = 1, 3
     sigmaxc (l, l) = - (etxc - vtxc) / omega
  enddo
  !
  !  xc contribution: add gradient corrections (non diagonal)
  !
  call stres_gradcorr ( rho%of_r, rho%of_g, rho_core, rhog_core, nspin, &
                        dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nnr, nl, &
                        ngm, g, alat, omega, sigmaxc)
  !
  ! core correction contribution
  !
  call stres_cc (sigmaxcc)
  !
  !  ewald contribution
  !
  call stres_ewa (alat, nat, ntyp, ityp, zv, at, bg, tau, omega, g, &
       gg, ngm, gstart, gamma_only, gcutm, sigmaewa)
  !
  !  semi-empirical dispersion contribution
  !
  sigmalon ( : , : ) = 0.d0
  !
  IF ( llondon ) &
    sigmalon = stres_london ( alat , nat , ityp , at , bg , tau , omega )

  ! xdm dispersion
  sigmaxdm = 0._dp
  if (lxdm) sigmaxdm = stress_xdm()
  !
  !  kinetic + nonlocal contribuition
  !
  call stres_knl (sigmanlc, sigmakin)
  !
  do l = 1, 3
     do m = 1, 3
        sigmabare (l, m) = sigmaloc (l, m) + sigmanlc (l, m)
     enddo
  enddo
  !
  !  Hubbard contribution
  !  (included by stres_knl if using beta as local projectors)
  !
  sigmah(:,:) = 0.d0
  IF ( lda_plus_u .AND. U_projection.NE.'pseudo' ) CALL stres_hub(sigmah)
  !
  !   Electric field contribution
  !
  sigmael(:,:)=0.d0
  sigmaion(:,:)=0.d0
  !the following is for calculating the improper stress tensor
!  call stress_bp_efield (sigmael )
!  call stress_ion_efield (sigmaion )

  sigma_ts = 0.0_DP
  if(ts_vdw) sigma_ts = -2.0_DP*alat*MATMUL( HtsvdW, transpose(at) )/omega
  !
  !   DFT-non_local contribution
  !
  sigma_nonloc_dft (:,:) = 0.d0
  call stres_nonloc_dft(rho%of_r, rho_core, nspin, sigma_nonloc_dft)
  !
  ! SUM
  !
  sigma(:,:) = sigmakin(:,:) + sigmaloc(:,:) + sigmahar(:,:) + &
               sigmaxc(:,:) + sigmaxcc(:,:) + sigmaewa(:,:) + &
               sigmanlc(:,:) + sigmah(:,:) + sigmael(:,:) +  &
               sigmaion(:,:) + sigmalon(:,:) + sigmaxdm(:,:) + &
               sigma_nonloc_dft(:,:) + sigma_ts(:,:)
  !
  IF (dft_is_hybrid()) THEN
     sigmaexx = exx_stress()
     CALL symmatrix ( sigmaexx )
     sigma(:,:) = sigma(:,:) + sigmaexx(:,:)
  ELSE
     sigmaexx = 0.d0
  ENDIF
  ! Resymmetrize the total stress. This should not be strictly necessary,
  ! but prevents loss of symmetry in long vc-bfgs runs

  CALL symmatrix ( sigma )
  !
  ! write results in Ryd/(a.u.)^3 and in kbar
  !
  WRITE( stdout, 9000) (sigma(1,1) + sigma(2,2) + sigma(3,3)) * ry_kbar/3d0, &
                  (sigma(l,1), sigma(l,2), sigma(l,3),                    &
            sigma(l,1)*ry_kbar, sigma(l,2)*ry_kbar, sigma(l,3)*ry_kbar, l=1,3)

  if ( iverbosity > 0 ) WRITE( stdout, 9005) &
     (sigmakin(l,1)*ry_kbar,sigmakin(l,2)*ry_kbar,sigmakin(l,3)*ry_kbar, l=1,3),&
     (sigmaloc(l,1)*ry_kbar,sigmaloc(l,2)*ry_kbar,sigmaloc(l,3)*ry_kbar, l=1,3),&
     (sigmanlc(l,1)*ry_kbar,sigmanlc(l,2)*ry_kbar,sigmanlc(l,3)*ry_kbar, l=1,3),&
     (sigmahar(l,1)*ry_kbar,sigmahar(l,2)*ry_kbar,sigmahar(l,3)*ry_kbar, l=1,3),&
     (sigmaxc (l,1)*ry_kbar,sigmaxc (l,2)*ry_kbar,sigmaxc (l,3)*ry_kbar, l=1,3),&
     (sigmaxcc(l,1)*ry_kbar,sigmaxcc(l,2)*ry_kbar,sigmaxcc(l,3)*ry_kbar, l=1,3),&
     (sigmaewa(l,1)*ry_kbar,sigmaewa(l,2)*ry_kbar,sigmaewa(l,3)*ry_kbar, l=1,3),&
     (sigmah  (l,1)*ry_kbar,sigmah  (l,2)*ry_kbar,sigmah  (l,3)*ry_kbar, l=1,3),&
     (sigmalon(l,1)*ry_kbar,sigmalon(l,2)*ry_kbar,sigmalon(l,3)*ry_kbar, l=1,3), &
     (sigmaxdm(l,1)*ry_kbar,sigmaxdm(l,2)*ry_kbar,sigmaxdm(l,3)*ry_kbar, l=1,3), &
     (sigma_nonloc_dft(l,1)*ry_kbar,sigma_nonloc_dft(l,2)*ry_kbar,sigma_nonloc_dft(l,3)*ry_kbar, l=1,3),&
     (sigma_ts(l,1)*ry_kbar,sigma_ts(l,2)*ry_kbar,sigma_ts(l,3)*ry_kbar, l=1,3)

  IF ( dft_is_hybrid() .AND. (iverbosity > 0) ) WRITE( stdout, 9006) &
     (sigmaexx(l,1)*ry_kbar,sigmaexx(l,2)*ry_kbar,sigmaexx(l,3)*ry_kbar, l=1,3)
9006 format (5x,'EXX     stress (kbar)',3f10.2/2(26x,3f10.2/)/ )

  if( lelfield .and. iverbosity > 0 ) then
     write(stdout,*) "Stress tensor electronic el field part:"
     write(stdout,*) (sigmael(l,1),sigmael(l,2),sigmael(l,3), l=1,3)
     write(stdout,*) "Stress tensor electronic el field part:"
     write(stdout,*) (sigmaion(l,1),sigmaion(l,2),sigmaion(l,3), l=1,3)
  endif

  call stop_clock ('stress')

  return
9000 format (10x,'total   stress  (Ry/bohr**3) ',18x,'(kbar)', &
             &5x,'P=',f8.2/3 (3f13.8,4x,3f10.2/))
9005 format &
         &  (5x,'kinetic stress (kbar)',3f10.2/2(26x,3f10.2/)/ &
         &   5x,'local   stress (kbar)',3f10.2/2(26x,3f10.2/)/ &
         &   5x,'nonloc. stress (kbar)',3f10.2/2(26x,3f10.2/)/ &
         &   5x,'hartree stress (kbar)',3f10.2/2(26x,3f10.2/)/ &
         &   5x,'exc-cor stress (kbar)',3f10.2/2(26x,3f10.2/)/ &
         &   5x,'corecor stress (kbar)',3f10.2/2(26x,3f10.2/)/ &
         &   5x,'ewald   stress (kbar)',3f10.2/2(26x,3f10.2/)/ &
         &   5x,'hubbard stress (kbar)',3f10.2/2(26x,3f10.2/)/ &
         &   5x,'london  stress (kbar)',3f10.2/2(26x,3f10.2/)/ &
         &   5x,'XDM     stress (kbar)',3f10.2/2(26x,3f10.2/)/ &
         &   5x,'dft-nl  stress (kbar)',3f10.2/2(26x,3f10.2/)/ &
         &   5x,'TS-vdW  stress (kbar)',3f10.2/2(26x,3f10.2/)/ )
end subroutine stress

