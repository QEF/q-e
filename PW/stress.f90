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
  USE constants,     ONLY : uakbar
  USE ener,          ONLY : etxc, vtxc
  USE gvect,         ONLY : ngm, gstart, nl, g, gg, gcutm
  USE fft_base,      ONLY : dfftp
  USE ldaU,          ONLY : lda_plus_u
  USE lsda_mod,      ONLY : nspin
  USE scf,           ONLY : rho, rho_core, rhog_core
  USE control_flags, ONLY : iverbosity, gamma_only, llondon
  USE noncollin_module, ONLY : noncolin
  USE funct,         ONLY : dft_is_meta, dft_is_gradient 
  USE symme,         ONLY : symmatrix
  USE bp,            ONLY : lelfield
  USE uspp,          ONLY : okvan
  USE london_module, ONLY : stres_london
#ifdef EXX
  USE exx,           ONLY : exx_stress
  USE funct,         ONLY : dft_is_hybrid, exx_is_active
#endif
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: sigma(3,3)
  !
  real(DP) :: sigmakin (3, 3), sigmaloc (3, 3), sigmahar (3, 3), &
       sigmaxc (3, 3), sigmaxcc (3, 3), sigmaewa (3, 3), sigmanlc (3, 3), &
       sigmabare (3, 3), sigmah (3, 3), sigmael( 3, 3), sigmaion(3, 3), &
       sigmalon ( 3 , 3 ), sigma_nonloc_dft (3 ,3)
#ifdef EXX
  real(DP) :: sigmaexx(3,3)
#endif
  integer :: l, m
  !
  WRITE( stdout, '(//5x,"entering subroutine stress ..."/)')

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
                        dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nnr, nl, ngm, g, alat, omega, sigmaxc)
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
  !
  sigmah(:,:) = 0.d0
  if (lda_plus_u) call stres_hub(sigmah)


  !
  !   Electric field contribution
  !
  sigmael(:,:)=0.d0
  sigmaion(:,:)=0.d0
  !the following is for calculating the improper stress tensor
!  call stress_bp_efield (sigmael )
!  call stress_ion_efield (sigmaion )


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
               sigmaion(:,:) + sigmalon(:,:) + sigma_nonloc_dft(:,:)
#ifdef EXX
  if (dft_is_hybrid()) then
    sigmaexx = exx_stress()
    CALL symmatrix ( sigmaexx )
    sigma(:,:) = sigma(:,:) + sigmaexx(:,:)
  else
    sigmaexx = 0.d0
  endif
#endif
  ! Resymmetrize the total stress. This should not be strictly necessary,
  ! but prevents loss of symmetry in long vc-bfgs runs

  CALL symmatrix ( sigma )
  !
  ! write results in Ryd/(a.u.)^3 and in kbar
  !
  WRITE( stdout, 9000) (sigma(1,1) + sigma(2,2) + sigma(3,3)) * uakbar / 3d0,  &
                  (sigma(l,1), sigma(l,2), sigma(l,3),                    &
            sigma(l,1)*uakbar, sigma(l,2)*uakbar, sigma(l,3)*uakbar, l=1,3)

  if (iverbosity.ge.1) WRITE( stdout, 9005) &
     (sigmakin(l,1)*uakbar,sigmakin(l,2)*uakbar,sigmakin(l,3)*uakbar, l=1,3),&
     (sigmaloc(l,1)*uakbar,sigmaloc(l,2)*uakbar,sigmaloc(l,3)*uakbar, l=1,3),&
     (sigmanlc(l,1)*uakbar,sigmanlc(l,2)*uakbar,sigmanlc(l,3)*uakbar, l=1,3),&
     (sigmahar(l,1)*uakbar,sigmahar(l,2)*uakbar,sigmahar(l,3)*uakbar, l=1,3),&
     (sigmaxc (l,1)*uakbar,sigmaxc (l,2)*uakbar,sigmaxc (l,3)*uakbar, l=1,3),&
     (sigmaxcc(l,1)*uakbar,sigmaxcc(l,2)*uakbar,sigmaxcc(l,3)*uakbar, l=1,3),&
     (sigmaewa(l,1)*uakbar,sigmaewa(l,2)*uakbar,sigmaewa(l,3)*uakbar, l=1,3),&
     (sigmah  (l,1)*uakbar,sigmah  (l,2)*uakbar,sigmah  (l,3)*uakbar, l=1,3),&
     (sigmalon(l,1)*uakbar,sigmalon(l,2)*uakbar,sigmalon(l,3)*uakbar, l=1,3), &
     (sigma_nonloc_dft(l,1)*uakbar,sigma_nonloc_dft(l,2)*uakbar,sigma_nonloc_dft(l,3)*uakbar, l=1,3)

#ifdef EXX
  if (iverbosity.ge.1) WRITE( stdout, 9006) &
     (sigmaexx(l,1)*uakbar,sigmaexx(l,2)*uakbar,sigmaexx(l,3)*uakbar, l=1,3)
9006 format (5x,'EXX     stress (kbar)',3f10.2/2(26x,3f10.2/)/ )
#endif

  if(lelfield .and. iverbosity >= 1) then
     write(stdout,*) "Stress tensor electronic el field part:"
     write(stdout,*) (sigmael(l,1),sigmael(l,2),sigmael(l,3), l=1,3)
     write(stdout,*) "Stress tensor electronic el field part:"
     write(stdout,*) (sigmaion(l,1),sigmaion(l,2),sigmaion(l,3), l=1,3)
  endif

  call stop_clock ('stress')

  return
9000 format (10x,'total   stress  (Ry/bohr**3) ',18x,'(kbar)', &
             &5x,'P=',f8.2/3 (3f13.8,4x,3f10.2/))
9001 format (5x,'Isostatic pressure: ',f8.2,' kbar')
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
         &   5x,'dft-nl  stress (kbar)',3f10.2/2(26x,3f10.2/)/ )
end subroutine stress

