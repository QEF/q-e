!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE A_h(npw,e,h,ah)
  !-----------------------------------------------------------------------
  USE kinds, ONLY: DP
  USE cell_base,ONLY : alat, omega, tpiba2
  USE uspp,     ONLY : vkb, nkb
  USE lsda_mod, ONLY : current_spin, nspin
  USE wvfct, ONLY: nbnd, npwx, g2kin
  USE wavefunctions_module,  ONLY: evc, psic
  USE scf,      ONLY : vrs, rho
  USE fft_base, ONLY : dffts, dfftp
  USE fft_interfaces, ONLY : fwfft, invfft
  USE gvect,    ONLY : gstart, nl, nlm, ngm, g, gg
  USE constants,  ONLY: degspin, e2, fpi
  USE becmod, ONLY: bec_type, becp, calbec
  USE cgcom
  USE funct, ONLY: dft_is_gradient
  !
  IMPLICIT NONE
  INTEGER :: npw, j, jkb, ibnd, na,nt,ih
  real(DP) :: e(nbnd)
  COMPLEX(DP) :: h(npwx,nbnd), ah(npwx,nbnd)
  !
  COMPLEX(DP) :: fp, fm
  COMPLEX(DP), POINTER :: dpsic(:), drhoc(:), dvxc(:)
  real(DP), POINTER  :: dv(:), drho(:)
  !
  CALL start_clock('a_h')
  !
  drho  => auxr
  dpsic => aux2
  drhoc => aux3
  !
  drho(:) = 0.d0
  !
  ! [(k+G)^2 - e ]psi
  DO ibnd = 1,nbnd
     ! set to zero the imaginary part of h at G=0
     ! needed for numerical stability
     IF (gstart==2) h(1,ibnd) = cmplx( dble(h(1,ibnd)),0.d0,kind=DP)
     DO j = 1,npw
        ah(j,ibnd) = (g2kin(j)-e(ibnd)) * h(j,ibnd)
     ENDDO
  ENDDO
  !     V_Loc psi
  DO ibnd = 1,nbnd, 2
     dpsic(:)= (0.d0, 0.d0)
     psic(:) = (0.d0, 0.d0)
     IF (ibnd<nbnd) THEN
        ! two ffts at the same time
        DO j = 1,npw
           psic (nl (j)) = evc(j,ibnd) + (0.0d0,1.d0)* evc(j,ibnd+1)
           dpsic(nl (j)) =   h(j,ibnd) + (0.0d0,1.d0)*   h(j,ibnd+1)
           psic (nlm(j))= conjg(evc(j,ibnd)-(0.0d0,1.d0)* evc(j,ibnd+1))
           dpsic(nlm(j))= conjg(  h(j,ibnd)-(0.0d0,1.d0)*   h(j,ibnd+1))
        ENDDO
     ELSE
        DO j = 1,npw
           psic (nl (j)) = evc(j,ibnd)
           dpsic(nl (j)) =   h(j,ibnd)
           psic (nlm(j)) = conjg( evc(j,ibnd))
           dpsic(nlm(j)) = conjg(   h(j,ibnd))
        ENDDO
     ENDIF
     CALL invfft ('Wave', psic, dffts)
     CALL invfft ('Wave',dpsic, dffts)
     DO j = 1,dfftp%nnr
        drho(j) = drho(j) - 2.0d0*degspin/omega *   &
              dble(psic(j)*conjg(dpsic(j)))
        dpsic(j) = dpsic(j) * vrs(j,current_spin)
     ENDDO
     CALL fwfft ('Wave',dpsic, dffts)
     IF (ibnd<nbnd) THEN
        ! two ffts at the same time
        DO j = 1,npw
           fp = (dpsic (nl(j)) + dpsic (nlm(j)))*0.5d0
           fm = (dpsic (nl(j)) - dpsic (nlm(j)))*0.5d0
           ah(j,ibnd  ) = ah(j,ibnd)  +cmplx( dble(fp), aimag(fm),kind=DP)
           ah(j,ibnd+1) = ah(j,ibnd+1)+cmplx(aimag(fp),- dble(fm),kind=DP)
        ENDDO
     ELSE
        DO j = 1,npw
           ah(j,ibnd) = ah(j,ibnd)  + dpsic (nl(j))
        ENDDO
     ENDIF
  ENDDO
  !
  NULLIFY(dpsic)
  ! V_NL psi
  CALL calbec ( npw, vkb, h, becp)
  IF (nkb > 0) CALL add_vuspsi (npwx, npw, nbnd, ah)
  !
  DO j = 1,dfftp%nnr
     drhoc(j) = cmplx(drho(j),0.d0,kind=DP)
  ENDDO
  CALL fwfft ('Dense', drhoc, dfftp)
  !
  ! drho is deltarho(r), drhoc is deltarho(g)
  !
  !  mu'(n(r)) psi(r) delta psi(r)
  !
  dvxc  => aux2
  DO j = 1,dfftp%nnr
     dvxc(j) = drho(j)*dmuxc(j)
  ENDDO
  !
  !  add gradient correction contribution (if any)
  !
  CALL start_clock('dgradcorr')
  IF (dft_is_gradient() ) CALL dgradcor1  &
       (rho%of_r, grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s,            &
        drho, drhoc, dfftp%nnr, nspin, nl, nlm, ngm, g, alat, omega, dvxc)
  CALL stop_clock('dgradcorr')
  NULLIFY (drho)
  !
  !  1/|r-r'| * psi(r') delta psi(r')
  !
  ! gstart is the first nonzero G vector (needed for parallel execution)
  !
  IF (gstart==2) drhoc(nl(1)) = 0.d0
  !
  DO j = gstart,ngm
     drhoc(nl (j)) = e2*fpi*drhoc(nl(j))/ (tpiba2*gg(j))
     drhoc(nlm(j)) = conjg(drhoc(nl (j)))
  ENDDO
  CALL invfft ('Dense', drhoc, dfftp)
  !
  ! drhoc now contains deltaV_hartree
  !
  dv => auxr
  DO j = 1,dfftp%nnr
     dv(j) = -  dble(dvxc(j)) - dble(drhoc(j))
  ENDDO
  !
  CALL vloc_psi_gamma(npwx, npw, nbnd, evc, dv, ah)
  !
  ! set to zero the imaginary part of ah at G=0
  ! needed for numerical stability
  IF (gstart==2) THEN
     DO ibnd = 1, nbnd
        ah(1,ibnd) = cmplx( dble(ah(1,ibnd)),0.d0,kind=DP)
     ENDDO
  ENDIF
  !
  CALL stop_clock('a_h')
  !
  RETURN
END SUBROUTINE A_h
